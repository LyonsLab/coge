#! /usr/bin/perl -w

#TODO -- FIX FTID SO THAT NO _ AT END OF ID NUM., EX. ID='12345_'

use strict;
use CGI;
use CGI::Ajax;
use CoGeX;
use CoGeX::Feature;
use Data::Dumper;
use CoGe::Accessory::LogUser;
use CoGe::Accessory::Web;
use HTML::Template;
use URI::Escape;
use Spreadsheet::WriteExcel;
use Benchmark;
use DBIxProfiler;


use vars qw($PAGE_NAME $TEMPDIR $USER $DATE $BASEFILE $coge $cogeweb $FORM);
$ENV{PATH}="/opt/apache/CoGe";
$DATE = sprintf( "%04d-%02d-%02d %02d:%02d:%02d",
		sub { ($_[5]+1900, $_[4]+1, $_[3]),$_[2],$_[1],$_[0] }->(localtime));
$PAGE_NAME = "FeatList.pl";
($USER) = CoGe::Accessory::LogUser->get_user();
$TEMPDIR = "/opt/apache/CoGe/tmp/";
$FORM = new CGI;
$coge = CoGeX->dbconnect();
#$coge->storage->debugobj(new DBIxProfiler());
#$coge->storage->debug(1);

my $pj = new CGI::Ajax(
    gen_html=>\&gen_html,
    gevo=>\&gevo,
    blast=>\&blast,
    get_fasta_seqs=>\&get_fasta_seqs,
    generate_excel_file=>\&generate_excel_file,
    codon_table=>\&codon_table,
    protein_table=>\&protein_table,
    gc_content=>\&gc_content,
    gen_data=>\&gen_data,
    send_to_featmap=>\&send_to_featmap,
    send_to_msa=>\&send_to_msa,
    send_to_featlist=>\&send_to_featlist,
    get_anno=>\&get_anno,
    get_gc=>\&get_gc,
    get_wobble_gc=>\&get_wobble_gc,
    save_FeatList_settings=>\&save_FeatList_settings,
    add_to_user_history=>\&add_to_user_history,
    );
$pj->js_encode_function('escape');
#my $t1 = new Benchmark;
print $pj->build_html($FORM, \&gen_html);
#print $FORM->header;print gen_html();
#my $t2 = new Benchmark;
#my $run_time = timestr(timediff($t2,$t1));
#print STDERR qq{
#Runtime:  $run_time
#};

sub gen_html
{
    my $html;    
    unless ($USER)
    {
	$html = login();
      }
    else
     {
       my ($body) = gen_body();
       my $template = HTML::Template->new(filename=>'/opt/apache/CoGe/tmpl/generic_page.tmpl');
       $template->param(TITLE=>'Feature List Viewer');
       $template->param(PAGE_TITLE=>'FeatList');
       $template->param(HELP=>'/wiki/index.php?title=FeatList');
       # print STDERR "user is: ",$USER,"\n";
       #add_to_user_history() unless $USER->user_name eq "public";
       my $name = $USER->user_name;
       $name = $USER->first_name if $USER->first_name;
       $name .= " ".$USER->last_name if $USER->first_name && $USER->last_name;
       $template->param(USER=>$name);
       $template->param(LOGO_PNG=>"FeatList-logo.png");
       $template->param(LOGON=>1) unless $USER->user_name eq "public";
       $template->param(DATE=>$DATE);
       $template->param(BOX_NAME=>'Feature List:');
       $template->param(BODY=>$body);
       $template->param(ADJUST_BOX=>1);
       $html .= $template->output;
   }
 }
 
sub gen_body
  {
    my $template = HTML::Template->new(filename=>'/opt/apache/CoGe/tmpl/FeatList.tmpl');
    my $form = $FORM;
    my $no_values;
    $BASEFILE = $form->param('basename');
    my $sort_by_type = $form->param('sort_type');
    my $sort_by_location = $form->param('sort_loc');
    my $prefs = load_settings(user=>$USER, page=>$PAGE_NAME);
    $prefs = {} unless $prefs;
    $prefs->{display}={
		       FeatNameD=>1,
		       TypeD=>1,
		       ChrD=>1,
		       StartD=>1,
		       StopD=>1,
		       StrandD=>1,
		       LengthD=>1,
		       OrgD=>1,
		       AnnoD=>1,
		      } unless $prefs->{display};
    foreach my $key (keys %{$prefs->{display}})
      {
	$template->param($key=>"checked");
      }

    $template->param(SAVE_DISPLAY=>1) unless $USER->user_name eq "public";
    #Don't show button to save list data unless valid user
    $template->param(SAVE_DATA=>1) unless $USER->user_name eq "public";

    my $feat_list = [];
    $feat_list = read_file() if $BASEFILE;#: $opts{feature_list};
    #feat ids may be in the format of <fid>_<gstid>
    foreach my $item ($form->param('fid'))
      {
	push @$feat_list, $item if $item =~ /^\d+_?\d*$/;
      }
    foreach my $item ($form->param('featid'))
      {
	push @$feat_list, $item if $item =~ /^\d+_?\d*$/;
      }
    my $dsid = $form->param('dsid') if $form->param('dsid');
    my $dsgid = $form->param('dsgid') if $form->param('dsgid');
    my $chr = $form->param('chr') if $form->param('chr');
    my $ftid = $form->param('ftid') if $form->param('ftid');
    my $start = $form->param('start') if $form->param('start');
    my $stop = $form->param('stop') if $form->param('stop');
    my $gstid = $form->param('gstid') if $form->param('gstid'); #genomic_sequence_type_id
    $template->param('GSTID'=>$gstid) if $gstid;
    push @$feat_list, @{get_fids(dsid=>$dsid, dsgid=>$dsgid, ftid=>$ftid, chr=>$chr, start=>$start, stop=>$stop)} if $dsid || $dsgid;
    my ($table, $feat_types, $count) = generate_table(feature_list=>$feat_list, ftid=>$ftid, gstid=>$gstid);
    $template->param('CDS_COUNT'=>$feat_types->{CDS});
    $template->param('FEAT_COUNT'=>$count);
    $template->param('SHOW_ALL_CODON_TABLES'=>$feat_types->{CDS}) if $feat_types->{CDS};

    my $type = qq{<SELECT ID="feature_type">};
    $type .= join ("\n", map {"<OPTION value=$_>".$_."</option>"} sort keys %$feat_types)."\n";
    $type .= "</select>";
    $template->param('FEAT_TYPES'=>$type);
    if ($table)
      {
	$template->param(INFO=>$table);
	return $template->output;
      }
    else
      {
	return "No feature ids were specified.";
      }
  }

sub add_to_user_history
  {
    my %opts = @_;
    if ($opts{archive})
      {
	$USER->add_to_works({
			     'name'=>$opts{work_name},
			     'archive'=>$opts{archive},
	    'page'=>$PAGE_NAME,
	    'parameter'=>$opts{url},
	    'description'=>$opts{description},
	    'note'=>$opts{note},
	    }); 
      }
    else{
      my $url = $ENV{'REQUEST_URI'};
      $USER->add_to_works({
			   'name'=>'FeatList-'.$DATE,
			   'archive'=>0,
			   'page'=>$PAGE_NAME,
			   'parameter'=>$url,
			   'description'=>'Feature List created on '.$DATE,
			  });
      
    }
  }

sub get_fids
  {
    my %opts = @_;
    my $dsid = $opts{dsid};
    my $dsgid = $opts{dsgid};
    my $ftid = $opts{ftid};
    my $chr = $opts{chr};
    my $start = $opts{start};
    my $stop = $opts{stop};
    my $search;
    my @dsids;
    push @dsids, $dsid if $dsid;
    if ($dsgid)
      {
	my $dsg = $coge->resultset('DatasetGroup')->find($dsgid);
	if ($dsg)
	  {
	    foreach my $ds ($dsg->datasets)
	      {
		push @dsids, $ds->id;
	      }
	  }
	else
	  {
	    warn "unable to create dsg object for id $dsgid\n";
	  }
      }
    if (@dsids)
      {
	$search->{-or}=[dataset_id=>[@dsids]];
      }
    $search->{feature_type_id}=$ftid if $ftid;
    $search->{chromosome}=$chr if $chr;
    my $join={};
    my @ids;
    if ($start)
      {
	@ids = map{$_->id} $coge->get_features_in_region(dataset=>$dsid, chr=>$chr, start=>$start, stop=>$stop);
      }
    else
      {
	@ids = map{$_->id} $coge->resultset('Feature')->search($search);
      }

    return \@ids;
  }

sub generate_table
  {
    my %opts = @_;
    my $feat_list = $opts{feature_list};
    my $ftid = $opts{ftid};
    my $gstid = $opts{gstid};
    $gstid = 1 unless defined $gstid;
    return unless @$feat_list;
    my @table;
    my %feat_types;
    my $count = 1;
    my %feats;
    foreach my $item (@$feat_list)
      {
	my ($fid, $gstidt);
	if ($item =~ /_/)
	  {
	    ($fid, $gstidt) = split/_/, $item;
	  }
	else
	  {
	    $fid = $item;
	    $gstidt = $gstid if $gstid;
	  }
	my ($feat) = $coge->resultset("Feature")->search(
							 {'me.feature_id'=>$fid},
							 {
							  join =>['feature_names', 'locations', 'feature_type',{'dataset'=>{'dataset_connectors'=>{dataset_group=>'organism'}}}],
							  prefetch=>['feature_names', 'locations', 'feature_type', {'dataset'=>{'dataset_connectors'=>{dataset_group=>'organism'}}}],
							  }
							);
	next unless $feat;
	$feats{$item} = {fid=>$fid,
			 feat=>$feat,
			 gstid=>$gstidt,
			};
      }
    foreach my $item(sort {$feats{$a}{feat}->organism->name cmp $feats{$b}{feat}->organism->name || $feats{$a}{feat}->type->name cmp $feats{$b}{feat}->type->name || $feats{$a}{feat}->chromosome cmp $feats{$b}{feat}->chromosome|| $feats{$a}{feat}->start <=> $feats{$b}{feat}->start} keys %feats)
    {
      my $feat = $feats{$item}{feat};
      unless ($feat)
	{
#	  warn "feature id $featid failed to return a valid feature object\n";
	  next;
	}
      if ($ftid) 
	{
	  next unless $feat->type->id eq $ftid;
	}
      $feat_types{$feat->type->name}++;
      my $featid = $feat->id;
      $item = $featid."_".$feats{$item}{gstid} unless $item =~/_/;
      my ($name) = $feat->names;
      my $hpp = qq{<div id=anno_$count class="link" onclick="get_anno(['args__fid','args__$item'],['anno_$count'])">}.qq{Get Annotation}.qq{</div>};
      my $other;
      my $cds_count = $feat_types{CDS};
      $other .= "<div class=link id=codon_usage$cds_count><DIV onclick=\" \$('#codon_usage$cds_count').removeClass('link'); codon_table(['args__fid','args__$item'],['codon_usage$cds_count'])\">"."Click for codon usage"."</DIV></DIV><input type=hidden id=CDS$cds_count value=$item>" if $feat->type->name eq "CDS";
      my $gc = qq{Get GC};
      my $at = qq{Get GC};
      my ($wat, $wgc);
      if ($feat->type->name eq "CDS")
	{
	  $wgc = qq{Get GC};
	  $wat = qq{Get GC};
	}
      push @table,{
		   COUNT=>$count,
		   FEATID=>$item,
		   NAME=>$name,
		   TYPE=>$feat->type->name,
		   CHR=>$feat->chr,
		   STRAND=>$feat->strand,
		   START=>$feat->start,
		   STOP=>$feat->stop,
		   ORG=>$feat->organism->name." (v".$feat->version.")",
		   HPP=>$hpp, 
		   LENGTH=>$feat->length(),
		   OTHER=>$other,
		   AT=>$at,
		   GC=>$gc,
		   WAT=>$wat,
		   WGC=>$wgc,

		  };
      $count++;
    }
    $count--;
   return \@table, \%feat_types, $count;
  }

  sub gen_data
  {
    my $message = shift;
    return qq{<font class="loading">$message. . .</font>};
  }

  sub read_file
  {
    my $file = "$TEMPDIR/$BASEFILE.featlist";
    my @featlist;
    unless (-r $file)
      {
	warn "unable to read file $file for feature ids\n";
	return \@featlist;
      }
    open (IN, $file) || die "can't open $file for reading: $!";
    while (<IN>)
      {
	chomp;
	push @featlist,$_;
      }
    close IN;
    return \@featlist;
  }
  
  sub gevo #Send to GEvo
  {
    my $accn_list = shift;
    $accn_list =~ s/^,//;
    $accn_list =~ s/,$//;
    my $url = "/CoGe/GEvo.pl?";
    my $count = 1;
    #print STDERR $url,"\n";
    foreach my $featid (split /,/,$accn_list)
      {
		$url .= "fid$count=$featid&";
		$count ++;
      }
    $count--;
    return ("alert",$count) if $count > 10;
    $url .= "num_seqs=$count";
    return $url;
  }
  
  sub send_to_featmap
  {
    my $accn_list = shift;
    $accn_list =~ s/^,//;
    $accn_list =~ s/,$//;
    my $url = "/CoGe/FeatMap.pl?";
    foreach my $featid (split /,/,$accn_list)
      {
		$url .= "fid=$featid&";
      }
    $url =~ s/&$//;
    return $url;
  }
  sub send_to_featlist
  {
    my $accn_list = shift;
    $accn_list =~ s/^,//;
    $accn_list =~ s/,$//;
    my $url = "/CoGe/FeatList.pl?";
    foreach my $featid (split /,/,$accn_list)
      {
		$url .= "fid=$featid&";
      }
    $url =~ s/&$//;
    return $url;
  }
  
    sub send_to_msa
  {
    my $accn_list = shift;
    $accn_list =~ s/^,//;
    $accn_list =~ s/,$//;
    my $url = "/CoGe/CoGeAlign.pl?";
    foreach my $featid (split /,/,$accn_list)
      {
		$url .= "fid=$featid&";
      }
$url =~ s/&$//;
    return $url;
  }
  
  
  sub blast #send to cogeblast
    {
    my $accn_list = shift;
    $accn_list =~ s/^,//;
    $accn_list =~ s/,$//;
    my $url = "/CoGe/CoGeBlast.pl?fid=$accn_list";
    return $url;
  }
  
  sub get_fasta_seqs
  {
    my $accn_list = shift;
    $accn_list =~ s/^,//;
    $accn_list =~ s/,$//;
    my $url = "FastaView.pl?";
    foreach my $featid (split /,/,$accn_list)
    {
		$url .= "fid=$featid&";
	}
	$url =~s/&$//;
	return $url;
  }
  
sub generate_excel_file
  {
  	my $accn_list = shift;
	$accn_list =~ s/^,//;
	$accn_list =~ s/,$//;
  	$cogeweb = initialize_basefile(prog=>"FeatList");
  	my $basename = $cogeweb->basefile;
  	my ($filename) = $basename =~ /FeatList\/(FeatList_.+)/;
  	my $workbook = Spreadsheet::WriteExcel->new("$TEMPDIR/Excel_$filename.xls");
    $workbook->set_tempdir("$TEMPDIR");
    my $worksheet = $workbook->add_worksheet();
    my $i = 1;
       	 
   	 $worksheet->write(0,0,"Feature Name");
   	 $worksheet->write(0,1,"Type");
   	 $worksheet->write(0,2,"Location");
   	 $worksheet->write(0,3,"Strand");
   	 $worksheet->write(0,4,"Chromosome");
   	 $worksheet->write(0,5,"Length");
   	 $worksheet->write(0,6,"Percent GC");
   	 $worksheet->write(0,7,"Percent AT");
   	 $worksheet->write(0,8,"Percent Wobble GC");
   	 $worksheet->write(0,9,"Percent Wobble AT");
   	 $worksheet->write(0,10,"Organism (version)");
   	 $worksheet->write(0,11,"More information");
   	 $worksheet->write(0,12,"Sequence");
   	
   	foreach my $item (split /,/,$accn_list)
	  {
	    my ($featid, $gstid) = split /_/, $item;
	    my ($feat) = $coge->resultset("Feature")->find($featid);

   	   next unless $feat;
   	   my ($name) = sort $feat->names;
   	   my $app = $feat->annotation_pretty_print_html();
   	   $app =~ s/(<\/?span(\s*class=\"\w+\")?\s*>)?//ig;
   	   my ($anno) = $app =~ /annotation:<td>(.+)?/i;
#	    print STDERR $app,"\n\n", $anno,"\n\n";
   	   ($anno) = split (/<BR/, $anno);
   	   $anno =~ s/;/;\n/g;
	   my ($at, $gc) = $feat->gc_content;
	   $at*=100;
	   $gc*=100;
	   my ($wat, $wgc) = $feat->wobble_content;
	   $wat*=100;
	   $wgc*=100;
	   $worksheet->write($i,0,"http://synteny.cnr.berkeley.edu/CoGe/FeatView.pl?accn=$name",$name);
	   $worksheet->write($i,1,$feat->type->name);
	   $worksheet->write($i,2,$feat->start."-".$feat->stop);
	   $worksheet->write($i,3,$feat->strand);
	   $worksheet->write($i,4,$feat->chr);
	   $worksheet->write($i,5,$feat->length);
	   $worksheet->write($i,6,$gc);
	   $worksheet->write($i,7,$at);
	   $worksheet->write($i,8,$wgc);
	   $worksheet->write($i,9,$wat);
	   $worksheet->write($i,10,$feat->organism->name."(v ".$feat->version.")");
	   $worksheet->write($i,11,$anno);
	    my $seq = $feat->genomic_sequence(gstid=>$gstid);
	   $worksheet->write($i,12,$seq);
	   $i++;
	 };
   	$workbook->close() or die "Error closing file: $!";
   	return "tmp/Excel_$filename.xls";
      }
  
  sub gc_content
  {
    my %args = @_;
    my $featid = $args{featid};
    return unless $featid;
    my ($feat) = $coge->resultset('Feature')->find($featid);
    my ($gc, $at) = $feat->gc_content;
    my $html = "GC:".(100*$gc)."%".", AT:".(100*$at)."%" ;
    return $html;
  }

sub codon_table
  {
    my %args = @_;
    my $featid = $args{fid};
    my $gstid = $args{gstid};
    ($featid, $gstid) = split/_/, $featid if $featid =~/_/;
    return unless $featid;
    my ($feat) = $coge->resultset('Feature')->find($featid);
    my ($codon, $code_type) = $feat->codon_frequency(counts=>1, gstid=>$gstid);
    my %aa;
    my ($code) = $feat->genetic_code;
    my $count = 0;
    foreach my $tri (keys %$code)
      {
	$aa{$code->{$tri}}+=$codon->{$tri};
	$count += $codon->{$tri};
      }
    my $html;
    $html .= "<table><tr valign=top><td>";
    $html .= "Codon Usage: $code_type<br>";
    my ($at, $gc) = $feat->gc_content(gstid=>$gstid);
    $at*=100;
    $gc*=100;
    my ($wat, $wgc) = $feat->wobble_content(gstid=>$gstid);
    $wat*=100;
    $wgc*=100;
    $html .= "Codon Count: $count".", GC: $at% $gc%".", Wobble GC: $wat% $wgc%";
    $html .="<td>";
    $html .= "Predicted amino acid usage";
    $html .= "<tr valign=top><td>";
    $html .= CoGe::Accessory::genetic_code->html_code_table(data=>$codon, code=>$code, counts=>1);
#    $html .= "</div>";
#    $html .= "Predicted amino acid usage for $code_type genetic code:";
    $html .= "<td>";
    $html .= CoGe::Accessory::genetic_code->html_aa(data=>\%aa, counts=>1, split=>1);
    $html .= "</table>";
    return $html;
  }

sub protein_table
  {
    my %args = @_;
    my $featid = $args{featid};
    my ($feat) = $coge->resultset('Feature')->find($featid);
    my $aa = $feat->aa_frequency(counts=>1);
    my $html = "Amino Acid Usage";
    $html .= CoGe::Accessory::genetic_code->html_aa(data=>$aa, counts=>1);
    return $html;
  }

sub get_anno
  {
    my %opts = @_;
    my $fid = $opts{fid};
    my $gstid=$opts{gstid};
    ($fid, $gstid) = split/_/, $fid if ($fid =~ /_/);
    return unless $fid;
    my ($feat) = $coge->resultset('Feature')->find($fid);
    return "No feature for id $fid" unless $feat;
    return $feat->annotation_pretty_print_html(gstid=>$gstid);
  }

sub get_gc
  {
    my %opts = @_;
    my $fid = $opts{fid};
    my $gstid = $opts{gstid};
    return unless $fid;
    ($fid, $gstid) = split/_/, $fid if ($fid =~ /_/);
    my ($feat) = $coge->resultset('Feature')->find($fid);
    return "No feature for id $fid" unless $feat;
    my ($gc, $at, $n) = $feat->gc_content(gstid=>$gstid);
    $at*=100;
    $gc*=100;
    return ($gc, $at);
  }

sub get_wobble_gc
  {
    my %opts = @_;
    my $fid = $opts{fid};
    my $gstid = $opts{gstid};
    ($fid, $gstid) = split/_/, $fid if ($fid =~ /_/);
    return unless $fid;
    my ($feat) = $coge->resultset('Feature')->find($fid);
    return unless $feat->type->name eq "CDS";
    return "No feature for id $fid" unless $feat;
    my ($wgc, $wat) = $feat->wobble_content(gstid=>$gstid);
    $wat*=100;
    $wgc*=100;
    return ($wgc, $wat);
  }


sub save_FeatList_settings
  {
    my %opts = @_;
    my $display = $opts{display};
    my %save;
    if ($display)
      {
	my %settings = (
			1=>'FeatNameD',
			2=>'TypeD',
			3=>'ChrD',
			4=>'StartD',
			5=>'StopD',
			6=>'StrandD',
			7=>'LengthD',
			8=>'GCD',
			9=>'ATD',
			10=>'WGCD',
			11=>'WATD',
			12=>'OrgD',
			13=>'AnnoD',
			14=>'OtherD',
		       );
	foreach my $index (split/,/,$display)
	  {
	    $save{display}{$settings{$index}}=1
	  }
      }
    save_settings(opts=>\%save, user=>$USER, page=>$PAGE_NAME);
  }

