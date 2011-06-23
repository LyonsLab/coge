#! /usr/bin/perl -w


use strict;
use CGI;
use CoGeX;
use DBI;

use Data::Dumper;
use CoGe::Accessory::LogUser;
use CoGe::Accessory::Web;
use HTML::Template;
use URI::Escape;
use Spreadsheet::WriteExcel;
use Benchmark;
use DBIxProfiler;
no warnings 'redefine';


use vars qw($P $PAGE_NAME $TEMPDIR $USER $DATE $BASEFILE $coge $cogeweb $FORM $URL);
$P = CoGe::Accessory::Web::get_defaults();
$ENV{PATH} = $P->{COGEDIR};
$URL = $P->{URL};
$DATE = sprintf( "%04d-%02d-%02d %02d:%02d:%02d",
		sub { ($_[5]+1900, $_[4]+1, $_[3]),$_[2],$_[1],$_[0] }->(localtime));
$PAGE_NAME = "FeatList.pl";
($USER) = CoGe::Accessory::LogUser->get_user();
$TEMPDIR = $P->{TEMPDIR};
$FORM = new CGI;
$coge = CoGeX->dbconnect();
#$coge->storage->debugobj(new DBIxProfiler());
#$coge->storage->debug(1);



$SIG{'__WARN__'} = sub { }; #silence warnings

my %FUNCTION = (
		gen_html=>\&gen_html,
		get_feature_counts=>\&get_feature_counts,
		get_gc=>\&get_gc,
		get_aa_usage=>\&get_aa_usage,
		get_codon_usage=>\&get_codon_usage,

		       gevo=>\&gevo,
		       blast=>\&blast,
		       get_fasta_seqs=>\&get_fasta_seqs,
		       generate_excel_file=>\&generate_excel_file,
		       protein_table=>\&protein_table,
		       gc_content=>\&gc_content,
		       gen_data=>\&gen_data,
		       send_to_featmap=>\&send_to_featmap,
		       send_to_msa=>\&send_to_msa,
		       send_to_featlist=>\&send_to_featlist,
		       send_to_SynFind=>\&send_to_SynFind,
		       get_anno=>\&get_anno,
		       get_wobble_gc=>\&get_wobble_gc,
		       save_FeatList_settings=>\&save_FeatList_settings,
		       add_to_user_history=>\&add_to_user_history,
		       export_CodeOn=>\&export_CodeOn,
    );
#my $pj = new CGI::Ajax(%FUNCTION);
#$pj->js_encode_function('escape');
#my $t1 = new Benchmark;
if ($FORM->param('jquery_ajax'))
  {
    dispatch();
  }
else
  {
    print $FORM->header,"\n",gen_html();
  }
sub dispatch
{
    my %args = $FORM->Vars;
    my $fname = $args{'fname'};
    if($fname)
    {
	#my %args = $cgi->Vars;
	#print STDERR Dumper \%args;
	if($args{args}){
	    my @args_list = split( /,/, $args{args} );
	    print $FORM->header, $FUNCTION{$fname}->(@args_list);
       	}
	else{
	    print $FORM->header, $FUNCTION{$fname}->(%args);
	}
    }
}
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
       my $template = HTML::Template->new(filename=>$P->{TMPLDIR}.'generic_page.tmpl');
       $template->param(PAGE_TITLE=>'GenomeList');
       $template->param(HELP=>'/wiki/index.php?title=GenomeList');
       # print STDERR "user is: ",$USER,"\n";
       #add_to_user_history() unless $USER->user_name eq "public";
       my $name = $USER->user_name;
       $name = $USER->first_name if $USER->first_name;
       $name .= " ".$USER->last_name if $USER->first_name && $USER->last_name;
       $template->param(USER=>$name);
       $template->param(LOGO_PNG=>"GenomeList-logo.png");
       $template->param(LOGON=>1) unless $USER->user_name eq "public";
       $template->param(DATE=>$DATE);
       my $list_name = $FORM->param('list_name') || $FORM->param('ln');
       my $box_name = "Genome List:";
       $box_name .= " $list_name" if $list_name;
       $template->param(BOX_NAME=>$box_name);
       $template->param(BODY=>$body);
       $template->param(ADJUST_BOX=>1);
       $html .= $template->output;
   }
 }
 
sub gen_body
  {
    my $template = HTML::Template->new(filename=>$P->{TMPLDIR}.'GenomeList.tmpl');
    my $form = $FORM;
    my $no_values;
    $BASEFILE = $form->param('basename');
    my $sort_by_type = $form->param('sort_type');
    my $sort_by_location = $form->param('sort_loc');
    my $prefs =CoGe::Accessory::Web::load_settings(user=>$USER, page=>$PAGE_NAME);
    $prefs = {} unless $prefs;
    $prefs->{display}={
		       NameD=>1,
		       DescD=>1,
		       VerD=>1,
		       ChrCountD=>1,
		       LengthD=>1,
		       GCD=>1,
		       ATD=>1,
		       AAD=>1,
		       ND=>1,
		       FeatD=>1,
		       CodonD=>1,
		       CDSGCHistD=>1,
		       CDSwGCHistD=>1,
		      } unless $prefs->{display};
    foreach my $key (keys %{$prefs->{display}})
      {
	$template->param($key=>"checked");
      }

    $template->param(SAVE_DISPLAY=>1) unless $USER->user_name eq "public";
    #Don't show button to save list data unless valid user
    $template->param(SAVE_DATA=>1) unless $USER->user_name eq "public";

    my $dsgids = [];
    $dsgids = read_file() if $BASEFILE;#: $opts{feature_list};
    foreach my $item ($form->param('dsgid'))
      {
	foreach my $item2 (split/(::)|(,)/, $item)
	  {
	    push @$dsgids, $item2 if $item2 =~ /^\d+_?\d*$/;
	  }
      }
    my ($table, $count) = generate_table(dsgids=>$dsgids);

    $template->param('GENOME_COUNT'=>$count);

    if ($table)
      {
	$template->param(INFO=>$table);
	return $template->output;
      }
    else
      {
	return "No dataset_group (genome) ids were specified.";
      }
  }

sub add_to_user_history
  {
    my %opts = @_;
#    my $url_params = $ENV{'QUERY_STRING'};
    my $url = $opts{url};
    if ($opts{archive})
      {
	$USER->add_to_works({
			     'name'=>$opts{work_name},
			     'archive'=>$opts{archive},
			     'page'=>$PAGE_NAME,
			     'parameter'=>$url,
			     'description'=>$opts{description},
			     'note'=>$opts{note},
			    }); 
      }
    else
      {
	$USER->add_to_works({
			     'name'=>'FeatList-'.$DATE,
			     'archive'=>0,
			     'page'=>$PAGE_NAME,
			     'parameter'=>$url,
			     'description'=>'Feature List created on '.$DATE,
			    });
      }
    return();
  }


sub generate_table
  {
    my %opts = @_;
    my $dsgids = $opts{dsgids};
    return unless @$dsgids;
    my @table;
    my $count = 1;
    foreach my $dsgid (@$dsgids)
      {
	my $dsg = $coge->resultset('DatasetGroup')->find($dsgid);
	my $name = $dsg->name ? $dsg->name : $dsg->organism->name;
	my $desc = $dsg->description ? $dsg->description : join("; ", map { qq{<span class=link onclick=window.open('OrganismView.pl?org_desc=$_')>$_</span>} } split /;\s*/, $dsg->organism->description);
	
	my @chr = $dsg->chromosomes;
	my $chr = scalar @chr;
	my $length = $dsg->length;
	push @table,{
		     COUNT=>$count,
		     DSGID=>$dsgid,
		     NAME=>$name,
		     DESC=>$desc,
		     VER=>$dsg->version,
		     CHR_COUNT=>commify($chr),
		     LENGTH=>commify($length),
		    };
	$count++;
      }
    $count--;
    return \@table, $count;
  }

sub get_feature_counts
  {
    my %opts = @_;
    my $dsgid = $opts{dsgid};
    return "No specified dataset group id" unless $dsgid;
    my $dsg = $coge->resultset('DatasetGroup')->find($dsgid);
    my $query = qq{
SELECT count(distinct(feature_id)), ft.name, ft.feature_type_id
  FROM feature
  JOIN feature_type ft using (feature_type_id)
  JOIN dataset_connector dc using (dataset_id)
 WHERE dataset_group_id = $dsgid
  GROUP BY ft.name

};
    my $coge = CoGeX->dbconnect();
    my $dbh = DBI->connect($coge->db_connection_string,$coge->db_name,$coge->db_passwd);
    my $sth = $dbh->prepare($query);
    $sth->execute;
    my $feats = {};
    while (my $row = $sth->fetchrow_arrayref)
      {
	my $name = $row->[1];
	$name =~ s/\s+/&nbsp/g;
	$feats->{$name} = {count=>$row->[0],
			   id=>$row->[2],
			   name=>$name,
			  };
      }
  my $feat_string .= qq{<table class="small ui-widget-content ui-corner-all">
<thead></thead><tbody>};
  $feat_string .= "<tr valign=top>". join ("\n<tr valign=top>",map {
      "<td valign=top><div id=$_>".$feats->{$_}{name}."</div>".
	    "<td valign=top align=right>".$feats->{$_}{count}
} sort {$a cmp $b} keys %$feats);
  $feat_string .= "</tbody></table>";
  $feat_string .= "None" unless keys %$feats;
  return $feat_string;
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
    my $url = $URL."GEvo.pl?";
    my $count = 1;
    #print STDERR $url,"\n";
    foreach my $featid (split /,/,$accn_list)
      {
		$url .= "fid$count=$featid&";
		$count ++;
      }
    $count--;
    return ("alert",$count) if $count > 20;
    $url .= "num_seqs=$count";
    return $url;
  }
  
sub send_to_SynFind
  {
    my $accn_list = shift;
    $accn_list =~ s/^,//;
    $accn_list =~ s/,$//;
    my ($fid) = split/,/, $accn_list;
    my $url = $URL."SynFind.pl?fid=$fid";
    return $url;
  }
  
  sub send_to_featmap
  {
    my $accn_list = shift;
    $accn_list =~ s/^,//;
    $accn_list =~ s/,$//;
    my $url = $URL."FeatMap.pl?";
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
    my $url = $URL."FeatList.pl?";
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
    my $url = $URL."CoGeAlign.pl?";
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
    my $url = $URL."CoGeBlast.pl?fid=$accn_list";
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

  sub export_CodeOn
    {
      my $accn_list = shift;
      $accn_list =~ s/^,//;
      $accn_list =~ s/,$//;
      my $url = "CodeOn.pl?fid=";
      my @list;
      foreach my $accn (split /,/,$accn_list)
	{
	  next if $accn =~ /no$/;
	  my ($featid, $hspnum, $dsgid) = $accn =~ m/^(\d+)_(\d+)?_?(\d+)?$/;
	  push @list,$featid;
	}
      my %seen = ();
      @list = grep {!$seen{$_}++} @list;
      $url .= join ("::", @list);
      $url =~s/&$//;
      return $url;
    }
  
  
sub generate_excel_file
  {
    my $accn_list = shift;
    $accn_list =~ s/^,//;
    $accn_list =~ s/,$//;
    $cogeweb = CoGe::Accessory::Web::initialize_basefile(prog=>"FeatList");
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
    $worksheet->write(0,12,"Link");
    $worksheet->write(0,13,"Sequence DNA");
    $worksheet->write(0,14,"Sequence Protein");
   	
   	foreach my $item (split /,/,$accn_list)
	  {
	    my ($featid, $gstid) = split /_/, $item;
	    my ($feat) = $coge->resultset("Feature")->find($featid);

   	   next unless $feat;
   	   my ($name) = sort $feat->names;
   	   my $app = $feat->annotation_pretty_print();
   	   $app =~ s/(<\/?span(\s*class=\"\w+\")?\s*>)?//ig;
#   	   my ($anno) = $app =~ /annotation:<td>(.+)?/i;
#   	   ($anno) = split (/<BR/, $anno);
#   	   $anno =~ s/;/;\n/g;
	   my ($at, $gc) = $feat->gc_content;
	   $at*=100;
	   $gc*=100;
	   my ($wat, $wgc) = $feat->wobble_content;
	   $wat*=100;
	   $wgc*=100;
	   $worksheet->write($i,0,$P->{SERVER}."FeatView.pl?accn=$name",$name);
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
	   $worksheet->write($i,11,$app);
	    if (my ($geneid) = $app =~ /geneid.*?(\d+)/i)
	      {
		my $link = "http://www.ncbi.nlm.nih.gov/sites/entrez?db=gene&cmd=Retrieve&dopt=full_report&list_uids=".$geneid;
		$worksheet->write($i,12, $link);
	      }
	    my $seq = $feat->genomic_sequence(gstid=>$gstid);
	   $worksheet->write($i,13,$seq);
	    if ($feat->type->name eq "CDS")
	      {
		$worksheet->write($i,14,$feat->protein_sequence());
	      }
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

sub get_codon_usage
  {
    my %opts = @_;
    my $dsid = $opts{dsid};
    my $chr = $opts{chr};
    my $dsgid = $opts{dsgid};
    my $gstid = $opts{gstid};
    return unless $dsid || $dsgid;

    my $search;
    $search = {"feature_type.name"=>"CDS"};
    $search->{"me.chromosome"}=$chr if defined $chr;

    my @dsids;
    push @dsids, $dsid if $dsid;
    if ($dsgid)
      {
	my $dsg = $coge->resultset('DatasetGroup')->find($dsgid);
	unless ($dsg)
	  {
	    my $error =  "unable to create dsg object using id $dsgid\n";
	    return $error;
	  }
	$gstid = $dsg->type->id;
	foreach my $ds ($dsg->datasets())
	  {
	    push @dsids, $ds->id;
	  }
      }
    my %codons;
    my $codon_total = 0;
    my $feat_count = 0;
    my ($code, $code_type);

    foreach my $dsidt (@dsids)
      {
	my $ds = $coge->resultset('Dataset')->find($dsidt);
	my %seqs; #let's prefetch the sequences with one call to genomic_sequence (slow for many seqs)
	if (defined $chr)
	  {
	    $seqs{$chr} = $ds->genomic_sequence(chr=>$chr, seq_type=>$gstid);
	  }
	else
	  {
	    %seqs= map {$_, $ds->genomic_sequence(chr=>$_, seq_type=>$gstid)} $ds->chromosomes;
	  }
	foreach my $feat ($ds->features($search,{join=>["feature_type",'locations', {'dataset'=>{'dataset_connectors'=>'dataset_group'}}],
						 prefetch=>['locations',{'dataset'=>{'dataset_connectors'=>'dataset_group'}}]}))
	  {
	    my $seq = substr($seqs{$feat->chromosome}, $feat->start-1, $feat->stop-$feat->start+1);
	    $feat->genomic_sequence(seq=>$seq);
	    $feat_count++;
	    ($code, $code_type) = $feat->genetic_code() unless $code;
	    my ($codon) = $feat->codon_frequency(counts=>1);
	    grep {$codon_total+=$_} values %$codon;
	    grep {$codons{$_}+=$codon->{$_}} keys %$codon;
	    print STDERR ".($feat_count)" if !$feat_count%10;
	  }
      }
    %codons = map {$_,$codons{$_}/$codon_total} keys %codons;

    #Josh put some stuff in here so he could get raw numbers instead of percentages for aa usage. He should either make this an option or delete this code when he is done. REMIND HIM ABOUT THIS IF YOU ARE EDITING ORGVIEW!
    my $html = "Codon Usage: $code_type";
    $html .= CoGe::Accessory::genetic_code->html_code_table(data=>\%codons, code=>$code);
    return $html
  }

sub get_aa_usage
  {
    my %opts = @_;
    my $dsgid = $opts{dsgid};
    my $table_name = $opts{table_name};
    return "no dsgid specified" unless $dsgid;
    my $dsg = $coge->resultset('DatasetGroup')->find($dsgid);
    return "Unable able to create DB object for $dsgid" unless $dsg;
    my %codons;
    my $codon_total = 0;
    my %aa;
    my $aa_total=0;
    my $feat_count = 0;
    my ($code, $code_type);
    my $trans_type = $dsg->translation_type;
    foreach my $feat ($dsg->features({feature_type_id=>3}))#,{join=>['locations'],
#					     prefetch=>['locations']}))
      {
	$feat->trans_type($trans_type) if $trans_type;
	$feat_count++;
	($code, $code_type) = $feat->genetic_code() unless $code;
	my $aa = $feat->aa_frequency(counts=>1);
	while (my ($k, $v) = each %$aa)
	  {
	    $aa{$k}+=$v;
	  }
	$trans_type = $feat->trans_type;
      }
    my $html .= $dsg->organism->name."<br>Predicted amino acid usage using $code_type";
    $html .= CoGe::Accessory::genetic_code->html_aa_new(data=>\%aa, counts=>1, table_name=>$table_name);
    return ($html);
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
    my $dsgid = $opts{dsgid};
    return "No dataset group id" unless $dsgid;
    my ($dsg) = $coge->resultset('DatasetGroup')->find($dsgid);
    return "No dsg object for id $dsgid" unless $dsg;
    my ($gc, $at, $n) = $dsg->percent_gc();
    $at*=100;
    $gc*=100;
    $n *=100;
    return ($gc."_".$at."_".$n);
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
   CoGe::Accessory::Web::save_settings(opts=>\%save, user=>$USER, page=>$PAGE_NAME);
  }

sub commify 
    {
      my $text = reverse $_[0];
      $text =~ s/(\d\d\d)(?=\d)(?!\d*\.)/$1,/g;
      return scalar reverse $text;
    }
