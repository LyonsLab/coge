#! /usr/bin/perl -w
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

use vars qw( $TEMPDIR $USER $DATE $BASEFILE $coge $cogeweb $FORM);

$DATE = sprintf( "%04d-%02d-%02d %02d:%02d:%02d",
		sub { ($_[5]+1900, $_[4]+1, $_[3]),$_[2],$_[1],$_[0] }->(localtime));
($USER) = CoGe::Accessory::LogUser->get_user();
$TEMPDIR = "/opt/apache/CoGe/tmp/";
$FORM = new CGI;

$coge = CoGeX->dbconnect();

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
			);
$pj->js_encode_function('escape');
print $pj->build_html($FORM, \&gen_html);
#print $FORM->header;
#print gen_html();

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
    $template->param(HELP=>'BLAST');
   # print STDERR "user is: ",$USER,"\n";
    my $name = $USER->user_name;
        $name = $USER->first_name if $USER->first_name;
        $name .= " ".$USER->last_name if $USER->first_name && $USER->last_name;
        $template->param(USER=>$name);
    $template->param(LOGO_PNG=>"FeatList-logo.png");
    $template->param(LOGON=>1) unless $USER->user_name eq "public";
    $template->param(DATE=>$DATE);
    $template->param(BOX_NAME=>'CoGe: Blast');
    $template->param(BODY=>$body);
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
    my $feat_list = [];
    $feat_list = read_file() if $BASEFILE;#: $opts{feature_list};
    foreach my $item ($form->param('fid'))
      {
	push @$feat_list, $item if $item =~ /^\d+$/;
      }
    my $dsid = $form->param('dsid') if $form->param('dsid');
    my $chr = $form->param('chr') if $form->param('chr');
    my $ftid = $form->param('ftid') if $form->param('ftid');
    my $start = $form->param('start') if $form->param('start');
    my $stop = $form->param('stop') if $form->param('stop');
    push @$feat_list, @{get_fids_from_dataset(dsid=>$dsid, ftid=>$ftid, chr=>$chr, start=>$start, stop=>$stop)} if $dsid;
    my ($table, $feat_types) = generate_table(feature_list=>$feat_list, ftid=>$ftid);
    $template->param('CDS_COUNT'=>$feat_types->{CDS});
    $template->param('SHOW_ALL_CODON_TABLES'=>$feat_types->{CDS}) if $feat_types->{CDS};
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

sub get_fids_from_dataset
  {
    my %opts = @_;
    my $dsid = $opts{dsid};
    my $ftid = $opts{ftid};
    my $chr = $opts{chr};
    my $start = $opts{start};
    my $stop = $opts{stop};
    my $search = {dataset_id=>$dsid};
    my $join={};
    if ($ftid)
      {
	$search->{feature_type_id}=$ftid;
      }
    if ($chr)
      {
	$search->{chromosome}=$chr;
      }
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
    return unless @$feat_list;
    my @table;
    my %feat_types;
    my $count = 1;
    $feat_list = [map {$coge->resultset("Feature")->find($_)} @$feat_list];
    $feat_list = [sort {$a->organism->name cmp $b->organism->name || $a->type->name cmp $b->type->name || $a->chromosome cmp $b->chromosome|| $a->start <=> $b->start}@$feat_list];
    foreach my $feat(@$feat_list)
    {
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
      my ($name) = $feat->names;
      my $hpp = $feat->annotation_pretty_print_html();
      #my $row_style = $count%2 ? "even" : "odd";
      my $other;
      my $cds_count = $feat_types{CDS};
      $other .= "<div class=link id=codon_usage$cds_count><DIV onclick=\" \$('#codon_usage$cds_count').removeClass('link'); gen_data(['args__loading'],['codon_usage$cds_count']); codon_table(['args__featid','args__$featid'],['codon_usage$cds_count'])\">"."Click for codon usage"."</DIV></DIV><input type=hidden id=CDS$cds_count value=$featid>" if $feat->type->name eq "CDS";
    my ($gc, $at) = $feat->gc_content;
    $at*=100;
    $gc*=100;
    my ($wgc, $wat) = $feat->wobble_content;
    $wat*=100 if $wat;
    $wgc*=100 if $wgc;

      push @table,{
		   FEATID=>$featid,
		   NAME=>$name,
		   TYPE=>$feat->type->name,
		   LOC=>"chr ".$feat->chr.": ".$feat->start."-".$feat->stop." (".$feat->strand.")",
#		   CHR=>$feat->chr,
#		   STRAND=>$feat->strand,
		   ORG=>$feat->organism->name." (v".$feat->version.")",
		   HPP=>$hpp, 
#		   TABLE_ROW=>$row_style,
		   LENGTH=>$feat->length(),
		   OTHER=>$other,
		   AT=>$at,
		   GC=>$gc,
		   WAT=>$wat,
		   WGC=>$wgc,

		  };
      $count++;
    }
   return \@table, \%feat_types;
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
    my $url = "/CoGe/CoGeBlast.pl?featid=$accn_list";
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
		$url .= "featid=$featid&";
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
   	
   	foreach my $featid (split /,/,$accn_list)
	  {
   	   my ($feat) = $coge->resultset("Feature")->find($featid);
   	   next unless $feat;
   	   my ($name) = sort $feat->names;
   	   my $app = $feat->annotation_pretty_print_html();
   	   $app =~ s/(<\/?span(\s*class=\"\w+\")?\s*>)?//ig;
   	   my ($anno) = $app =~ /annotation: (.+)?/i;
   	   ($anno) = split (/<BR/, $anno);
   	   $anno =~ s/;/;\n/g;
	   my ($at, $gc) = $feat->gc_content;
	   $at*=100;
	   $gc*=100;
	   my ($wat, $wgc) = $feat->wobble_content;
	   $wat*=100;
	   $wgc*=100;
	   $worksheet->write($i,0,"http://toxic.berkeley.edu/CoGe/FeatView.pl?accn=$name",$name);
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
	   $worksheet->write($i,12,$feat->genomic_sequence);
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
    my $featid = $args{featid};

    return unless $featid;
    my ($feat) = $coge->resultset('Feature')->find($featid);
    my ($codon, $code_type) = $feat->codon_frequency(counts=>1);
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
    my ($at, $gc) = $feat->gc_content;
    $at*=100;
    $gc*=100;
    my ($wat, $wgc) = $feat->wobble_content;
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

