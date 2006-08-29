#! /usr/bin/perl -w
use strict;
use CGI;
use CGI::Carp 'fatalsToBrowser';
use CGI::Ajax;
use CoGe::Accessory::LogUser;
use HTML::Template;
use Data::Dumper;
use CoGe::Genome;
use CoGe::Graphics::Chromosome;
use CoGe::Graphics::Feature;
use CoGe::Graphics::Feature::Gene;
use CoGe::Graphics::Feature::NucTide;
use POSIX;

$ENV{PATH} = "/opt/apache2/CoGe/";

use vars qw( $DATE $DEBUG $TEMPDIR $TEMPURL $USER $DB $FORM $ACCN $FID);

# set this to 1 to print verbose messages to logs
$DEBUG = 0;
$TEMPDIR = "/opt/apache/CoGe/tmp";
$TEMPURL = "/CoGe/tmp";
$| = 1; # turn off buffering
$DATE = sprintf( "%04d-%02d-%02d %02d:%02d:%02d",
		 sub { ($_[5]+1900, $_[4]+1, $_[3]),$_[2],$_[1],$_[0] }->(localtime));

$FORM = new CGI;
$ACCN = $FORM->param('accn');
($USER) = CoGe::Accessory::LogUser->get_user();
$DB = new CoGe::Genome;
my $pj = new CGI::Ajax(
		       db_lookup=>\&db_lookup,
		       source_search=>\&get_data_source_info_for_accn,
		       get_types=>\&get_types,
		       cogesearch=>\&cogesearch,
		       clear_div=>\&clear_div,
		       get_anno=>\&get_anno,
		       show_location=>\&show_location,
		       show_express=>\&show_express,
		       gen_data=>\&gen_data,
		       get_dna_seq_for_feat => \&get_dna_seq_for_feat,
		       get_prot_seq_for_feat => \&get_prot_seq_for_feat,
		      );
$pj->JSDEBUG(0);
$pj->DEBUG(0);
$pj->js_encode_function('escape');
print $pj->build_html($FORM, \&gen_html);
#print gen_html();



sub get_prot_seq_for_feat
  {
    my $featid = shift;
    my ($feat) = $DB->get_feat_obj->search(feature_id=>$featid);
    my ($seq) = $DB->get_protein_sequence_for_feature($feat);
    $seq = "No sequence available" unless $seq;
    return $seq;
  }

sub get_dna_seq_for_feat
  {
    my $featid = shift;
    my ($feat) = $DB->get_feat_obj->search(feature_id=>$featid);
    my $seq = $DB->get_genomic_sequence_for_feature($feat);
    $seq = "No sequence available" unless $seq;
    return $seq;
  }

sub gen_data
  {
    my $message = shift;
    return qq{<font class="loading">$message. . .</font>};
  }

sub clear_div
  {
    return;
  }

sub get_types
  {
    my ($accn, $dsid) = @_;
    my $html;
    my $blank = qq{<input type="hidden" id="Type_name">};
    my %seen;
    my @opts = sort map {"<OPTION>$_</OPTION>"} grep {! $seen{$_}++} map {$_->type->name} $DB->get_features_by_name_and_dataset_id(name=>$accn, id=>$dsid);
    $html .= "<font class=small>Type count: ".scalar @opts."</font>\n<BR>\n";
    $html .= qq{<SELECT id="Type_name" SIZE="5" MULTIPLE onChange="get_anno(['accn_select','Type_name', 'dsid'],['anno'])" >\n};
    $html .= join ("\n", @opts);
    $html .= "\n</SELECT>\n";
    $html =~ s/OPTION/OPTION SELECTED/;
    return $blank unless $html =~ /OPTION/;
    return ($html, 1);
  }


sub cogesearch
  {
    my $accn = shift;
    my $type = shift;
    my $org = shift;
    my $blank = qq{<input type="hidden" id="accn_select">};
#    print STDERR Dumper @_;
    return $blank unless length($accn) > 2 || $type || $org;
    my $html;
    print STDERR "type: $type\n";
    my %seen;
    my @opts = sort map {"<OPTION>$_</OPTION>"} grep {! $seen{$_}++} map {uc($_)} $DB->get_feature_name_obj->power_search(accn=>$accn."%", type=>$type, org=>$org);
    if (@opts > 5000)
      {
	return $blank."Search results over 1000, please refine your search.\n";
      }
    $html .= "<font class=small>Name count: ".scalar @opts."</font>\n<BR>\n";
    $html .= qq{<SELECT id="accn_select" SIZE="5" MULTIPLE onChange="source_search_chain(); " >\n};
    $html .= join ("\n", @opts);
    $html .= "\n</SELECT>\n";
    $html =~ s/<OPTION/<OPTION SELECTED/;
    return $blank."No results found.\n" unless $html =~ /OPTION/;
    return $html;
  }

sub get_anno
  {
    my $accn = shift;
    return unless $accn;
    my $type = shift;
    my $dataset_id = shift;
    my @feats;
    foreach my $feat ($DB->get_features_by_name_and_dataset_id(name=>$accn, id=>$dataset_id))
      {
	push @feats, $feat if ($feat->type->name eq $type);
      }
    my $anno;
    $anno .= "<font class=small>Annotation count: ".scalar @feats."</font>\n<BR>\n" if scalar @feats;
    my $i = 0;
    foreach my $feat (@feats)
      {
	$i++;
	my $chr = $feat->chr;
	my $di = $feat->info->id;
	my $x = $feat->begin_location;
	my $z = 10;
	$anno .= join "\n<BR><HR><BR>\n", $feat->annotation_pretty_print_html();
	$anno .= qq{<DIV id="loc$i"><input type="button" value = "Click for chromosomal view" onClick="window.open('GeLo.pl?chr=$chr&di=$di&INITIAL_CENTER=$x,0&z=$z');"></DIV>};
	$anno .= qq{<DIV id="exp$i"><input type="button" value = "Click for expression tree" onClick="gen_data(['args__Generating expression view image'],['exp$i']);show_express(['args__}.$accn.qq{','args__}.'1'.qq{','args__}.$i.qq{'],['exp$i']);"></DIV>};
	$anno .= qq{<DIV id="dnaseq$i"><input type="button" value = "Click for DNA sequence" onClick="gen_data(['args__retrieving sequence'],['dnaseq$i']);get_dna_seq_for_feat(['args__}.$feat->id.qq{'],['dnaseq$i']);"></DIV>};
	$anno .= qq{<DIV id="protseq$i"><input type="button" value = "Click for protein sequence" onClick="gen_data(['args__retrieving sequence'],['protseq$i']);get_prot_seq_for_feat(['args__}.$feat->id.qq{'],['protseq$i']);"></DIV>\n\n};

	$anno = "<font class=\"annotation\">No annotations for this entry</font>" unless $anno;
      }
    return ($anno);
  }

sub show_location
  {
#    my %opts = @_;
    my $start = shift;#$opts{start};
    my $stop = shift; #$opts{stop};
    my $chr = shift; # $opts{chr};
    my $info_id = shift; # = $opts{info_id};
    my $loc = shift;
   
    my $z = 7;#ceil (log((1000+$stop-$start)/10)/log(2));
#    my $link = qq{<img src="tiler.pl?start=$start&chr=$chr&di=$info_id&iw=2000&z=$z">\n};
    my $html  = qq{
<form action='javascript:;' >
<input id='iw' type ='hidden' value=256 />
<input id='chr' type ='hidden' value=$chr />
<input id='di' type ='hidden' value=$info_id />
<input id='z' type ='hidden' value='$z' />
<input id='x' type ='hidden' value='$start' />
<button id='zi' name='zi' class="tileClass" > Zoom In </button>
<button name='zo' id='zo' class="tileClass" > Zoom Out </button>

<div id='loc' style="width:800px;height:200px;border:1px solid black" > </div>

<div id='container' class="tileClass" style='position:relative;left:6px;top:10px;width:1000px;height:350px;border:1px solid black' > </div>
</form>
};

#    print STDERR $link;
    return $html, $loc;
  }


sub show_express
  {
    my %opts = @_;
    my $accn = shift;
    my $log = shift;
    my $div = shift;
    $log = 1 unless defined $log;
    $accn =~ s/\..*//;
    my $link = qq{<img src="expressiontree.pl?locus=$accn&label=1&rw=80&rh=8&name=1&legend=1&mean=1&log_trans=$log">\n};
    $log = $log ? 0 : 1;
    my $type = $log ? "log transformed" : "normal";
    print STDERR $link;
    $link .= qq{<br><input type="button" value = "Click for $type expression tree" onClick="gen_image([],['exp$div']);show_express(['args__}.$accn.qq{','args__}.$log.qq{','args__}.$div.qq{'],['exp$div']);">};
    return $link;
  }


sub gen_html
  {
    my $template = HTML::Template->new(filename=>'/opt/apache/CoGe/tmpl/generic_page.tmpl');
    $template->param(LOGO_PNG=>"FeatView-logo.png");
    $template->param(TITLE=>'CoGe: Feature Viewer');
    $template->param(USER=>$USER);
    $template->param(DATE=>$DATE);
    $template->param(bOX_NAME=>"Feature Selection");
    my $body = gen_body();
    $template->param(BODY=>$body);
    my $html;
    $html .= $template->output;
    return $html;
  }

sub gen_body
  {
    my $template = HTML::Template->new(filename=>'/opt/apache/CoGe/tmpl/FeatView.tmpl');

    $template->param(ACCN=>$ACCN);
    $template->param(TYPE_LOOP=> [{TYPE=>"<OPTION VALUE=0>All</OPTION>"},map {{TYPE=>"<OPTION value=\"".$_->id."\">".$_->name."</OPTION>"}} sort {$a->name cmp $b->name} $DB->all_feature_types]);
    $template->param(ORG_LOOP=> [{ORG=>"<OPTION VALUE=0>All</OPTION>"},map {{ORG=>"<OPTION value=\"".$_->id."\">".$_->name."</OPTION>"}} sort {$a->name cmp $b->name} $DB->all_orgs]);
    my $html = $template->output;
#    $html =~ s/(>gene<\/OPTION)/ SELECTED$1/; 
    return $html;
  }

sub get_data_source_info_for_accn
  {
    my $accn = shift;
    my $blank = qq{<input type="hidden" id="dsid">};
    return $blank unless $accn;
    my @feats = $DB->get_feats_by_name($accn);
    my %sources;
    foreach my $feat (@feats)
      {
	my $val = $feat->data_info;
	my $name = $val->name;
	my $ver = $val->version;
	my $desc = $val->description;
	my $sname = $val->data_source->name;
	my $ds_name = $val->name;
	my $org = $val->org->name;
	my $title = "$org: $ds_name ($sname, v$ver)";
	$sources{$title} = $val->id;
      }
    my $html;
    $html .= qq{
<SELECT name = "dsid" id="dsid" MULTIPLE SIZE="5" onChange="get_types_chain();">
};
    my $count = 0;
    foreach my $title (reverse sort keys %sources)
      {
	my $id = $sources{$title};
	$html .= qq{  <option value="$id" >$title\n};
	$html =~ s/option/option selected/ unless $count;
	$count++;
      }
    $html .= qq{</SELECT>\n};
    return ("<font class=small>Dataset count: ".$count ."</font>\n<BR>\n".$html, 1);
  }

