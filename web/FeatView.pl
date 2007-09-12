#! /usr/bin/perl -w
use strict;
use CGI;
use CGI::Carp 'fatalsToBrowser';
use CGI::Ajax;
use CoGe::Accessory::LogUser;
use CoGe::Accessory::Web;
use CoGe::Accessory::Restricted_orgs;
use HTML::Template;
use Data::Dumper;
use CoGe::Graphics::Chromosome;
use CoGe::Graphics::Feature;
use CoGe::Graphics::Feature::Gene;
use CoGe::Graphics::Feature::NucTide;
use CoGeX;
use POSIX;

$ENV{PATH} = "/opt/apache2/CoGe/";

use vars qw( $DATE $DEBUG $TEMPDIR $TEMPURL $USER $FORM $ACCN $FID $coge);

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
my $connstr = 'dbi:mysql:dbname=genomes;host=biocon;port=3306';
$coge = CoGeX->connect($connstr, 'cnssys', 'CnS' );

my $pj = new CGI::Ajax(
		       db_lookup=>\&db_lookup,
		       source_search=>\&get_data_source_info_for_accn,
		       get_types=>\&get_types,
		       cogesearch=>\&cogesearch,
		       get_anno=>\&get_anno,
		       show_express=>\&show_express,
		       gen_data=>\&gen_data,
		      );
$pj->JSDEBUG(0);
$pj->DEBUG(0);
$pj->js_encode_function('escape');
print $pj->build_html($FORM, \&gen_html);
#print gen_html();



sub gen_data
  {
    my $message = shift;
    return qq{<font class="loading">$message. . .</font>};
  }

sub get_types
  {
    my ($accn, $dsid) = @_;
    my $html;
    my $blank = qq{<input type="hidden" id="Type_name">};
    my %seen;
    my @opts = sort map {"<OPTION>$_</OPTION>"} grep {! $seen{$_}++} map {$_->type->name} $coge->resultset('Feature')->search({
															       'feature_names.name'=>$accn,
															       dataset_id=>$dsid
															      },{join=>'feature_names'});

    $html .= "<font class=small>Type count: ".scalar @opts."</font>\n<BR>\n";
    $html .= qq{<SELECT id="Type_name" SIZE="5" MULTIPLE onChange="get_anno(['accn_select','Type_name', 'dsid'],[show_anno])" >\n};
    $html .= join ("\n", @opts);
    $html .= "\n</SELECT>\n";
    $html =~ s/OPTION/OPTION SELECTED/;
    return $blank unless $html =~ /OPTION/;
    return ($html, 1);
  }


sub cogesearch
  {
    my %opts = @_;
    my $accn = $opts{accn};
    my $type = $opts{type};
    my $org = $opts{org};
    my $blank = qq{<input type="hidden" id="accn_select">};
#    print STDERR "cogesearch: $accn\n";
#    print STDERR Dumper @_;
    return $blank unless length($accn) > 2 || $type || $org;
    ($USER) = CoGe::Accessory::LogUser->get_user();
    my $restricted_orgs = restricted_orgs(user=>$USER);
    my $html;
    my %seen;
    my @opts;
    my $search = {'me.name'=>{like=>$accn."%"}};
    $search->{feature_type_id}=$type if $type;
    $search->{organism_id}=$org if $org;
    foreach my $name ($coge->resultset('FeatureName')->search(
							      $search,
							      {join=>{'feature'=>'dataset'},
							       prefetch=>'feature',
							       order_by=>'name ASC',
							      },
							     ))
      {
	my $item = $name->name;
	next if $seen{uc($item)};
	if (%{$restricted_orgs})
	  {
	    next if $restricted_orgs->{$name->feature->dataset->organism->name};
	  }
	$seen{uc($item)}++;
	
	push @opts, "<OPTION>$item</OPTION>"
      }
    if (@opts > 5000)
      {
	return $blank."Search results over 5000, please refine your search.\n";
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
    foreach my $feat ($coge->resultset('Feature')->search({
							   'feature_names.name'=>$accn,
							   dataset_id=>$dataset_id
							  },{join=>'feature_names'}))
      {
	push @feats, $feat if ($feat->type->name eq $type);
      }
    my $anno;
    $anno .= "<font class=small>Annotation count: ".scalar @feats."</font>\n<BR>\n" if scalar @feats;
    my $i = 0;
    foreach my $feat (@feats)
      {
	$i++;
	my $featid = $feat->id;
	my $chr = $feat->chr;
	my $rc = 0;
	my $pro = 0;
	my $ds = $feat->dataset->id;
	my $x = $feat->start;
	my $z = 10;
	$anno .= join "\n<BR><HR><BR>\n", $feat->annotation_pretty_print_html(loc_link=>0);
	$anno .= qq{<DIV id="loc$i"><input type="button" value = "Click for Genome view" onClick="window.open('GenomeView.pl?chr=$chr&ds=$ds&x=$x&z=$z');"></DIV>};
#	$anno .= qq{<DIV id="exp$i"><input type="button" value = "Click for expression tree" onClick="gen_data(['args__Generating expression view image'],['exp$i']);show_express(['args__}.$accn.qq{','args__}.'1'.qq{','args__}.$i.qq{'],['exp$i']);"></DIV>};
	$anno .= qq{<DIV id="dnaseq$i"><input type="button" value = "Click for Sequence" onClick="window.open('SeqView.pl?featid=$featid&dsid=$ds&chr=$chr&featname=$accn');"></DIV>};
	$anno = "<font class=\"annotation\">No annotations for this entry</font>" unless $anno;
      }
    return ($anno);
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
#    print STDERR $link;
    $link .= qq{<br><input type="button" value = "Click for $type expression tree" onClick="gen_image([],['exp$div']);show_express(['args__}.$accn.qq{','args__}.$log.qq{','args__}.$div.qq{'],['exp$div']);">};
    return $link;
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
	my $template = HTML::Template->new(filename=>'/opt/apache/CoGe/tmpl/generic_page.tmpl');
	$template->param(LOGO_PNG=>"FeatView-logo.png");
	$template->param(TITLE=>'Feature Viewer');
	$template->param(USER=>$USER);
	$template->param(DATE=>$DATE);
	$template->param(bOX_NAME=>"Feature Selection");
	my $body = gen_body();
	$template->param(BODY=>$body);
	
	$html .= $template->output;
      }
    return $html;
  }

sub gen_body
  {
    my $template = HTML::Template->new(filename=>'/opt/apache/CoGe/tmpl/FeatView.tmpl');

    $template->param(ACCN=>$ACCN);
    $template->param(TYPE_LOOP=> [{TYPE=>"<OPTION VALUE=0>All</OPTION>"},map {{TYPE=>"<OPTION value=\"".$_->id."\">".$_->name."</OPTION>"}} sort {uc($a->name) cmp uc($b->name)} $coge->resultset('FeatureType')->all]);
    my @orgs;
    ($USER) = CoGe::Accessory::LogUser->get_user();
    my $restricted_orgs = restricted_orgs(user=>$USER);
    foreach my $org ($coge->resultset('Organism')->all)
      {
	push @orgs, $org unless $restricted_orgs->{$org->name};
      }
    $template->param(ORG_LOOP=> [{ORG=>"<OPTION VALUE=0>All</OPTION>"},map {{ORG=>"<OPTION value=\"".$_->id."\">".$_->name."</OPTION>"}} sort {uc($a->name) cmp uc($b->name)} @orgs]);
    my $html = $template->output;
#    $html =~ s/(>gene<\/OPTION)/ SELECTED$1/; 
    return $html;
  }

sub get_data_source_info_for_accn
  {
    my $accn = shift;
    my $blank = qq{<input type="hidden" id="dsid">};
    return $blank unless $accn;
    my @feats = $coge->resultset('Feature')->search({'feature_names.name'=>$accn},{join=>'feature_names'});
    my %sources;
    ($USER) = CoGe::Accessory::LogUser->get_user();
    my $restricted_orgs = restricted_orgs(user=>$USER);
    foreach my $feat (@feats)
      {
	my $val = $feat->dataset;
	unless ($val)
	  {
	    warn "error with feature: ".$feat->id ." for name $accn\n";
	    next;
	  }
	my $name = $val->name;
	my $ver = $val->version;
	my $desc = $val->description;
	my $sname = $val->datasource->name if $val->datasource;
	my $ds_name = $val->name;
	my $org = $val->organism->name if $val->organism;
	my $title = "$org: $ds_name ($sname, v$ver)";
	next if $restricted_orgs->{$org};
	$sources{$title} = $val->id;
      }
    my $html;
    $html .= qq{
<SELECT name = "dsid" id="dsid" MULTIPLE SIZE="5" onChange="get_types_chain();">
};
    my $count = 0;
    foreach my $title (sort {$b cmp $a} keys %sources)
      {
	my $id = $sources{$title};
	$html .= qq{  <option value="$id" >$title\n};
	$html =~ s/option/option selected/ unless $count;
	$count++;
      }
    $html .= qq{</SELECT>\n};
    return ("<font class=small>Dataset count: ".$count ."</font>\n<BR>\n".$html, 1);
  }

