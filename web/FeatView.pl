#! /usr/bin/perl -w
use strict;
use CGI;
use CGI::Carp 'fatalsToBrowser';
use CGI::Ajax;
use GS::LogUser;
use HTML::Template;
use Data::Dumper;
use CoGe::Genome;

$ENV{PATH} = "/opt/apache2/CoGe/";

use vars qw( $DATE $DEBUG $TEMPDIR $TEMPURL $USER $DB $FORM $ACCN);

# set this to 1 to print verbose messages to logs
$DEBUG = 0;
$TEMPDIR = "/opt/apache/CoGe/tmp";
$TEMPURL = "/CoGe/tmp";
$| = 1; # turn off buffering
$DATE = sprintf( "%04d-%02d-%02d %02d:%02d:%02d",
		 sub { ($_[5]+1900, $_[4]+1, $_[3]),$_[2],$_[1],$_[0] }->(localtime));

$FORM = new CGI;
$ACCN = $FORM->param('accn');
$USER = GS::LogUser->get_user();
$DB = new CoGe::Genome;
my $pj = new CGI::Ajax(
		       db_lookup=>\&db_lookup,
		       source_search=>\&get_data_source_info_for_accn,
		       get_types=>\&get_types,
		       accn_search=>\&accn_search,
		       clear_div=>\&clear_div,
		       get_anno=>\&get_anno,
		      );
$pj->JSDEBUG(0);
$pj->DEBUG(0);
print $pj->build_html($FORM, \&gen_html);


sub clear_div
  {
    return;
  }

sub get_types
  {
    my ($accn, $info_ver) = @_;
    my $html;
    my $blank = qq{<input type="hidden" id="Type_name">};
    my %seen;
    my @opts = sort map {"<OPTION>$_</OPTION>"} grep {! $seen{$_}++} map {$_->type->name} $DB->get_features_by_name_and_version(name=>$accn, version=>$info_ver);
    $html .= "<font class=small>Type count: ".scalar @opts."</font>\n<BR>\n";
    $html .= qq{<SELECT id="Type_name" SIZE="5" MULTIPLE onChange="get_anno(['accn_select','infoid','Type_name'],['anno'])" >\n};
    $html .= join ("\n", @opts);
    $html .= "\n</SELECT>\n";
    $html =~ s/OPTION/OPTION SELECTED/;
    return $blank unless $html =~ /OPTION/;
    return ($html, 1);
  }

sub accn_search
  {
    my $accn = shift;

    my $blank = qq{<input type="hidden" id="accn_select">};
    return $blank unless length($accn) > 2;
    my $html;

    my %seen;
    my @opts = sort map {"<OPTION>$_</OPTION>"} grep {! $seen{$_}++} map {uc($_)} $DB->get_feature_name_obj->search_name($accn."%");
    $html .= "<font class=small>Name count: ".scalar @opts."</font>\n<BR>\n";
    $html .= qq{<SELECT id="accn_select" SIZE="5" MULTIPLE onChange="source_search_chain(); " >\n};
    $html .= join ("\n", @opts);
    $html .= "\n</SELECT>\n";
    $html =~ s/<OPTION/<OPTION SELECTED/;
    return $blank."Accession not found.\n" unless $html =~ /OPTION/;
    return $html;
  }

sub get_anno
  {
    my $accn = shift;
    return unless $accn;
    my $info_ver = shift;
    my $type = shift;
    my @feats;
    foreach my $feat ($DB->get_features_by_name_and_version(name=>$accn, version=>$info_ver))
      {
	push @feats, $feat if ($feat->type->name eq $type);
      }
    my $anno;
    $anno .= "<font class=small>Annotation count: ".scalar @feats."</font>\n<BR>\n" if scalar @feats;
    $anno .= join "\n<BR><HR><BR>\n", map {$_->annotation_pretty_print_html} @feats;
    $anno = "<font class=\"annotation\">No annotations for this entry</font>" unless $anno;
    return $anno;
  }

sub gen_html
  {
    my $template = HTML::Template->new(filename=>'/opt/apache/CoGe/tmpl/FeatView.tmpl');
    $template->param(TITLE=>'Genome Feature Viewer');
    $template->param(USER=>$USER);
    $template->param(DATE=>$DATE);
    $template->param(TABLE_NAME1=>"Feature Selection");
    $template->param(ACCN=>$ACCN);	
#    $template->param(ADDIIONAL_JS=>'accn_search_chain(0);');
    if ($ACCN)
      {

      }
    my $html;# =  "Content-Type: text/html\n\n";
    $html .= $template->output;
    return $html;
  }

sub get_data_source_info_for_accn
  {
    my $accn = shift;
    my $blank = qq{<input type="hidden" id="infoid">};
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
	my $org = $val->org->name;
	my $title = "$org: $sname, version $ver";
	$sources{$title} = $ver
      }
    my $html;
    $html .= qq{
<SELECT name = "infoid" id="infoid" MULTIPLE SIZE="5" onChange="get_types_chain();">
};
    my $count = 0;
    foreach my $title (reverse sort keys %sources)
      {
	my $ver = $sources{$title};
	$html .= qq{  <option value="$ver" >$title\n};
	$html =~ s/option/option selected/ unless $count;
	$count++;
      }
    $html .= qq{</SELECT>\n};
    return ("<font class=small>Version count: ".$count ."</font>\n<BR>\n".$html, 1);
  }

sub get_data_source_info_for_accn_old
  {
    my $accn = shift;
    my $blank = qq{<input type="hidden" id="infoid">};
    return $blank unless $accn;
    my @feats = $DB->get_feats_by_name($accn);
    my %sources;
    foreach my $feat (@feats)
      {
	$sources{$feat->data_info->id} = $feat->data_info;
      }
    my $html;
    $html .= "<font class=small>Data count: ".scalar keys (%sources) ."</font>\n<BR>\n";
    $html .= qq{
<SELECT name = "infoid" id="infoid" MULTIPLE SIZE="5" onChange="get_types_chain();">
};
#<SELECT name = "infoid" id="infoid" MULTIPLE SIZE="5" onChange="get_types(['accn_select', 'infoid'],['FeatType']);">
    my $count = 0;
    foreach my $id (sort {$b <=> $a} keys %sources)
      {
	my $val = $sources{$id};
	my $name = $val->name;
	my $ver = $val->version;
	my $desc = $val->description;
	my $sname = $val->data_source->name;
	my $org = $val->org->name;
	$html .= qq{  <option value="$id" >$org: $sname, version $ver\n};
	$html =~ s/option/option selected/ unless $count;
	$count++;
      }
    $html .= qq{</SELECT>\n};
    return ($html);
  }

