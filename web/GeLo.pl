#! /usr/bin/perl -w
use strict;
use CGI;
use CGI::Carp 'fatalsToBrowser';
use CoGe::Accessory::LogUser;
use HTML::Template;
use Data::Dumper;
use CGI::Ajax;
use CoGe::Genome;

$ENV{PATH} = "/opt/apache2/CoGe/";

use vars qw( $DATE $DEBUG $TEMPDIR $TEMPURL $USER $FORM $NEATO $DOT $DB);

# set this to 1 to print verbose messages to logs
$DEBUG = 0;
$TEMPDIR = "/opt/apache/CoGe/tmp";
$TEMPURL = "/CoGe/tmp";
$NEATO = "/usr/bin/neato";
$DOT = "/usr/bin/dot";
$| = 1; # turn off buffering
$DATE = sprintf( "%04d-%02d-%02d %02d:%02d:%02d",
		 sub { ($_[5]+1900, $_[4]+1, $_[3]),$_[2],$_[1],$_[0] }->(localtime));

$FORM = new CGI;
$USER = CoGe::Accessory::LogUser->get_user();
$DB = new CoGe::Genome;

my $pj = new CGI::Ajax(
		       get_data_info => \&get_data_info,
		       get_data_info_info => \&get_data_info_info,
		       gen_data => \&gen_data,
		      );
$pj->JSDEBUG(0);
$pj->DEBUG(0);
print $pj->build_html($FORM, \&gen_html);
#print "Content-Type: text/html\n\n";
#print gen_html($FORM);

sub gen_html
  {
    my ($body, $seq_names, $seqs) = gen_body();
    my $template = HTML::Template->new(filename=>'/opt/apache/CoGe/tmpl/generic_page.tmpl');

    $template->param(TITLE=>'CoGe: Genome Location Viewer');
    $template->param(HEAD=>qq{
<link rel="stylesheet" type="text/css" href="css/tiler.css">
<script src=js/Dom.js> </script>
<script src=js/common.js> </script>
<script src=js/extendEvent.js> </script>
<script src=js/panner.js> </script>
<script src=js/tiler.js> </script>
<script src=js/tilerConfig.js> </script>
<script src=js/kaj.stable.js> </script>
});
    $template->param(USER=>$USER);
    $template->param(DATE=>$DATE);
    $template->param(LOGO_PNG=>"GeLo-logo.png");
    $template->param(BODY=>$body);
#    $template->param(HEAD=>'<SCRIPT language="JavaScript" type="text/javascript" src="./js/kaj.stable.js"></SCRIPT>');
    my $html;
    $html .= $template->output;
    return $html;
  }

sub gen_body
  {
    my $form = shift || $FORM;
    my $chr = $form->param('chr');# || 1;
    my $di = $form->param('di');# || 6;
    my $z = $form->param('z');# || 7;
    my $x = $form->param('x');# || 1;
    my $template = HTML::Template->new(filename=>'/opt/apache/CoGe/tmpl/GeLo.tmpl');
    if ($chr && $di)
      {
	$template->param(CHR=>$chr);
	$template->param(DI=>$di);
	$template->param(Z=>$z);
	$template->param(X=>$x);
      }
    else
      {
	$template->param(FRONT_PAGE => 1);
	$template->param(ORG_LIST=>get_orgs());
      }
    return $template->output;
  }

sub get_orgs
  {
    my $oid = shift;
    my @opts = sort map {"<OPTION value=\"$_->id\">".$_->name."</OPTION>"} $DB->get_org_obj->retrieve_all();
    my $html;
    $html .= qq{<FONT CLASS ="small">Organism count: }.scalar @opts.qq{</FONT>\n<BR>\n};
    $html .= qq{<SELECT id="org_id" SIZE="5" MULTIPLE onChange="dataset_chain()" >\n};
    $html .= join ("\n", @opts);
    $html .= "\n</SELECT>\n";
    $html =~ s/OPTION/OPTION SELECTED/;
    return $html;
  }

sub get_data_info
  {
    my $oid = shift;
    return unless $oid;
    my $org = $DB->get_org_obj->retrieve($oid);
    return unless $org;
    my @opts = map {"<OPTION value=\"".$_->id."\">".$_->name. " (v".$_->version.")</OPTION>"} sort {$b->version cmp $a->version || $a->name cmp $b->name} $org->data_info;
    my $html;
    $html .= qq{<FONT CLASS ="small">Dataset count: }.scalar @opts.qq{</FONT>\n<BR>\n};
    $html .= qq{<SELECT id="di_id" SIZE="5" MULTIPLE onChange="gen_data(['args__loading. . .'],['di_info']); get_data_info_info(['di_id'],['di_info','viewer'])" >\n};
    $html .= join ("\n", @opts);
    $html .= "\n</SELECT>\n";
    $html =~ s/OPTION/OPTION SELECTED/;
    return $html;
  }

sub get_data_info_info
  {
    my $did = shift;
    return unless $did;
    
    my $di = $DB->get_data_info_obj->retrieve($did);
    return unless $di;
    my $html .= "<table>";;
    my $ds = $di->data_source->name ." ". $di->data_source->description;
    my $link = $di->data_source->link;
    $link = "http://".$link if ($link && $link !~ /http/);
    $ds = "<a href =\"".$link."\">".$ds."</a>" if $di->data_source->link;
    $html .= qq{<TR><TD>Data Source:<TD>$ds}."\n";
    $html .= qq{<tr><td>Version:<td>}.$di->version."\n";
    my @chr = $DB->get_genomic_seq_obj->get_chromosome_for_data_information($di);
    my $viewer;

    if (@chr)
      {
#	$viewer .= "<form>";
	$viewer .= qq{<input type="hidden" id="di" value="$did">};
	$viewer .= qq{<input type="hidden" id="chr" value="$chr[0]">};
	$viewer .= qq{<input type="hidden" id="z" value="7">};
	$viewer .= "<br><br>";
	$viewer .= "<font class=\"oblique\">GeLo Viewer Launcher</font><br>";
	$viewer .= "<table>";
	$viewer .= "<tr><td class = \"ital\">Starting location: ";
	$viewer .= qq{<td><input type="text" size=10 value="1000" id="x">};
#	$viewer .= qq{<tr><td class = \"ital\">Zoom level:};
#	$viewer .= qq{<td><select id="z">};
#	my @opts = map {"<OPTION>".$_."</OPTION>"} (5..10);
#	$viewer .= join ("\n", @opts);
#	$viewer =~ s/OPTION/OPTION SELECTED/;
#	$viewer .= qq{</select>};
	$viewer .= "</table>";
	$viewer .= qq{<input type="submit" value = "Launch!" onClick="launch_viewer($did, $chr[0])">};
#	$viewer .= "</form>";
      }

    push @chr, "no genomic sequence" unless @chr;
    $html .= qq{<tr><td>Chromosome};
    $html .= "s" if scalar @chr > 1;
    $html .= ":";
    $html .= "<td>".join (" ", sort {$a cmp $b} @chr)."\n";
    my $length = commify( $DB->get_genomic_seq_obj->get_last_position($di) );
    $html .= qq{<tr><td>Nucleotides:<td>$length} if defined $length;
    my $feats = $di->get_feature_type_count;
    $html .= qq{<tr><td valign=top>Features:<td><table>};
    $html .= join ("\n<tr>",map {"<td>$_<td>".$feats->{$_} } sort {$feats->{$b}<=> $feats->{$a}} keys %$feats);
    $html .= "</table>";
    return $html, $viewer;
  }

sub gen_data
  {
    my $message = shift;
    return qq{<font class="loading">$message. . .</font>};
  }

sub commify 
    {
      my $text = reverse $_[0];
      $text =~ s/(\d\d\d)(?=\d)(?!\d*\.)/$1,/g;
      return scalar reverse $text;
    }
1;
