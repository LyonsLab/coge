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
($USER) = CoGe::Accessory::LogUser->get_user();
$DB = new CoGe::Genome;

my $pj = new CGI::Ajax(
		       get_data_info => \&get_data_info,
		       get_data_info_info => \&get_data_info_info,
		       get_data_info_chr_info => \&get_data_info_chr_info,
		       gen_data => \&gen_data,
		       get_genomic_seq => \&get_genomic_seq,
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
<script src="js/Base.js"> </script>
<script src="js/Async.js"> </script>
<script src="js/Iter.js"> </script>
<script src="js/DOM.js"> </script>
<script src="js/Logging.js"> </script>
<script src="js/Signal.js"> </script>
<script src="js/Style.js"> </script>
<script src="js/util.js"> </script>
<script src="js/panner.js"> </script>
<script src="js/tiler.js"> </script>
<script src=js/tilerConfig.js> </script>
<script src=js/kaj.stable.js> </script>
});
    $template->param(USER=>$USER);
    $template->param(DATE=>$DATE);
    $template->param(LOGO_PNG=>"GeLo-logo.png");
    $template->param(BODY=>$body);
    my $html;
    $html .= $template->output;
    return $html;
  }

sub gen_body
  {
    my $form = shift || $FORM;
    my $chr = $form->param('chr');# || 1;
    my $di = $form->param('di');# || 6;
    my $dio = $DB->get_dataset_obj->retrieve($di);
    my $z = $form->param('z');# || 7;
    my $x = $form->param('x');# || 1;
    my $template = HTML::Template->new(filename=>'/opt/apache/CoGe/tmpl/GeLo.tmpl');
    if ($chr && $di)
      {
	$template->param(CHR=>$chr);
	$template->param(DI=>$di);
	$template->param(Z=>$z);
	$template->param(X=>$x);
	my $org = $dio->org->name." (id".$dio->org->id.")";
	$template->param(ORG=>$org);
	my $ds = $dio->name."(v".$dio->version.", id".$dio->id.")";
	$template->param(DS=>$ds);
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
    my @opts = map {"<OPTION value=\"$_->id\">".$_->name."</OPTION>"} sort {$a->name cmp $b->name} $DB->get_org_obj->retrieve_all();
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
    my @opts = map {"<OPTION value=\"".$_->id."\">".$_->name. " (v".$_->version.", id",$_->id,")</OPTION>"} sort {$b->version cmp $a->version || $a->name cmp $b->name} $org->data_info;
    my $html;
    $html .= qq{<FONT CLASS ="small">Dataset count: }.scalar @opts.qq{</FONT>\n<BR>\n};
    $html .= qq{<SELECT id="di_id" SIZE="5" MULTIPLE onChange="gen_data(['args__loading. . .'],['di_info']); get_data_info_info(['di_id'],[dataset_info_chr_chain])" >\n};
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
    #push @chr, "no genomic sequence" unless @chr;
    if (@chr)
      {
	$html .= qq{<tr><td>Chromosome};
	$html .= "s" if scalar @chr > 1;
	$html .= ":";
	$html .= "<td>";
    my $select;
	$select .= qq{<SELECT id="chr" onChange="gen_data(['args__searching for features. . .'],['chr_info']); gen_data(['args__waiting. . .'],['viewer']); get_data_info_chr_info(['di_id', 'chr'],['chr_info','viewer'])" >\n};
	$select .= join ("\n", map {"<OPTION value=\"$_\">".$_."</OPTION>"} sort {$a cmp $b} @chr)."\n";
	$select =~ s/OPTION/OPTION SELECTED/;
	$select .= "\n</SELECT>\n";
    $html .= $select;
      }
    else {
      $html .= qq{<input type="hidden" id="chr" value="">};
      $html .= "No genomic sequence";
    }
#    my $length = commify( $DB->get_genomic_seq_obj->get_last_position($di) );
#    $html .= qq{<tr><td>Nucleotides:<td>$length} if $length;
#    my $feats = $di->get_feature_type_count;
#    $html .= qq{<tr><td valign=top>Features:<td><table>};
#    $html .= join ("\n<tr>",map {"<td>$_<td>".$feats->{$_} } sort {$feats->{$b}<=> $feats->{$a}} keys %$feats);
#    $html .= "</table>";
    return $html;
  }

sub get_data_info_chr_info
  {
    my $did = shift;
    my $chr = shift;
    return unless $did;    
    my $di = $DB->get_data_info_obj->retrieve($did);
    return unless $di;
    my $html .= "<table>";;
    my $length = commify( $DB->get_genomic_seq_obj->get_last_position(di=>$di, chr=>$chr) );
    $html .= qq{<tr><td>Nucleotides:<td>$length} if $length;
    my $feats = $di->get_feature_type_count(chr=>$chr);
    $html .= qq{<tr><td valign=top>Features:<td valign=top><table>};
    my $feat_string = join ("\n<tr>",map {"<td>$_<td>".$feats->{$_} } sort {$feats->{$b}<=> $feats->{$a}} keys %$feats);
    $feat_string = "None" unless $feat_string;
    $html .= "$feat_string</table>";

    my $viewer;
    if ($chr)
      {
	$viewer .= "<br><br>";
	$viewer .= "<font class=\"oblique\">GeLo Viewer Launcher</font><br>";
	$viewer .= "<table>";
	$viewer .= "<tr><td class = \"ital\">Starting location: ";
	$viewer .= qq{<td><input type="text" size=10 value="1000" id="x">};
	my $zoom;
   	$zoom .= qq{<tr><td class = "ital">Zoom level:};
   	$zoom .= qq{<td><SELECT name = "z" id="z" size = 1>};
   	my @opts = map {"<OPTION value=\"$_\">".$_."</OPTION>"} (5..15);
   	$zoom .= join ("\n", @opts);
   	$zoom =~ s/OPTION value="10"/OPTION SELECTED value="10"/;
   	$zoom .= qq{</SELECT>};
	$viewer .= qq{<tr><td class = "ital">Zoom level:<td><input type = "text" size=10 value ="10" id = "z">};
	#    $viewer .= $zoom;
	$viewer .= "</table>";
	#$viewer .= qq{<input type="hidden" id="z" value="7">};
	$viewer .= qq{<input type="submit" value = "Launch GeLo Viewer!" onClick="launch_viewer($did, '$chr')">};
      }
    my $seq_grab;
    if ($chr)
      {
	$seq_grab .= "<br><br>";
	$seq_grab .= qq{<font class="oblique">Genomic Sequence Retrieval</font><br>};
	$seq_grab .= qq{<table>};
	$seq_grab .= "<tr><td class = \"ital\">Start position: ";
	$seq_grab .= qq{<td><input type="text" size=10 value="1000" id="start">};
	$seq_grab .= "<tr><td class = \"ital\">End position: ";
	$seq_grab .= qq{<td><input type="text" size=10 value="100000" id="stop">};
	$seq_grab .= qq{</table>};
	$seq_grab .= qq{<input type="submit" value = "Get Genomic Sequence!" onClick="get_genomic_seq(['args__$did', 'args__$chr', 'start', 'stop'], ['gseq'])">};
	$seq_grab .= qq{<div id="gseq"></div>};
      }
    return $html, $viewer, $seq_grab;
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

sub get_genomic_seq
  {
    my $dsid = shift;
    my $chr = shift;
    my $start = shift;
    my $stop = shift;
    my $seq;
    my $ds = $DB->get_dataset_obj->retrieve($dsid);
    $seq .= ">".$ds->org->name.":".$chr.":".$start."-".$stop."\n";
    $seq .= $DB->get_genomic_seq_obj->get_sequence(start=>$start, stop => $stop, chr=>$chr, info_id=>$dsid);
    $seq =~ s/(.{80})/$1\n/g;
    return "<pre>".$seq."</pre>";
  }
1;
