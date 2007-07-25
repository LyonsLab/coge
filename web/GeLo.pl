#! /usr/bin/perl -w
use strict;
use CGI;
use CGI::Carp 'fatalsToBrowser';
use CoGe::Accessory::LogUser;
use CoGe::Accessory::Web;
use HTML::Template;
use Data::Dumper;
use CGI::Ajax;
use CoGe::Genome;
use Benchmark;

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
		       get_dataset => \&get_dataset,
		       get_dataset_info => \&get_dataset_info,
		       get_dataset_chr_info => \&get_dataset_chr_info,
		       gen_data => \&gen_data,
		       get_genomic_seq => \&get_genomic_seq,
		       get_orgs => \&get_orgs,
		       get_start_stop=>\&get_start_stop,
		      );
$pj->JSDEBUG(0);
$pj->DEBUG(0);
print $pj->build_html($FORM, \&gen_html);
#print "Content-Type: text/html\n\n";
#print gen_html($FORM);

sub gen_html
  {
    my $html;
    unless ($USER)
      {
	$html = login();
      }
    else
      {
	my ($body, $seq_names, $seqs) = gen_body();
	my $template = HTML::Template->new(filename=>'/opt/apache/CoGe/tmpl/generic_page.tmpl');
	
	
	if ($FORM->param("ds"))
	  {
	    $template->param(HEADER_LINK=>'/CoGe/GeLo.pl');
	  }
	$template->param(TITLE=>'CoGe: Genome Location Viewer');
	$template->param(HEAD=>qq{});
	$template->param(USER=>$USER);
	$template->param(DATE=>$DATE);
	$template->param(LOGO_PNG=>"GeLo-logo.png");
	$template->param(BODY=>$body);
	$html .= $template->output;
      }
    return $html;
  }

sub gen_body
  {
    my $form = shift || $FORM;
    my $template = HTML::Template->new(filename=>'/opt/apache/CoGe/tmpl/GeLo.tmpl');
    $template->param(ORG_LIST=>get_orgs());
    return $template->output;
  }

sub get_orgs
  {
    my $name = shift;
    my @db = $name ? $DB->get_org_obj->search_like({name=>"%".$name."%"}) :$DB->get_org_obj->retrieve_all();
    my %restricted;
    if (!$USER || $USER =~ /public/i)
      {
	$restricted{papaya} = 1;
      }
    my @opts;
    foreach my $item (sort {uc($a->name) cmp uc($b->name)} @db)
      {
	next if $restricted{$item->name};
	push @opts, "<OPTION value=\"".$item->id."\">".$item->name."</OPTION>";
      }
    my $html;
    $html .= qq{<FONT CLASS ="small">Organism count: }.scalar @opts.qq{</FONT>\n<BR>\n};
    unless (@opts) 
      {
	$html .=  qq{<input type = hidden name="org_id" id="org_id">};
	return $html;
      }

    $html .= qq{<SELECT id="org_id" SIZE="5" MULTIPLE onChange="dataset_chain()" >\n};
    $html .= join ("\n", @opts);
    $html .= "\n</SELECT>\n";
    $html =~ s/OPTION/OPTION SELECTED/;
    return $html;
  }

sub get_dataset
  {
    my $oid = shift;
    my $html = "";
    my $org = $DB->get_org_obj->retrieve($oid);
    my @opts = map {"<OPTION value=\"".$_->id."\">".$_->name. " (v".$_->version.", id".$_->id.")</OPTION>"} sort {$b->version cmp $a->version || $a->name cmp $b->name} $org->datasets if $org;
    $html = qq{<FONT CLASS ="small">Dataset count: }.scalar (@opts).qq{</FONT>\n<BR>\n};
    unless (@opts) 
      {
	$html .=  qq{<input type = hidden name="ds_id" id="ds_id">};
	return $html;
      }

    $html .= qq{<SELECT id="ds_id" SIZE="5" MULTIPLE onChange="gen_data(['args__loading. . .'],['ds_info']); dataset_info_chain()" >\n};
    $html .= join ("\n", @opts);
    $html .= "\n</SELECT>\n";
    $html =~ s/OPTION/OPTION SELECTED/;
    return $html;
  }

sub get_dataset_info
  {
    my $dsd = shift;
    my $html = "";
    unless ($dsd) # error flag for empty dataset
    {
    	$html .= qq{<input type="hidden" id="chr" value="">};
    	$html .= "No genomic sequence";
    	return $html;
    }

    my $ds = $DB->get_dataset_obj->retrieve($dsd);

    return $html unless $ds;
    $html = "<table>";
    my $ds_name = $ds->data_source->name ." ". $ds->data_source->description;
    my $link = $ds->data_source->link;

    $link = "http://".$link if ($link && $link !~ /http/);
    $ds_name = "<a href =\"".$link."\">".$ds_name."</a>" if $ds->data_source->link;
    $html .= qq{<TR><TD>Data Source:<TD>$ds_name}."\n";
    $html .= qq{<tr><td>Version:<td>}.$ds->version."\n";

    my %chr;
    map{$chr{$_}++} ($ds->chromosomes, $DB->get_genomic_seq_obj->get_chromosome_for_dataset($ds));
    my @chr = sort keys %chr;
    if (@chr)
      {
	$html .= qq{<tr><td>Chromosome};
	$html .= "s". " (".scalar @chr.")" if scalar @chr > 1;
	$html .= ":";
	$html .= "<td>";
    my $select;
	$select .= qq{<SELECT id="chr" onChange="gen_data(['args__searching for features. . .'],['chr_info']); gen_data(['args__waiting. . .'],['viewer']); gen_data(['args__'],['get_seq']); get_dataset_chr_info(['ds_id', 'chr'],['chr_info','viewer', 'get_seq'])" >\n};
	$select .= join ("\n", map {"<OPTION value=\"$_\">".$_."</OPTION>"} @chr)."\n";
	$select =~ s/OPTION/OPTION SELECTED/;
	$select .= "\n</SELECT>\n";
    $html .= $select;
      }
    else {
      $html .= qq{<input type="hidden" id="chr" value="">};
      $html .= "<tr><td>No genomic sequence";
    }
    return $html;
  }

sub get_dataset_chr_info
  {
    my $dsd = shift;
    my $chr = shift;
	unless ($dsd) # error flag for empty dataset
	{
		return "", "", "";
	}
    my $start = "'start'";
    my $stop = "'stop'";	
    my $html .= "<table>";
    return $html unless $dsd;
    my $ds = $DB->get_dataset_obj->retrieve($dsd);
    return $html unless $ds;
    my $length = commify( $DB->get_genomic_seq_obj->get_last_position(ds=>$ds, chr=>$chr) );
    $html .= qq{<tr><td>Nucleotides:<td>$length} if $length;
    my $feats = $ds->get_feature_type_count(chr=>$chr);
    $html .= qq{<tr><td valign=top>Features:<td valign=top><table>};
    my $feat_string = join ("\n<tr valing=top>",map {"<td valign=top>$_<tdvalign=top>".$feats->{$_} } sort {$feats->{$b}<=> $feats->{$a}} keys %$feats);
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
   	$zoom =~ s/OPTION value="7"/OPTION SELECTED value="7"/;
   	$zoom .= qq{</SELECT>};
	$viewer .= qq{<tr><td class = "ital">Zoom level:<td><input type = "text" size=10 value ="10" id = "z">};
	#    $viewer .= $zoom;
	$viewer .= "</table>";
	#$viewer .= qq{<input type="hidden" id="z" value="7">};
	$viewer .= qq{<input type="submit" value = "Launch GeLo Viewer!" onClick="launch_viewer($dsd, '$chr')">};
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
	$seq_grab .= qq{<input type="submit" value = "Get Genomic Sequence!" onClick="launch_seqview($dsd, $chr)">};
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


  sub get_organism_name
  {
  	my $blah = $_[0];
    print $blah;

	my $html .= "<table>";
	$html .= "</tr></table>";
#alert(this.value)
    return $html;
  }

1;
