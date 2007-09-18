#! /usr/bin/perl -w
use strict;
use CGI;
use CGI::Carp 'fatalsToBrowser';
use CoGe::Accessory::LogUser;
use CoGe::Accessory::Web;
use CoGe::Accessory::Restricted_orgs;
use HTML::Template;
use Data::Dumper;
use CGI::Ajax;
use CoGeX;
use Benchmark;

$ENV{PATH} = "/opt/apache2/CoGe/";

use vars qw( $DATE $DEBUG $TEMPDIR $TEMPURL $USER $FORM $NEATO $DOT $coge $connstr);

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
$connstr = 'dbi:mysql:dbname=genomes;host=biocon;port=3306';
$coge = CoGeX->connect($connstr, 'cnssys', 'CnS' );
#$coge->storage->debugobj(new DBIxProfiler());
#$coge->storage->debug(1);

my $pj = new CGI::Ajax(
		       get_dataset => \&get_dataset,
		       get_dataset_info => \&get_dataset_info,
		       get_dataset_chr_info => \&get_dataset_chr_info,
		       gen_data => \&gen_data,
		       get_orgs => \&get_orgs,
		       get_start_stop=>\&get_start_stop,
		       get_feature_counts => \&get_feature_counts,
		      );
$pj->JSDEBUG(0);
$pj->DEBUG(0);
print $pj->build_html($FORM, \&gen_html);
#print "Content-Type: text/html\n\n";print gen_html($FORM);

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
    my @db = $name ? $coge->resultset("Organism")->search({name=>{like=>"%".$name."%"}}) : $coge->resultset("Organism")->all;
    ($USER) = CoGe::Accessory::LogUser->get_user();
    my $restricted_orgs = restricted_orgs(user=>$USER);
    my @opts;
    foreach my $item (sort {uc($a->name) cmp uc($b->name)} @db)
      {
	next if $restricted_orgs->{$item->name};
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
    my $org = $coge->resultset("Organism")->find($oid);
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

    my $ds = $coge->resultset("Dataset")->find($dsd);

    return $html unless $ds;
    $html = "<table>";
    my $dataset = $ds->name.": ". $ds->description;
#    $dataset .= " <a href= ".$ds->link.">(".$ds->link.")" if $ds->link;
    my $source_name = $ds->datasource->name .": ". $ds->datasource->description;
    my $link = $ds->datasource->link;

    $link = "http://".$link if ($link && $link !~ /http/);
    $source_name = "<a href =\"".$link."\">".$source_name."</a>" if $ds->datasource->link;
    $html .= qq{<tr><td>Name: <td>$dataset}."\n";
    $html .= qq{<TR><TD>Data Source: <TD>$source_name}."\n";
    $html .= qq{<tr><td>Version: <td>}.$ds->version."\n";
    $html .= qq{<tr><td>Date deposited: <td>}.$ds->date."\n";
    my $org = $ds->organism->name;
    $org .= ": ".$ds->organism->description if $ds->organism->description;
    $org .="\n" ;
    $html .= qq{<tr><td>Organism: <td>$org\n};
    my %chr;
#    map{$chr{$_}++} ($ds->chromosomes, $DB->get_genomic_seq_obj->get_chromosome_for_dataset($ds));
    map{$chr{$_}++} ($ds->get_chromosomes);
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
    my $ds = $coge->resultset("Dataset")->find($dsd);
    return $html unless $ds;
    my $length = commify( $ds->last_chromosome_position($chr) );
    $html .= qq{<tr><td>Nucleotides:<td>$length} if $length;
#    my $feat_string = get_feature_counts($dsd, $chr);
    my $feat_string = qq{
<div id=feature_count onclick="gen_data(['args__loading. . .'],['feature_count_data']);\$('#feature_count').hide(0);get_feature_counts(['ds_id', 'chr'],['feature_count_data']);" >Click here for feature count information</div><div id=feature_count_data></div>
};
    $html .= "</table>$feat_string";

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
	$viewer .= qq{<tr><td class = "ital">Zoom level:<td><input type = "text" size=10 value ="8" id = "z">};
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

sub get_feature_counts
  {
    my $dsd = shift;
    my $chr = shift;
    my $query = qq{
SELECT count(distinct(feature_id)), ft.name 
  FROM feature
  JOIN feature_type ft using (feature_type_id)
  JOIN location l USING (feature_id)
 WHERE dataset_id = $dsd
   AND l.chromosome = $chr
 GROUP BY ft.name
};
    my $dbh = DBI->connect($connstr, 'cnssys', 'CnS' );
    my $sth = $dbh->prepare($query);
    $sth->execute;
    my $feats = {};
    while (my $row = $sth->fetchrow_arrayref)
      {
	$feats->{$row->[1]} = $row->[0];
      }
#    my $feats = $ds->get_feature_type_count(chr=>$chr);
    my $feat_string .= qq{<tr><td valign=top>Features:<td valign=top><table>};
    $feat_string .= join ("\n<tr valing=top>",map {"<td valign=top> ".$_."<td valign=top> ".$feats->{$_} } sort {($feats->{$b})<=>($feats->{$a})} keys %$feats);
    $feat_string .= "None" unless $feat_string;
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
