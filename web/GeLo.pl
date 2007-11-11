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
		       gen_gc_for_chromosome=> \&gen_gc_for_chromosome,
		       gen_gc_for_chromosome_and_type =>\&gen_gc_for_chromosome_and_type,
		       get_codon_usage_for_chromosome=>\&get_codon_usage_for_chromosome,
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
    my $name = $form->param('org_name');
    $name = "Search" unless $name;
    $template->param(ORG_NAME=>$name);
    $name = "" if $name =~ /Search/;
    $template->param(ORG_LIST=>get_orgs($name));
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
	push @opts, "<OPTION value=\"".$item->id."\">".$item->name." (id".$item->id.")</OPTION>";
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
    my %opts = @_;
    my $oid = $opts{oid};
    my $dsname = $opts{dsname};
    my $html; 
    my $html2;
    my @opts;
    if ($oid)
      {
	my $org = $coge->resultset("Organism")->find($oid);
	@opts = map {"<OPTION value=\"".$_->id."\">".$_->name. " (v".$_->version.", id".$_->id.")</OPTION>"} sort {$b->version cmp $a->version || $a->name cmp $b->name} $org->datasets if $org;

	}
    elsif ($dsname)
      {
	my @ds = $coge->resultset("Dataset")->search({name=>{like=>"%".$dsname."%"}});
	($USER) = CoGe::Accessory::LogUser->get_user();
	my $restricted_orgs = restricted_orgs(user=>$USER);
	my @orgs;
	foreach my $item (sort {uc($a->name) cmp uc($b->name)} @ds)
	  {
	    next if $restricted_orgs->{$item->organism->name};
	    push @opts, "<OPTION value=\"".$item->id."\">".$item->name."(v".$item->version.", id".$item->id.")</OPTION>";
	    push @orgs, "<OPTION value=\"".$item->organism->id."\">".$item->organism->name." (id".$item->organism->id.")</OPTION>";
	  }
	$html2 .= qq{<FONT CLASS ="small">Organism count: }.scalar @orgs.qq{</FONT>\n<BR>\n};
	$html2 .= qq{<SELECT id="org_id" SIZE="5" MULTIPLE onChange="dataset_chain()" >\n};
	$html2 .= join ("\n", @orgs);
	$html2 .= "\n</SELECT>\n";
	$html2 =~ s/OPTION/OPTION SELECTED/;
      }
    if (@opts) 
      {
	$html = qq{<FONT CLASS ="small">Dataset count: }.scalar (@opts).qq{</FONT>\n<BR>\n};
	$html .= qq{<SELECT id="ds_id" SIZE="5" MULTIPLE onChange="gen_data(['args__loading. . .'],['ds_info']); dataset_info_chain()" >\n};
	$html .= join ("\n", @opts);
	$html .= "\n</SELECT>\n";
	$html =~ s/OPTION/OPTION SELECTED/;
      }
    else
      {
	$html .=  qq{<input type = hidden name="ds_id" id="ds_id">};
      }
    return $html, $html2;
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
    $html .= qq{<tr><td>Nucleotides:<td>$length<td><div id=chromosome_gc class="link" onclick="\$('#chromosome_gc').removeClass('link'); gen_gc_for_chromosome(['args__dsid','ds_id','args__chr','chr'],['chromosome_gc']);">Click for percent GC content</div>} if $length;
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
   AND l.chromosome = "$chr"
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
    my $feat_string .= qq{<tr><td valign=top>Features:<td valign=top><span class=small>(click feature name for percent gc)</span><table>};
    $feat_string .= join ("\n<tr valing=top>",map {"<td valign=top><div id=$_ class=\"link\" onclick=\" \$('#$_').removeClass('link'); gen_data(['args__loading. . .'],['".$_."_type']); gen_gc_for_chromosome_and_type(['args__dsid','ds_id','args__chr','chr','args__type','args__$_'],['".$_."_type'])\">".$_."</div><td valign=top> ".$feats->{$_}."<td><div id=".$_."_type></div>" } sort {($feats->{$b})<=>($feats->{$a})} keys %$feats);
    $feat_string .= "</table>";
    if ($feats->{CDS})
      {
	$feat_string .= "<div class=link id=codon_usage onclick=\" \$('#codon_usage').removeClass('link'); gen_data(['args__loading. . .'],['codon_usage']); get_codon_usage_for_chromosome(['args__dsid','ds_id','args__chr','chr'],['codon_usage'])\">"."Click for codon usage"."</div>";
      }
    $feat_string .= "None" unless keys %$feats;
    return $feat_string;
  }

sub gen_data
  {
    my $message = shift;
    return qq{<font class="loading">$message. . .</font>};
  }

sub gen_gc_for_chromosome_and_type
  {
    my %args = @_;
    my $dsid = $args{dsid};
    my $chr = $args{chr};
    my $type = $args{type};
    return unless $dsid;
    my $ds = $coge->resultset('Dataset')->find($dsid);
    my $gc = 0;
    my $at = 0;
    my $count = 0;
    foreach my $feat ($ds->features({"feature_type.name"=>$type},{join=>"feature_type"}))
      {
	my @gc = $feat->gc_content(counts=>1);
	$gc+=$gc[0] if $gc[0] =~ /^\d+$/;
	$at+=$gc[1] if $gc[1] =~ /^\d+$/;
      }
    my $total = $gc+$at;
    return "error" unless $total;
    return "GC: ".sprintf("%.2f",100*$gc/($total))."%  AT: ".sprintf("%.2f",100*$at/($total))."%";
  }

sub gen_gc_for_chromosome
  {
    my %args = @_;
    my $dsid = $args{dsid};
    my $chr = $args{chr};
    return unless $dsid;
    my $ds = $coge->resultset('Dataset')->find($dsid);
    my $gc=$ds->percent_gc(chr=>$chr);
    return "GC: ".(100*$gc)."%  AT: ".(100*(1-$gc))."%";
  }

sub get_codon_usage_for_chromosome
  {
    my %args = @_;
    my $dsid = $args{dsid};
    my $chr = $args{chr};
    return unless $dsid;
    my $ds = $coge->resultset('Dataset')->find($dsid);
    my %codons;
    my ($code, $code_type);
    my $codon_total = 0;
    my %aa;
    my $aa_total=0;
    my $feat_count = 0;
    foreach my $feat ($ds->features({"feature_type.name"=>"CDS"},{join=>"feature_type"}))
      {
	$feat_count++;
#	print STDERR "~$feat_count~";
	($code, $code_type) = $feat->genetic_code() unless $code;
	my ($codon) = $feat->codon_frequency(counts=>1);
	grep {$codon_total+=$_} values %$codon;
	grep {$codons{$_}+=$codon->{$_}} keys %$codon;
	foreach (keys %$codon)
	  {
	    next if length ($_) eq 3;
#	    print STDERR join ("\t", $_, $feat->id, $feat->names),"\n"
	  }
	foreach my $tri (keys %$code)
	  {
	    $aa{$code->{$tri}}+=$codon->{$tri};
	    $aa_total+=$codon->{$tri};
	  }
	print STDERR ".($feat_count)" if !$feat_count%10;
      }
#    foreach my $tri (keys %code)
#      {
#	$codons{$tri} = 0 unless $codons{$tri};
#      }
    my $count = 0;
    my $html = "Codon Usage: $code_type<table>";
     my ($max_codon) = sort {$b<=>$a} values %codons;
    $max_codon = sprintf("%.2f",100*$max_codon/$codon_total);
    foreach (sort { sort_nt2(substr($a, 0, 1)) <=> sort_nt2(substr($b,0, 1)) || sort_nt2(substr($a,1,1)) <=> sort_nt2(substr($b,1,1)) || sort_nt(substr($a,2,1)) <=> sort_nt(substr($b,2,1)) } keys %$code)
      {
	my $current_val = sprintf("%.2f",100*$codons{$_}/$codon_total);
	 my $color = color_by_usage($max_codon, $current_val);
	my $str = "<tr style=\"background-color: rgb($color,255,$color)\" ><td>".$_."(".$code->{$_}.")<td>".$codons{$_}."<td>(".$current_val."%)";
	delete $codons{$_};
	$html .= "</table><tr><td><table>" unless $count;
	if ($count%4)
	  {
	  }
	else
	  {
	    $html .= "</table><td nospan><table>"; 
	  }

	$html .= $str;
	$count++;
	$count = 0 if $count == 16;
      }
     $count = 0;
     foreach (sort keys %codons)
       {
	 my $current_val = sprintf("%.2f",100*$codons{$_}/$codon_total);
	 my $color = color_by_usage($max_codon, $current_val);
 	my $str = "<tr style=\"background-color: rgb($color,255,$color)\" ><td>".$_."<td>".$codons{$_}."<td>(".$current_val."%)";
	$html .= "</table><tr><td><table>" unless $count;
 	if ($count%4)
 	  {
 	  }
 	else
 	  {
 	    $html .= "</table><td nospan><table>"; 
 	  }
 	$html .= $str;
 	$count++;
 	$count = 0 if $count == 16;
      }
    $html .="</table></table>";
     my ($max_aa) = sort {$b<=>$a} values %aa;
    $max_aa = sprintf("%.2f",100*$max_aa/$aa_total);
    $html .= "Predicted amino acid usage using $code_type<table>";
    my %aa_sort = map {$_,{}} keys %aa;
    foreach my $codon (keys %$code)
      {
	my $gc = $codon =~ tr/GCgc/GCgc/;
	my $at = $codon =~ tr/ATat/ATat/;
	$aa_sort{$code->{$codon}}{AT}+=$at;
	$aa_sort{$code->{$codon}}{GC}+=$gc;
      }
    %aa_sort = map {$_,($aa_sort{$_}{GC}/($aa_sort{$_}{AT}+$aa_sort{$_}{GC}))} keys %aa_sort;

    foreach (sort {$aa_sort{$b} <=> $aa_sort{$a} }keys %aa)
      {	
	my $current_val = sprintf("%.2f",100*$aa{$_}/$aa_total);
	my $color = color_by_usage($max_aa,$current_val);
	$html .= "<tr style=\"background-color: rgb($color,255,$color)\"><td>$_ (GC:".sprintf("%.0f",100*$aa_sort{$_})."%)<td>".$aa{$_}."<td>(".$current_val."%)";
  }
    $html .= "</table>";
    return $html;
  }


sub sort_nt
  {
    my $chr = uc(shift);

    $chr = substr($chr, -1,1) if length($chr)>1;
    my $val = 0;
    if ($chr eq "C")
      {
	$val = 1;
      }
    elsif ($chr eq "G")
      {
	$val = 2;
      }
    elsif ($chr eq "A")
      {
	$val = 3;
      }
    return $val;
  }

sub sort_nt2
  {
    my $chr = uc(shift);

    $chr = substr($chr, -1,1) if length($chr)>1;
    my $val = 0;
    if ($chr eq "A")
      {
	$val = 1;
      }
    elsif ($chr eq "G")
      {
	$val = 2;
      }
    elsif ($chr eq "C")
      {
	$val = 3;
      }
    return $val;
  }

sub commify
    {
      my $text = reverse $_[0];
      $text =~ s/(\d\d\d)(?=\d)(?!\d*\.)/$1,/g;
      return scalar reverse $text;
    }
sub color_by_usage
	{
		my ($max,$value) = @_;
		my $g = 255*(($max - $value) / $max);
		return int($g + .5);
	}


1;
