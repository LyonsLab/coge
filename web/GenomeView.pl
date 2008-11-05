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

use vars qw( $DATE $DEBUG $TEMPDIR $TEMPURL $USER $FORM $NEATO $DOT $coge);

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
$coge = CoGeX->dbconnect();
$coge->storage->debugobj(new DBIxProfiler());
$coge->storage->debug(1);

my $pj = new CGI::Ajax(
		       get_dataset => \&get_dataset,
		       get_dataset_info => \&get_dataset_info,
		       get_dataset_chr_info => \&get_dataset_chr_info,
		       gen_data => \&gen_data,
		       get_orgs => \&get_orgs,
		       get_recent_orgs=>\&get_recent_orgs,
		       get_start_stop=>\&get_start_stop,
		       get_feature_counts => \&get_feature_counts,
		       gen_gc_for_chromosome=> \&gen_gc_for_chromosome,
		       gen_gc_for_noncoding=> \&gen_gc_for_noncoding,
		       gen_gc_for_chromosome_and_type =>\&gen_gc_for_chromosome_and_type,
		       get_codon_usage_for_chromosome=>\&get_codon_usage_for_chromosome,
		       get_total_length_for_ds=>\&get_total_length_for_ds,
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
	$template->param(TITLE=>'Organism Overview');
	$template->param(HEAD=>qq{});
	my $name = $USER->user_name;
        $name = $USER->first_name if $USER->first_name;
        $name .= " ".$USER->last_name if $USER->first_name && $USER->last_name;
        $template->param(USER=>$name);

	$template->param(LOGON=>1) unless $USER->user_name eq "public";
	$template->param(DATE=>$DATE);
	$template->param(LOGO_PNG=>"OrganismView-logo.png");
	$template->param(BODY=>$body);
	$html .= $template->output;
      }
    return $html;
  }

sub gen_body
  {
    my $form = shift || $FORM;
    my $template = HTML::Template->new(filename=>'/opt/apache/CoGe/tmpl/GenomeView.tmpl');
    my $name = $form->param('org_name');
    my $desc = $form->param('org_desc');
    my $oid = $form->param('oid');
    my $org = $coge->resultset('Organism')->resolve($oid) if $oid;
    $name = $org->name if $org;
    my $dsname = $form->param('dsname');
    my $dsid = $form->param('dsid');
    $name = "Search" unless $name;
    $template->param(ORG_NAME=>$name) if $name;
    $desc = "Search" unless $desc;
    $template->param(ORG_DESC=>$desc) if $desc;
    $name = "" if $name =~ /Search/;
    $template->param(ORG_LIST=>get_orgs(name=>$name));
    #$template->param(RECENT=>get_recent_orgs());
    my ($ds) = $coge->resultset('Dataset')->resolve($dsid) if $dsid;
    $dsname = $ds->name if $ds;
    $dsname = "Search" unless $dsname;
    $template->param(DS_NAME=>$dsname);
    $dsname = "" if $dsname =~ /Search/;
    my ($dslist,$orglist) = get_dataset(dsname=>$dsname, dsid=>$dsid) if $dsname;
    $template->param(DS_LIST=>$dslist) if $dslist;
    $template->param(ORG_LIST=>$orglist) if $orglist;
    return $template->output;
  }

sub get_recent_orgs
  {
    my %opts = @_;
    my $limit = $opts{limit} || 100;
    my @db = $coge->resultset("Dataset")->search({},
						 {
						  distinct=>"organism.name",
						  join=>"organism",
						  order_by=>"me.date desc",
						  rows=>$limit}
						);
    ($USER) = CoGe::Accessory::LogUser->get_user();
    my $restricted_orgs = restricted_orgs(user=>$USER);
    my @opts;
    my %org_names;
    foreach my $item (@db)
      {
	my $date = $item->date;
	$date =~ s/\s.*//;
	next if $restricted_orgs->{$item->organism->name};
	next if $org_names{$item->organism->name};
	$org_names{$item->organism->name}=1;
	push @opts, "<OPTION value=\"".$item->organism->id."\">".$date." ".$item->organism->name." (id".$item->organism->id.") "."</OPTION>";
      }
    my $html;
#    $html .= qq{<FONT CLASS ="small">Organism count: }.scalar @opts.qq{</FONT>\n<BR>\n};
    unless (@opts) 
      {
	$html .=  qq{<input type = hidden name="org_id" id="org_id">};
	return $html;
      }

    $html .= qq{<SELECT id="recent_org_id" SIZE="5" MULTIPLE onChange="recent_dataset_chain()" >\n};
    $html .= join ("\n", @opts);
    $html .= "\n</SELECT>\n";
    $html =~ s/OPTION/OPTION SELECTED/;
    return $html;
  }


sub get_orgs
  {
    my %opts = @_;
    my $name = $opts{name};
    my $desc = $opts{desc};
    my @db;
    if ($name) 
      {
	@db = $coge->resultset("Organism")->search({name=>{like=>"%".$name."%"}});
      }
    elsif($desc)
      {
	@db = $coge->resultset("Organism")->search({description=>{like=>"%".$desc."%"}});
      }
    else
      {
	@db = $coge->resultset("Organism")->all;
      }
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
    my $dsid = $opts{dsid};
    if ($dsid)
      {
	my ($ds) = $coge->resultset('Dataset')->resolve($dsid);
	$dsname = $ds->name;
      }
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
	my %orgs;
	foreach my $item (sort {$b->version cmp $a->version || uc($a->name) cmp uc($b->name)} @ds)
	  {
	    next if $restricted_orgs->{$item->organism->name};
	    my $option = "<OPTION value=\"".$item->id."\">".$item->name."(v".$item->version.", id".$item->id.")</OPTION>";
	    if ($dsid && $dsid == $item->id)
	      {
		$option =~ s/(<OPTION)/$1 selected/;
	      }
	    push @opts, $option;
	    $orgs{$item->organism->id}=$item->organism;
	  }
	my @orgs = map	{"<OPTION value=\"".$_->id."\">".$_->name." (id".$_->id.")</OPTION>"} values %orgs;
	$html2 .= qq{<FONT CLASS ="small">Organism count: }.scalar @orgs.qq{</FONT>\n<BR>\n};
	$html2 .= qq{<SELECT id="org_id" SIZE="5" MULTIPLE onChange="dataset_chain()" >\n};
	$html2 .= join ("\n", @orgs);
	$html2 .= "\n</SELECT>\n";
	$html2 =~ s/OPTION/OPTION SELECTED/;
      }
    if (@opts) 
      {
	$html = qq{<FONT CLASS ="small">Dataset count: }.scalar (@opts).qq{</FONT>\n<BR>\n};
	$html .= qq{<SELECT id="ds_id" SIZE="5" MULTIPLE onChange="gen_data(['args__loading'],['ds_info']); dataset_info_chain()" >\n};
	$html .= join ("\n", @opts);
	$html .= "\n</SELECT>\n";
	$html =~ s/OPTION/OPTION SELECTED/ unless $dsid;
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
    $dataset = " <a href=\"".$ds->link."\">".$dataset."</a>" if $ds->link;
    my $source_name = $ds->datasource->name .": ". $ds->datasource->description;
    my $link = $ds->datasource->link;

    $link = "http://".$link if ($link && $link !~ /http/);
    $source_name = "<a href =\"".$link."\">".$source_name."</a>" if $ds->datasource->link;
    $html .= qq{<tr><td>Name: <td>$dataset}."\n";
    $html .= qq{<TR><TD>Data Source: <TD>$source_name (id}.$ds->datasource->id.qq{)}."\n";
    $html .= qq{<tr><td>Version: <td>}.$ds->version."\n";
    $html .= qq{<tr><td>Date deposited: <td>}.$ds->date."\n";
    my $org = $ds->organism->name;
    $org .= ": ".$ds->organism->description if $ds->organism->description;
    $org .="\n" ;
    $html .= qq{<tr><td>Organism: <td>$org\n};
    my %chr;
#    map{$chr{$_}++} ($ds->chromosomes, $DB->get_genomic_seq_obj->get_chromosome_for_dataset($ds));
    map{$chr{$_}++} ($ds->get_chromosomes);
    my $count = 100000;
    foreach my $item (sort keys %chr)
      {
	my ($num) = $item=~/(\d+)/;
	$num = $count unless $num;
	$chr{$item} = $num;
	$count++;
      }
    my @chr = sort {$chr{$a} <=> $chr{$b}}keys %chr;
    
    if (@chr)
      {
	my $size = scalar @chr;
	$size = 5 if $size > 5;
	$html .= qq{<tr><td>Chromosome};
	$html .= "s". " (".scalar @chr.")" if scalar @chr > 1;
	$html .= ":";
	$html .= "<td>";
    my $select;
	$select .= qq{<SELECT id="chr" size =$size onChange="gen_data(['args__searching for features. . .'],['chr_info']); gen_data(['args__waiting. . .'],['viewer']); gen_data(['args__'],['get_seq']); get_dataset_chr_info(['ds_id', 'chr'],['chr_info','viewer', 'get_seq'])" >\n};
	$select .= join ("\n", map {"<OPTION value=\"$_\">".$_."</OPTION>"} @chr)."\n";
	$select =~ s/OPTION/OPTION SELECTED/;
	$select .= "\n</SELECT>\n";
	$select .= qq{<tr><td>Total length of all chromosomes: <td><div id=ds_total_length class="link" onclick="\$('#ds_total_length').html('<span class=loading>loading. . .</span>');\$('#ds_total_length').removeClass('link'); get_total_length_for_ds(['args__dsid','ds_id'],['ds_total_length']);">Click for total length</div>} if scalar @chr>1;
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
    my $length = $ds->last_chromosome_position($chr);
    my $type = $ds->sequence_type;
    my $type_html = $type->name if $type;
    $type_html .= ": ".$type->description if $type && $type->description;
    $type_html = "Unknown" unless $type_html;
    my $gc = $length < 10000000? gen_gc_for_chromosome(dsid=>$ds->id, chr=>$chr): 0;
    $length = commify($length);
    $gc = $gc ? $gc : qq{<div id=chromosome_gc class="link" onclick="\$('#chromosome_gc').removeClass('link'); gen_gc_for_chromosome(['args__dsid','ds_id','args__chr','chr'],['chromosome_gc']);">Click for percent GC content</div>};
    $html .= qq{
<tr><td class = oblique>Specifics for chromosome $chr:
<tr><td>Nucleotides:<td>$length<td>$gc
<tr><td>Sequence Type:<td colspan=2>$type_html
<tr><td>Noncoding sequence:<td colspan=2><div id=noncoding_gc class="link" onclick = "gen_data(['args__loading'],['noncoding_gc']);\$('#noncoding_gc').removeClass('link');  gen_gc_for_noncoding(['args__dsid','ds_id','args__chr','chr'],['noncoding_gc']);">Click for percent GC content</div>
} if $length;
#    my $feat_string = get_feature_counts($dsd, $chr);
    my $feat_string = qq{
<div id=feature_count onclick="gen_data(['args__loading'],['feature_count_data']);\$('#feature_count').hide(0);get_feature_counts(['ds_id', 'chr'],['feature_count_data']);" >Click here for feature count information</div><div id=feature_count_data></div>
};
    $html .= "</table>$feat_string";

    my $viewer;
    if ($chr)
     {
	$viewer .= "<br><br>";
	$viewer .= "<font class=\"oblique\">Genome Viewer Launcher</font><br>";
	$viewer .= "<table>";
	$viewer .= "<tr><td nowrap class = \"ital\">Starting location: ";
	$viewer .= qq{<td><input type="text" size=10 value="10000" id="x">};
	my $zoom;
   	$zoom .= qq{<tr><td class = "ital">Zoom level:};
   	$zoom .= qq{<td><SELECT name = "z" id="z" size = 1>};
   	my @opts = map {"<OPTION value=\"$_\">".$_."</OPTION>"} (5..15);
   	$zoom .= join ("\n", @opts);
   	$zoom =~ s/OPTION value="4"/OPTION SELECTED value="4"/;
   	$zoom .= qq{</SELECT>};
	$viewer .= qq{<tr><td class = "ital">Zoom level:<td><input type = "text" size=10 value ="4" id = "z">};
	#    $viewer .= $zoom;
	$viewer .= "</table>";
	#$viewer .= qq{<input type="hidden" id="z" value="7">};
	$viewer .= qq{<input type="submit" value = "Launch Genome Viewer!" onClick="launch_viewer($dsd, '$chr')">};
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
	$seq_grab .= qq{<input type="submit" value = "Get Genomic Sequence!" onClick="launch_seqview($dsd, '$chr')">};
	$seq_grab .= qq{<div id="gseq"></div>};
      }
    return $html, $viewer, $seq_grab;
  }

sub get_feature_counts
  {
    my $dsd = shift;
    my $chr = shift;
    my $query = qq{
SELECT count(distinct(feature_id)), ft.name, ft.feature_type_id
  FROM feature
  JOIN feature_type ft using (feature_type_id)
  JOIN location l USING (feature_id)
 WHERE dataset_id = $dsd
   AND l.chromosome = "$chr"
 GROUP BY ft.name
};
    my $coge = CoGeX->dbconnect();
    my $dbh = DBI->connect($coge->db_connection_string,$coge->db_name,$coge->db_passwd);
    my $sth = $dbh->prepare($query);
    $sth->execute;
    my $feats = {};
    while (my $row = $sth->fetchrow_arrayref)
      {
	$feats->{$row->[1]} = {count=>$row->[0],
			       id=>$row->[2]};
      }
#    my $feats = $ds->get_feature_type_count(chr=>$chr);
    my $feat_string .= qq{<tr><td valign=top>Features:<td valign=top><table>};
    $feat_string .= join ("\n<tr valing=top>",map {
"<td valign=top><div id=$_  >".$_."</div>".
"<td valign=top>".$feats->{$_}{count}.
"<td><div id=".$_."_type class=\"link small\" 
  onclick=\" \$('#".$_."_type').removeClass('link'); 
  \$('#".$_."_type').removeClass('small'); 
  gen_data(['args__loading'],['".$_."_type']); 
  gen_gc_for_chromosome_and_type(['args__dsid','ds_id','args__chr','chr','args__type','args__$_'],['".$_."_type'])\">".'show %GC?</div>'.
"<td class='small link' onclick=\"window.open('FeatList.pl?dsid=$dsd&ftid=".$feats->{$_}{id}.";chr=$chr')\">Feature List?"
} sort {$a cmp $b} keys %$feats);
    $feat_string .= "</table>";
    if ($feats->{CDS})
      {
	$feat_string .= "<div class=link id=codon_usage onclick=\" \$('#codon_usage').removeClass('link'); gen_data(['args__loading'],['codon_usage']); get_codon_usage_for_chromosome(['args__dsid','ds_id','args__chr','chr'],['codon_usage'])\">"."Click for codon usage"."</div>";
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
    return "Total length: ".commify($total).", GC: ".sprintf("%.2f",100*$gc/($total))."%  AT: ".sprintf("%.2f",100*$at/($total))."%";
  }

sub gen_gc_for_chromosome
  {
    my %args = @_;
    my $dsid = $args{dsid};
    my $chr = $args{chr};
    return unless $dsid;
    my $ds = $coge->resultset('Dataset')->find($dsid);
    my ($gc, $at)=$ds->percent_gc(chr=>$chr);
    return "GC: ".(100*$gc)."%  AT: ".(100*($at))."%";
  }

sub gen_gc_for_noncoding
  {
    my %args = @_;
    my $dsid = $args{dsid};
    my $chr = $args{chr};
    return "no dsid specified" unless $dsid;
    my $ds = $coge->resultset('Dataset')->find($dsid);
    my $seq = $ds->get_genomic_sequence(chr=>$chr);
    foreach my $feat ($ds->features)
      {
	next unless $feat->type->name eq "CDS" ||
	  $feat->type->name =~ "RNA";
	foreach my $loc ($feat->locations)
	  {
	    substr($seq, $loc->start-1,($loc->stop-$loc->start+1)) = "N"x($loc->stop-$loc->start+1);
	  }
      }
    my $gc = $seq=~tr/GCgc/GCgc/;
    my $at = $seq=~tr/ATat/ATat/;
    my $total = $gc+$at;
    return "error" unless $total;
    return "Total length: ".commify($total).", GC: ".sprintf("%.2f",100*$gc/($total))."%  AT: ".sprintf("%.2f",100*$at/($total))."%";

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
    %codons = map {$_,$codons{$_}/$codon_total} keys %codons;
    %aa = map {$_,$aa{$_}/$aa_total} keys %aa;
    my $html = "Codon Usage: $code_type";
    $html .= CoGe::Accessory::genetic_code->html_code_table(data=>\%codons, code=>$code);

    $html .= "Predicted amino acid usage using $code_type";
    $html .= CoGe::Accessory::genetic_code->html_aa(data=>\%aa);
    return $html;
  }

sub get_total_length_for_ds
  {
    my %opts = @_;
    my $dsid = $opts{dsid};
    my $ds = $coge->resultset('Dataset')->find($dsid);
    my $length = 0;
    map {$length+=$ds->last_chromosome_position($_)} $ds->get_chromosomes();
    return commify($length);
  }

sub commify
    {
      my $text = reverse $_[0];
      $text =~ s/(\d\d\d)(?=\d)(?!\d*\.)/$1,/g;
      return scalar reverse $text;
    }

1;
