#! /usr/bin/perl -w
use strict;
use CGI;
use CGI::Carp 'fatalsToBrowser';
use CGI::Ajax;
use CoGe::Accessory::LogUser;
use CoGe::Accessory::Web;
use CoGe::Accessory::genetic_code;
use CoGe::Accessory::Restricted_orgs;
use HTML::Template;
use Data::Dumper;
use CoGe::Graphics::Chromosome;
use CoGe::Graphics::Feature;
use CoGe::Graphics::Feature::Gene;
use CoGe::Graphics::Feature::NucTide;
use CoGeX;
use POSIX;
use DBIxProfiler;

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

$coge = CoGeX->dbconnect();
#$coge->storage->debugobj(new DBIxProfiler());
#$coge->storage->debug(1);
my $pj = new CGI::Ajax(
		       db_lookup=>\&db_lookup,
		       source_search=>\&get_data_source_info_for_accn,
		       get_types=>\&get_types,
		       cogesearch=>\&cogesearch,
		       get_anno=>\&get_anno,
		       show_express=>\&show_express,
		       gen_data=>\&gen_data,
		       get_feature_types=>\&get_feature_types,
		       codon_table=>\&codon_table,
		       protein_table=>\&protein_table,
		       gc_content=>\&gc_content,
		       update_featlist=>\&update_featlist,
		       parse_for_FeatList=>\&parse_for_FeatList,
		       get_orgs=>\&get_orgs,
		       codon_aa_alignment=>\&codon_aa_alignment,
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
    my $blank = qq{<input type="hidden" id="Type_name">-------};
    my %seen;
    my @opts = sort map {"<OPTION>$_</OPTION>"} grep {! $seen{$_}++} map {$_->type->name} $coge->resultset('Feature')->search({
															       'feature_names.name'=>$accn,
															       dataset_id=>$dsid
															       },{join=>'feature_names'});

    $html .= "<font class=small>Type count: ".scalar @opts."</font>\n<BR>\n";
    $html .= qq{<SELECT id="Type_name" SIZE="10" MULTIPLE onChange="get_anno(['args__accn','accn_select','args__type','Type_name', 'args__dsid','dsid'],[show_anno])" >\n};
    $html .= join ("\n", @opts);
    $html .= "\n</SELECT>\n";
    $html =~ s/OPTION/OPTION SELECTED/;
    return $blank unless $html =~ /OPTION/;
    return ($html, 1);
  }

sub update_featlist
  {
    my $accn = shift;
    return unless $accn;
    my $type = shift;
    my $featid = shift;
    $accn .= " ($type)";
    return $accn,$featid;
  }
  
sub parse_for_FeatList
  {
	my $featlist = shift;
	my $url = "/CoGe/FeatList.pl?";
	foreach my $featid (split /,/,$featlist)
	{
		next unless $featid;
		$url .="fid=$featid&";
	}
	$url =~ s/&$//;
	return $url;
  }

sub cogesearch
  {
    my %opts = @_;
    my $accn = $opts{accn};
    $accn =~ s/^\s+//;
    $accn =~ s/\s+$//;
    my $anno = $opts{anno};
    $anno =~ s/^\s+//;
    $anno =~ s/\s+$//;
    my $type = $opts{type};
    my $org_id = $opts{org_id};
    my $org_name = $opts{org_name};
    my $org_desc = $opts{org_desc};
    my @org_ids;
    $org_id = "all" unless $org_id;
    if ($org_id eq "all")
      {
	my ($type, $search) = ("name", $org_name) if $org_name && $org_name ne "Search";
	($type, $search) = ("desc", $org_desc) if $org_desc && $org_desc ne "Search";
	@org_ids = get_orgs(id_only=>1, type=>$type, search=>$search);
      }
    else
      {
	push @org_ids, $org_id;
      }
    my $feat_accn_wild = $opts{feat_name_wild};
    my $feat_anno_wild = $opts{feat_anno_wild};
    my $blank = qq{<input type="hidden" id="accn_select">};
    my $weak_query = "Query needs to be better defined.";
    if (!$accn && !$anno)
      {
	return $weak_query.$blank unless $org_id && $type;
      }
    ($USER) = CoGe::Accessory::LogUser->get_user();
    my $restricted_orgs = restricted_orgs(user=>$USER);
    my $html;
    my %seen;
    my @opts;
    $accn = "%".$accn if $accn && ($feat_accn_wild eq "both" || $feat_accn_wild eq "left");
    $accn = $accn."%" if $accn && ($feat_accn_wild eq "both" || $feat_accn_wild eq "right");
    $anno = "%".$anno if $anno && ($feat_anno_wild eq "both" || $feat_anno_wild eq "left");
    $anno = $anno."%" if $anno && ($feat_anno_wild eq "both" || $feat_anno_wild eq "right");
    my $search = {'me.name'=>{like=>$accn}} if $accn;
    $search->{annotation}={like=>$anno} if $anno;
    $search->{feature_type_id}=$type if $type;
    $search->{organism_id}{ -not_in}=[values %$restricted_orgs] if values %$restricted_orgs;
    $search->{organism_id}{ -in}=[@org_ids] if @org_ids;
    my $join = {'feature'=>'dataset'};
    $join->{'feature'} = ['dataset','annotations'] if $anno;
    foreach my $name ($coge->resultset('FeatureName')->search(
							      $search,
							      {join=>$join,
							       order_by=>'name ASC',
							      },
							     ))
      {
	my $item = $name->name;
	next if $seen{uc($item)};
#	if (%{$restricted_orgs})
#	  {
#	    next if $restricted_orgs->{$name->feature->dataset->organism->name};
#	  }
	$seen{uc($item)}++;
	push @opts, "<OPTION>$item</OPTION>"
      }
    if (@opts > 5000)
      {
	return $blank."Search results over 5000, please refine your search.\n";
      }
    $html .= "<font class=small>Name count: ".scalar @opts."</font>\n<BR>\n";
    $html .= qq{<SELECT id="accn_select" SIZE="10" MULTIPLE onChange="source_search_chain(); " >\n};
    $html .= join ("\n", @opts);
    $html .= "\n</SELECT>\n";
    $html =~ s/<OPTION/<OPTION SELECTED/;
    return $blank."No results found\n" unless $html =~ /OPTION/;
    return $html;
  }

sub get_anno
  {
    my %opts = @_;
    my $accn = $opts{accn};
    my $fid = $opts{fid};
    return unless $accn || $fid;
    my $type = $opts{type};
    my $dataset_id = $opts{dsid};
    my @feats;
    if ($accn)
      {
	foreach my $feat ($coge->resultset('Feature')->search({
							       'feature_names.name'=>$accn,
							       dataset_id=>$dataset_id
							      },{join=>'feature_names'}))
	  {
	    push @feats, $feat if ($feat->type->name eq $type);
	  }
      }
    else
      {
	push @feats,$coge->resultset('Feature')->find($fid);
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
	my $z = 4;
	$anno .= join "\n<BR><HR><BR>\n", $feat->annotation_pretty_print_html(loc_link=>0);
	$anno .= "<table>";
	$anno .= "<tr>";
	$anno .= qq{<td><DIV id="dnaseq$i"><input type="button" value = "SeqView" onClick="window.open('SeqView.pl?featid=$featid&dsid=$ds&chr=$chr&featname=$accn');"></DIV>};
	$anno .= qq{<td><DIV id="dnaseq$i"><input type="button" value = "Blast" onClick="window.open('CoGeBlast.pl?featid=$featid');"></DIV>};
	$anno .= qq{<td><DIV id="loc$i"><input type="button" value = "GenomeView" onClick="window.open('GeLo.pl?chr=$chr&ds=$ds&x=$x&z=$z');"></DIV>};
#	$anno .= qq{<DIV id="exp$i"><input type="button" value = "Click for expression tree" onClick="gen_data(['args__Generating expression view image'],['exp$i']);show_express(['args__}.$accn.qq{','args__}.'1'.qq{','args__}.$i.qq{'],['exp$i']);"></DIV>};
	$anno .= qq{<td><DIV id="addfeat$featid"><input type="button" value = "Add to List" onClick="\$('#addfeat$featid').html('<i>$accn ($type) has been added to Feature List</i><br><br>');update_featlist(['args__$accn','args__$type','args__$featid'],[add_to_featlist]);"></DIV>};
	$anno .= "</table>";
	$anno .= qq{<DIV id="gc_info$i"><input type="button" value = "GC content" onClick="gc_content(['args__featid','args__$featid'],['gc_info$i'])"></DIV>};
	$anno .= qq{<DIV id="codon_info$i"><input type="button" value = "Codon usage" onClick="codon_table(['args__featid','args__$featid'],['codon_info$i'])"></DIV>} if $feat->type->name eq "CDS";
	$anno .= qq{<DIV id="protein_info$i"><input type="button" value = "Amino acid usage" onClick="protein_table(['args__featid','args__$featid'],['protein_info$i'])"></DIV>} if $feat->protein_sequence;
	$anno .= qq{<DIV id="codon_aa_align$i"><input type="button" value = "Codon/AA alignment" onClick="codon_aa_alignment(['args__featid','args__$featid'],['codon_aa_align$i'])"></DIV>} if $feat->type->name eq "CDS";


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
	my $name = $USER->user_name;
        $name = $USER->first_name if $USER->first_name;
        $name .= " ".$USER->last_name if $USER->first_name && $USER->last_name;
        $template->param(USER=>$name);

	$template->param(LOGON=>1) unless $USER->user_name eq "public";
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
    $template->param(FEAT_TYPE=> get_feature_types());
    $template->param(ORG_LIST=>get_orgs(type=>"none"));
    my $html = $template->output;
    return $html;
  }


sub get_feature_types
  {
    my %args = @_;
    my $orgid = $args{orgid};
    my $html=qq{<select  name="type" id="type" />
<OPTION VALUE=0>All</OPTION>
};
    if ($orgid)
      {
	my @tmp = $coge->resultset('FeatureType')->search (
																				   {
																				    organism_id=>$orgid,
																				   },
																				   {
																				    distinct=>['me.name'],
																				    join =>{'features'=>'dataset'}
																				   }
																				  );
	$html.= join("\n",map {"<OPTION value=\"".$_->id."\">".$_->name."</OPTION>"} sort {uc($a->name) cmp uc($b->name)} @tmp);
      }
    else
      {
	$html.= join ("\n",map {"<OPTION value=\"".$_->id."\">".$_->name."</OPTION>"} sort {uc($a->name) cmp uc($b->name)} $coge->resultset('FeatureType')->all);
      }
    return $html;
  }

sub get_data_source_info_for_accn
  {
    my %opts = @_;
    my $accn = $opts{accn};
    my $org_id = $opts{org_id};
    my $org_name = $opts{org_name};
    my $org_desc = $opts{org_desc};
    $org_id = "all" unless $org_id;
    my %org_ids;
    if ($org_id eq "all")
      {
	my ($type, $search);
	($type, $search) = ("name", $org_name) if $org_name && $org_name ne "Search";
	($type, $search) = ("desc", $org_desc) if $org_desc && $org_desc ne "Search";
	%org_ids = map {$_=>1} get_orgs(id_only=>1, type=>$type, search=>$search) if $type && $search;
      }
    else
      {
	$org_ids{$org_id}=1;
      }
    my $blank = qq{<input type="hidden" id="dsid">------------};
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
	my $org = $val->organism->name if $val->organism;
	next if $restricted_orgs->{$org};
	if (keys %org_ids) {next unless $org_ids{$val->organism->id};}
	my $name = $val->name;
	my $ver = $val->version;
	my $desc = $val->description;
	my $sname = $val->datasource->name if $val->datasource;
	my $ds_name = $val->name;
	my $title = "$org: $ds_name ($sname, v$ver)";
	$sources{$title}{id} = $val->id;
	$sources{$title}{v} = $ver;
      }
    my $html;
    $html .= qq{
<SELECT name = "dsid" id="dsid" MULTIPLE SIZE="10" onChange="get_types_chain();">
};
    my $count = 0;
    foreach my $title (sort { $a cmp $b || $sources{$b}{v} <=> $sources{$a}{v}} keys %sources)
      {
	my $id = $sources{$title}{id};
	$html .= qq{  <option value="$id" >$title\n};
	$html =~ s/option/option selected/ unless $count;
	$count++;
      }
    $html .= qq{</SELECT>\n};
    return ("<font class=small>Dataset count: ".$count ."</font>\n<BR>\n".$html, 1);
  }

sub get_orgs
  {
    my %opts = @_;
    my $type = $opts{type};
    my $search = $opts{search};
    my $id_only = $opts{id_only};
    my @db;
    if ($type && $type eq "name")
      {
	@db = $coge->resultset("Organism")->search({name=>{like=>"%".$search."%"}});
      }
    elsif($type && $type eq "desc")
      {
	@db = $coge->resultset("Organism")->search({description=>{like=>"%".$search."%"}});
      }
    else
      {
	@db = $coge->resultset("Organism")->all;
      }
    return map {$_->id} @db if $id_only;
    #my @db = $name ? $coge->resultset('Organism')->search({name=>{like=>"%".$name."%"}})
    #  : $coge->resultset('Organism')->all();
    ($USER) = CoGe::Accessory::LogUser->get_user();
    my $restricted_orgs = restricted_orgs(user=>$USER);
    my @opts;
    foreach my $item (sort {uc($a->name) cmp uc($b->name)} @db)
      {
	next if $restricted_orgs->{$item->name};
	push @opts, "<OPTION value=\"".$item->id."\" id=\"o".$item->id."\">".$item->name."</OPTION>";
      }
    my $html;
    $html .= qq{<FONT CLASS ="small" id="org_count">Organism count: }.scalar @opts.qq{</FONT>\n<BR>\n};
    unless (@opts) 
      {
	$html .=  qq{<input type = hidden name="org_id" id="org_id"><br>};
	$html .= "No results";
	return $html;
      }
	unshift(@opts,"<OPTION value=\"all\" id=\"all\">All Listed Organisms</OPTION>");
    $html .= qq{<SELECT id="org_id" SIZE="8" MULTIPLE >\n};
    $html .= join ("\n", @opts);
    $html .= "\n</SELECT>\n";
    $html =~ s/OPTION/OPTION SELECTED/;
    return $html;
  }


sub gc_content
  {
    my %opts = @_;
    my $featid = $opts{featid};
    return unless $featid;
    my ($feat) = $coge->resultset('Feature')->find($featid);
    my ($gc, $at) = $feat->gc_content;
    my $html = "GC:".(100*$gc)."%".", AT:".(100*$at)."%" ;
    return $html;
  }

sub codon_table
  {
    my %opts = @_;
    my $featid = $opts{featid};

    return unless $featid;
    my ($feat) = $coge->resultset('Feature')->find($featid);
    my ($codon, $code_type) = $feat->codon_frequency(counts=>1);
    my %aa;
    my ($code) = $feat->genetic_code;
    foreach my $tri (keys %$code)
      {
	$aa{$code->{$tri}}+=$codon->{$tri};
      }
    my $html = "Codon Usage: $code_type";
    $html .= CoGe::Accessory::genetic_code->html_code_table(data=>$codon, code=>$code, counts=>1);
    $html .= "Predicted amino acid usage for $code_type genetic code:";
    $html .= CoGe::Accessory::genetic_code->html_aa(data=>\%aa, counts=>1);
    return $html;
  }

sub protein_table
  {
    my %opts = @_;
    my $featid = $opts{featid};
    my ($feat) = $coge->resultset('Feature')->find($featid);
    my $aa = $feat->aa_frequency(counts=>1);
    my $html = "Amino Acid Usage";
    $html .= CoGe::Accessory::genetic_code->html_aa(data=>$aa, counts=>1);
    return $html;
  }

sub codon_aa_alignment
  {
    my %opts = @_;
    my $featid = $opts{featid};
    my ($feat) = $coge->resultset('Feature')->find($featid);
    my $seq = join (" ", $feat->genomic_sequence()=~ /(...)/g);
    my $aa = join ("   ",split //,$feat->protein_sequence());
    my @seq = $seq =~ /(.{1,80})/g;
    my @aa = $aa =~ /(.{1,80})/g;

    my $aln = "<pre>";
    for (my $i = 0; $i < @seq; $i++)
      {
	$aln .= $aa[$i]."\n".$seq[$i]."\n\n"
      }
    $aln .="</pre>";
    return $aln;
  }
