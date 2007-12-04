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
my $connstr = 'dbi:mysql:dbname=genomes;host=biocon;port=3306';
$coge = CoGeX->connect($connstr, 'cnssys', 'CnS' );
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
    $html .= qq{<SELECT id="Type_name" SIZE="10" MULTIPLE onChange="get_anno(['accn_select','Type_name', 'dsid'],[show_anno])" >\n};
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
    my $org_id = $opts{org_id} if $opts{org_id} && $opts{org_id} ne "all";
    my $feat_accn_wild = $opts{feat_name_wild};
    my $feat_anno_wild = $opts{feat_anno_wild};
#    print STDERR Dumper \%opts;
    my $blank = qq{<input type="hidden" id="accn_select">};
#    print STDERR "cogesearch: $accn\n";
#    print STDERR Dumper @_;
    my $weak_query = "Query needs to be better defined.";
#    return $weak_query.$blank unless length($accn) > 2 || $type || $org || length($anno) > 5;
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
    $search->{organism_id}={ -not_in=>[values %$restricted_orgs]} if values %$restricted_orgs;
    $search->{organism_id}=$org_id if $org_id;
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
	my $z = 7;
	$anno .= join "\n<BR><HR><BR>\n", $feat->annotation_pretty_print_html(loc_link=>0);
	$anno .= qq{<DIV id="loc$i"><input type="button" value = "Click for Genome view" onClick="window.open('GeLo.pl?chr=$chr&ds=$ds&x=$x&z=$z');"></DIV>};
#	$anno .= qq{<DIV id="exp$i"><input type="button" value = "Click for expression tree" onClick="gen_data(['args__Generating expression view image'],['exp$i']);show_express(['args__}.$accn.qq{','args__}.'1'.qq{','args__}.$i.qq{'],['exp$i']);"></DIV>};
	$anno .= qq{<DIV id="dnaseq$i"><input type="button" value = "Click for Sequence" onClick="window.open('SeqView.pl?featid=$featid&dsid=$ds&chr=$chr&featname=$accn');"></DIV>};
	$anno .= qq{<DIV id="gc_info$i"><input type="button" value = "Click for GC content" onClick="gc_content(['args__featid','args__$featid'],['gc_info$i'])"></DIV>};
	$anno .= qq{<DIV id="codon_info$i"><input type="button" value = "Click for codon usage" onClick="codon_table(['args__featid','args__$featid'],['codon_info$i'])"></DIV>} if $feat->type->name eq "CDS";
	$anno .= qq{<DIV id="protein_info$i"><input type="button" value = "Click for amino acid usage" onClick="protein_table(['args__featid','args__$featid'],['protein_info$i'])"></DIV>} if $feat->protein_sequence;
	$anno .= qq{<DIV id="addfeat$featid"><input type="button" value = "Add Feature to List" onClick="\$('#addfeat$featid').html('<i>$accn ($type) has been added to Feature List</i><br><br>');update_featlist(['args__$accn','args__$type','args__$featid'],[add_to_featlist]);"></DIV>};

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
    $template->param(FEAT_TYPE=> get_feature_types());
#    my @orgs;
#    ($USER) = CoGe::Accessory::LogUser->get_user();
#    my $restricted_orgs = restricted_orgs(user=>$USER);
#    foreach my $org ($coge->resultset('Organism')->all)
#      {
#	push @orgs, $org unless $restricted_orgs->{$org->name};
#      }
#    $template->param(ORG_LOOP=> [{ORG=>"<OPTION VALUE=0>All</OPTION>"},map {{ORG=>"<OPTION value=\"".$_->id."\">".$_->name."</OPTION>"}} sort {uc($a->name) cmp uc($b->name)} @orgs]);
    $template->param(ORG_LIST=>get_orgs());
    my $html = $template->output;
#    $html =~ s/(>gene<\/OPTION)/ SELECTED$1/; 
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
    my $accn = shift;
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
<SELECT name = "dsid" id="dsid" MULTIPLE SIZE="10" onChange="get_types_chain();">
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

sub get_orgs
  {
    my $name = shift;
    my @db = $name ? $coge->resultset('Organism')->search({name=>{like=>"%".$name."%"}})
      : $coge->resultset('Organism')->all();
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
	unshift(@opts,"<OPTION value=\"all\" id=\"all\">All Organisms</OPTION>") if !$name and @opts;
    $html .= qq{<SELECT id="org_id" SIZE="8" MULTIPLE >\n};
    $html .= join ("\n", @opts);
    $html .= "\n</SELECT>\n";
    $html =~ s/OPTION/OPTION SELECTED/;
    return $html;
  }


sub gc_content
  {
    my %args = @_;
    my $featid = $args{featid};
    return unless $featid;
    my ($feat) = $coge->resultset('Feature')->find($featid);
    my ($gc, $at) = $feat->gc_content;
    my $html = "GC:".(100*$gc)."%"."AT:".(100*$at)."%" ;
    return $html;
  }

sub codon_table
  {
    my %args = @_;
    my $featid = $args{featid};

    return unless $featid;
    my ($feat) = $coge->resultset('Feature')->find($featid);
    my ($codon, $code_type) = $feat->codon_frequency(counts=>1);
    my $codon_total = 0;
    grep {$codon_total+=$_} values %$codon;
    my $count = 0;
    my %aa;
    my $aa_total=0;
    my ($code) = $feat->genetic_code;
    foreach my $tri (keys %$code)
      {
	$aa{$code->{$tri}}+=$codon->{$tri};
	$aa_total+=$codon->{$tri};
      }
    my $html = "Codon Usage: $code_type<table>";
    foreach (sort { sort_nt2(substr($a, 0, 1)) <=> sort_nt2(substr($b,
    0, 1)) || sort_nt2(substr($a,1,1)) <=> sort_nt2(substr($b,1,1)) ||
    sort_nt(substr($a,2,1)) <=> sort_nt(substr($b,2,1)) } keys %$code)
      {
	my $str = "<tr><td>".$_."(".$code->{$_}.")<td>".$codon->{$_}."<td>(".sprintf("%.2f",100*$codon->{$_}/$codon_total)."%)";
	delete $codon->{$_};
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
     foreach (sort keys %$codon)
       {
 	my $str = "<tr><td>".$_."<td>".$codon->{$_}."<td>(".sprintf("%.2f",100*$codon->{$_}/$codon_total)."%)";
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
    $html .= "(Predicted amino acid usage for $code_type genetic code)<table>";
    $html .= join ("<tr>",map  {"<td>$_<td>".$aa{$_}."<td>(".sprintf("%.2f",100*$aa{$_}/$aa_total)."%)"} sort keys %aa);
    $html .= "</table>";
    return $html;
  }

sub protein_table
  {
    my %args = @_;
    my $featid = $args{featid};
    my ($feat) = $coge->resultset('Feature')->find($featid);
    my $aa = $feat->aa_frequency(counts=>1);
    my $total = 0;
    grep {$total+=$_} values %$aa;
    my $html = "Amino Acid Usage"."<table>";
    $html .= join ("<tr>",map  {"<td>$_<td>".$aa->{$_}."<td>(".sprintf("%.2f",100*$aa->{$_}/$total)."%)"} sort keys %$aa);
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
