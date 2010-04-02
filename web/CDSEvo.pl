#! /usr/bin/perl -w

use strict;
use CGI;
use CGI::Ajax;
use CoGeX;
use CoGeX::Feature;
use Data::Dumper;
use CoGe::Accessory::LogUser;
use CoGe::Accessory::Web;
use HTML::Template;
use URI::Escape;
use Spreadsheet::WriteExcel;
use Benchmark;
use DBIxProfiler;
use CoGe::Accessory::genetic_code;
use Statistics::Basic::Mean;
use POSIX;


use vars qw($PAGE_NAME $TEMPDIR $USER $DATE $BASEFILE $coge $cogeweb $FORM);
$ENV{PATH}="/opt/apache/CoGe";
$DATE = sprintf( "%04d-%02d-%02d %02d:%02d:%02d",
		sub { ($_[5]+1900, $_[4]+1, $_[3]),$_[2],$_[1],$_[0] }->(localtime));
$PAGE_NAME = "FeatAAEvo.pl";
($USER) = CoGe::Accessory::LogUser->get_user();
$TEMPDIR = "/opt/apache/CoGe/tmp/";
$FORM = new CGI;
$coge = CoGeX->dbconnect();
#$coge->storage->debugobj(new DBIxProfiler());
#$coge->storage->debug(1);

my $pj = new CGI::Ajax(
		       gen_html=>\&gen_html,
		       get_orgs=>\&get_orgs,
		       go=>\&go,
		      );
$pj->js_encode_function('escape');
print $pj->build_html($FORM, \&gen_html);
#print $FORM->header,gen_html();


sub gen_html
{
    my $html;    
    my ($body) = gen_body();
    my $template = HTML::Template->new(filename=>'/opt/apache/CoGe/tmpl/generic_page.tmpl');
    $template->param(TITLE=>'Coding Sequence Evolution');
    $template->param(PAGE_TITLE=>'CDSEvo');
    $template->param(HELP=>'/wiki/index.php?title=CDSEvo');
    my $name = $USER->user_name;
    $name = $USER->first_name if $USER->first_name;
    $name .= " ".$USER->last_name if $USER->first_name && $USER->last_name;
    $template->param(USER=>$name);
    $template->param(LOGO_PNG=>"CDSEvo-logo.png");
    $template->param(LOGON=>1) unless $USER->user_name eq "public";
    $template->param(DATE=>$DATE);
    $template->param(BOX_NAME=>'');
    $template->param(BODY=>$body);
    $template->param(ADJUST_BOX=>0);
    $html .= $template->output;

    return $html;
 }
 
sub gen_body
  {
    my $template = HTML::Template->new(filename=>'/opt/apache/CoGe/tmpl/CDSEvo.tmpl');
    my $form = $FORM;
    $template->param(INITIALIZE=>1);
    my $html = $template->output;
    return $html;
  }
  
sub get_orgs
  {
    my %opts = @_;
    my $type = $opts{type};
    my $search = $opts{search};
    my $id_only = $opts{id_only};
    my @db;
    my $count;
    if ($type && $type eq "name")
      {
        @db = $coge->resultset("Organism")->search({name=>{like=>"%".$search."%"}});
        $count = scalar @db;
      }
    elsif($type && $type eq "desc")
      {
        @db = $coge->resultset("Organism")->search({description=>{like=>"%".$search."%"}});
        $count = scalar @db;
      }
    else
      {
        $count = $coge->resultset("Organism")->count();
#       @db = $coge->resultset("Organism")->all;
      }
    return map {$_->id} @db if $id_only;
    #my @db = $name ? $coge->resultset('Organism')->search({name=>{like=>"%".$name."%"}})
    #  : $coge->resultset('Organism')->all();
    ($USER) = CoGe::Accessory::LogUser->get_user();
    my @opts;
    foreach my $item (sort {uc($a->name) cmp uc($b->name)} @db)
      {
        next if $USER->user_name =~ /public/i && $item->restricted;
        push @opts, "<OPTION value=\"".$item->id."\" id=\"o".$item->id."\">".$item->name."</OPTION>";
      }
    my $html;
    $html .= qq{<FONT CLASS ="small" id="org_count">Organism count: }.$count.qq{</FONT>\n<BR>\n};
    if ($search && !@opts) 
      {
        $html .=  qq{<input type = hidden name="org_id" id="org_id"><br>};
        $html .= "No results";
        return $html;
      }
        unshift(@opts,"<OPTION value=\"all\" id=\"all\">All Listed Organisms</OPTION>");
    my $size = scalar @opts;
    $size = 8 if $size > 8;
    $html .= qq{<SELECT id="org_id" SIZE="$size" MULTIPLE >\n};
    $html .= join ("\n", @opts);
    $html .= "\n</SELECT>\n";
    $html =~ s/OPTION/OPTION SELECTED/;
    return $html;
  }

sub go{
  my %opts = @_;
  my $accn = $opts{accn};
  my $anno = $opts{anno};
  my $org_id = $opts{org_id};
  my $org_name = $opts{org_name};
  my $org_desc = $opts{org_desc};
  my ($data, $feats) = get_features(accn=>$accn, anno=>$anno, org_id=>$org_id, org_name=>$org_name, org_desc=>$org_desc);
  return $data unless ref ($data) =~ /hash/i;

  my $template = HTML::Template->new(filename=>'/opt/apache/CoGe/tmpl/CDSEvo.tmpl');
  $template->param(RESULTS=>1);
  my $aa_sort = CoGe::Accessory::genetic_code->sort_aa_by_gc();
  my $table_head = "<th>".join ("<th>", "GC% (org count)", map {$_."% (".$data->{$_}{bin_count}.")" } sort keys %$data);
  $template->param(GC_HEAD=>$table_head);
  my $max_aa = 0;
  my $min_aa = 100;
  foreach my $aa (keys %$aa_sort)
    {
      next if $aa eq "*";
      my (@tmp) = map {$data->{$_}{data}{$aa}} sort {$data->{$b}{data}{$aa} <=> $data->{$a}{data}{$aa}} keys %$data;
      $max_aa = $tmp[0] if $tmp[0] > $max_aa;
      $min_aa = $tmp[-1] if $tmp[-1] < $min_aa;
    }
  my @rows;
  foreach my $aa (sort {$aa_sort->{$b} <=> $aa_sort->{$a} || $a cmp $b}keys %$aa_sort)
    {
      my @row;
      next if $aa eq "*";
      push @row, "<td>$aa (GC:".sprintf("%.0f",100*$aa_sort->{$aa})."%)";
      foreach my $bin (sort keys %$data)
	{
	  my $aa_val = sprintf("%.2f",100*$data->{$bin}{data}{$aa});
	  my $rel_val = ($data->{$bin}{data}{$aa}-$min_aa) / ($max_aa-$min_aa);
	  my $color = get_color(val=>$rel_val);
	  push @row, "<td style=\"color: black; background-color: rgb(".join(",",@$color).")\">".$aa_val."%";
	}
      push @rows, {RESULTS_ROW=>join ("\n",@row)};
    }
  my $featlistlink = "/CoGe/FeatList.pl?fid=".join ("::",keys %$feats);
  $template->param(FEAT_COUNT=>scalar(keys %$feats));
  $template->param(INFO=>\@rows);
  $template->param(FEATLISTLINK=>$featlistlink);
  return $template->output;
}

sub get_features
  {
    my %opts = @_;
    my $accn = $opts{accn};
    my $anno = $opts{anno};
    my $org_id = $opts{org_id};
    my $org_name = $opts{org_name};
    my $org_desc = $opts{org_desc};
    my $fid = $opts{fid};
    $org_name = undef if $org_name =~ /search/i;
    $org_desc = undef if $org_desc =~ /search/i;
   my $weak_query = "Query needs to be better defined.";
    if (!$accn && !$anno && !$fid)
      {
        return $weak_query unless $org_id;
      }
    ($USER) = CoGe::Accessory::LogUser->get_user();

    my $search ={};
    $search->{feature_type_id}=3;
    $search->{"organism.restricted"}=0 if $USER->name =~ /public/i;
    $search->{'dataset_group.organism_id'}=$org_id if $org_id;
    unless ($org_id)
      {
	$search->{'organism.name'}={like=>"%".$org_name."%"} if $org_name;
	$search->{'organism.description'}={like=>"%".$org_desc."%"} if $org_desc;
      }
    my $join = {join=>[{feature=>{'dataset'=>{'dataset_connectors'=>{'dataset_group'=>'organism'}}}}]};
#    push @{$join->{join}}, 'annotations' if $anno;#=>['annotation',]};
#    push @{$join->{join}}, 'feature_names' if $accn;#=>['annotation',]};


    my %feats;
    if($accn)
      {
        map {$feats{$_->feature->id}= $_->feature} $coge->resultset('FeatureName')->search($search,$join)->search_literal('MATCH(me.name) AGAINST (?)',$accn);
      }
    if ($anno)
      {
	map {$feats{$_->feature->id}= $_->feature} $coge->resultset('Annotation')->search($search,$join)->search_literal('MATCH(annotation) AGAINST (?)',$anno);
      }
    unless ($accn || $anno)
      { #org only
	my @org_ids;
	$org_id = "all" unless $org_id;
	if ($org_id eq "all")
	  {
	    my ($otype, $search) = ("name", $org_name) if $org_name && $org_name ne "Search";
	    ($otype, $search) = ("desc", $org_desc) if $org_desc && $org_desc ne "Search";
	    @org_ids = get_orgs(id_only=>1, type=>$otype, search=>$search);
	  }
	else
	  {
	    push @org_ids, $org_id;
	  }
	foreach my $oid (@org_ids)
	  {
	    my $org = $coge->resultset('Organism')->find($oid);
	    next unless $org;
	    foreach my $dsg ($org->dataset_groups)
	      {
		foreach my $ds ($dsg->datasets)
		  {
		    map{$feats{$_->id}=$_} $ds->features({feature_type_id=>3});
		  }
	      }
	  }
      }

    my %data;
    my $count = 0;
    my @aa;

    my %codes; #store genetic codes for reuse;
    foreach my $feat (values %feats)
      {
	my ($code) = $codes{$feat->dataset->id} ? $codes{$feat->dataset->id} : $feat->genetic_code;
	$codes{$feat->dataset->id} = $code;
	my ($gc, $at) = $feat->gc_content();
	next unless $gc && $at;
	$gc *= 100;
	my $div = floor($gc/5);
	my $mod = $gc%5;
#	print join ("\t", $gc, $div, $mod),"\n";
	$div ++ if $mod >2;
	my $aa_usage = $feat->aa_frequency(counts=>1, code=>$code);
	@aa = keys %$aa_usage unless @aa;
	foreach my $aa (keys %$aa_usage)
	  {
	    push @{$data{5*$div}{$aa}}, $aa_usage->{$aa};
	  }
	$count++;
#	push @feats, {f=>$feat, gc=>$gc, at=>$at};
#	last if $count > 50;
      }
    my %return_data;
    foreach my $bin (sort keys %data)
      {
	$return_data{$bin}{bin_count}=scalar @{$data{$bin}{"A"}};
	foreach my $aa (@aa)
	  {
	    my $ave = sprintf("%.4f",Statistics::Basic::Mean->new($data{$bin}{$aa})->query);	    
	    $return_data{$bin}{data}{$aa}=$ave;
	  }
	$return_data{$bin}{data}{"*"} = 0 unless $return_data{$bin}{data}{"*"};
      }
    foreach my $bin (keys %return_data)
      {
	my $hash = $return_data{$bin}{data};
	my $total = 0;
	map {$total+=$hash->{$_}} keys %$hash;
	unless ($total)
	  {
	    delete $return_data{$bin};
	    next;
	  }
	map {$hash->{$_}=$hash->{$_}/$total} keys %$hash;
	map {$total+=$hash->{$_}} keys %$hash;
      }
  
    return \%return_data, \%feats;
  }

sub commify
    {
      my $text = reverse $_[0];
      $text =~ s/(\d\d\d)(?=\d)(?!\d*\.)/$1,/g;
      return scalar reverse $text;
    }

sub get_color
  {
    my %opts = @_;
    my $val = $opts{val};
    return [0,0,0] unless defined $val;
    return [125,125,125] if $val < 0;
    my @colors = (
                  [255,0,0], #red
                  [255,126,0], #orange
                  [255,255,0], #yellow
                  [0,255,0], # green
                  [0,255,255], # cyan
                  [0,0,255], # blue
#                 [255,0,255], #magenta
                  [126,0,126], #purple
                 );
    @colors = reverse @colors;
    my ($index1, $index2) = ((floor((scalar(@colors)-1)*$val)), ceil((scalar(@colors)-1)*$val));

    my $color=[];
    my $step = 1/(scalar (@colors)-1);
    my $scale = $index1*$step;
    my $dist = ($val-$scale)/($step);
    for (my $i=0; $i<=2; $i++)
      {
        my $diff = ($colors[$index1][$i]-$colors[$index2][$i])*$dist;
        push @$color, sprintf("%.0f", $colors[$index1][$i]-$diff);
      }
    return $color;
  }
