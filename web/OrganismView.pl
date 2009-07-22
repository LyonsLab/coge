#! /usr/bin/perl -w
use strict;
use CGI;
use CGI::Carp 'fatalsToBrowser';
use CoGe::Accessory::LogUser;
use CoGe::Accessory::Web;
use HTML::Template;
use Data::Dumper;
use CGI::Ajax;
use CoGeX;
use Benchmark;
use File::Path;
use Benchmark qw(:all);

$ENV{PATH} = "/opt/apache2/CoGe/";

use vars qw( $DATE $DEBUG $TEMPDIR $TEMPURL $USER $FORM $NEATO $DOT $coge $HISTOGRAM);

# set this to 1 to print verbose messages to logs
$DEBUG = 0;
$TEMPDIR = "/opt/apache/CoGe/tmp/OrganismView";
mkpath ($TEMPDIR, 0,0777) unless -d $TEMPDIR;

$TEMPURL = "/CoGe/tmp/OrganismView";
$NEATO = "/usr/bin/neato";
$DOT = "/usr/bin/dot";
$HISTOGRAM = "/opt/apache/CoGe/bin/histogram.pl";
$| = 1; # turn off buffering
$DATE = sprintf( "%04d-%02d-%02d %02d:%02d:%02d",
		 sub { ($_[5]+1900, $_[4]+1, $_[3]),$_[2],$_[1],$_[0] }->(localtime));

$FORM = new CGI;
($USER) = CoGe::Accessory::LogUser->get_user();
$coge = CoGeX->dbconnect();
#$coge->storage->debugobj(new DBIxProfiler());
#$coge->storage->debug(1);



my $pj = new CGI::Ajax(
		       get_dataset_groups=>\&get_dataset_groups,
		       get_dataset_group_info=>\&get_dataset_group_info,
		       get_dataset => \&get_dataset,
		       get_dataset_info => \&get_dataset_info,
		       get_dataset_chr_info => \&get_dataset_chr_info,
		       gen_data => \&gen_data,
		       get_orgs => \&get_orgs,
		       get_org_info=>\&get_org_info,
		       get_recent_orgs=>\&get_recent_orgs,
		       get_start_stop=>\&get_start_stop,
		       get_feature_counts => \&get_feature_counts,
		       gen_gc_for_chromosome=> \&gen_gc_for_chromosome,
		       get_gc_for_noncoding=> \&get_gc_for_noncoding,
		       gen_gc_for_feature_type =>\&gen_gc_for_feature_type,
		       get_codon_usage=>\&get_codon_usage,
		       get_aa_usage=>\&get_aa_usage,
		       get_wobble_gc=>\&get_wobble_gc,
		       get_wobble_gc_diff=>\&get_wobble_gc_diff,
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
	$template->param(PAGE_TITLE=>'OrgView');
	$template->param(HEAD=>qq{});
	my $name = $USER->user_name;
        $name = $USER->first_name if $USER->first_name;
        $name .= " ".$USER->last_name if $USER->first_name && $USER->last_name;
        $template->param(USER=>$name);

	$template->param(LOGON=>1) unless $USER->user_name eq "public";
	$template->param(DATE=>$DATE);
	$template->param(LOGO_PNG=>"OrganismView-logo.png");
	$template->param(BODY=>$body);
#	$template->param(ADJUST_BOX=>1);
	$html .= $template->output;
      }
    return $html;
  }

sub gen_body
  {
    my $form = shift || $FORM;
    my $template = HTML::Template->new(filename=>'/opt/apache/CoGe/tmpl/OrganismView.tmpl');
    my $org_name = $form->param('org_name');
    my $desc = $form->param('org_desc');
    my $oid = $form->param('oid');
    my $org = $coge->resultset('Organism')->resolve($oid) if $oid;
    my $dsname = $form->param('dsname');
    my $dsid = $form->param('dsid');
    my $dsgid = $form->param('dsgid');
    my ($dsg) = $coge->resultset('DatasetGroup')->find($dsgid) if $dsgid;
    $org = $dsg->organism if $dsg;

    $org_name = $org->name if $org;
    $org_name = "Search" unless $org_name;
    $template->param(ORG_NAME=>$org_name) if $org_name;
    $desc = "Search" unless $desc;
    $template->param(ORG_DESC=>$desc) if $desc;
    $org_name = "" if $org_name =~ /Search/;
    my ($org_list, $org_count) = get_orgs(name=>$org_name, oid=>$oid, dsgid=>$dsgid);
    $template->param(ORG_LIST=>$org_list);
    $template->param(ORG_COUNT=>$org_count);
    #$template->param(RECENT=>get_recent_orgs());
    my ($ds) = $coge->resultset('Dataset')->resolve($dsid) if $dsid;
    $dsname = $ds->name if $ds;
    $dsname = "Search" unless $dsname;
    $template->param(DS_NAME=>$dsname);
    $dsname = "" if $dsname =~ /Search/;
    my ($dslist,$dscount) = get_dataset(dsname=>$dsname, dsid=>$dsid) if $dsname;
    $template->param(DS_LIST=>$dslist) if $dslist;
    $template->param(DS_COUNT=>$dscount) if $dscount;
    my $dsginfo = "<input type=hidden id=gstid>";
    $dsginfo .= $dsgid ? "<input type=hidden id=dsg_id value=$dsgid>" : "<input type=hidden id=dsg_id>";
    $template->param(DSG_INFO=>$dsginfo);
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
    my @opts;
    my %org_names;
    foreach my $item (@db)
      {
	my $date = $item->date;
	$date =~ s/\s.*//;
	next if $USER->user_name =~ /public/i && $item->organism->restricted;
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
    my $oid = $opts{oid};
    my $dsgid = $opts{dsgid};
    my $dsg = $coge->resultset('DatasetGroup')->find($dsgid) if $dsgid;
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
    my @opts;
    foreach my $item (sort {uc($a->name) cmp uc($b->name)} @db)
      {
	next if $USER->user_name =~ /public/i && $item->restricted;
	my $option = "<OPTION value=\"".$item->id."\"";
	$option .= " SELECTED" if $oid && $item->id == $oid;
	$option .= " SELECTED" if $dsg && $item->id == $dsg->organism->id;
	my $name = $item->name;
	if (length($name) > 50)
	  {
	    $name = substr($name, 0,50)."...";
	  }
	$option .= ">".$name." (id".$item->id.")</OPTION>";
	push @opts, $option;
      }
    my $html;
    unless (@opts) 
      {
	$html .=  qq{<input type = hidden name="org_id" id="org_id">No organisms found};
	return $html,0;
      }

    $html .= qq{<SELECT id="org_id" SIZE="5" MULTIPLE onChange="get_org_info_chain()" >\n};
    $html .= join ("\n", @opts);
    $html .= "\n</SELECT>\n";
    $html =~ s/OPTION/OPTION SELECTED/ unless $html =~ /SELECTED/;
    return $html, scalar @opts;
  }

sub get_org_info
  {
    my %opts = @_;
    my $oid = $opts{oid};
    return " " unless $oid;
    my $org = $coge->resultset("Organism")->find($oid);
    return "Unable to find an organism for id: $oid\n" unless $org;
    my $html;# = qq{<div class="backbox small">};
    $html .= $org->name."<br>";
    $html .= $org->description if $org->description;
    $html.= "<br><span class=alert>Restricted Organism!  Authorized Use Only!</span>" if $org->restricted;
#    $html .= "</div>";
    return $html;
  }

sub get_dataset_groups
    {
      my %opts = @_;
      my $oid = $opts{oid};
      my $dsgid = $opts{dsgid};
      my $org = $coge->resultset("Organism")->find($oid);
      my @opts;
      my %selected;
      $selected{$dsgid} = "SELECTED" if $dsgid;
      if ($org)
	{
	  my @dsg = $org->dataset_groups;
	  foreach my $dsg (@dsg)
	    {
	      $dsg->name($org->name) unless $dsg->name;
	      $selected{$dsg->id} = " " unless $selected{$dsg->id};
	    }
	  @opts = map {"<OPTION value=\"".$_->id."\" ".$selected{$_->id} .">".$_->name." (v".$_->version.", dsgid".$_->id. "): ". $_->genomic_sequence_type->name."</OPTION>"} sort {$b->version <=> $a->version || $a->name cmp $b->name} @dsg;
	}
      my $html;
      if (@opts) 
      {
#	$html = qq{<FONT CLASS ="small">Dataset group count: }.scalar (@opts).qq{</FONT>\n<BR>\n};
	$html .= qq{<SELECT id="dsg_id" SIZE="5" MULTIPLE onChange="get_dataset_group_info(['args__dsgid','dsg_id'],[dataset_chain]);" >\n};
	$html .= join ("\n", @opts);
	$html .= "\n</SELECT>\n";
	$html =~ s/OPTION/OPTION SELECTED/ unless $html =~ /SELECTED/i;
      }
    else
      {
	$html .=  qq{<input type = hidden name="dsg_id" id="dsg_id">};
      }
      return $html, scalar @opts;
    }

sub get_dataset_group_info
  {
    my %opts = @_;
    my $dsgid = $opts{dsgid};
    return " " unless $dsgid;
    my $dsg = $coge->resultset("DatasetGroup")->find($dsgid);
    return "Unable to create dataset_group object for id: $dsgid" unless $dsg;
    my $html;# = qq{<div style="overflow:auto; max-height:78px">};
    my $total_length;
    my @gs = sort {$a->chromosome cmp $b->chromosome} $dsg->genomic_sequences;
    my $chr_num = scalar @gs;
    
#    $html .=qq{<table class=small ><tr class=small valign=top>};
    my $count =0;
    foreach my $gs (@gs)
      {
	my $chr = $gs->chromosome;
	my $length = $gs->sequence_length;
	$total_length += $length;
	$length = commify($length);
	$html .= qq{$chr:  $length bp<br>};
	$count++;
      }
    $html = qq{<div style="overflow:auto; max-height:78px">};
    $html .= "<table class='small annotation_table'>";
    $html .= qq{<tr><td>Chromosome count: <td>$chr_num<br>}. qq{<div style="float:left;">};
    $html .= qq{<tr><td>Sequence Type: <td>}.$dsg->genomic_sequence_type->name.qq{<input type=hidden id=gstid value=}.$dsg->genomic_sequence_type->id.qq{>};
    $html .= qq{<tr><td>Total length: };
    $html .= qq{<td>}.commify($total_length)." bp";
    my $gc = $total_length < 10000000? gen_gc_for_chromosome(dsgid=>$dsgid): 0;
    $gc = $gc ? " -- ".$gc : qq{  </div><div style="float: left; text-indent: 1em;" id=datasetgroup_gc class="link" onclick="\$('#datasetgroup_gc').removeClass('link'); gen_gc_for_chromosome(['args__dsgid','dsg_id','args__gstid', 'gstid'],['datasetgroup_gc']);">  Click for percent GC content</div>};
    $html .= "<td>$gc";


    $html .= qq{
<tr><td>Noncoding sequence:<td colspan=2><div id=dsg_noncoding_gc class="link" onclick = "gen_data(['args__loading...'],['dsg_noncoding_gc']);\$('#dsg_noncoding_gc').removeClass('link');  get_gc_for_noncoding(['args__dsgid','dsg_id','args__gstid', 'gstid'],['dsg_noncoding_gc']);">Click for percent GC content</div>
} if $total_length;

    my $feat_string = qq{
<tr><td><div id=dsg_feature_count class="small link" onclick="get_feature_counts(['args__dsgid','dsg_id', 'args__gstid','gstid'],['feature_count_data']);" >Click for feature counts</div>};
    $html .= $feat_string;
    $html .= "</table>";
    return $html;
  }
  
sub get_dataset
  {
    my %opts = @_;
    my $dsgid = $opts{dsgid};
    my $dsname = $opts{dsname};
    my $dsid = $opts{dsid};
    return "<hidden id='ds_id'>",0 unless  $dsid || $dsname|| $dsgid;
    if ($dsid)
      {
	my ($ds) = $coge->resultset('Dataset')->resolve($dsid);
	$dsname = $ds->name;
      }
    my $html; 
    my @opts;
    if ($dsgid)
      {
	my $dsg = $coge->resultset("DatasetGroup")->find($dsgid);
	@opts = map {"<OPTION value=\"".$_->id."\">".$_->name. " (v".$_->version.", dsid".$_->id.")</OPTION>"} sort {$b->version <=> $a->version || $a->name cmp $b->name} $dsg->datasets if $dsg;

	}
    elsif ($dsname)
      {
	my @ds = $coge->resultset("Dataset")->search({name=>{like=>"%".$dsname."%"}});
	($USER) = CoGe::Accessory::LogUser->get_user();
	my %orgs;
	foreach my $item (sort {$b->version <=> $a->version || uc($a->name) cmp uc($b->name)} @ds)
	  {
	    next if $USER->user_name =~ /public/i && $item->organism->restricted;
	    my $option = "<OPTION value=\"".$item->id."\">".$item->name."(v".$item->version.", id".$item->id.")</OPTION>";
	    if ($dsid && $dsid == $item->id)
	      {
		$option =~ s/(<OPTION)/$1 selected/;
	      }
	    push @opts, $option;
	    $orgs{$item->organism->id}=$item->organism;
	  }
      }
    if (@opts) 
      {
	$html .= qq{<SELECT id="ds_id" SIZE="5" MULTIPLE onChange="dataset_info_chain()" >\n};
	$html .= join ("\n", @opts);
	$html .= "\n</SELECT>\n";
	$html =~ s/OPTION/OPTION SELECTED/ unless $dsid;
      }
    else
      {
	$html .=  qq{<input type = hidden name="ds_id" id="ds_id">};
      }
    return $html, scalar @opts;
  }

sub get_dataset_info
  {
    my $dsd = shift;
    return qq{<input type="hidden" id="chr" value="">}, " ",0 unless ($dsd); # error flag for empty dataset

    my $ds = $coge->resultset("Dataset")->find($dsd);
    my $html = "";
    return "unable to find dataset object for id: $dsd"  unless $ds;
    $html = "<table class=\"small annotation_table\">";
    my $dataset = $ds->name;
    $dataset .= ": ". $ds->description if $ds->description;
    $dataset = " <a href=\"".$ds->link."\" target=_new\>".$dataset."</a>" if $ds->link;
    my $source_name = $ds->data_source->name ;
    $source_name.=": ". $ds->data_source->description if $ds->data_source->description;
    my $link = $ds->data_source->link;

    $link = "http://".$link if ($link && $link !~ /http/);
    $source_name = "<a href =\"".$link."\" target=_new\>".$source_name."</a>" if $ds->data_source->link;
    $html .= qq{<tr><td>Name: <td>$dataset}."\n";
    $html .= qq{<TR><TD>Data Source: <TD>$source_name (id}.$ds->data_source->id.qq{)}."\n";
    $html .= qq{<tr><td>Version: <td>}.$ds->version."\n";
    $html .= qq{<tr><td>Organism:<td class="link"><a href="OrganismView.pl?oid=}.$ds->organism->id.qq{" target=_new>}.$ds->organism->name."</a>\n";
    $html .= qq{<tr><td>Date deposited: <td>}.$ds->date."\n";


    my $html2;
    my %chr;
    map{$chr{$_->chromosome}={length=>$_->stop}} ($ds->get_chromosomes(ftid=>301, length=>1)); #the chromosome feature type in coge is 301
    my $count = 100000;
    foreach my $item (sort keys %chr)
      {
	my ($num) = $item=~/(\d+)/;
	$num = $count unless $num;
	$chr{$item}{num} = $num;
	$count++;
      }
    my @chr = scalar keys %chr > 1000 ? sort {$chr{$b}{length} <=> $chr{$a}{length}} keys %chr
      : sort {$chr{$a}{num} <=> $chr{$b}{num} || $a cmp $b}keys %chr;
    my $length =0;
    if (@chr)
      {
	my $size = scalar @chr;
	$size = 5 if $size > 5;
	my $select;
	$select .= qq{<SELECT id="chr" size =$size onChange="dataset_chr_info_chain()" >\n};
	if (scalar @chr > 1000)
	  {
	    my @tmp = @chr[0..999];
	    	    $select .= join ("\n", map {"<OPTION value=\"$_\">".$_." (".commify($chr{$_}{length})." bp)</OPTION>"} @tmp)."\n";
	  }
	else
	  {
	    $select .= join ("\n", map {"<OPTION value=\"$_\">".$_." (".commify($chr{$_}{length})." bp)</OPTION>"} @chr)."\n";
	  }
	$select =~ s/OPTION/OPTION SELECTED/;
	$select .= "\n</SELECT>\n";

	$html2 .= $select;
	map{$length += $chr{$_}{length}} @chr;
      }
    else {
      $html2 .= qq{<input type="hidden" id="chr" value="">};
      $html2 .= "<tr><td>No chromosomes";
    }
    $html .= "<tr><td>Total length:<td><div style=\"float: left;\">".commify($length);
    my $gc = $length < 10000000? gen_gc_for_chromosome(dsid=>$ds->id): 0;
    $gc = $gc ? $gc : qq{  </div><div style="float: left; text-indent: 1em;" id=dataset_gc class="link" onclick="\$('#dataset_gc').removeClass('link'); gen_gc_for_chromosome(['args__dsid','ds_id','args__gstid', 'gstid'],['dataset_gc']);">  Click for percent GC content</div>} if $length;
    $html .= " -- ".$gc if $gc;
    $html .= qq{</table>};
    my $feat_string = qq{
<div id=ds_feature_count class="small link" onclick="get_feature_counts(['args__dsid','ds_id','args__gstid', 'gstid'],['feature_count_data']);" >Click for feature counts</div>};
    $html .= $feat_string;
    my $chr_count = scalar (@chr);
    $chr_count .= " <span class=small>Only 1000 largest shown</span>" if ($chr_count >1000); 
    return $html, $html2, $chr_count;
  }

sub get_dataset_chr_info
  {
    my $dsid = shift;
    my $chr = shift;
    my $dsgid = shift;
    $dsgid = 0 unless defined $dsgid;
    $dsid = 0 unless $dsid;
    unless ($dsid && $chr) # error flag for empty dataset
	{
		return "", "", "";
	}
    my $start = "'start'";
    my $stop = "'stop'";	
    my $html .= "<table class=\"small annotation_table\">";
    my $ds = $coge->resultset("Dataset")->find($dsid);
    return $html unless $ds;
    my $length = 0;
    $length = $ds->last_chromosome_position($chr) if $chr;
    my $gc = $length < 10000000? gen_gc_for_chromosome(dsid=>$ds->id, chr=>$chr): 0;
    $gc = $gc ? " -- ".$gc : qq{<div id=chromosome_gc class="link" onclick="\$('#chromosome_gc').removeClass('link'); gen_gc_for_chromosome(['args__dsid','ds_id','args__chr','chr','args__gstid', 'gstid'],['chromosome_gc']);">Click for percent GC content</div>};
    $length = commify($length)." bp";
    $html .= qq{
<tr><td class = oblique>Specifics for chromosome $chr:
<tr><td>Nucleotides:<td>$length<td>$gc
};

    $html .= qq{
<tr><td>Noncoding sequence:<td colspan=2><div id=noncoding_gc class="link" onclick = "gen_data(['args__loading..'],['noncoding_gc']);\$('#noncoding_gc').removeClass('link');  get_gc_for_noncoding(['args__dsid','ds_id','args__chr','chr','args__gstid', 'gstid'],['noncoding_gc']);">Click for percent GC content</div>
} if $length;

    $html .= "</table>";
    my $feat_string = qq{
<div class=small id=feature_count onclick="get_feature_counts(['args__dsid','ds_id','args__chr','chr','args__gstid', 'gstid'],['feature_count_data']);" >Click for feature counts</div>};
    $html .= "$feat_string";

    my $viewer;
    if ($chr)
     {
	$viewer .= "<font class=\"oblique\">Genome Viewer</font><br>";
	$viewer .= "<table class=\"small backbox\">";
	$viewer .= "<tr><td nowrap class = \"ital\">Starting location: ";
	$viewer .= qq{<td><input type="text" size=10 value="20000" id="x">};
	$viewer .= qq{<tr><td class = "ital">Zoom level:<td><input type = "text" size=10 value ="5" id = "z">};
	$viewer .= "</table>";
	$viewer .= qq{<span class='ui-button ui-button-icon-left ui-state-default ui-corner-all' onClick="launch_viewer('$dsgid', '$chr', '$dsid')"><span class="ui-icon ui-icon-newwin"></span>Launch Genome Viewer</span>};
      }
    my $seq_grab;
    if ($chr)
      {
	$seq_grab .= qq{<font class="oblique">Genomic Sequence Retrieval</font><br>};
	$seq_grab .= qq{<table class=\"small backbox\">};
	$seq_grab .= "<tr><td class = \"ital\">Start position: ";
	$seq_grab .= qq{<td><input type="text" size=10 value="1" id="start">};
	$seq_grab .= "<tr><td class = \"ital\">End position: ";
	$seq_grab .= qq{<td><input type="text" size=10 value="100000" id="stop">};
	$seq_grab .= qq{</table>};
	$seq_grab .= qq{<span class='ui-button ui-button-icon-left ui-state-default ui-corner-all' onClick="launch_seqview('$dsgid', '$chr','$dsid')"><span class="ui-icon ui-icon-newwin"></span>Get Sequence</span>};
	$seq_grab .= qq{<div id="gseq"></div>};
      }
    return $html, $viewer, $seq_grab;
  }

sub get_feature_counts
  {
    my %opts = @_;
    my $dsid = $opts{dsid};
    my $dsgid = $opts{dsgid};
    my $gstid=$opts{gstid};
    my $chr = $opts{chr};
    my $query;
    my $name;
    if ($dsid)
      {
	my $ds = $coge->resultset('Dataset')->find($dsid);
	$name = "dataset ". $ds->name;
	$query = qq{
SELECT count(distinct(feature_id)), ft.name, ft.feature_type_id
  FROM feature
  JOIN feature_type ft using (feature_type_id)
 WHERE dataset_id = $dsid
};
	$query .= qq{AND chromosome = '$chr'} if $chr;
	$query.= qq{
  GROUP BY ft.name
};
	$name .= " chromosome $chr" if $chr;
      }
    elsif ($dsgid)
      {
	my $dsg = $coge->resultset('DatasetGroup')->find($dsgid);
	$name = "dataset group ";
	$name .= $dsg->name ? $dsg->name : $dsg->organism->name;
	$query = qq{
SELECT count(distinct(feature_id)), ft.name, ft.feature_type_id
  FROM feature
  JOIN feature_type ft using (feature_type_id)
  JOIN dataset_connector dc using (dataset_id)
 WHERE dataset_group_id = $dsgid
  GROUP BY ft.name

};
      }

    my $coge = CoGeX->dbconnect();
    my $dbh = DBI->connect($coge->db_connection_string,$coge->db_name,$coge->db_passwd);
    my $sth = $dbh->prepare($query);
    $sth->execute;
    my $feats = {};
    while (my $row = $sth->fetchrow_arrayref)
      {
	my $name = $row->[1];
	$name =~ s/\s+/_/g;
	$feats->{$name} = {count=>$row->[0],
			   id=>$row->[2],
			   name=>$row->[1],
			  };
      }
    #    my $feats = $ds->get_feature_type_count(chr=>$chr);

    my $gc_args;
    $gc_args .= qq{'args__dsid','ds_id',} if $dsid;
    $gc_args .= qq{'args__dsgid','dsg_id',} if $dsgid;
    $gc_args .= qq{'args__chr','chr',} if $chr;
    my $feat_list_string = $dsid ? "dsid=$dsid" : "dsgid=$dsgid";
    $feat_list_string .= ";chr=$chr" if $chr;
    my $feat_string .= qq{<div class=oblique>Features for $name</div>};
    $feat_string .= qq{<table class = " backbox small">};
    $feat_string .= "<tr valign=top>". join ("\n<tr valign=top>",map {
      "<td valign=top><div id=$_  >".$feats->{$_}{name}."</div>".
	    "<td valign=top>".$feats->{$_}{count}.
	      "<td><div id=".$_."_type class=\"link small\" 
  onclick=\"
  \$('#gc_histogram').dialog('option','title', 'Histogram of GC content for ".$feats->{$_}{name}."s');
  \$('#gc_histogram').html('loading...');
  \$('#gc_histogram').dialog('open');
  gen_gc_for_feature_type([$gc_args 'args__gstid','gstid','args__typeid','args__".$feats->{$_}{id}."'],['gc_histogram'])\">".'show %GC?</div>'.
    "<td class='small link' onclick=\"window.open('FeatList.pl?$feat_list_string"."&ftid=".$feats->{$_}{id}.";gstid=$gstid')\">Feature List?"
  } sort {$a cmp $b} keys %$feats);
	$feat_string .= "</table>";
    
    if ($feats->{CDS})
      {
	my $args;
	$args .= "'args__dsid','ds_id'," if $dsid;
	$args .= "'args__dsgid','dsg_id'," if $dsgid;
	$args .= "'args__chr','chr'," if $chr;
	$feat_string .= "<div class=\"small link\" id=wobble_gc onclick=\"\$('#wobble_gc_histogram').html('loading...');\$('#wobble_gc_histogram').dialog('open');get_wobble_gc([$args],['wobble_gc_histogram'])\">"."Click for codon wobble GC content"."</div>";
	$feat_string .= "<div class=\"small link\" id=wobble_gc_diff onclick=\"\$('#wobble_gc_diff_histogram').html('loading...');\$('#wobble_gc_diff_histogram').dialog('open');get_wobble_gc_diff([$args],['wobble_gc_diff_histogram'])\">"."Click for diff(CDS GC vs. codon wobble GC) content"."</div>";
	$feat_string .= "<div class=\"small link\" id=codon_usage onclick=\"
        \$('#codon_usage_table').html('loading...');
        \$('#codon_usage_table').dialog('open');
        get_codon_usage([$args],['codon_usage_table']); 
        \">".
        "Click for codon usage table"."</div>";
	$feat_string .= "<div class=\"small link\" id=aa_usage onclick=\"
        \$('#aa_usage_table').html('loading...');
        \$('#aa_usage_table').dialog('open');
        get_aa_usage([$args],['aa_usage_table']); 
        \">".
        "Click for amino acid usage table"."</div>";

      }
    $feat_string .= "None" unless keys %$feats;
    return $feat_string;
  }

sub gen_data
  {
    my $message = shift;
    return qq{<font class="small alert">$message. . .</font>};
  }

sub gen_gc_for_feature_type
  {
    my %args = @_;
    my $dsid = $args{dsid};
    my $dsgid = $args{dsgid};
    my $chr = $args{chr};
    my $typeid = $args{typeid};
    my $gstid = $args{gstid}; #genomic sequence type id
    return unless $dsid || $dsgid;
    my $gc = 0;
    my $at = 0;
    my $n = 0;
    my $type = $coge->resultset('FeatureType')->find($typeid);
    my $search;
    $search = {"feature_type_id"=>$typeid};
    $search->{"me.chromosome"}=$chr if defined $chr;
    my @data;
    my @dsids;
    push @dsids, $dsid if $dsid;
    if ($dsgid)
      {
	my $dsg = $coge->resultset('DatasetGroup')->find($dsgid);
	unless ($dsg)
	  {
	    my $error =  "unable to create dsg object using id $dsgid\n";
	    return $error;
	  }
	foreach my $ds ($dsg->datasets())
	  {
	    push @dsids, $ds->id;
	  }
      }
    foreach my $dsidt (@dsids)
      {
	my $ds = $coge->resultset('Dataset')->find($dsidt);
	unless ($ds)
	  {
	    warn "no dataset object found for id $dsidt\n";
	    next;
	  }
	my $t1 = new Benchmark;
	my %seqs; #let's prefetch the sequences with one call to genomic_sequence (slow for many seqs)
	if ($chr)
	  {
	    $seqs{$chr} = $ds->genomic_sequence(chr=>$chr, seq_type=>$gstid);
	  }
	else
	  {
	    %seqs= map {$_, $ds->genomic_sequence(chr=>$_, seq_type=>$gstid)} $ds->chromosomes;
	  }
	my $t2 = new Benchmark;
	my @feats = $ds->features($search,{join=>['locations', {'dataset'=>{'dataset_connectors'=>'dataset_group'}}],
					   prefetch=>['locations',{'dataset'=>{'dataset_connectors'=>'dataset_group'}}],
					  });
#	my $unpack_str;# = join (" ", map{"x".($_->start-1)." A".($_->stop-$_->start+1)} @feats);
#	my $pos = 0;
#	print STDERR "Mark!\n";
#	foreach my $feat (@feats)
# 	  {
# 	    $pos = $feat->start-1-$pos;
# 	    $unpack_str .= " x$pos";
# 	    my $dist = $feat->stop-$feat->start+1;
# 	    $unpack_str .=" A$dist";
# 	    print STDERR $pos,"\t";
# 	    $pos = $feat->stop;
# 	    print STDERR join ("\t", $pos, $feat->start, $dist),"\n";
# 	  }
# 	$unpack_str =~ s/^\s+//;
# #	print STDERR "!".$unpack_str,"!\n";

# 	my $i=0;
# 	foreach my $seq (unpack($unpack_str, $seqs{$feats[0]->chromosome}))
# 	  {
# 	    $feats[$i]->genomic_sequence($seq);
# 	    $i++;
# 	  }
 	foreach my $feat (@feats)
 	  {

#	    my $str = $seqs{$feat->chromosome};
#	    my $offset = "x".($feat->start-1);
#	    my $length = "A".($feat->stop-$feat->start+1);
#	    my $seq = unpack("x".($feat->start-1)." A".($feat->stop-$feat->start+1), $seqs{$feat->chromosome});
	    #reason this is slow is due to substr handling a large sequence.  Unpack is the same, but perhaps using multiple retrieves in one go will help.
	    my $seq = substr($seqs{$feat->chromosome}, $feat->start-1, $feat->stop-$feat->start+1);
#	    next;
	    $feat->genomic_sequence(seq=>$seq);
	    my @gc = $feat->gc_content(counts=>1);
	    $gc+=$gc[0] if $gc[0] =~ /^\d+$/;
	    $at+=$gc[1] if $gc[1] =~ /^\d+$/;
	    $n+=$gc[2] if $gc[2] =~ /^\d+$/;
	    my $total = 0;
	    $total += $gc[0] if $gc[0];
	    $total += $gc[1] if $gc[1];
	    $total += $gc[2] if $gc[2];
	    push @data, sprintf("%.2f",100*$gc[0]/$total) if $total;
	  }
	my $t3 = new Benchmark;
	my $get_seq_time = timestr(timediff($t2,$t1));
	my $process_seq_time = timestr(timediff($t3,$t2));
# 	print STDERR qq{

# -----------------------------------------
# Organism View: sub gen_gc_for_feature_type

#  Time to get sequence:        $get_seq_time
#  Time to process sequence:    $process_seq_time
# -----------------------------------------

# };
       }
    my $total = $gc+$at+$n;
    return "error" unless $total;

    my $file = $TEMPDIR."/".join ("_",@dsids);
    $file .= "_".$chr."_" if defined $chr;
    $file .= "_".$type->name."_gc.txt";
    open(OUT, ">".$file);
    print OUT "#wobble gc for dataset ids: ".join (" ", @dsids),"\n";
    print OUT join ("\n", @data),"\n";
    close OUT;
    my $cmd = $HISTOGRAM;
    $cmd .= " -f $file";
    my $out = $file;
    $out =~ s/txt$/png/;
    $cmd .= " -o $out";
    $cmd .= " -t \"".$type->name." gc content\"";
    $cmd .= " -min 0";
    $cmd .= " -max 100";
    `$cmd`;
    

    my $info= "<div class = small>Total length: ".commify($total)." bp, GC: ".sprintf("%.2f",100*$gc/($total))."%  AT: ".sprintf("%.2f",100*$at/($total))."%  N: ".sprintf("%.2f",100*($n)/($total))."%</div>";
    $out =~ s/$TEMPDIR/$TEMPURL/;
    my $hist_img = "<img src=\"$out\">";
    return $info."<br>". $hist_img;
  }

sub gen_gc_for_chromosome
  {
    my %args = @_;
    my $dsid = $args{dsid};
    my $chr = $args{chr};
    my $gstid = $args{gstid};
    my $dsgid = $args{dsgid};
    my @ds;
    if ($dsid)
      {
	my $ds = $coge->resultset('Dataset')->find($dsid);
	push @ds, $ds if $ds;
      }
    if ($dsgid)
      {
	my $dsg = $coge->resultset('DatasetGroup')->find($dsgid);
	$gstid = $dsg->type->id;
	map {push @ds,$_} $dsg->datasets;
      }
    return unless @ds;
    my ($gc, $at, $n) = (0,0,0);
    my %chr;
    foreach my $ds (@ds)
      {
	if ($chr)
	  {
	    $chr{$chr}=1;
	  }
	else
	  {
	    map {$chr{$_}=1} $ds->chromosomes;
	  }
	foreach my $chr(keys %chr)
	  {
	    my @gc =$ds->percent_gc(chr=>$chr, seq_type=>$gstid, counts=>1);
	    $gc+= $gc[0] if $gc[0];
	    $at+= $gc[1] if $gc[1];
	    $n+= $gc[2] if $gc[2];
	  }
      }
    my $total = $gc+$at+$n;
    my $results = "GC: ".sprintf("%.2f",100*$gc/$total)."%  AT: ".sprintf("%.2f",100*$at/$total)."%  N: ".sprintf("%.2f",100*$n/$total)."%" if $total;
    return $results;
  }

sub get_gc_for_noncoding
  {
    my %args = @_;
    my $dsid = $args{dsid};
    my $dsgid = $args{dsgid};
    my $chr = $args{chr};
    my $gstid = $args{gstid}; #genomic sequence type id
    return "error" unless $dsid || $dsgid;
    my $gc = 0;
    my $at = 0;
    my $n = 0;
    my $search;
    $search = {"feature_type_id"=>3};
    $search->{"me.chromosome"}=$chr if $chr;
    my @data;
    my @dsids;
    push @dsids, $dsid if $dsid;
    if ($dsgid)
      {
	my $dsg = $coge->resultset('DatasetGroup')->find($dsgid);
	unless ($dsg)
	  {
	    my $error =  "unable to create dsg object using id $dsgid\n";
	    return $error;
	  }
	$gstid = $dsg->type->id;
	foreach my $ds ($dsg->datasets())
	  {
	    push @dsids, $ds->id;
	  }
      }
    my %seqs; #let's prefetch the sequences with one call to genomic_sequence (slow for many seqs)
    foreach my $dsidt (@dsids)
      {
	my $ds = $coge->resultset('Dataset')->find($dsidt);
	unless ($ds)
	  {
	    warn "no dataset object found for id $dsidt\n";
	    next;
	  }
	
	if ($chr)
	  {
	    $seqs{$chr} = $ds->genomic_sequence(chr=>$chr, seq_type=>$gstid);
	  }
	else
	  {
	    map {$seqs{$_}= $ds->genomic_sequence(chr=>$_, seq_type=>$gstid)} $ds->chromosomes;
	  }
	foreach my $feat ($ds->features($search,{join=>['locations', {'dataset'=>{'dataset_connectors'=>'dataset_group'}}],
						 prefetch=>['locations',{'dataset'=>{'dataset_connectors'=>'dataset_group'}}]}))
	  {
	    foreach my $loc ($feat->locations)
	      {
		if ($loc->stop > length ($seqs{$feat->chromosome}))
		  {
		    print STDERR "feature ".$feat->id ." stop exceeds sequence length: ".$loc->stop." :: ".length($seqs{$feat->chromosome}),"\n";
		  }
		substr($seqs{$feat->chromosome}, $loc->start-1,($loc->stop-$loc->start+1)) = "-"x($loc->stop-$loc->start+1);
	      }
#	    push @data, sprintf("%.2f",100*$gc[0]/$total) if $total;
	  }
      }
    foreach my $seq (values %seqs)
      {
	$gc += $seq=~tr/GCgc/GCgc/;
	$at += $seq=~tr/ATat/ATat/;
	$n += $seq =~ tr/xnXN/xnXN/;
      }
    my $total = $gc+$at+$n;
    return "error" unless $total;
    return commify($total)." bp -- GC: ".sprintf("%.2f",100*$gc/($total))."%  AT: ".sprintf("%.2f",100*$at/($total))."% N: ".sprintf("%.2f",100*($n)/($total))."%";





    my $file = $TEMPDIR."/".join ("_",@dsids)."_wobble_gc.txt";
    open(OUT, ">".$file);
    print OUT "#wobble gc for dataset ids: ".join (" ", @dsids),"\n";
    print OUT join ("\n", @data),"\n";
    close OUT;
    my $cmd = $HISTOGRAM;
    $cmd .= " -f $file";
    my $out = $file;
    $out =~ s/txt$/png/;
    $cmd .= " -o $out";
    $cmd .= " -t \"CDS wobble gc content\"";
    $cmd .= " -min 0";
    $cmd .= " -max 100";
    `$cmd`;
    my $info =  "<div class = small>Total: ".commify($total)." codons, GC: ".sprintf("%.2f",100*$gc/($total))."%  AT: ".sprintf("%.2f",100*$at/($total))."%  N: ".sprintf("%.2f",100*($n)/($total))."%</div>";
    $out =~ s/$TEMPDIR/$TEMPURL/;
    my $hist_img = "<img src=\"$out\">";

    return $info, $hist_img;
  }

sub get_codon_usage
  {
    my %args = @_;
    my $dsid = $args{dsid};
    my $chr = $args{chr};
    my $dsgid = $args{dsgid};
    my $gstid = $args{gstid};
    return unless $dsid || $dsgid;

    my $search;
    $search = {"feature_type.name"=>"CDS"};
    $search->{"me.chromosome"}=$chr if $chr;

    my @dsids;
    push @dsids, $dsid if $dsid;
    if ($dsgid)
      {
	my $dsg = $coge->resultset('DatasetGroup')->find($dsgid);
	unless ($dsg)
	  {
	    my $error =  "unable to create dsg object using id $dsgid\n";
	    return $error;
	  }
	$gstid = $dsg->type->id;
	foreach my $ds ($dsg->datasets())
	  {
	    push @dsids, $ds->id;
	  }
      }
    my %codons;
    my $codon_total = 0;
    my $feat_count = 0;
    my ($code, $code_type);

    foreach my $dsidt (@dsids)
      {
	my $ds = $coge->resultset('Dataset')->find($dsidt);
	my %seqs; #let's prefetch the sequences with one call to genomic_sequence (slow for many seqs)
	if ($chr)
	  {
	    $seqs{$chr} = $ds->genomic_sequence(chr=>$chr, seq_type=>$gstid);
	  }
	else
	  {
	    %seqs= map {$_, $ds->genomic_sequence(chr=>$_, seq_type=>$gstid)} $ds->chromosomes;
	  }
	foreach my $feat ($ds->features($search,{join=>["feature_type",'locations', {'dataset'=>{'dataset_connectors'=>'dataset_group'}}],
						 prefetch=>['locations',{'dataset'=>{'dataset_connectors'=>'dataset_group'}}]}))
	  {
	    my $seq = substr($seqs{$feat->chromosome}, $feat->start-1, $feat->stop-$feat->start+1);
	    $feat->genomic_sequence(seq=>$seq);
	    $feat_count++;
	    ($code, $code_type) = $feat->genetic_code() unless $code;
	    my ($codon) = $feat->codon_frequency(counts=>1);
	    grep {$codon_total+=$_} values %$codon;
	    grep {$codons{$_}+=$codon->{$_}} keys %$codon;
	    print STDERR ".($feat_count)" if !$feat_count%10;
	  }
      }
    %codons = map {$_,$codons{$_}/$codon_total} keys %codons;

    #Josh put some stuff in here so he could get raw numbers instead of percentages for aa usage. He should either make this an option or delete this code when he is done. REMIND HIM ABOUT THIS IF YOU ARE EDITING ORGVIEW!
    my $html = "Codon Usage: $code_type";
    $html .= CoGe::Accessory::genetic_code->html_code_table(data=>\%codons, code=>$code);
    return $html
  }

sub get_aa_usage
  {
    my %args = @_;
    my $dsid = $args{dsid};
    my $chr = $args{chr};
    my $dsgid = $args{dsgid};
    my $gstid = $args{gstid};
    return unless $dsid || $dsgid;

    my $search;
    $search = {"feature_type.name"=>"CDS"};
    $search->{"me.chromosome"}=$chr if $chr;

    my @dsids;
    push @dsids, $dsid if $dsid;
    if ($dsgid)
      {
	my $dsg = $coge->resultset('DatasetGroup')->find($dsgid);
	unless ($dsg)
	  {
	    my $error =  "unable to create dsg object using id $dsgid\n";
	    return $error;
	  }
	$gstid = $dsg->type->id;
	foreach my $ds ($dsg->datasets())
	  {
	    push @dsids, $ds->id;
	  }
      }
    my %codons;
    my $codon_total = 0;
    my %aa;
    my $aa_total=0;
    my $feat_count = 0;
    my ($code, $code_type);

    foreach my $dsidt (@dsids)
      {
	my $ds = $coge->resultset('Dataset')->find($dsidt);
	my %seqs; #let's prefetch the sequences with one call to genomic_sequence (slow for many seqs)
	if ($chr)
	  {
	    $seqs{$chr} = $ds->genomic_sequence(chr=>$chr, seq_type=>$gstid);
	  }
	else
	  {
	    %seqs= map {$_, $ds->genomic_sequence(chr=>$_, seq_type=>$gstid)} $ds->chromosomes;
	  }
	foreach my $feat ($ds->features($search,{join=>["feature_type",'locations', {'dataset'=>{'dataset_connectors'=>'dataset_group'}}],
						 prefetch=>['locations',{'dataset'=>{'dataset_connectors'=>'dataset_group'}}]}))
	  {
	    my $seq = substr($seqs{$feat->chromosome}, $feat->start-1, $feat->stop-$feat->start+1);
	    $feat->genomic_sequence(seq=>$seq);
	    $feat_count++;
	    ($code, $code_type) = $feat->genetic_code() unless $code;
	    my ($codon) = $feat->codon_frequency(counts=>1);
	    grep {$codon_total+=$_} values %$codon;
	    grep {$codons{$_}+=$codon->{$_}} keys %$codon;
	    foreach my $tri (keys %$code)
	      {
		$aa{$code->{$tri}}+=$codon->{$tri};
		$aa_total+=$codon->{$tri};
	      }
	    print STDERR ".($feat_count)" if !$feat_count%10;
	  }
      }
    %codons = map {$_,$codons{$_}/$codon_total} keys %codons;

    #Josh put some stuff in here so he could get raw numbers instead of percentages for aa usage. He should either make this an option or delete this code when he is done. REMIND HIM ABOUT THIS IF YOU ARE EDITING ORGVIEW!
    %aa = $USER->user_name =~ /jkane/i ? map {$_,$aa{$_}} keys %aa : map {$_,$aa{$_}/$aa_total} keys %aa;
    
#    my $html1 = "Codon Usage: $code_type";
#    $html1 .= CoGe::Accessory::genetic_code->html_code_table(data=>\%codons, code=>$code);
    
    my $html2 .= "Predicted amino acid usage using $code_type";
    $html2 .= "<br/>Total Amino Acids: $aa_total" if $USER->user_name =~ /jkane/i;
    $html2 .= CoGe::Accessory::genetic_code->html_aa(data=>\%aa);
    $html2 =~ s/00.00%//g if $USER->user_name =~ /jkane/i;
    return $html2;
#    return $html1, $html2;
  }

sub get_wobble_gc
  {
    my %args = @_;
    my $dsid = $args{dsid};
    my $dsgid = $args{dsgid};
    my $chr = $args{chr};
    my $gstid = $args{gstid}; #genomic sequence type id
    return "error" unless $dsid || $dsgid;
    my $gc = 0;
    my $at = 0;
    my $n = 0;
    my $search;
    $search = {"feature_type_id"=>3};
    $search->{"me.chromosome"}=$chr if $chr;
    my @data;
    my @dsids;
    push @dsids, $dsid if $dsid;
    if ($dsgid)
      {
	my $dsg = $coge->resultset('DatasetGroup')->find($dsgid);
	unless ($dsg)
	  {
	    my $error =  "unable to create dsg object using id $dsgid\n";
	    return $error;
	  }
	$gstid = $dsg->type->id;
	foreach my $ds ($dsg->datasets())
	  {
	    push @dsids, $ds->id;
	  }
      }
    foreach my $dsidt (@dsids)
      {
	my $ds = $coge->resultset('Dataset')->find($dsidt);
	unless ($ds)
	  {
	    warn "no dataset object found for id $dsidt\n";
	    next;
	  }
	foreach my $feat ($ds->features($search,{join=>['locations', {'dataset'=>{'dataset_connectors'=>'dataset_group'}}],
						 prefetch=>['locations',{'dataset'=>{'dataset_connectors'=>'dataset_group'}}],
						}))
	  {
	    my @gc = $feat->wobble_content(counts=>1);
	    $gc+=$gc[0] if $gc[0] && $gc[0] =~ /^\d+$/;
	    $at+=$gc[1] if $gc[1] && $gc[1] =~ /^\d+$/;
	    $n+=$gc[2] if $gc[2] && $gc[2] =~ /^\d+$/;
	    my $total = 0;
	    $total += $gc[0] if $gc[0];
	    $total += $gc[1] if $gc[1];
	    $total += $gc[2] if $gc[2];
	    push @data, sprintf("%.2f",100*$gc[0]/$total) if $total;
	  }
      }
    my $total = $gc+$at+$n;
    return "error" unless $total;
    
    my $file = $TEMPDIR."/".join ("_",@dsids)."_wobble_gc.txt";
    open(OUT, ">".$file);
    print OUT "#wobble gc for dataset ids: ".join (" ", @dsids),"\n";
    print OUT join ("\n", @data),"\n";
    close OUT;
    my $cmd = $HISTOGRAM;
    $cmd .= " -f $file";
    my $out = $file;
    $out =~ s/txt$/png/;
    $cmd .= " -o $out";
    $cmd .= " -t \"CDS wobble gc content\"";
    $cmd .= " -min 0";
    $cmd .= " -max 100";
    `$cmd`;
    my $info =  "<div class = small>Total: ".commify($total)." codons, GC: ".sprintf("%.2f",100*$gc/($total))."%  AT: ".sprintf("%.2f",100*$at/($total))."%  N: ".sprintf("%.2f",100*($n)/($total))."%</div>";
    $out =~ s/$TEMPDIR/$TEMPURL/;
    my $hist_img = "<img src=\"$out\">";
    return $info."<br>". $hist_img;
  }

sub get_wobble_gc_diff
  {
    my %args = @_;
    my $dsid = $args{dsid};
    my $dsgid = $args{dsgid};
    my $chr = $args{chr};
    my $gstid = $args{gstid}; #genomic sequence type id
    return "error"," " unless $dsid || $dsgid;
    my $search;
    $search = {"feature_type_id"=>3};
    $search->{"me.chromosome"}=$chr if $chr;
    my @data;
    my @dsids;
    push @dsids, $dsid if $dsid;
    if ($dsgid)
      {
	my $dsg = $coge->resultset('DatasetGroup')->find($dsgid);
	unless ($dsg)
	  {
	    my $error =  "unable to create dsg object using id $dsgid\n";
	    return $error;
	  }
	$gstid = $dsg->type->id;
	foreach my $ds ($dsg->datasets())
	  {
	    push @dsids, $ds->id;
	  }
      }
    foreach my $dsidt (@dsids)
      {
	my $ds = $coge->resultset('Dataset')->find($dsidt);
	unless ($ds)
	  {
	    warn "no dataset object found for id $dsidt\n";
	    next;
	  }
	foreach my $feat ($ds->features($search,{join=>['locations', {'dataset'=>{'dataset_connectors'=>'dataset_group'}}],
						 prefetch=>['locations',{'dataset'=>{'dataset_connectors'=>'dataset_group'}}]}
				       ))
	  {
	    my @wgc = $feat->wobble_content();
	    my @gc = $feat->gc_content();
	    my $diff = $gc[0]-$wgc[0] if defined $gc[0] && defined $wgc[0];
	    push @data, sprintf("%.2f", 100*$diff) if $diff;
	  }
      }
    return "error"," " unless @data;
    my $file = $TEMPDIR."/".join ("_",@dsids)."_wobble_gc_diff.txt";
    open(OUT, ">".$file);
    print OUT "#wobble gc for dataset ids: ".join (" ", @dsids),"\n";
    print OUT join ("\n", @data),"\n";
    close OUT;
    my $cmd = $HISTOGRAM;
    $cmd .= " -f $file";
    my $out = $file;
    $out =~ s/txt$/png/;
    $cmd .= " -o $out";
    $cmd .= " -t \"CDS GC - wobble gc content\"";
    `$cmd`;
    my $sum=0;
    map {$sum+=$_}@data;
    my $mean = sprintf ("%.2f", $sum/scalar @data);
    my $info = "Mean $mean%";
    $info .= " (";
    $info .= $mean > 0 ? "CDS" : "wobble";
    $info .= " is more GC rich)";
    $out =~ s/$TEMPDIR/$TEMPURL/;
    my $hist_img = "<img src=\"$out\">";
    return $info."<br>". $hist_img;
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
