#! /usr/bin/perl -w


use strict;
use CGI;
use CoGeX;
use DBI;

use Data::Dumper;
use CoGe::Accessory::LogUser;
use CoGe::Accessory::Web;
use HTML::Template;
use URI::Escape;
use Spreadsheet::WriteExcel;
use Benchmark;
use DBIxProfiler;
use File::Path;
no warnings 'redefine';


use vars qw($P $DBNAME $DBHOST $DBPORT $DBUSER $DBPASS $connstr $PAGE_NAME $TEMPDIR $USER $DATE $BASEFILE $coge $cogeweb $FORM $URL $HISTOGRAM $TEMPURL);
$P = CoGe::Accessory::Web::get_defaults($ENV{HOME}.'coge.conf');
$ENV{PATH} = $P->{COGEDIR};
$URL = $P->{URL};
$DATE = sprintf( "%04d-%02d-%02d %02d:%02d:%02d",
		sub { ($_[5]+1900, $_[4]+1, $_[3]),$_[2],$_[1],$_[0] }->(localtime));
$PAGE_NAME = "GenomeList.pl";
$HISTOGRAM = $P->{HISTOGRAM};	
($USER) = CoGe::Accessory::LogUser->get_user();
$TEMPDIR = $P->{TEMPDIR}."GenomeList/";
mkpath ($TEMPDIR, 0,0777) unless -d $TEMPDIR;
$TEMPURL = $P->{TEMPURL}."GenomeList/";
$FORM = new CGI;
$DBNAME = $P->{DBNAME};
$DBHOST = $P->{DBHOST};
$DBPORT = $P->{DBPORT};
$DBUSER = $P->{DBUSER};
$DBPASS = $P->{DBPASS};
$connstr = "dbi:mysql:dbname=".$DBNAME.";host=".$DBHOST.";port=".$DBPORT;
$coge = CoGeX->connect($connstr, $DBUSER, $DBPASS );
#$coge->storage->debugobj(new DBIxProfiler());
#$coge->storage->debug(1);



$SIG{'__WARN__'} = sub { }; #silence warnings

my %FUNCTION = (
		gen_html=>\&gen_html,
		get_feature_counts=>\&get_feature_counts,
		get_gc=>\&get_gc,
		get_aa_usage=>\&get_aa_usage,
		get_codon_usage=>\&get_codon_usage,
	 	get_fasta_seqs=>\&get_fasta_seqs,
	 	generate_excel_file=>\&generate_excel_file,
	 	generate_csv_file=>\&generate_csv_file,
	 	gc_content=>\&gc_content,
	 	gen_data=>\&gen_data,
	 	send_to_msa=>\&send_to_msa,
	 	send_to_SynFind=>\&send_to_SynFind,
	 	send_to_CoGeBlast=>\&send_to_CoGeBlast,
	 	send_to_GenomeList=>\&send_to_GenomeList,
	 	send_to_SynFind=>\&send_to_SynFind,
	 	get_wobble_gc=>\&get_wobble_gc,
	 	save_FeatList_settings=>\&save_FeatList_settings,
	 	add_to_user_history=>\&add_to_user_history,
		cds_wgc_hist=>\&cds_wgc_hist,
		get_gc_for_feature_type=>\&get_gc_for_feature_type,
    );
#my $pj = new CGI::Ajax(%FUNCTION);
#$pj->js_encode_function('escape');
#my $t1 = new Benchmark;
if ($FORM->param('jquery_ajax'))
  {
    dispatch();
  }
else
  {
    print $FORM->header,"\n",gen_html();
  }
sub dispatch
{
    my %args = $FORM->Vars;
    my $fname = $args{'fname'};
    if($fname)
    {
	#my %args = $cgi->Vars;
	#print STDERR Dumper \%args;
	if($args{args}){
	    my @args_list = split( /,/, $args{args} );
	    print $FORM->header, $FUNCTION{$fname}->(@args_list);
       	}
	else{
	    print $FORM->header, $FUNCTION{$fname}->(%args);
	}
    }
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
       my ($body) = gen_body();
       my $template = HTML::Template->new(filename=>$P->{TMPLDIR}.'generic_page.tmpl');
       $template->param(PAGE_TITLE=>'GenomeList');
       $template->param(HELP=>'/wiki/index.php?title=GenomeList');
       # print STDERR "user is: ",$USER,"\n";
       #add_to_user_history() unless $USER->user_name eq "public";
       my $name = $USER->user_name;
       $name = $USER->first_name if $USER->first_name;
       $name .= " ".$USER->last_name if $USER->first_name && $USER->last_name;
       $template->param(USER=>$name);
       $template->param(LOGO_PNG=>"GenomeList-logo.png");
       $template->param(LOGON=>1) unless $USER->user_name eq "public";
       $template->param(DATE=>$DATE);
       my $list_name = $FORM->param('list_name') || $FORM->param('ln');
       my $box_name = "Genome List:";
       $box_name .= " $list_name" if $list_name;
       $template->param(BOX_NAME=>$box_name);
       $template->param(BODY=>$body);
       $template->param(ADJUST_BOX=>1);
       $html .= $template->output;
   }
 }
 




sub cds_wgc_hist
  {
    my %opts = @_;
    my $dsid = $opts{dsid};
    my $dsgid = $opts{dsgid};
    my $chr = $opts{chr};
    my $gstid = $opts{gstid}; #genomic sequence type id
    my $min = $opts{min}; #limit results with gc values greater than $min;
    my $max = $opts{max}; #limit results with gc values smaller than $max;
    my $hist_type = $opts{hist_type};
    return "Error: No dsid or dsgid passed." unless $dsid || $dsgid;
    my $gc = 0;
    my $at = 0;
    my $n = 0;
    my $search;
    $search = {"feature_type_id"=>3};
    $search->{"me.chromosome"}=$chr if defined $chr;
    my @data;
    my @fids;
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
	    my $perc_gc = 100*$gc[0]/$total if $total;
	    next unless $perc_gc; #skip if no values
	    next if defined $min && $min =~/\d+/ && $perc_gc < $min; #check for limits
	    next if defined $max && $max =~/\d+/ && $perc_gc > $max; #check for limits
	    push @data, sprintf("%.2f",$perc_gc);
	    push @fids, $feat->id."_".$gstid;
	    #push @data, sprintf("%.2f",100*$gc[0]/$total) if $total;
	  }
      }
    my $total = $gc+$at+$n;
    return "error" unless $total;
    
    my $file = $TEMPDIR."/".join ("_",@dsids);#."_wobble_gc.txt";
    ($min) = $min =~ /(.*)/ if defined $min;
    ($max) = $max =~ /(.*)/ if defined $max;
    ($chr) = $chr =~ /(.*)/ if defined $chr;
    $file .= "_".$chr."_" if defined $chr;
    $file .= "_min".$min if defined $min;
    $file .= "_max".$max if defined $max;
    $file .= "_$hist_type" if $hist_type;
    $file .= "_wobble_gc.txt";
    my $out = $file;
    $out =~ s/txt$/png/;
    unless (-r $out)
      {
	open(OUT, ">".$file);
	print OUT "#wobble gc for dataset ids: ".join (" ", @dsids),"\n";
	print OUT join ("\n", @data),"\n";
	close OUT;
	my $cmd = $HISTOGRAM;
	$cmd .= " -f $file";
	$cmd .= " -o $out";
	$cmd .= " -t \"CDS wobble gc content\"";
	$cmd .= " -min 0";
	$cmd .= " -max 100";
	$cmd .= " -ht $hist_type" if $hist_type;
	`$cmd`;
      }
    $min = 0 unless defined $min && $min =~/\d+/;
    $max = 100 unless defined $max && $max =~/\d+/;
    my $info;
    $info .= qq{<div class="small">
Min: <input type="text" size="3" id="wobble_gc_min" value="$min">
Max: <input type=text size=3 id=wobble_gc_max value=$max>
Type: <select id=wobble_hist_type>
<option value ="counts">Counts</option>
<option value = "percentage">Percentage</option>
</select>
};
    $info =~ s/>Per/ selected>Per/ if $hist_type =~/per/;
    my $args;
    $args .= "'args__dsid','ds_id'," if $dsid;
    $args .= "'args__dsgid','dsg_id'," if $dsgid;
    $args .= "'args__chr','chr'," if defined $chr;
    $args .= "'args__min','wobble_gc_min',";
    $args .= "'args__max','wobble_gc_max',";
    $args .= "'args__max','wobble_gc_max',";
    $args .= "'args__hist_type', 'wobble_hist_type',";
    $info .= qq{<span class="link" onclick="get_wobble_gc([$args],['wobble_gc_histogram']);\$('#wobble_gc_histogram').html('loading...');">Regenerate histogram</span>};
    $info .= "</div>";

    $info .=  "<div class = small>Total: ".commify($total)." codons.  Mean GC: ".sprintf("%.2f",100*$gc/($total))."%  AT: ".sprintf("%.2f",100*$at/($total))."%  N: ".sprintf("%.2f",100*($n)/($total))."%</div>";
    if ($min || $max)
      {
	$min = 0 unless defined $min;
	$max = 100 unless defined $max;
	$info .= qq{<div class=small style="color: red;">Limits set:  MIN: $min  MAX: $max</div>
} 
  }
    my $stuff = join "::",@fids;
    $info .= qq{<div class="link small" onclick="window.open('FeatList.pl?fid=$stuff')">Open FeatList of Features</div>};
    $out =~ s/$TEMPDIR/$TEMPURL/;
    my $hist_img = "<img src=\"$out\">";
    return $info."<br>". $hist_img;
  }

sub get_gc_for_feature_type
  {
    my %opts = @_;
    my $dsid = $opts{dsid};
    my $dsgid = $opts{dsgid};
    my $chr = $opts{chr};
    my $typeid = $opts{typeid};
    my $gstid = $opts{gstid};#genomic sequence type id
    my $min = $opts{min}; #limit results with gc values greater than $min;
    my $max = $opts{max}; #limit results with gc values smaller than $max;
    my $hist_type = $opts{hist_type};
    $hist_type = "counts" unless $hist_type;
    $min = undef if $min && $min eq "undefined";
    $max = undef if $max && $max eq "undefined";
    $chr = undef if $chr && $chr eq "undefined";
    $dsid = undef if $dsid && $dsid eq "undefined";
    $hist_type = undef if $hist_type && $hist_type eq "undefined";
    $typeid = 1 if $typeid eq "undefined";
    return unless $dsid || $dsgid;
    my $gc = 0;
    my $at = 0;
    my $n = 0;
    my $type = $coge->resultset('FeatureType')->find($typeid);
    my @data;
    my @fids; #storage for fids that passed.  To be sent to FeatList
    my @dsids;
    push @dsids, $1 if $dsid && $dsid =~ /(\d+)/;
    if ($dsgid)
      {
	my $dsg = $coge->resultset('DatasetGroup')->find($dsgid);
	unless ($dsg)
	  {
	    my $error =  "unable to create dsg object using id $dsgid\n";
	    return $error;
	  }
	$gstid = $dsg->type->id;
	if (!$dsid)
	  {
	    foreach my $ds ($dsg->datasets())
	      {
		push @dsids, $ds->id;
	      }
	  }
      }
    my $search;
    $search = {"feature_type_id"=>$typeid};
    $search->{"me.chromosome"}=$chr if defined $chr;
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
	if (defined $chr)
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
 	foreach my $feat (@feats)
 	  {
	    my $seq = substr($seqs{$feat->chromosome}, $feat->start-1, $feat->stop-$feat->start+1);

	    $feat->genomic_sequence(seq=>$seq);
	    my @gc = $feat->gc_content(counts=>1);

	    $gc+=$gc[0] if $gc[0] =~ /^\d+$/;
	    $at+=$gc[1] if $gc[1] =~ /^\d+$/;
	    $n+=$gc[2] if $gc[2] =~ /^\d+$/;
	    my $total = 0;
	    $total += $gc[0] if $gc[0];
	    $total += $gc[1] if $gc[1];
	    $total += $gc[2] if $gc[2];
	    my $perc_gc = 100*$gc[0]/$total if $total;
	    next unless $perc_gc; #skip if no values
	    next if defined $min && $min =~/\d+/ && $perc_gc < $min; #check for limits
	    next if defined $max && $max =~/\d+/ && $perc_gc > $max; #check for limits
	    push @data, sprintf("%.2f",$perc_gc);
	    push @fids, $feat->id."_".$gstid;
	  }
	my $t3 = new Benchmark;
	my $get_seq_time = timestr(timediff($t2,$t1));
	my $process_seq_time = timestr(timediff($t3,$t2));
       }
    my $total = $gc+$at+$n;
    return "error" unless $total;

    my $file = $TEMPDIR."/".join ("_",@dsids);
    #perl -T flag
    ($min) = $min =~ /(.*)/ if defined $min;
    ($max) = $max =~ /(.*)/ if defined $max;
    ($chr) = $chr =~ /(.*)/ if defined $chr;
    $file .= "_".$chr."_" if defined $chr;
    $file .= "_min".$min if defined $min;
    $file .= "_max".$max if defined $max;
    $file .= "_$hist_type" if $hist_type;
    $file .= "_".$type->name."_gc.txt";
    my $out = $file;
    $out =~ s/txt$/png/;
    unless (-r $out)
      {
	open(OUT, ">".$file);
	print OUT "#wobble gc for dataset ids: ".join (" ", @dsids),"\n";
	print OUT join ("\n", @data),"\n";
	close OUT;
	my $cmd = $HISTOGRAM;
	$cmd .= " -f $file";
	$cmd .= " -o $out";
	$cmd .= " -t \"".$type->name." gc content\"";
	$cmd .= " -min 0";
	$cmd .= " -max 100";
	$cmd .= " -ht $hist_type" if $hist_type;
	`$cmd`;
      }
	
    $min = 0 unless defined $min && $min =~/\d+/;
    $max = 100 unless defined $max && $max =~/\d+/;
    my $info;
    $info .= qq{<div class="small">
Min: <input type="text" size="3" id="feat_gc_min" value="$min">
Max: <input type=text size=3 id=feat_gc_max value=$max>
Type: <select id=feat_hist_type>
<option value ="counts">Counts</option>
<option value = "percentage">Percentage</option>
</select>
};
    $info =~ s/>Per/ selected>Per/ if $hist_type =~/per/;
    my $gc_args;
    $gc_args = "chr: '$chr'," if defined $chr;
    $gc_args .= "dsid: $dsid," if $dsid; #set a var so that histograms are only calculated for the dataset and not hte dataset_group
    $gc_args .= "typeid: '$typeid'";
    $info .= qq{<span class="link" onclick="get_feat_gc({$gc_args})">Regenerate histogram</span>};
    $info .= "</div>";
    $info .= "<div class = small>Total length: ".commify($total)." bp, GC: ".sprintf("%.2f",100*$gc/($total))."%  AT: ".sprintf("%.2f",100*$at/($total))."%  N: ".sprintf("%.2f",100*($n)/($total))."%</div>";
    if ($min || $max)
      {
	$min = 0 unless defined $min;
	$max = 100 unless defined $max;
	$info .= qq{<div class=small style="color: red;">Limits set:  MIN: $min  MAX: $max</div>
} 
      }
    my $stuff = join "::",@fids;
    $info .= qq{<div class="link small" onclick="window.open('FeatList.pl?fid=$stuff')">Open FeatList of Features</div>};
    
    $out =~ s/$TEMPDIR/$TEMPURL/;
    $info .= "<br><img src=\"$out\">";
    return $info;
  }


sub gen_body
  {
    my $template = HTML::Template->new(filename=>$P->{TMPLDIR}.'GenomeList.tmpl');
    my $form = $FORM;
    my $no_values;
    $BASEFILE = $form->param('basename');
    my $sort_by_type = $form->param('sort_type');
    my $sort_by_location = $form->param('sort_loc');
    my $prefs =CoGe::Accessory::Web::load_settings(user=>$USER, page=>$PAGE_NAME);
    $prefs = {} unless $prefs;
    $prefs->{display}={
		       NameD=>1,
		       DescD=>1,
		       TypeD=>1,
		       SourceD=>1,
		       ProvenanceD=>1,
		       VerD=>1,
		       ChrCountD=>1,
		       LengthD=>1,
		       GCD=>1,
		       ATD=>1,
		       AAD=>1,
		       ND=>1,
		       FeatD=>1,
		       CodonD=>1,
		       CDSGCHistD=>1,
		       CDSwGCHistD=>1,
		      } unless $prefs->{display};
    foreach my $key (keys %{$prefs->{display}})
      {
	$template->param($key=>"checked");
      }

    $template->param(SAVE_DISPLAY=>1) unless $USER->user_name eq "public";
    #Don't show button to save list data unless valid user
    $template->param(SAVE_DATA=>1) unless $USER->user_name eq "public";

    my $dsgids = [];
    $dsgids = read_file() if $BASEFILE;#: $opts{feature_list};
    foreach my $item ($form->param('dsgid'))
      {
	foreach my $item2 (split/(::)|(,)/, $item)
	  {
	    push @$dsgids, $item2 if $item2 =~ /^\d+_?\d*$/;
	  }
      }
    my ($table, $count) = generate_table(dsgids=>$dsgids);

    $template->param('GENOME_COUNT'=>$count);

    if ($table)
      {
	$template->param(INFO=>$table);
	return $template->output;
      }
    else
      {
	return "No dataset_group (genome) ids were specified.";
      }
  }

sub add_to_user_history
  {
    my %opts = @_;
#    my $url_params = $ENV{'QUERY_STRING'};
    my $url = $opts{url};
    if ($opts{archive})
      {
	$USER->add_to_works({
			     'name'=>$opts{work_name},
			     'archive'=>$opts{archive},
			     'page'=>$PAGE_NAME,
			     'parameter'=>$url,
			     'description'=>$opts{description},
			     'note'=>$opts{note},
			    }); 
      }
    else
      {
	$USER->add_to_works({
			     'name'=>'FeatList-'.$DATE,
			     'archive'=>0,
			     'page'=>$PAGE_NAME,
			     'parameter'=>$url,
			     'description'=>'Feature List created on '.$DATE,
			    });
      }
    return $opts{work_name};
  }


sub generate_table
  {
    my %opts = @_;
    my $dsgids = $opts{dsgids};
    return unless @$dsgids;
    my @table;
    my $count = 1;
    foreach my $dsgid (@$dsgids)
      {
	my $dsg = $coge->resultset('DatasetGroup')->find($dsgid);
	my $name = $dsg->name ? $dsg->name : $dsg->organism->name;
	my $desc = join("; ", map { qq{<span class=link onclick=window.open('OrganismView.pl?org_desc=$_')>$_</span>} } split /;\s*/, $dsg->organism->description);
	
	my $chr_count = $dsg->chromosome_count;
	my $length = $dsg->length;
	my $type = $dsg->type->name;
	my ($ds_source) = $dsg->source;
	my $source = $ds_source->name;
	my $source_link = $ds_source->link;
	$source_link = "http://".$source_link unless !$source_link || $source_link =~ /http/;
	$source = qq{<span class=link onclick="window.open('}.$source_link.qq{')">$source</span>}if $source_link;
	my $provenance;
	foreach my $ds ($dsg->datasets)
	  {
	    my $item = $ds->name;
	    my $link = $ds->link;
	    $link = "http://".$link unless !$link || $link =~ /http/;
	    $item = qq{<span class=link onclick="window.open('}.$link.qq{')">$item</span>}if $link;
	    $provenance .= $item."<br>";
	  }
	$provenance =~ s/<br>$//;
	push @table,{
		     COUNT=>$count,
		     DSGID=>$dsgid,
		     NAME=>$name,
		     DESC=>$desc,
		     SOURCE=>$source,
		     PROVENANCE=>$provenance,
		     VER=>$dsg->version,
		     TYPE=>$type,
		     CHR_COUNT=>commify($chr_count),
		     LENGTH=>commify($length),
		
		    };
	$count++;
      }
    $count--;
    return \@table, $count;
  }

sub get_feature_counts
  {
    my %opts = @_;
    my $dsgid = $opts{dsgid};
    return "No specified dataset group id" unless $dsgid;
    my $dsg = $coge->resultset('DatasetGroup')->find($dsgid);
    my $query = qq{
SELECT count(distinct(feature_id)), ft.name, ft.feature_type_id
  FROM feature
  JOIN feature_type ft using (feature_type_id)
  JOIN dataset_connector dc using (dataset_id)
 WHERE dataset_group_id = $dsgid
  GROUP BY ft.name

};
    my $dbh = DBI->connect($connstr,$DBUSER,$DBPASS);
    my $sth = $dbh->prepare($query);
    $sth->execute;
    my $feats = {};
    while (my $row = $sth->fetchrow_arrayref)
      {
	my $name = $row->[1];
	$name =~ s/\s+/&nbsp/g;
	$feats->{$name} = {count=>$row->[0],
			   id=>$row->[2],
			   name=>$name,
			  };
      }
  my $feat_string .= qq{<table class="small ui-widget-content ui-corner-all">
<thead></thead><tbody>};
  $feat_string .= "<tr valign=top>". join ("\n<tr valign=top>",map {
      "<td valign=top><div id=$_>".$feats->{$_}{name}."</div>".
	    "<td valign=top align=right>".$feats->{$_}{count}
} sort {$a cmp $b} keys %$feats);
  $feat_string .= "</tbody></table>";
  $feat_string .= "None" unless keys %$feats;
  return $feat_string;
}


  sub gen_data
  {
    my $message = shift;
    return qq{<font class="loading">$message. . .</font>};
  }

  sub read_file
  {
    my $file = "$TEMPDIR/$BASEFILE.featlist";
    my @featlist;
    unless (-r $file)
      {
	warn "unable to read file $file for feature ids\n";
	return \@featlist;
      }
    open (IN, $file) || die "can't open $file for reading: $!";
    while (<IN>)
      {
	chomp;
	push @featlist,$_;
      }
    close IN;
    return \@featlist;
  }
  
  sub gevo #Send to GEvo
  {
    my $accn_list = shift;
    $accn_list =~ s/^,//;
    $accn_list =~ s/,$//;
    my $url = $URL."GEvo.pl?";
    my $count = 1;
    #print STDERR $url,"\n";
    foreach my $featid (split /,/,$accn_list)
      {
		$url .= "fid$count=$featid&";
		$count ++;
      }
    $count--;
    return ("alert",$count) if $count > 20;
    $url .= "num_seqs=$count";
    return $url;
  }
  
sub send_to_SynFind
  {
    my $accn_list = shift;
    $accn_list =~ s/^,//;
    $accn_list =~ s/,$//;
    my ($fid) = split/,/, $accn_list;
    my $url = $URL."SynFind.pl?fid=$fid";
    return $url;
  }
  
  
sub send_to_msa
  {
    my %opts = @_;
    my $accn_list = $opts{accn};;
    $accn_list =~ s/^,//;
    $accn_list =~ s/,$//;
    my $url = $URL."GenomeAlign.pl?dsgid=$accn_list";
    $url =~ s/&$//;
    return $url;
  }
  
  
sub send_to_CoGeBlast #send to cogeblast
  {
    my %opts = @_;
    my $accn_list = $opts{accn};
    $accn_list =~ s/^,//;
    $accn_list =~ s/,$//;
    my $url = $URL."CoGeBlast.pl?dsgid=$accn_list";
    return $url;
  }
  
sub send_to_GenomeList #send to GenomeList
  {
    my %opts = @_;
    my $accn_list = $opts{accn};
    $accn_list =~ s/^,//;
    $accn_list =~ s/,$//;
    my $url = $URL."GenomeList.pl?dsgid=$accn_list";
    return $url;
  }
  
sub send_to_SynFind #send to SynFind
  {
    my %opts = @_;
    my $accn_list = $opts{accn};
    $accn_list =~ s/^,//;
    $accn_list =~ s/,$//;
    my $url = $URL."SynFind.pl?dsgid=$accn_list";
    return $url;
  }
  
  sub get_fasta_seqs
  {
    my %opts = @_;
    my $accn_list = $opts{accn};
    
    $accn_list =~ s/^,//;
    $accn_list =~ s/,$//;
    $cogeweb = CoGe::Accessory::Web::initialize_basefile(prog=>"GenomeList");
    my $basename = $cogeweb->basefilename;
    my $file = $TEMPDIR."$basename.faa";
    open (OUT, ">$file");
    foreach my $dsgid (split /,/,$accn_list)
      {
	next unless $dsgid;
	my ($dsg) = $coge->resultset('DatasetGroup')->find($dsgid);
	next unless $dsg;
	print OUT $dsg->fasta;
      }
    close OUT;
    $file =~ s/$TEMPDIR/$TEMPURL/;
    return $file;
  }

sub generate_excel_file
  {
    my %args = @_;
    my $accn_list = $args{accn};
    $accn_list =~ s/^,//;
    $accn_list =~ s/,$//;
    $cogeweb = CoGe::Accessory::Web::initialize_basefile(prog=>"GenomeList");
    my $basename = $cogeweb->basefilename;
    my $file = "$TEMPDIR/Excel_$basename.xls";
    my $workbook = Spreadsheet::WriteExcel->new($file);
    $workbook->set_tempdir("$TEMPDIR");
    my $worksheet = $workbook->add_worksheet();
    my $i = 1;
    
    $worksheet->write(0,0,"Name");
    $worksheet->write(0,1,"Description");
    $worksheet->write(0,2,"Source");
    $worksheet->write(0,3,"Provenance");
    $worksheet->write(0,4,"Sequence Type");
    $worksheet->write(0,5,"Chr Count");
    $worksheet->write(0,6,"Length (bp)");
    $worksheet->write(0,7,"Percent GC");
    $worksheet->write(0,8,"Percent AT");
    $worksheet->write(0,9,"Percent N|X");
    $worksheet->write(0,10,"OrganismView Link");
	
    
   	
    foreach my $dsgid (split /,/,$accn_list)
      {
	
	my ($dsg) = $coge->resultset("DatasetGroup")->find($dsgid);
	
	next unless $dsg;
	my $name = $dsg->name ? $dsg->name : $dsg->organism->name;
	my $desc = $dsg->description ? $dsg->description :  $dsg->organism->description;
	my ($ds_source) = $dsg->source;
	my $source = $ds_source->name;
	my $provenance = join (" ", map {$_->name} $dsg->datasets);
	my $length = $dsg->length;
	my $chr_count = $dsg->chromosome_count;
	my $type = $dsg->type->name;
	my ($gc, $at, $n) = $dsg->percent_gc();
	$at*=100;
	$gc*=100;
	#my ($wgc, $wat) = $dsg->wobble_content;
	#$wat*=100;
	#$wgc*=100;
	
	$worksheet->write($i,0,$name);
	$worksheet->write($i,1,$desc);
	$worksheet->write($i,2,$source);
	$worksheet->write($i,3,$provenance);
	$worksheet->write($i,4,$type);
	$worksheet->write($i,5,$chr_count);
	$worksheet->write($i,6,$length);
	$worksheet->write($i,7,$gc.'%');
	$worksheet->write($i,8,$at.'%');
	$worksheet->write($i,9,$n.'%');
	$worksheet->write($i,10,$P->{SERVER}.'OrganismView.pl?dsgid='.$dsgid);
	
	$i++;
      };
    $workbook->close() or die "Error closing file: $!";
    $file =~ s/$TEMPDIR/$TEMPURL/;
    return $file;
  }
  sub generate_csv_file
  {
    my %args = @_;
    my $accn_list = $args{accn};
    $accn_list =~ s/^,//;
    $accn_list =~ s/,$//;
    $cogeweb = CoGe::Accessory::Web::initialize_basefile(prog=>"GenomeList");
    my $basename = $cogeweb->basefilename;
    my $file = "$TEMPDIR/$basename.csv";
    open (OUT, ">$file");
    print OUT join ("\t", "COGE_DSG_ID", "Name", "Description","Source","Provenance","Sequence Type","Chr Count","Length (bp)","Percent GC","Percent AT","Percent N|X","OrganismView Link"),"\n";
    foreach my $dsgid (split /,/,$accn_list)
      {
	next unless $dsgid;
	my ($dsg) = $coge->resultset("DatasetGroup")->find($dsgid);	
	next unless $dsg;
	my $name = $dsg->name ? $dsg->name : $dsg->organism->name;
	my $desc = $dsg->description ? $dsg->description :  $dsg->organism->description;
	my ($ds_source) = $dsg->source;
	my $source = $ds_source->name;
	my $provenance = join ("||", map {$_->name} $dsg->datasets);
	my $chr_count = $dsg->chromosome_count;
	my $length = $dsg->length;
	my $type = $dsg->type->name;
	my ($gc, $at, $n) = $dsg->percent_gc();
	$at*=100;
	$gc*=100;
	print OUT join ("\t", $dsgid, $name,$desc,$source, $provenance,$type,$chr_count,$length,$gc,$at,$n,$P->{SERVER}.'OrganismView.pl?dsgid='.$dsgid),"\n";
      };
    close OUT;
    $file =~ s/$TEMPDIR/$TEMPURL/;
    return $file;
  }
  
sub gc_content
  {
    my %args = @_;
    my $featid = $args{featid};
    return unless $featid;
    my ($feat) = $coge->resultset('Feature')->find($featid);
    my ($gc, $at) = $feat->gc_content;
    my $html = "GC:".(100*$gc)."%".", AT:".(100*$at)."%" ;
    return $html;
  }

sub get_codon_usage
  {
    my %opts = @_;
    my $dsid = $opts{dsid};
    my $chr = $opts{chr};
    my $dsgid = $opts{dsgid};
    my $gstid = $opts{gstid};
    return unless $dsid || $dsgid;

    my $search;
    $search = {"feature_type.name"=>"CDS"};
    $search->{"me.chromosome"}=$chr if defined $chr;

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
	if (defined $chr)
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
    my %opts = @_;
    my $dsgid = $opts{dsgid};
    my $table_name = $opts{table_name};
    return "no dsgid specified" unless $dsgid;
    my $dsg = $coge->resultset('DatasetGroup')->find($dsgid);
    return "Unable able to create DB object for $dsgid" unless $dsg;
    my %codons;
    my $codon_total = 0;
    my %aa;
    my $aa_total=0;
    my $feat_count = 0;
    my ($code, $code_type);
    my $trans_type = $dsg->translation_type;
    foreach my $feat ($dsg->features({feature_type_id=>3}))#,{join=>['locations'],
#					     prefetch=>['locations']}))
      {
	$feat->trans_type($trans_type) if $trans_type;
	$feat_count++;
	($code, $code_type) = $feat->genetic_code() unless $code;
	my $aa = $feat->aa_frequency(counts=>1);
	while (my ($k, $v) = each %$aa)
	  {
	    $aa{$k}+=$v;
	  }
	$trans_type = $feat->trans_type;
      }
    my $html .= $dsg->organism->name."<br>Predicted amino acid usage using $code_type";
    $html .= CoGe::Accessory::genetic_code->html_aa_new(data=>\%aa, counts=>1, table_name=>$table_name);
    return ($html);
  }



sub codon_table
  {
    my %args = @_;
    my $featid = $args{fid};
    my $gstid = $args{gstid};
    ($featid, $gstid) = split/_/, $featid if $featid =~/_/;
    return unless $featid;
    my ($feat) = $coge->resultset('Feature')->find($featid);
    my ($codon, $code_type) = $feat->codon_frequency(counts=>1, gstid=>$gstid);
    my %aa;
    my ($code) = $feat->genetic_code;
    my $count = 0;
    foreach my $tri (keys %$code)
      {
	$aa{$code->{$tri}}+=$codon->{$tri};
	$count += $codon->{$tri};
      }
    my $html;
    $html .= "<table><tr valign=top><td>";
    $html .= "Codon Usage: $code_type<br>";
    my ($at, $gc) = $feat->gc_content(gstid=>$gstid);
    $at*=100;
    $gc*=100;
    my ($wat, $wgc) = $feat->wobble_content(gstid=>$gstid);
    $wat*=100;
    $wgc*=100;
    $html .= "Codon Count: $count".", GC: $at% $gc%".", Wobble GC: $wat% $wgc%";
    $html .="<td>";
    $html .= "Predicted amino acid usage";
    $html .= "<tr valign=top><td>";
    $html .= CoGe::Accessory::genetic_code->html_code_table(data=>$codon, code=>$code, counts=>1);
#    $html .= "</div>";
#    $html .= "Predicted amino acid usage for $code_type genetic code:";
    $html .= "<td>";
    $html .= CoGe::Accessory::genetic_code->html_aa(data=>\%aa, counts=>1, split=>1);
    $html .= "</table>";
    return $html;
  }



sub get_gc
  {
    my %opts = @_;
    my $dsgid = $opts{dsgid};
    return "No dataset group id" unless $dsgid;
    my ($dsg) = $coge->resultset('DatasetGroup')->find($dsgid);
    return "No dsg object for id $dsgid" unless $dsg;
    my ($gc, $at, $n) = $dsg->percent_gc();
    $at*=100;
    $gc*=100;
    $n *=100;
    return ($gc."_".$at."_".$n);
  }

sub get_wobble_gc
  {
    my %opts = @_;
    my $fid = $opts{fid};
    my $gstid = $opts{gstid};
    ($fid, $gstid) = split/_/, $fid if ($fid =~ /_/);
    return unless $fid;
    my ($feat) = $coge->resultset('Feature')->find($fid);
    return unless $feat->type->name eq "CDS";
    return "No feature for id $fid" unless $feat;
    my ($wgc, $wat) = $feat->wobble_content(gstid=>$gstid);
    $wat*=100;
    $wgc*=100;
    return ($wgc, $wat);
  }


sub save_FeatList_settings
  {
    my %opts = @_;
    my $display = $opts{display};
    my %save;
    if ($display)
      {
	my %settings = (
			1=>'FeatNameD',
			2=>'TypeD',
			3=>'ChrD',
			4=>'StartD',
			5=>'StopD',
			6=>'StrandD',
			7=>'LengthD',
			8=>'GCD',
			9=>'ATD',
			10=>'WGCD',
			11=>'WATD',
			12=>'OrgD',
			13=>'AnnoD',
			14=>'OtherD',
		       );
	foreach my $index (split/,/,$display)
	  {
	    $save{display}{$settings{$index}}=1
	  }
      }
   CoGe::Accessory::Web::save_settings(opts=>\%save, user=>$USER, page=>$PAGE_NAME);
  }

sub commify 
    {
      my $text = reverse $_[0];
      $text =~ s/(\d\d\d)(?=\d)(?!\d*\.)/$1,/g;
      return scalar reverse $text;
    }
