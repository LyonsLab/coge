#! /usr/bin/perl -w
use strict;
use CGI;
use CGI::Ajax;
use CoGe::Accessory::LogUser;
use CoGe::Accessory::Web;
use HTML::Template;
use Text::Wrap qw($columns &wrap);
use Data::Dumper;
use LWP::Simple;
use LWP::Simple::Post qw(post post_xml);
use URI::Escape;
use CoGeX;
use POSIX;
use Digest::MD5 qw(md5_hex);
use DBIxProfiler;
use File::Temp;
use File::Basename;
use CoGe::Accessory::blast_report;
use CoGe::Accessory::blastz_report;
use CoGe::Accessory::Restricted_orgs;
use CoGe::Graphics::GenomeView;
use CoGe::Graphics;
use CoGe::Graphics::Chromosome;
use CoGe::Graphics::Feature::HSP;
use Spreadsheet::WriteExcel;
use Benchmark qw(:all);
use Parallel::ForkManager;

$ENV{PATH} = "/opt/apache/CoGe/";
$ENV{BLASTDB}="/opt/apache/CoGe/data/blast/db/";
$ENV{BLASTMAT}="/opt/apache/CoGe/data/blast/matrix/";
use vars qw( $PAGE_NAME $TEMPDIR $TEMPURL $DATADIR $FASTADIR $BLASTDBDIR $FORMATDB $BLAST $BLASTZ $FORM $USER $DATE $coge $cogeweb $RESULTSLIMIT $MAX_PROC $connstr);
#refresh again?
$PAGE_NAME = "CoGeBlast.pl";
$TEMPDIR = "/opt/apache/CoGe/tmp/CoGeBlast";
$DATADIR = "/opt/apache/CoGe/data/";
$FASTADIR = $DATADIR.'/fasta/';
$BLASTDBDIR = $DATADIR.'/blast/db/';
$TEMPURL = "/CoGe/tmp/CoGeBlast";
$FORMATDB = "/usr/bin/formatdb";
$BLAST = "/usr/bin/blast -a 8 -K 100";
$BLASTZ = "/usr/bin/blastz";
$RESULTSLIMIT=500;
$MAX_PROC=8;

$DATE = sprintf( "%04d-%02d-%02d %02d:%02d:%02d",
		sub { ($_[5]+1900, $_[4]+1, $_[3]),$_[2],$_[1],$_[0] }->(localtime));
($USER) = CoGe::Accessory::LogUser->get_user();
$FORM = new CGI;
my %ajax = CoGe::Accessory::Web::ajax_func();

$coge = CoGeX->dbconnect();
#$coge->storage->debugobj(new DBIxProfiler());
#$coge->storage->debug(1);


my $pj = new CGI::Ajax(
		       gen_html=>\&gen_html,
		       source_search=>\&get_data_source_info_for_accn,
		       get_types=>\&get_types,
		       get_sequence=>\&get_sequence,
		       get_url=>\&get_url,
		       check_seq=>\&check_seq,
		       set_seq=>\&set_seq,
		       blast_param=>\&blast_param,
		       database_param=>\&database_param,
		       get_orgs=>\&get_orgs,
		       get_seq_types=>\&get_seq_types,
		       get_from_id=>\&get_from_id,
		       blast_search=>\&blast_search,
		       generate_feat_info=>\&generate_feat_info,
		       get_hsp_info=>\&get_hsp_info,
		       generate_overview_image=>\&generate_overview_image,
		       overlap_feats_parse=>\&overlap_feats_parse,
		       get_nearby_feats=>\&get_nearby_feats,
		       export_fasta_file=>\&export_fasta_file,
		       export_to_excel=>\&export_to_excel,
		       generate_tab_deliminated=>\&generate_tab_deliminated,
		       generate_feat_list=>\&generate_feat_list,
		       dataset_description_for_org=>\&dataset_description_for_org,
		       export_hsp_info=>\&export_hsp_info,
		       export_hsp_query_fasta=>\&export_hsp_query_fasta,
		       export_hsp_subject_fasta=>\&export_hsp_subject_fasta,
		       export_alignment_file=>\&export_alignment_file,
		       save_settings_cogeblast=>\&save_settings_cogeblast,
		       generate_basefile=>\&generate_basefile,
		       %ajax,
		      );
$pj->js_encode_function('escape');
print $pj->build_html($FORM, \&gen_html);
#print $FORM->header;
#print gen_html();


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
    my $template = HTML::Template->new(filename=>'/opt/apache/CoGe/tmpl/generic_page.tmpl');
    $template->param(TITLE=>'CoGe BLAST Analysis');
    $template->param(HELP=>'BLAST');
    my $name = $USER->user_name;
        $name = $USER->first_name if $USER->first_name;
        $name .= " ".$USER->last_name if $USER->first_name && $USER->last_name;
        $template->param(USER=>$name);

    $template->param(LOGON=>1) unless $USER->user_name eq "public";
    $template->param(DATE=>$DATE);
    $template->param(LOGO_PNG=>"CoGeBlast-logo.png");
    $template->param(BOX_NAME=>'CoGeBlast Settings');
    $template->param(ADJUST_BOX=>1);
    $template->param(BODY=>$body);
    my $prebox = HTML::Template->new(filename=>'/opt/apache/CoGe/tmpl/CoGeBlast.tmpl');
	$prebox->param(RESULTS_DIV=>1);
	$template->param(PREBOX=>$prebox->output);
    $html .= $template->output;
    }
  }

sub gen_body
  {
    my $template = HTML::Template->new(filename=>'/opt/apache/CoGe/tmpl/CoGeBlast.tmpl');
    my $form = $FORM;
    my $featid = join (",",$form->param('featid')) || 0;
    my $chr = $form->param('chr') || 0;
    my $upstream = $form->param('upstream') || 0;
    my $downstream = $form->param('downstream') || 0;
    my $dsid = $form->param('dsid') || 0;
    #my $feat_name = $form->param('featname');
    my $rc = $form->param('rc') || 0;
    my $seq = $form->param('seq');
    my $prefs = load_settings(user=>$USER, page=>$PAGE_NAME);
    $template->param(JAVASCRIPT=>1);
    $template->param(BLAST_FRONT_PAGE=>1);
    $template->param(UPSTREAM=>$upstream);
    $template->param(DOWNSTREAM=>$downstream);
    $template->param(DSID=>$dsid);
    $template->param(ORG_LIST=>get_orgs());
    $template->param(RC=>$rc);
    $template->param(FEATID=>$featid);
    $template->param(CHR=>$chr);
    if ($featid)
    {
    	$template->param(SEQVIEW=>1);
#    	my $seq = get_sequence($featid, $dsid, 0, 1, $upstream, $downstream,$rc);
#        $template->param(SEQUENCE=>$seq);
    }
    elsif ($chr)
    {
    	$template->param(SEQVIEW=>2);
#    	my $seq = get_sequence($chr, $dsid, 0, 2, $upstream, $downstream,$rc);
#    	$template->param(SEQUENCE=>$seq);
    }
    elsif($seq)
      {
	$template->param(SEQVIEW=>0);
        $template->param(SEQUENCE=>$seq);
    
      }
    else{
        $template->param(SEQVIEW=>0);
        $template->param(SEQUENCE=>'Enter a fasta sequence here');
    }
    $template->param(USER_NAME=>$USER->user_name);
    #$template->param(DEFAULT_PARAM=>$param);
    $template->param(REST=>1);
    #populate user specified default values
    my $db_list;
    if ($prefs->{orgids})
      {
	foreach my $orgid (split /,/,$prefs->{orgids})
	  {
	    my ($id, $org) = get_from_id($orgid);
	    next unless ($id);
	    $db_list .= qq{
        add_to_list('$id', '$org');
}
	  }
      }
    $template->param(document_ready=>$db_list) if $db_list;
    my $resultslimit = 200;
    $resultslimit = $prefs->{'resultslimit'} if $prefs->{'resultslimit'};
    $template->param(RESULTSLIMIT=>$resultslimit);
    $template->param(SAVE_ORG_LIST=>1) unless $USER->user_name eq "public";
    return $template->output;
  }
  
sub get_sequence
  {
    my %opts = @_;
    my $fids = $opts{fid};
    my $dsid = $opts{dsid};
    my $chr = $opts{chr};
    my $start = $opts{start};
    my $stop = $opts{stop};
    my $blast_type = $opts{blast_type};
    my $upstream = $opts{upstream};
    my $downstream = $opts{downstream};
    my $rc = $opts{rc};
    my $fasta;
    my $prot = $blast_type =~ /blast_type_p/ ? 1 : 0;
    if ($fids)
      {
	foreach my $fid (split /,/, $fids)
	  {
	    my $feat = $coge->resultset('Feature')->find($fid);
	    $fasta .= ref($feat) =~ /Feature/i ?
	      $feat->fasta(
			   prot=>$prot,
			   rc=>$rc,
			   upstream=>$upstream,
			   downstream=>$downstream,
			  )
		:
		  ">Unable to retrieve Feature object for id: $fid\n";
	  }
      }
    else
      {
	my $ds = $coge->resultset('Dataset')->find($dsid);
	$fasta = ref ($ds) =~ /dataset/i ? 
	  $ds->fasta
	    (
	     start=>$start,
	     stop=>$stop,
	     chr=>$chr,
	     prot=>$prot,
	     rc=>$rc,
	    )
	      :
		">Unable to retrieve dataset object for id: $dsid";
      }
    return $fasta
  }

sub get_url
  {
    my $url = shift;
    if ($url eq "blastn") {
      $url = "http://www.ncbi.nlm.nih.gov/blast/Blast.cgi?PAGE=Nucleotides&PROGRAM=blastn&MEGABLAST=on&BLAST_PROGRAMS=megaBlast&PAGE_TYPE=BlastSearch&SHOW_DEFAULTS=on";}
    elsif ($url eq "blastp") {
      $url = "http://www.ncbi.nlm.nih.gov/BLAST/Blast.cgi?PAGE=Proteins&PROGRAM=blastp&BLAST_PROGRAMS=blastp&PAGE_TYPE=BlastSearch&SHOW_DEFAULTS=on";}
      elsif ($url eq "blastx") {
        $url = "http://www.ncbi.nlm.nih.gov/BLAST/Blast.cgi?PAGE=Translations&PROGRAM=blastx&BLAST_PROGRAMS=blastx&PAGE_TYPE=BlastSearch&SHOW_DEFAULTS=on";}
    elsif ($url eq "tblastn") {
      $url = "http://www.ncbi.nlm.nih.gov/BLAST/Blast.cgi?PAGE=Translations&PROGRAM=tblastn&BLAST_PROGRAMS=tblastn&PAGE_TYPE=BlastSearch&SHOW_DEFAULTS=on";}
    elsif ($url eq "tblastx") {
      $url = "http://www.ncbi.nlm.nih.gov/BLAST/Blast.cgi?PAGE=Translations&PROGRAM=tblastx&BLAST_PROGRAMS=tblastx&PAGE_TYPE=BlastSearch&SHOW_DEFAULTS=on";}
    else {
      $url = 1;}
    return $url;
  }

sub blast_param
{
    my %opts = @_;
    my $seq_type = $opts{blast_type} || "blast_type_n";
    my $translate = $opts{translate};
    my $version = $opts{version};
    my $pro;
    my $template = HTML::Template->new(filename=>'/opt/apache/CoGe/tmpl/CoGeBlast.tmpl');
    if ($seq_type =~ "blast_type_n") {
      if($version && $version =~ /coge_radio/) {$template->param(BLAST_NU=>1);}
      else {$template->param(NCBI_BLAST_NU=>1);}        
    }
    else {
      $pro = 1;
      if($version && $version =~ /coge_radio/) {$template->param(BLAST_PRO=>1);}
      else {$template->param(NCBI_BLAST_PRO=>1);}
	  
	  unless ($translate)
	  {
	     if($version =~ /coge_radio/) {$template->param(BLAST_PRO_COMP=>1);}
	     else {$template->param(NCBI_BLAST_PRO_COMP=>1);}	  
	  }
    }
    my $html = $template->output;
    return $html,$version,$pro;
}

sub database_param
  {
    my $program = shift || "blastn";
    my $template = HTML::Template->new(filename=>'/opt/apache/CoGe/tmpl/CoGeBlast.tmpl');
    if ($program eq "blastn")
      {$template->param(NU_DB=>1);}
    elsif (($program eq "blastp") || ($program eq "blastx"))
      {$template->param(PRO_DB=>1);}
    else
      {$template->param(T_DB=>1);}
    my $html = $template->output;
    return $html;
  }
  
sub get_orgs
  {
    my %opts = @_;
    my $name = $opts{name};
    my $desc = $opts{desc};
    $name = "" if $name && $name =~ /^Search$/i;
    $desc = "" if $desc && $desc =~ /^Search$/i;
    my $html;
    my @db;
    if ($name)
      {
	@db = $coge->resultset('Organism')->search({name=>{like=>"%".$name."%"}});
      }
    elsif($desc)
      {
	@db = $coge->resultset("Organism")->search({description=>{like=>"%".$desc."%"}});
      }
    else
      {
	$html .= qq{<FONT CLASS ="small" id="org_count">Organism count: }.$coge->resultset('Organism')->count().qq{</FONT>\n<BR>\n};
	return $html;
      }
    ($USER) = CoGe::Accessory::LogUser->get_user();
    my $restricted_orgs = restricted_orgs(user=>$USER);
    my @opts;
    foreach my $item (sort {uc($a->name) cmp uc($b->name)} @db)
      {
	next if $restricted_orgs->{$item->name};
	push @opts, "<OPTION value=\"".$item->id."\" id=\"o".$item->id."\">".$item->name."</OPTION>";
      }
    
    $html .= qq{<FONT CLASS ="small" id="org_count">Organism count: }.scalar @opts.qq{</FONT>\n<BR>\n};
    unless (@opts) 
      {
	$html .=  qq{<input type = hidden name="org_id" id="org_id"><br>};
	$html .= "No results";
	return $html;
      }
    $html .= qq{<SELECT id="org_id" SIZE="8" MULTIPLE onClick="\$('#remove').hide(0);\$('#add').show(0);get_seq_types(['org_id'],['org_seq_types']);" ondblclick="get_from_id(['org_id', 'seq_type_id'],[add_to_list]);">\n};
    $html .= join ("\n", @opts);
    $html .= "\n</SELECT>\n";
    $html =~ s/OPTION/OPTION SELECTED/;
    return $html;
  }

sub get_seq_types
  {
    my $org_id = shift;
    return unless $org_id;
    my $org = $coge->resultset('Organism')->find($org_id);
    my %types;
    foreach my $ds ($org->datasets)
      {
	my $type = $ds->sequence_type;
	next unless $type;
	$types{$type->id}=$type;
      }
    my $html = "Sequence Types<br>";
    $html .= qq{<SELECT id="seq_type_id" size=}.scalar keys (%types).qq{ multiple onSelect="\$('#remove').hide(0);\$('#add').show(0);get_seq_types(['org_id'],['org_seq_types']);" ondblclick="get_from_id(['org_id', 'seq_type_id'],[add_to_list]);">};
    $html .= join ("\n", map{qq{<option value="}.$types{$_}->id.qq{" id="}.$types{$_}->id.qq{">}.$types{$_}->name.qq{</OPTION>}} sort keys %types);
    $html .= "</SELECT>";
    $html =~ s/OPTION/OPTION SELECTED/i;
    return $html,"\n";
  }

sub get_from_id
  {
    my $id = shift;
    my $seq_type_id = shift || 1;
    ($id, $seq_type_id) = split /_/, $id if $id =~ /_/;
    my ($obj) = $coge->resultset('Organism')->find($id);
    return unless $obj;
    my $seq_type = $coge->resultset('GenomicSequenceType')->find($seq_type_id);
    my $org = $obj->name;
    $org .= " (".$seq_type->name.")";
    return ($id."_".$seq_type_id,$org);
  }

sub generate_basefile
{
	$cogeweb = initialize_basefile(prog=>"CoGeBlast");
	return $cogeweb->basefilename;
}

sub blast_search
  {
    my %opts = @_;
    my $color_hsps_by_query_seq =  $opts{color_hsps_by_query_seq};
    my $program = $opts{program};
    my $expect = $opts{expect};
    my $job_title = $opts{job_title};
    my $wordsize = $opts{wordsize};
    my $comp = $opts{comp};
    my $matrix = $opts{matrix};
    my $gapcost = $opts{gapcost};
    my $match_score = $opts{matchscore};
    my $filter_query = $opts{filter_query};
    my $resultslimit = $opts{resultslimit} || $RESULTSLIMIT;
    my $basename = $opts{basename};
	$cogeweb = initialize_basefile(basename=>$basename, prog=>"CoGeBlast");
	
    #blastz params
    my $zwordsize = $opts{zwordsize};
    my $zgap_start = $opts{zgap_start};
    my $zgap_extension = $opts{zgap_extension};
    my $zchaining = $opts{zchaining};
    my $zthreshold = $opts{zthreshold};
    my $zmask = $opts{zmask};


    my $seq = $opts{seq};
    my $blastable = $opts{blastable};
    my $width = $opts{width};

   # exit;

    my $t1 = new Benchmark;
    my @org_ids = split(/,/,$blastable);
    my ($fasta_file, $query_seqs_info) = create_fasta_file($seq);
    my $opts;
    my $pre_command;
    if ($program eq "blastz")
      {
	$pre_command = $BLASTZ;
	$pre_command .= " $fasta_file";

	$opts .= " W=" .$zwordsize if defined $zwordsize;
	$opts .= " C=" .$zchaining if defined $zchaining;
	$opts .= " K=" .$zthreshold if defined $zthreshold;
	$opts .= " M=" .$zmask if defined $zmask;
	$opts .= " O=" .$zgap_start if defined $zgap_start;
	$opts .= " E=" .$zgap_extension if defined $zgap_extension;
	my $tmp;
	($tmp, $opts) = check_taint($opts);
      }
    else
      {
	my ($nuc_penalty,$nuc_reward,$exist,$extent);
	if ($gapcost =~/^(\d+)\s+(\d+)/) {($exist,$extent) = ($1,$2);}
	
	if ($match_score=~/^(\d+)\,(-\d+)/) {($nuc_penalty,$nuc_reward) = ($2,$1);}
	$pre_command = "$BLAST -p $program -i $fasta_file";
	if ($program =~ /^blastn$/i)
	  {
	    $pre_command .= " -q $nuc_penalty -r $nuc_reward";
	  }
	else
	  {
	    $pre_command .= " -M $matrix";
	  }
	$pre_command .=" -W $wordsize";
	$pre_command .= " -G $exist -E $extent" if $exist && $extent;
	$pre_command .= " -e $expect";
	$pre_command .= " -C $comp" if $program =~ /tblastn/i;
	$pre_command .= " -F F " unless $filter_query;
      }
    my $x;
    ($x, $pre_command) = check_taint($pre_command);
    my @results;
    my $count =1;
    my $t2 = new Benchmark;
    foreach my $orgid (@org_ids)
      {
	my ($org, $db, $fasta_file) = get_blast_db($orgid);
	next unless $db;
	my $command;
	my $outfile;
	my $report;
	if ($program eq "blastz")
	  {
	    $command = $pre_command." $fasta_file $opts";
	    $outfile = $cogeweb->basefile."-$count.blastz";
	  }
	else
	  {
	    $command = $pre_command." -d $db";
	    $outfile = $cogeweb->basefile."-$count.blast";
	  }
	push @results, {
			command=>$command,
			file=>$outfile,
			organism=>$org,
		       };
	$count++;
      }
    my $pm = new Parallel::ForkManager($MAX_PROC);
    foreach my $item (@results)
      {
	$pm->start and next;
	my $command = $item->{command};
	my $organism_name = $item->{organism};
	my $outfile = $item->{file};
	write_log("*$organism_name* running $command" ,$cogeweb->logfile);
	`$command > $outfile`;
	write_log("*$organism_name* blast analysis complete",$cogeweb->logfile);
	$pm->finish;
      }
    $pm->wait_all_children;
    foreach my $item (@results)
      {
	my $command = $item->{command};
	my $outfile = $item->{file};
	my $ta = new Benchmark;
	my $report = $outfile =~ /blastz/ ? new CoGe::Accessory::blastz_report({file=>$outfile}) : new CoGe::Accessory::blast_report({file=>$outfile});
#	my $report = $outfile =~ /blastz/ ? new CoGe::Accessory::blastz_report({file=>$outfile, limit=>$resultslimit}) : new CoGe::Accessory::blast_report({file=>$outfile, limit=>$resultslimit});
	my $tb = new Benchmark;
	my $itime = timestr(timediff($tb,$ta));
	$item->{report}= $report;
	my $file = $report->file();
	$file =~ s/$TEMPDIR//;
	$file = $TEMPURL."/".$file;
	$item->{link}=$file;
      }
    my $t3 = new Benchmark;
    write_log("Initializing sqlite database",$cogeweb->logfile);
    initialize_sqlite();
    my $t4 = new Benchmark;
    write_log("Generating Results",$cogeweb->logfile);
    my ($html,$click_all_links) = gen_results_page(results=>\@results,width=>$width,resultslimit=>$resultslimit,prog=>$program, color_hsps_by_query_seq =>  $color_hsps_by_query_seq, query_seqs_info=>$query_seqs_info);
    my $t5 = new Benchmark;
    my $init_time = timestr(timediff($t2,$t1));
    my $blast_time = timestr(timediff($t3,$t2));
    my $dbinit_time = timestr(timediff($t4,$t3));
    my $resultpage_time = timestr(timediff($t5,$t4));
    my $benchmark = qq{
Time to initialize:              $init_time
Time to blast:                   $blast_time
Time to initialize sqlite:       $dbinit_time
Time to generate results page:   $resultpage_time
};
      write_log("$benchmark" ,$cogeweb->logfile);

    write_log("Finished!", $cogeweb->logfile);
    return $html, $click_all_links;
  }
 
 
sub gen_results_page
   {
     my %opts = @_;
     my $results = $opts{results};
     my $width = $opts{width};
     my $resultslimit = $opts{resultslimit};
     my $color_hsps_by_query_seq =  $opts{color_hsps_by_query_seq};
     my $prog=$opts{prog};
     my $query_seqs_info = $opts{query_seqs_info};
     my $click_all_links;
     my $null;
     my %hsp_count;
     my $length;
     my @check;
     my @hsp;
     my $t0 = new Benchmark;
     my %query_hit_count;
     foreach my $set (@$results)
       {
	 $hsp_count{$set->{organism}}=0 unless $hsp_count{$set->{organism}};
	 if (@{$set->{report}->hsps()})
	   {
	     foreach my $hsp (sort {$a->number <=> $b->number} @{$set->{report}->hsps()})
	       {
		 my ($dsid) = $hsp->subject_name =~ /id: (\d+)/;
		 my ($chr) = $hsp->subject_name =~ /chromosome: (.*?),/;
		 my ($org) = $set->{organism};#$hsp->subject_name =~ /^\s*(.*?)\s*\(/;
		 next unless $dsid && defined $chr;
		 $hsp_count{$org}++;
		 last if ($hsp_count{$org} > $resultslimit);
		 
		 populate_sqlite($hsp,$dsid,$org);
		 my $id = $hsp->number."_".$dsid;
		 $click_all_links .= $id.",";
		 my $feat_link = qq{<span class="link" onclick="fill_nearby_feats('$id','true')">Click for Closest Feature</span>};
		 my $qname = $hsp->query_name =~ /Name: (\S*)/ ? $1 : $hsp->query_name;
		 $qname =~ s/,$//;
		 $query_hit_count{$qname}{$org}++;
		 my $coverage = $query_seqs_info->{$hsp->query_name} ? sprintf("%.1f", $hsp->length/$query_seqs_info->{$hsp->query_name}*100)."%" : "Error";
		 push @hsp, {
			     CHECKBOX=>$id."_".$chr."_".$hsp->subject_start."no",
			     ID=>$id,
			     QUERY_SEQ=>$qname,
			     HSP_ORG=>$org,
			     HSP=>qq{<span class="link" title="Click for HSP information" onclick="update_hsp_info('table_row$id');\$('#middle_column_button_hide').show(0);\$('#middle_column_button_show').hide(0);">}.$hsp->number."</span>",
			     HSP_EVAL=>$hsp->pval,
			     HSP_LENGTH=>$hsp->length,
			     COVERAGE=>$coverage,
			     HSP_PID=>$hsp->percent_id."%",
			     HSP_SCORE=>$hsp->score,
			     HSP_POS=>($hsp->subject_start),
			     HSP_CHR=>$chr,
			     HSP_LINK=>$feat_link};
	       }
	   }
	 $hsp_count{$set->{organism}} = scalar @{$set->{report}->hsps()}
       }
     my $t1 = new Benchmark;
     my ($chromosome_data, $chromosome_data_large) = generate_chromosome_images(results=>$results,large_width=>$width,, resultslimit=>$resultslimit, color_hsps_by_query_seq =>  $color_hsps_by_query_seq);
     my $t2 = new Benchmark;
     unless (@hsp) 
       {
	 $null = "null";
       }
     my $template = HTML::Template->new(filename=>'/opt/apache/CoGe/tmpl/CoGeBlast.tmpl');
	 $template->param(RESULT_TABLE=>1);
     # ERIC, i added this so it wouldnt fail
     $template->param(NULLIFY=>$null) if $null;
     my $hsp_limit_flag =0;
     my $hsp_count;
     $hsp_count .= qq{<table class="small resultborder"><tr><th>Query Seq<th>}.join "<th>", sort keys %hsp_count;
     my $class = "even";
     foreach my $query (sort keys %query_hit_count)
       {
	 $hsp_count.= qq{<tr class="$class"><td>$query. <span class = species>Length: }.$query_seqs_info->{$query}."</span>";
	 foreach my $org (sort keys %hsp_count)
	   {
	     my $count = $query_hit_count{$query}{$org} ? $query_hit_count{$query}{$org} : 0;
	     $hsp_count .= qq{<td align=center>$count};
	   }
	 $class = $class eq "even" ? "odd" : "even";
       }
     $hsp_count.= qq{<tr class="$class"><th>Total};
     foreach my $org (sort keys %hsp_count)
       {
	 my $count = $hsp_count{$org};
	 $count = "<span class=alert>$count</span>" if $count > $resultslimit;
	 $hsp_count .= qq{<td align=center class=species>$count</td>};
	 
       }
     $hsp_count .= "</table>";
     foreach my $org (keys %hsp_count)
       {
	 if ($hsp_count{$org} > $resultslimit)
	   {
	     $hsp_count .= "<span class=\"small alert\">Only top $resultslimit HSPs shown for $org.</span><br>";
	     $hsp_limit_flag = 1;
	   }
       }
     $hsp_count .= "<span class=\"small alert\">All results are in the blast report.</span>" if $hsp_limit_flag;

     $template->param(HSP_COUNT=>$hsp_count);

     if (@hsp)
     {
       @hsp = sort {$a->{HSP_ORG} cmp $b->{HSP_ORG} || $a->{HSP} cmp $b->{HSP}} @hsp if @hsp;
       $template->param(HSP_TABLE=>1);
       $template->param(HSPS=>\@hsp);
     }
     
     my $hsp_results = $template->output;
     $template->param(HSP_TABLE=>0);
     $template->param(RESULT_TABLE=>0);
     $template->param(HSP_RESULTS=>$hsp_results);
     $template->param(CHROMOSOMES_IF=>1);
     $template->param(CHROMOSOME_LOOP=>$chromosome_data);
     my $chromosome_element = $template->output;
     $template->param(CHROMOSOMES_IF=>0);
     if (scalar(@$chromosome_data_large) > 0)
     {
       $template->param(CHROMOSOMES_LARGE_IF=>1);
       $template->param(CHROMOSOME_LOOP_LARGE=>$chromosome_data_large);
       my $chr_large_element = $template->output;
       $template->param(CHROMOSOMES_LARGE_IF=>0);
       $template->param(CHR_LARGE=>$chr_large_element);
     }
     $template->param(CHROMOSOMES=>$chromosome_element);
     $template->param(BLAST_RESULTS=>1);
     $template->param(DATA_FILES=>gen_data_file_summary(prog=>$prog, results=>$results));
     my $html = $template->output;
     my $box_template = HTML::Template->new(filename=>'/opt/apache/CoGe/tmpl/box.tmpl');
     $box_template->param(BOX_NAME=>"CoGeBlast Results");
     $box_template->param(BODY=>$html);
     my $outhtml = $box_template->output;
     my $t3 = new Benchmark;
    my $table_time = timestr(timediff($t1,$t0));
    my $figure_time = timestr(timediff($t2,$t1));
    my $render_time = timestr(timediff($t3,$t2));
    my $benchmark = qq{
Time to gen tables:              $table_time
Time to gen images:              $figure_time
Time to gen results:             $render_time
};
     write_log($benchmark, $cogeweb->logfile);
     return $outhtml, $click_all_links;
   }

sub gen_data_file_summary
  {
    my %opts = @_;
    my $prog = $opts{prog};
    my $results = $opts{results};
    my $html = "<table><tr>";
    $html .= qq{<td class = small valign="top">Data Download};
    $html .= "<div class=xsmall><A HREF=\"#\" onClick=\"get_all_hsp_data();\">HSP Data</A></DIV>\n";
    $html .= "<div class=xsmall><A HREF=\"#\" onClick=\"get_query_fasta();\">Query HSP FASTA File</A></DIV>\n";
    if ($prog eq "tblastn")
      {
	$html .= "<div class=xsmall><A HREF=\"#\" onClick=\"get_subject_fasta();\">Subject HSP Protein FASTA File</A></DIV>\n";
	$html .= "<div class=xsmall><A HREF=\"#\" onClick=\"get_subject_fasta('dna');\">Subject HSP DNA FASTA File</A></DIV>\n";
      }
    else
      {
	$html .= "<div class=xsmall><A HREF=\"#\" onClick=\"get_subject_fasta();\">Subject HSP FASTA File</A></DIV>\n";
      }
    $html .= "<div class=xsmall><A HREF=\"#\" onClick=\"get_alignment_file();\">Alignment File</A></DIV>\n";
    $html .= qq{<td class = small valign="top">Analysis Files};
    #my $dbname = $TEMPURL."/".basename($cogeweb->sqlitefile);
    my $dbname = $TEMPURL."/".$cogeweb->basefilename.".sqlite";
    $html .= "<div class=xsmall><A HREF=\"$dbname\" target=_new>SQLite DB file</A></DIV>\n";
    foreach my $item (@$results)
      {
	my $blast_file = $item->{link};
	my $org = $item->{organism};
	$html .= qq{<div class=xsmall><a href = "$blast_file" target=_new>Blast file for $org</div>\n};
      }
#    my $blast = $TEMPURL."/".$cogeweb->basefilename.".sqlite";
#    $html .= "<div class=xsmall><A HREF=\"$dbname\" target=_new>SQLite DB file</A></DIV>\n";
    $html .= qq{<td class = small valign="top">Log File};
    #my $logfile = $TEMPURL."/".basename($cogeweb->logfile);
    my $logfile = $TEMPURL."/".$cogeweb->basefilename.".log";
    $html .= "<div class=xsmall><A HREF=\"$logfile\" target=_new>Log</A></DIV>\n";
    $html .= qq{</table>};
  }

sub generate_chromosome_images
  {
     my %opts = @_;
     my $results = $opts{results};
     my $hsp_type = $opts{hsp_type} || "NA";
     my $width = $opts{width} || 400;
     my $large_width = $opts{large_width} || 3*$width;
     my $resultslimit = $opts{resultslimit};
     my $imagefile_name = $opts{filename} || "null";
     my $color_hsps_by_query_seq =  $opts{color_hsps_by_query_seq};
     $color_hsps_by_query_seq = 0 if $color_hsps_by_query_seq =~ /false/;
     my $height = ($width / 16);
     my $large_height = ($large_width / 16) <= 64 ? ($large_width / 16) : 64;
     my $scale = $opts{scale} || 'linear';
     my %data;
     my $filename;
     my (@data, @large_data,@no_data);
     my %hsp_count;
     my %query_seqs;
     if ($color_hsps_by_query_seq)
       {
	 foreach my $set (@$results)
	   {
	     map {$query_seqs{$_->query_name}=1} @{$set->{report}->hsps()};
	   }
	 my $colors = color_pallet(num_seqs=>scalar keys %query_seqs);
	 my $count = 0;
	 foreach my $key (sort keys %query_seqs)
	   {
	     $query_seqs{$key} = $colors->[$count];
	     $count++;
	   }
       }

     foreach my $set (@$results)
    {
	my $org = $set->{organism};
	$data{$org}{file}=$set->{link};
	$filename = $imagefile_name."_*.png";
	if ($imagefile_name ne "null")
	{
	  if (`ls $filename 2>>/dev/null`)
	  {
	    $data{$org}{image} = $imagefile_name;
	    $data{$org}{skip} = 1;
	    next;
	  }
	}
	$filename = $set->{link};
	if (@{$set->{report}->hsps()})
	  {
	    foreach my $hsp (@{$set->{report}->hsps()})
	      {
		#first, initialize graphic
		$org =~ s/\s+$//;
		$hsp_count{$org}++;
		next if $hsp_count{$org}> $resultslimit;
		my ($chr) = $hsp->subject_name =~ /chromosome: (.*?),/;
		$data{$org}{image} =new CoGe::Graphics::GenomeView({color_band_flag=>1, image_width=>$width, chromosome_height=>$height}) unless $data{$org}{image};
		$data{$org}{large_image} =new CoGe::Graphics::GenomeView({color_band_flag=>1, image_width=>$large_width, chromosome_height=>$large_height}) unless $data{$org}{large_image};
		my ($dsid) = $hsp->subject_name =~ /id: (\d+)/;
		#add chromosome to graphic
		unless ($data{$org}{chr}{$chr})
		  {
		    next unless $dsid;
		    my $ds = $coge->resultset('Dataset')->find($dsid);
		    my $last_pos = $ds->last_chromosome_position($chr);
		    $data{$org}{image}->add_chromosome(name=>"Chr: $chr",
						       end=>$last_pos,
						      );
		    $data{$org}{chr}{$chr}=1;
		  }
		my $num = $hsp->number."_".$dsid;
		my $up = $hsp->strand eq "++" ? 1 : 0;
		my $color = $color_hsps_by_query_seq ? $query_seqs{$hsp->query_name} : [0,200,0];
  		$data{$org}{image}->add_feature(name=>$hsp->number,
  						start=>$hsp->sstart,
  						stop=>$hsp->sstop,
  						chr=>"Chr: $chr",
  						imagemap=>qq/class="imagemaplink" title="HSP No. /.$hsp->number.qq/" onclick="hide_big_picture();show_hsp_div();loading('image_info','Information');loading('query_image','Image');loading('subject_image','Image');get_hsp_info(['args__blastfile','args__/.$cogeweb->basefile.qq/','args__num','args__$num'],['image_info','query_image','subject_image']);"/,
  						up=>$up,
  						color=>$color,
  					       );
	      }
	  }
      }
    my $count = 1;
    foreach my $org (sort keys %data)
    {
	if ($data{$org}{image})
	  {
	    my $image_file;
	    my $image_map;
	    my $large_image_file;
	    my $image_map_large;
	    unless ($data{$org}{skip})	    
	      {
		my $x;
		$large_image_file = $cogeweb->basefile."_".$hsp_type."_$count"."_large.png";
		($x, $large_image_file) = check_taint($large_image_file);
		$image_file = $cogeweb->basefile."_".$hsp_type."_$count.png";
		($x, $image_file) = check_taint($image_file);
		$data{$org}{image}->generate_png(filename=>$image_file);
		$image_map = $data{$org}{image}->generate_imagemap(mapname=>$cogeweb->basefilename."_".$count);
		my $map_file = $cogeweb->basefile."_$count.$hsp_type.map";
		($x, $map_file) = check_taint($map_file);
		open (MAP, ">$map_file");
		print MAP $image_map;
		close MAP;
		$data{$org}{image}->image_width($large_width);
		$data{$org}{image}->chromosome_height($large_height);
		$data{$org}{image}->show_count(1);
		$data{$org}{image}->generate_png(filename=>$large_image_file);
		$image_map_large = $data{$org}{image}->generate_imagemap(mapname=>$cogeweb->basefilename."_".$count."_large");
		$map_file = $cogeweb->basefile."_$count.$hsp_type.large.map";
		($x, $map_file) = check_taint($map_file);
		open (MAP, ">$map_file");
		print MAP $image_map_large;
		close MAP;

	    }
	    else
	      {
		my $x;
		$image_file = $data{$org}{image}."_$count.png";
		$image_map = get_map($cogeweb->basefile."_$count.$hsp_type.map");
		$large_image_file = $data{$org}{image}."_$count"."_large.png";
		$image_map_large = get_map($cogeweb->basefile."_$count.$hsp_type.large.map");
		($x, $image_file) = check_taint($image_file);
		($x, $large_image_file) = check_taint($large_image_file);
	      }
	    
	    $image_file =~ s/$TEMPDIR/$TEMPURL/;
	    $large_image_file =~ s/$TEMPDIR/$TEMPURL/;
	    
	    push @large_data,  {DB_NAME_LARGE=>"<a href=".$data{$org}{file}. " target=_new>$org <span class='small link'>Blast Report</span></a><br>", CHR_IMAGE_LARGE=>"<img src=$large_image_file ismap usemap='#".$cogeweb->basefilename."_"."$count"."_large' border=0>$image_map_large",IMAGE_ID_LARGE=>$count,};
	    push @data,  {DB_NAME=>"<a href=".$data{$org}{file}. " target=_new>$org <span class='small link'>Blast Report</span></a><br>", CHR_IMAGE=>"<img src=$image_file ismap usemap='#".$cogeweb->basefilename."_"."$count' border=0>$image_map",HIT=>1,IMAGE_ID=>$count,};
	    $count++;
	  }
	else
	  {
	    push @no_data,  {DB_NAME=>"No Hits: <a href=".$data{$org}{file}. " target=_new>$org</a>"};
	  }
      }


    return [@data,@no_data], \@large_data;
  }

sub get_map
  {
    my $file = shift;
    my $map;
    open (IN, $file) || die "$!";
    while (<IN>)
      {
	$map .= $_;
      }
    close IN;
    return $map;
  }

sub create_fasta_file
  {
    my $seq = shift;
    my %seqs; #names and lengths
    if ($seq =~ />/)
      {
	foreach (split />/, $seq)
	  {
	    next unless $_;
	    my ($name, $tmp) = split/\n/,$_,2;
	    $tmp =~ s/\n//g;
	    $seqs{$name}=length($tmp);
	  }
      }
    write_log("creating user's fasta file",$cogeweb->logfile);
    open(NEW,"> ".$cogeweb->basefile.".fasta");
    print NEW $seq;
    close NEW;
    return $cogeweb->basefile.".fasta", \%seqs;
  }
    

sub get_blast_db
  {
    my $id = shift;
    my ($orgid, $seqtypeid) = $id =~ /_/ ? split/_/, $id : ($id,1);
    my $org = $coge->resultset('Organism')->find($orgid);
    my $seqtype = $coge->resultset('GenomicSequenceType')->find($seqtypeid);
    my @ds = $org->current_datasets(type=>$seqtypeid);
    @ds = sort {$a->id <=> $b->id }@ds;
    return unless @ds;
    my $org_name = $ds[0]->organism->name;
    $org_name .= " (".$seqtype->name.")" if $seqtype;
    my $title = $org_name ." (";
    my %vers = map {$_->version,1} @ds;
    if (keys %vers > 1)
      {
	my @chrs;
	foreach my $ds (@ds)
	  {
	    push @chrs, join (", ", map {"chr:".$_." v:".$ds->version." ds:".$ds->id} $ds->get_chromosomes);
	  }
	$title .= join (", ", @chrs);
      }
    else
      {
	$title .= "v".join ("",keys %vers)." ds:".$ds[0]->id();
      }
    $title .= ")";
    $title =~ s/(`|')//g;
    my $md5 = md5_hex($title);
    my $file = $FASTADIR."/$md5.fasta";
    my $res;
    while (-e "$file.running")
      {
	sleep 60;
      }
    if (-r $file)
      {
	write_log("*$org_name* fasta file ($md5) exists", $cogeweb->logfile);
	$res = 1;
      }
    else
      {
	system "touch $file.running"; #track that a blast anlaysis is running for this
	$res = generate_fasta(dslist=>\@ds, file=>$file) unless -r $file;
	system "rm $file.running" if -r "$file.running"; #remove track file
      }    
    my $blastdb = "$BLASTDBDIR/$md5";
     while (-e "$blastdb.running")
      {
	sleep 60;
      }
    if (-r $blastdb.".nsq")
      {
	write_log("*$org_name* blastdb file ($md5) exists", $cogeweb->logfile);
	$res = 1;
      }
    else
      {
	system "touch $blastdb.running"; #track that a blast anlaysis is running for this
	$res = generate_blast_db(fasta=>$file, blastdb=>$blastdb, org=>$org_name);
	system "rm $blastdb.running" if -r "$blastdb.running"; #remove track file
      }
    return $org_name, $blastdb, $file if $res;
    return 0;
    
  }


sub generate_fasta
  {
    my %opts = @_;
    my $dslist = $opts{dslist};
    my $file = $opts{file};
    $file = $FASTADIR."/$file" unless $file =~ /$FASTADIR/;
    write_log("creating fasta file.", $cogeweb->logfile);
    open (OUT, ">$file") || die "Can't open $file for writing: $!";;
    foreach my $ds (@$dslist)
      {
	foreach my $chr (sort $ds->get_chromosomes)
	  {
	    my $title =  $ds->organism->name." (v". $ds->version.") "."chromosome: $chr".", CoGe database id: ".$ds->id;
	    $title =~ s/^>+/>/;
	    write_log("adding sequence $title", $cogeweb->logfile);
	    print OUT ">".$title."\n";
	    print OUT $ds->get_genomic_sequence(chr=>$chr),"\n";
	  }
      }
    close OUT;
    write_log("Completed fasta creation", $cogeweb->logfile);
    return 1 if -r $file;
    write_log("Error with fasta file creation", $cogeweb->logfile);
    return 0;
  }

sub generate_blast_db
  {
    my %opts = @_;
    my $fasta = $opts{fasta};
    my $blastdb = $opts{blastdb};
#    my $title= $opts{title};
    my $org= $opts{org};
    my $command = $FORMATDB." -p F";
    $command .= " -i '$fasta'";
    $command .= " -t '$org'";
    $command .= " -n '$blastdb'";
    write_log("creating blastdb for $org ($blastdb)",$cogeweb->logfile);
    `$command`;
    return 1 if -r "$blastdb.nsq";
    write_log("error creating blastdb for $org ($blastdb)",$cogeweb->logfile);
    return 0;
  }

sub generate_feat_info 
  {
    my $featid = shift;
    my $checkbox = shift;
    $featid =~ s/^table_row//;
    #$featid =~ s/^no_feat//;
    $featid =~ s/_\d+_\d+$//;
    my ($feat) = $coge->resultset("Feature")->find($featid);
    unless (ref($feat) =~ /Feature/i)
    {
      return "Unable to retrieve Feature object for id: $featid";
    }
    my $html = qq{<a href="#" onClick="\$('#overlap_box').slideToggle(pageObj.speed);" style="float: right;"><img src='/CoGe/picts/delete.png' width='16' height='16' border='0'></a>} unless $checkbox;
    $html .= $feat->annotation_pretty_print_html();
    return $html;
  }

sub get_hsp_info
  {
    my %opts = @_;
    my $hsp_id = $opts{num};
    my $filename = $opts{blastfile};
    $filename =~ s/$TEMPDIR//;
    $cogeweb = initialize_basefile(basename=>$filename, prog=>"CoGeBlast");
    my $dbfile = $cogeweb->sqlitefile;
    my $dbh = DBI->connect("dbi:SQLite:dbname=$dbfile","","");
    
    $hsp_id =~ s/^table_row// if $hsp_id =~ /table_row/;
    $hsp_id =~ s/^\d+_// if $hsp_id =~ tr/_/_/ > 1;
    
    my ($hsp_num, $pval, $pid,$psim, $score, $qgap, $sgap, $match,$qmismatch, $smismatch, $strand, $length,$qstart, $qstop,$sstart, $sstop,$qalign,$salign,$align,$qname,$sname);
    
    my $sth = $dbh->prepare(qq{SELECT * FROM hsp_data WHERE name = ?});
    
    ($hsp_num) = $hsp_id =~ /^(\d+)_\d+$/;
    $sth->execute($hsp_id) || die "unable to execute";
    while (my $info = $sth->fetchrow_hashref())
      
      {
	$pval = $info->{eval};
	$pid = $info->{pid};
	$psim = $info->{psim};
	$score = $info->{score};
	$qgap = $info->{qgap};
	$sgap = $info->{sgap};
	$match = $info->{match};
	$qmismatch = $info->{qmismatch};
	$smismatch = $info->{smismatch};
	$strand = $info->{strand};
	$length = $info->{length};
	$qstart = $info->{qstart};
	$qstop =  $info->{qstop};
	$sstart = $info->{sstart};
	$sstop =  $info->{sstop};
	$qalign = $info->{qalign};
	$salign = $info->{salign};
	$align = $info->{align};
	$qname = $info->{qname};
	$sname = $info->{sname};
      }
    #$sth->execute($name, $pval, $pid,$psim, $score, $qgap, $sgap,$match,$qmismatch, $smismatch, $strand, $length,$qposition,$sposition,$qalign,$salign,$align);
    my $qlength = $qstop - $qstart;
    
    my ($sub_chr) = $sname =~ /chromosome: (.*?),/;
    my $query_name = "<pre>".$qname."</pre>";
    $query_name = wrap('','',$query_name);
    $query_name =~ s/\n/<br>/g;
    
    my $subject_name = "<pre>".$sname."</pre>";
    $subject_name = wrap('','',$subject_name);
    $subject_name =~ s/\n/<br>/g;
#     
    my @table1 = ({HSP_PID_QUERY=>$pid,
		   HSP_PSIM_QUERY=>$psim,
		   HSP_GAP_QUERY=>$qgap,
		   HSP_PID_SUB=>$pid,
		   HSP_PSIM_SUB=>$psim,
		   HSP_GAP_SUB=>$sgap,
		   HSP_MATCH_QUERY=>$match,
		   HSP_MISMATCH_QUERY=>$qmismatch,
		   HSP_MATCH_SUB=>$match,
		   HSP_MISMATCH_SUB=>$smismatch,
		   HSP_POSITION_QUERY=>$qstart."-".$qstop,
		   HSP_POSITION_SUB=>$sstart."-".$sstop,
		  });
    my @table2 = ({HSP_STRAND=>$strand,
		   HSP_EVAL=>$pval,
		   HSP_SCORE=>$score,
		   HSP_LENGTH=>$length,
		   HSP_CHR=>$sub_chr,
    		  });
     my $template = HTML::Template->new(filename=>'/opt/apache/CoGe/tmpl/CoGeBlast.tmpl');
     $template->param(HSP_IF=>1);
     $template->param(HSP_NUM=>$hsp_num);
     $template->param(HSP_QS=>\@table1);
     $template->param(HSP_HSP=>\@table2);

     my ($sub_dsid) = $sname =~ /id: (\d+)/;
     my $align_str = "";
     my $query_seq = $qalign;
     $query_seq = wrap('','',$query_seq);
     my @query = split(/\n/,$query_seq);
     $query_seq =~ s/[^atgc]//ig;
     $query_seq = wrap('','',$query_seq);
     $query_seq =~ s/\n/<br>/g;
     $query_seq =~ tr/atgc/ATGC/;
     $query_seq = qq{<pre>$query_seq</pre>};
     
     my $sub_seq = $salign;
     $sub_seq = wrap('','',$sub_seq);
     my @sub = split(/\n/, $sub_seq);
     $sub_seq =~ s/[^atgc]//ig;
     $sub_seq = wrap('','',$sub_seq);
     $sub_seq =~ s/\n/<br>/g;
     $sub_seq =~ tr/atgc/ATGC/;
     $sub_seq = qq{<pre>$sub_seq</pre>};
     

    
     my $alignment = $align;
    $alignment =~ s/ /\./g;
     $alignment = wrap('','',$alignment);
     my @align = split(/\n/,$alignment);
     
     for(my $i=0;$i<scalar(@sub);$i++)
     {
       $align_str .= $query[$i]."<br>".$align[$i]."<br>".$sub[$i]."<br>";
     }
    $align_str =~ s/<br>$//;
    $align_str =~ s/\./ /g;
     $align_str = "<pre>$align_str</pre>";
     
     $template->param(QUERY_SEQ=>qq{<a href="#" class="small" onclick="show_seq('$query_seq','$query_name',1,'seqObj','seqObj','}.$qstart."','".$qstop.qq{')">Query Sequence Hit</a>});
     my $rc = $strand =~ /-/ ? 1 : 0;
#     $template->param(SUB_SEQ=>qq{<a href="#" onclick="show_seq('$sub_seq','$subject_name',2,'$sub_dsid','$sub_chr','}.$sstart."','".$sstop.qq{','$rc')">Click for Subject Sequence</a>});
     $template->param(SUB_SEQ=>qq{<a class="small" href="SeqView.pl?dsid=$sub_dsid;chr=$sub_chr;start=$sstart;stop=$sstop;rc=$rc" target=_new>Subject Sequence Hit</a>});
     
     $template->param(ALIGNMENT=>qq{<a href="#" class="small" onclick="show_seq('$align_str','HSP No. $hsp_num',0,0,0,0)">Alignment</a>});
     
     my $html = $template->output;
     $template->param(HSP_IF=>0);
     
    #get query sequence total length
    my $query = $dbh->selectall_arrayref(qq{SELECT * FROM sequence_info WHERE type = "query" AND name = "$qname"});
    my $subject = $dbh->selectall_arrayref(qq{SELECT * FROM sequence_info WHERE type = "subject" AND name = "$sname"});
    my ($query_image, $subject_image) = generate_hit_image(hsp_num=>$hsp_num, hspdb=>$dbh, hsp_name=> $hsp_id);
    my $query_link = qq{
<div class=small>Query: $qname</div>
<img src=$query_image border=0>
};
    $query_link =~ s/$TEMPDIR/$TEMPURL/;

    my ($dsid) = $sname =~ /id: (\d+)/;
    my ($chr) = $sname =~ /chromosome: (.*?),/;
    
    my $subject_link = qq{
<div class=small>Subject: $sname</div>
<a href = 'GeLo.pl?chr=$chr&ds=$dsid&x=$sstart&z=5' target=_new border=0><img src=$subject_image border=0></a>
};
    $subject_link =~ s/$TEMPDIR/$TEMPURL/;
     return $html, $query_link, $subject_link;
   }
  
sub generate_overview_image
  {
     my %opts = @_;
     my $basename = $opts{basename};
     my $type = $opts{type};
     my $image_width = $opts{width};
     my @set = split/\n/, `ls $TEMPDIR/$basename*.blast`;
     my @reports;
     my $count = 1;
     $cogeweb = initialize_basefile(basename=>$basename,prog=>"CoGeBlast");
     foreach my $blast (@set){
       my $report = new CoGe::Accessory::blast_report({file=>$blast});
       my ($org_name) = $report->hsps->[$count-1]->subject_name =~ /^\s*(.*?)\s*\(/;
       push @reports,{report=>$report,organism=>$org_name,link=>$TEMPURL."/".$basename."-".$count.".blast",};
       $count++;
     }
     my $image_filename = $cogeweb->basefile."_".$type;
     my ($chromosome_data, $chromosome_data_large) = generate_chromosome_images(results=>\@reports,hsp_type=>$type,large_width=>$image_width,filename=>$image_filename);
     my $template = HTML::Template->new(filename=>'/opt/apache/CoGe/tmpl/CoGeBlast.tmpl');
     $template->param(CHROMOSOMES_IF=>1);
     $template->param(CHROMOSOME_LOOP=>$chromosome_data);
     my $chromosome_element = $template->output;
     $template->param(CHROMOSOMES_IF=>0);
     my $chr_large_element;
     if (scalar(@$chromosome_data_large) > 0)
     {
       $template->param(CHROMOSOMES_LARGE_IF=>1);
       $template->param(CHROMOSOME_LOOP_LARGE=>$chromosome_data_large);
       $chr_large_element = $template->output;
       $template->param(CHROMOSOMES_LARGE_IF=>0);
     }
     my $html = $chromosome_element.$chr_large_element;
     return $html;
  }
     

sub generate_hit_image
  {
    my %opts = @_;
    my $hsp_name = $opts{hsp_name};
    my $width = $opts{width} || 400;
    my $dbh = $opts{dbh};
    my $dbfile = $cogeweb->sqlitefile;
    $dbh = DBI->connect("dbi:SQLite:dbname=$dbfile","","") unless $dbh;
    my $sth = $dbh->prepare(qq{SELECT * FROM hsp_data WHERE name = ?});
    $sth->execute($hsp_name) || die "unable to execute";
    my $hsp = $sth->fetchrow_hashref();
    my $query = $dbh->selectall_arrayref(qq{SELECT * FROM sequence_info WHERE type = "query" AND name = "}.$hsp->{qname}.qq{"});
    my $subject = $dbh->selectall_arrayref(qq{SELECT * FROM sequence_info WHERE type = "subject" AND name = "}.$hsp->{sname}.qq{"});
    my ($hsp_num) = $hsp_name =~ /^(\d+)_\d+$/;
    #generate_query_image
    my $cq = new CoGe::Graphics::Chromosome ();
    $cq->iw($width);
    $cq->draw_chromosome(1);
    $cq->draw_ruler(1);
    $cq->draw_chr_end(0);
    $cq->minor_tick_labels(0);
    $cq->major_tick_labels(1);
    $cq->draw_hi_qual(0);
    $cq->padding(2);
    $cq->set_region(start=>1, stop=>$query->[0][3]);
    $cq->feature_height(10);
    $cq->feature_labels(1);
    my $strand = $hsp->{strand} =~ /-/ ? "-1" : 1;
    my ($qstart, $qstop) = ($hsp->{qstart}, $hsp->{qstop});
    ($qstart,$qstop) = ($qstop,$qstart) if $qstart > $qstop;
    my $feat = CoGe::Graphics::Feature::HSP->new({start=>$qstart, stop=>$qstop, strand=>$strand, label=>$hsp_num, type=>"HSP"});
    $feat->color([255,200,0]);
    $cq->add_feature($feat);
    my ($dsid) = $hsp->{sname} =~ /id: (\d+)/;
    my ($ds) = $coge->resultset('Dataset')->resolve($dsid);
    my ($chr) = $hsp->{sname} =~ /chromosome: (.*?),/;
    my $len = $hsp->{sstop} - $hsp->{sstart}+1;
    my $start = $hsp->{sstart}-5000;
    $start = 1 if $start < 1;
    my $stop = $hsp->{sstop}+5000;    
    my $cs = new CoGe::Graphics::Chromosome ();
#    $cs->chr_length($subject->[0][3]);
    $cs->iw($width);
    $cs->draw_chromosome(1);
    $cs->draw_ruler(1);
    $cs->draw_chr_end(0);
    $cs->minor_tick_labels(0);
    $cs->major_tick_labels(1);
    $cs->draw_hi_qual(0);
    $cs->padding(2);
    $cs->set_region(start=>$start, stop=>$stop);
    $cs->feature_height(10);
    $cs->overlap_adjustment(0);
    $cs->feature_labels(1);
    my ($sstart, $sstop) = ($hsp->{sstart}, $hsp->{sstop});
    ($sstart,$sstop) = ($sstop,$sstart) if $sstart > $sstop;
    $feat = CoGe::Graphics::Feature::HSP->new({start=>$sstart, stop=>$sstop, strand=>$strand, order=>2, label=>$hsp_num, type=>"HSP"});
    $feat->color([255,200,0]);
    $cs->add_feature($feat);
    
    my $graphics = new CoGe::Graphics;
    $graphics->process_features(c=>$cs, layers=>{features=>{gene=>1, cds=>1, mrna=>1, rna=>1, cns=>1}}, ds=>$ds, chr=>$chr, coge=>$coge);
    $cs->overlap_adjustment(1);
    $cq->overlap_adjustment(1);

    #find neighboring hsps
    my $sname = $hsp->{sname};
    my $qname = $hsp->{qname};
    my $statement = qq{
select * from hsp_data where
((sstart <= $stop AND
sstart >= $start) OR
(sstop <= $stop AND
sstop >= $start)) AND
sname = "$sname" AND
qname = "$qname"
};
    $sth = $dbh->prepare($statement);
    $sth->execute;
    while (my $data = $sth->fetchrow_hashref)
      {
	next if $data->{name} eq $hsp_name;
	my ($label) = $data->{name} =~ /^(\d+)_\d+$/;
	$feat = CoGe::Graphics::Feature::HSP->new({start=>$data->{sstart}, stop=>$data->{sstop}, strand=>$data->{strand}, order=>2, label=>$label, type=>"HSP"});
	$feat->color([255,0,0]);
	$cs->add_feature($feat);
	$feat = CoGe::Graphics::Feature::HSP->new({start=>$data->{qstart}, stop=>$data->{qstop}, strand=>$data->{strand}, order=>1, label=>$label, type=>"HSP"});
	$feat->color([255,0,0]);
	$cq->add_feature($feat);
      }

    my $query_file = $cogeweb->basefile.".q.".$hsp_name.".png";
    $cq->generate_png(file=>$query_file);
    my $sub_file = $cogeweb->basefile.".s.".$hsp_name.".png";
    $cs->generate_png(file=>$sub_file);
    return $query_file, $sub_file;
  }

sub overlap_feats_parse #Send to GEvo
  {
    my $accn_list = shift;
    $accn_list =~ s/^,//;
    $accn_list =~ s/,$//;
    my @list;
     my @no_feats;
    my $url = "/CoGe/GEvo.pl?";
    my ($chr,$dsid,$loc);
    my $count = 1;
    foreach my $featid (split /,/,$accn_list)
    {
		if($featid=~/no/)
		{
			($dsid,$chr,$loc) = $featid =~/^\d+_(\d+)_(.+?)_(\d+)no$/;
			push @no_feats,{dsid=>$dsid,chr=>$chr,loc=>$loc};
		}
		else{
			$featid =~ s/_.*$//;
			push @list, $featid;			
		}
    }
    my %seen = ();
    @list = grep {!$seen{$_}++} @list;
    foreach my $featid( @list)
      {
      	 #my ($feat) = $coge->resultset("Feature")->find($featid);
		#my ($feat_name) = sort $feat->names;#something
		$url .= "fid$count=$featid&";
		$count ++;
      }
	foreach my $no_feat (@no_feats)
	{
		$url .= "dsid$count=$no_feat->{dsid}&chr$count=$no_feat->{chr}&x$count=$no_feat->{loc}&";
		$count++
	}
    $count--;
    return ("alert",$count) if $count > 20;
    $url .= "num_seqs=$count";
    return $url;
  }
	
sub initialize_sqlite
  {
    my $dbfile = $cogeweb->sqlitefile;
#    $dbfile = $TEMPDIR."/".$dbfile unless $dbfile =~ /$TEMPDIR/;
    return if -r $dbfile;
    my $dbh = DBI->connect("dbi:SQLite:dbname=$dbfile","","")  || die "cant connect to db";
    my $create = qq{
CREATE TABLE hsp_data
(
id INTEGER PRIMARY KEY AUTOINCREMENT,
name varchar(50),
eval varchar(10),
pid varchar(10),
psim varchar(10),
score varchar(10),
qgap integer(10),
sgap integer(10),
match integer(10),
qmismatch integer(10),
smismatch integer(10),
strand varchar(2),
length integer(10),
qstart integer(50),
qstop integer(50),
sstart integer(50),
sstop integer(50),
qalign text,
salign text,
align text,
qname text,
sname text,
hsp_num integer,
org text
)
};
    $dbh->do($create);
     my $index = qq{
  ALTER TABLE 'hsp_data' ADD AUTO_INCREMENT 'id';
 };
#     $dbh->do($index);
     $index = qq{
 CREATE INDEX name ON hsp_data (name)
 };
     $dbh->do($index);
     $index = qq{
 CREATE INDEX qname ON hsp_data (qname)
 };
     $dbh->do($index);
     $index = qq{
 CREATE INDEX sname ON hsp_data (sname)
 };
     $dbh->do($index);
    $create = qq{
 CREATE TABLE sequence_info
 (
 id INTEGER PRIMARY KEY AUTOINCREMENT,
 name varchar(255),
 type varchar(10),
 length integer(1024)
 )
 };
     $dbh->do($create);
      $index = qq{
  CREATE INDEX seqname ON sequence_info (name)
  };
      $dbh->do($index);
      $index = qq{
  CREATE INDEX type ON sequence_info (type)
  };
      $dbh->do($index);
    system "chmod +rw $dbfile";
  }
  
sub populate_sqlite
  {
     my ($hsp,$dsid,$org) = @_;
     my $pval = $hsp->pval;
     my $pid = $hsp->percent_id;
     my $qgap = $hsp->query_gaps;
     my $sgap = $hsp->subject_gaps;
     my $score = $hsp->score;
     my $psim = $hsp->percent_sim;
     my $length = $hsp->length;
     my $strand = $hsp->strand;
     my $match = $hsp->match;
     my $qalign = $hsp->query_alignment;
     my $salign = $hsp->subject_alignment;
     my $align = $hsp->alignment;
     my $qname = $hsp->query_name;
     my $sname = $hsp->subject_name;
     my $name = $hsp->number."_".$dsid;
     my $hsp_num = $hsp->number;
    
     my $dbh = DBI->connect("dbi:SQLite:dbname=".$cogeweb->sqlitefile,"","");
     my $query_length = $hsp->query_stop > $hsp->query_start ? (($hsp->query_stop) - ($hsp->query_start) + 1) : (($hsp->query_start) - ($hsp->query_stop) + 1);
     my ($qstart, $qstop) = $hsp->query_stop > $hsp->query_start ? ($hsp->query_start, $hsp->query_stop) : ($hsp->query_stop, $hsp->query_start);
#     ($qstart, $qstop) = ($qstop, $qstart) if $qstop < $qstart;
     my $qmismatch = $query_length - $hsp->match;
     
     my $subject_length = $hsp->subject_stop > $hsp->subject_start ? (($hsp->subject_stop) - ($hsp->subject_start) + 1) : (($hsp->subject_start) - ($hsp->subject_stop) + 1);
     my ($sstart, $sstop) = ($hsp->subject_start,$hsp->subject_stop);
#     ($sstart, $sstop) = ($sstop, $sstart) if $sstop < $sstart;
     my $smismatch = $subject_length - $hsp->match;
     
     my $statement = qq{
       INSERT INTO hsp_data (name, eval, pid, psim, score, qgap, sgap,match,qmismatch,smismatch, strand, length, qstart,qstop, sstart, sstop,qalign,salign,align,qname,sname,hsp_num,org) values ("$name", "$pval", "$pid","$psim", "$score", $qgap, $sgap, $match,$qmismatch, $smismatch, "$strand",$length,$qstart, $qstop, $sstart, $sstop,"$qalign","$salign","$align","$qname","$sname",$hsp_num,"$org") 
     };
     print STDERR $statement unless $dbh->do($statement);

     #populate sequence_info table

     $statement = "SELECT name FROM sequence_info where name = '$qname'";
     my $val = $dbh->selectall_arrayref($statement);
     unless ($val->[0][0])
       {
	 my $qlength = $hsp->query_length;
	 $statement = qq{
INSERT INTO sequence_info (name, type, length) values ("$qname","query","$qlength")
};
	 print STDERR $statement unless $dbh->do($statement);
       }
     $statement = "SELECT name FROM sequence_info where name = '$sname'";
     $val = $dbh->selectall_arrayref($statement);
     unless ($val->[0][0])
       {
	 my $slength = $hsp->subject_length;
	 $statement = qq{
INSERT INTO sequence_info (name, type, length) values ("$sname","subject","$slength")
};
	 print STDERR $statement unless $dbh->do($statement);
       }
   }
  
sub get_nearby_feats
  {
    my %opts = @_;
    my $hsp_id = $opts{num};
    my $filename = $opts{basefile};
    $cogeweb = initialize_basefile(basename=>$filename, prog=>"CoGeBlast");
    my $dbh = DBI->connect("dbi:SQLite:dbname=".$cogeweb->sqlitefile,"","");
    $hsp_id =~ s/^table_row// if $hsp_id =~ /table_row/;
    $hsp_id =~ s/^\d+_// if $hsp_id =~ tr/_/_/ > 1;
    
    my $name = "None";
#    my $fid;
    #my $checkbox = " ";
    my $distance = ">250";
    my $sth = $dbh->prepare(qq{SELECT * FROM hsp_data WHERE name = ?});
    my ($sstart, $sstop,$sname);
    my ($hsp_num,$dsid) = $hsp_id =~ /^(\d+)_(\d+)$/;
    $sth->execute($hsp_id) || die "unable to execute";
    while (my $info = $sth->fetchrow_hashref())
      {
	$sstart=$info->{sstart};
	$sstop = $info->{sstop};
	$sname = $info->{sname};
      }
    my ($start,$stop) = ($sstart,$sstop);
    my ($chr) = $sname =~ /chromosome: (.*?),/;
    my @feat;
    my $count = 0;
    my $mid = ($stop+$start)/2;
    my $coge = CoGeX->dbconnect();
    my $cogedb = DBI->connect($coge->db_connection_string,$coge->db_name,$coge->db_passwd);
    my $query =qq!

select * from (
  (SELECT * FROM ((SELECT * FROM feature where start<=$mid and dataset_id = $dsid and chromosome = '$chr' ORDER BY start DESC  LIMIT 1) 
   UNION (SELECT * FROM feature where start>=$mid and dataset_id = $dsid and chromosome = '$chr' ORDER BY start LIMIT 1)) as u)
  UNION
  (SELECT * FROM ((SELECT * FROM feature where stop<=$mid and dataset_id = $dsid and chromosome = '$chr' ORDER BY stop   DESC  LIMIT 1) 
   UNION (SELECT * FROM feature where stop>=$mid and dataset_id = $dsid and chromosome = '$chr' ORDER BY stop LIMIT 1)) as v)
   ) as w
order by abs((start + stop)/2 - $mid) LIMIT 10

!;
    my $handle = $cogedb->prepare($query);
    $handle->execute();
    my $feat;
    my $new_checkbox_info;
    my $min_dist;
    while (my $res = $handle->fetchrow_arrayref())
      {
	my $fid = $res->[0];
	my ($tmpfeat) = $coge->resultset('Feature')->find($fid);
	next unless $tmpfeat->type->name =~ /gene/i || $tmpfeat->type->name =~ /rna/i || $tmpfeat->type->name =~ /cds/i;
	$feat = $tmpfeat unless $feat;
	$min_dist = abs($tmpfeat->start-$mid) unless defined $min_dist;
	my $newmin = abs($tmpfeat->start-$mid) < abs($tmpfeat->stop-$mid) ?  abs($tmpfeat->start-$mid) : abs($tmpfeat->stop-$mid);
	$feat = $tmpfeat if $newmin < $min_dist;
	$min_dist = $newmin if $newmin < $min_dist;
      }
    if ($feat)
      {
	if (($start >= $feat->start && $start <= $feat->stop) || ($stop >= $feat->start && $stop <= $feat->stop) )
	  {
	    $distance = "overlapping";
	  }
	else
	  {
	    $distance = abs(($stop+$start)/2-($feat->stop+$feat->start)/2);
	  }
	($name) = $feat->names;
	$name = qq{<a href="#" title="Click for Feature Information" onclick=update_info_box('}.$feat->id."_".$hsp_num."_".$dsid."')>$name</a>";
	$new_checkbox_info = #$distance eq "overlapping" ? 
	  $hsp_id."_".$chr."_".$sstart."no,".$feat->id."_".$hsp_id;# : 
#	  $hsp_id."_".$chr."_".$sstart."no,".$hsp_id."_".$chr."_".$sstart."_".$feat->id."_".$distance;
      }
    else
      {
	$distance = "No neighboring features found";
      }
    
    return $name,$distance,$hsp_id,$new_checkbox_info;
    
  }


sub export_fasta_file
  {
    my $accn_list = shift;
    $accn_list =~ s/^,//;
    $accn_list =~ s/,$//;
    my $url = "FastaView.pl?";
    my @list;
    foreach my $accn (split /,/,$accn_list)
    {
		next if $accn =~ /no$/;
		my ($featid) = $accn =~ m/^(\d+)_\d+_\d+$/;
		push @list,$featid;
	}
	my %seen = ();
    @list = grep {!$seen{$_}++} @list;
    foreach my $featid( @list)
    {
    	$url .= "featid=$featid&";
    }
	$url =~s/&$//;
	return $url;
  }
  
sub export_to_excel
  {
    my $accn_list = shift;
    my $filename = shift;
    $cogeweb = initialize_basefile(basename=>$filename, prog=>"CoGeBlast");
    my $dbh = DBI->connect("dbi:SQLite:dbname=".$cogeweb->sqlitefile,"","");
    my $sth = $dbh->prepare(qq{SELECT * FROM hsp_data WHERE name = ?});

    my $workbook = Spreadsheet::WriteExcel->new("$TEMPDIR/Excel_$filename.xls");
    $workbook->set_tempdir("$TEMPDIR");
    my $worksheet = $workbook->add_worksheet();
    $accn_list =~ s/^,//;
    $accn_list =~ s/,$//;
    my $i = 1;
    if ($accn_list =~/no/)
    {  	 
   	 $worksheet->write(0,0,"Organism");
   	 $worksheet->write(0,1,"Chr");
   	 $worksheet->write(0,2,"Position");
   	 $worksheet->write(0,3,"HSP No.");
   	 $worksheet->write(0,4,"E-value");
   	 $worksheet->write(0,5,"Percent ID");
   	 $worksheet->write(0,6,"Score");
    }
    else
    {
     $worksheet->write(0,0,"Organism");
   	 $worksheet->write(0,1,"Chr");
   	 $worksheet->write(0,2,"Position");
   	 $worksheet->write(0,3,"HSP No.");
   	 $worksheet->write(0,4,"E-value");
   	 $worksheet->write(0,5,"Percent ID");
   	 $worksheet->write(0,6,"Score");
   	 $worksheet->write(0,7,"Closest Feature");
   	 $worksheet->write(0,8,"Distance");
    }
    
    my ($org,$chr,$pos,$hsp_no,$eval,$pid,$score,$distance,$featid,$dsid);
    foreach my $accn (split /,/,$accn_list)
    {
      if($accn =~/no/)
      {
      	 my $accn_with_commas = $accn;
        $accn_with_commas =~ tr/_/,/;
      	($hsp_no,$dsid,$chr,$pos) = $accn_with_commas =~ /(\d+),(\d+),(\w*_?\d+),(\d+)no/;
      	my $ds = $coge->resultset("Dataset")->find($dsid);
      	$org = $ds->organism->name;
      	$sth->execute($hsp_no."_".$dsid) || die "unable to execute";
      	while (my $info = $sth->fetchrow_hashref())
	      {
 	        $eval = $info->{eval};
 	        $pid = $info->{pid};
 	        $score = $info->{score};
	      }
      	
      	$worksheet->write($i,0,$org);
   		$worksheet->write($i,1,$chr);
   	 	$worksheet->write($i,2,$pos);
   	 	$worksheet->write($i,3,$hsp_no);
   	 	$worksheet->write($i,4,$eval);
   	 	$worksheet->write($i,5,$pid);
   	 	$worksheet->write($i,6,$score);
      	
      }
      else
      {
      	if ($accn =~ tr/_/_/ > 2)
      	{
      		my $accn_with_commas = $accn;
      		$accn_with_commas =~ tr/_/,/;
      		($hsp_no,$dsid,$chr,$pos,$featid,$distance) = $accn_with_commas =~ /(\d+),(\d+),(\w*_?\d+),(\d+),(\d+),(\d+.?\d*)/;
      	}
      	else
      	{
      		($featid,$hsp_no,$dsid) = $accn =~ /(\d+)_(\d+)_(\d+)/;
      		$distance = "overlapping";
      	}
      	
      	    my $ds = $coge->resultset("Dataset")->find($dsid);
      		$org = $ds->organism->name;
      		$sth->execute($hsp_no."_".$dsid) || die "unable to execute";
      		while (my $info = $sth->fetchrow_hashref())
	     	 {
 	       	  $eval = $info->{eval};
 	      	  $pid = $info->{pid};
 	      	  $score = $info->{score};
 	      	  $pos = $info->{sstart} unless $pos;
	     	 }
           my ($feat) = $coge->resultset("Feature")->find($featid);
   	       my ($name) = sort $feat->names;
   	       $chr = $feat->chr if $accn =~ tr/_/_/ < 3;
   	    
   	    $worksheet->write($i,0,$org);
   		$worksheet->write($i,1,$chr);
   	 	$worksheet->write($i,2,$pos);
   	 	$worksheet->write($i,3,$hsp_no);
   	 	$worksheet->write($i,4,$eval);
   	 	$worksheet->write($i,5,$pid);
   	 	$worksheet->write($i,6,$score);
   	 	$worksheet->write($i,7,"http://synteny.cnr.berkeley.edu/CoGe/FeatView.pl?accn=$name",$name);
   	 	$worksheet->write($i,8,$distance);
      
      }
   	 
   	 $i++;
   	}
   	
   	$workbook->close() or die "Error closing file: $!";
   	 
   	return "$TEMPURL/Excel_$filename.xls";
  }
   	 
sub generate_tab_deliminated
  {
    my $accn_list = shift;
    my $filename = shift;
    my $cogeweb = initialize_basefile(basename=>$filename, prog=>"CoGeBlast");
    my $dbh = DBI->connect("dbi:SQLite:dbname=".$cogeweb->sqlitefile,"","");
    
    my $sth = $dbh->prepare(qq{SELECT * FROM hsp_data WHERE name = ?});

    $accn_list =~ s/^,//;
    $accn_list =~ s/,$//;
    
    my $str = "Name\tHSP No.\tE-value\tPerc ID\tScore\tOrganism\n";
    
    foreach my $accn (split /,/,$accn_list)
      {
	next if $accn =~ /no$/;
	my ($featid,$hsp_num,$dsid) = $accn =~ m/^(\d+)_(\d+)_(\d+)$/;
	my $ds = $coge->resultset("Dataset")->find($dsid);
	my ($feat) = $coge->resultset("Feature")->find($featid);
	
	my ($name) = sort $feat->names;
	my $org = $ds->organism->name;
	
	$sth->execute($hsp_num."_".$dsid) || die "unable to execute";
	my ($pval,$pid,$score);
	while (my $info = $sth->fetchrow_hashref())
	  {
	    $pval = $info->{eval};
	    $pid = $info->{pid};
	    $score = $info->{score};
	  }
	$str .= "$name\t$hsp_num\t$pval\t$pid\t$score\t$org\n";
      }
    $str =~ s/\n$//;
	
    open(NEW,"> $TEMPDIR/tab_delim$filename.tabbed");
    print  NEW $str;
    close NEW;
    return "$TEMPURL/tab_delim$filename.tabbed";
  }
   	 
sub generate_feat_list
  {
    my $accn_list = shift;
    my $filename = shift;
    
    $accn_list =~ s/^,//;
    $accn_list =~ s/,$//;
    
    my $url = "FeatList.pl?";
    my @list;
    foreach my $accn (split /,/,$accn_list)
    {
		next if $accn =~ /no$/;
		my ($featid) = $accn =~ m/^(\d+)_\d+_\d+$/;
		push @list,$featid;
	}
	my %seen = ();
    @list = grep {!$seen{$_}++} @list;
    foreach my $featid( @list)
    {
    	$url .= "fid=$featid&";
    }
	$url =~s/&$//;
	return $url;
  }

sub dataset_description_for_org
  {
    my $org_id = shift;
    my $seq_type_id = shift;
    my ($org) = $coge->resultset('Organism')->find($org_id);
    my $html = "Current datasets for ".$org->name;
    $html .= ": ".$org->description if $org->description;
    $html .= "<table>";
    my $i = 0;
    my $total_length=0;
    my $type;
    foreach my $ds (sort {$a->name cmp $b->name} $org->current_datasets(type=>$seq_type_id))
      {
	$type = $ds->sequence_type unless $type;
	my $name = $ds->name;
	$name .= " (v ".$ds->version.")" if $ds->version;
        $name = "<a href=GenomeView.pl?dsid=".$ds->id." target=_new>".$name."</a>";
	my $length = 0;
	my $chr_count =0;
	my $plasmid = 0;
	my $contig =0;
	my @chrs = $ds->get_chromosomes();
	foreach my $chr (@chrs)
	  {
	    $plasmid = 1 if $chr =~ /plasmid/i;
	    $contig = 1 if $chr =~ /contig/i;
	    $length += $ds->last_chromosome_position($chr);
	    $chr_count++;
	  }
	$total_length += $length;
	$length = commify ($length);
	$html .= "<tr";
	$html .= " class='even'" if $i % 2 == 0;
	$html .= "><td>";#.join ("<td>", $name, $source)."\n";
	$html .= join (", length: ", $name, $length. "bp");
	$html .= "; chr count: $chr_count";
	$html .= " (".join (", ", @chrs).")" if scalar @chrs < 5;
	if ($plasmid || $contig)
	  {
	    $html .= " <span class=small>(";
	    $html .= "plastmid" if $plasmid;
	    $html .= " " if $plasmid && $contig;
	    $html .= "contig" if $contig;
	    $html .= ")</span>";
	  }
	$i++;
      }
    $html.="</table>";
    $total_length = commify($total_length);
    $html .= "<span class=species>Total Length: $total_length</span>";
    if ($type)
      {
	$html .= "<br><span class=small>Type: ".$type->name;
	$html .= ", ".$type->description if $type->description;
	$html .= "</span>";
      }
    return $html;
  }
      
sub commify {
        my $input = shift;
        $input = reverse $input;
        $input =~ s<(\d\d\d)(?=\d)(?!\d*\.)><$1,>g;
        return scalar reverse $input;
}

sub export_hsp_info
{
	my $accn = shift;
	my $filename = shift;
	my $cogeweb = initialize_basefile(basename=>$filename, prog=>"CoGeBlast");
    my $dbh = DBI->connect("dbi:SQLite:dbname=".$cogeweb->sqlitefile,"","");

    my $sth = $dbh->prepare(qq{SELECT * FROM hsp_data ORDER BY org,hsp_num});
    $sth->execute();
    
    my $str = qq{Org\tChr\tPosition\tStrand\tHSP No.\tPercent ID\tAligment length\tE-value\tScore\tMatch\tQuery Mismatch\tQuery Gap\tQuery Length\tSubject Mismatch\tSubject Gap\tSubject Length\tQuery HSP Sequence\tSubject HSP Sequence\n};
    
    my ($org,$chr,$pos,$hsp_num,$pid,$align_length,$eval,$score,$match,$strand,$query_mismatch,$query_gap,$query_length,$subject_mismatch,$subject_gap,$subject_length,$query_seq,$subject_seq,$align_seq,$sstart,$sstop,$sname,$align);
    while (my $row = $sth->fetchrow_hashref())
      {
        ($hsp_num) = $row->{name}=~/(\d+)_\d+/;
        $org = $row->{org};
        $eval = $row->{eval};
 	    $pid = $row->{pid};
 	    $score = $row->{score};
 	    $query_gap = $row->{qgap};
 	    $subject_gap = $row->{sgap};
 	    $match = $row->{match};
 	    $query_mismatch = $row->{qmismatch};
 	    $subject_mismatch = $row->{smismatch};
 	    $strand = $row->{strand};
 	    $align_length = $row->{length};
		$sstart = $row->{sstart};
		$sstop =  $row->{sstop};
 	    $query_seq = $row->{qalign};
  	    $subject_seq = $row->{salign};
	    $align = $row->{align};
	    $sname = $row->{sname};
	    
	    ($chr) = $sname =~ /chromosome: (.*?),/;
	    $pos = $sstart." - ".$sstop;
	    $align_seq = $query_seq."<br>".$align."<br>".$subject_seq;
	    
	    $query_seq =~ s/-//g;
	    $query_length = length $query_seq;
	    
	    $subject_seq =~ s/-//g;
	    $subject_length = length $subject_seq;
	    
	    $str .= "$org\t$chr\t$pos\t$strand\t$hsp_num\t$pid\t$align_length\t$eval\t$score\t$match\t$query_mismatch\t$query_gap\t$query_length\t$subject_mismatch\t$subject_gap\t$subject_length\t$query_seq\t$subject_seq\n";
      }
    open(NEW,"> $TEMPDIR/tab_delim_$filename.txt");
    print  NEW $str;
    close NEW;
    return "$TEMPURL/tab_delim_$filename.txt";
    
    }
    
sub export_hsp_query_fasta
{
	my $filename = shift;
	my $cogeweb = initialize_basefile(basename=>$filename, prog=>"CoGeBlast");
    my $dbh = DBI->connect("dbi:SQLite:dbname=".$cogeweb->sqlitefile,"","");

    my $sth = $dbh->prepare(qq{SELECT * FROM hsp_data});
    $sth->execute();
    my ($fasta,$query_seq,$name,$qstart,$qstop, $hsp_num);
    while (my $row = $sth->fetchrow_hashref())
      {
      	$query_seq = $row->{qalign};
      	$name = $row->{qname};
      	$name =~ s/^\s+//;
      	$name =~ s/\s+$//;
      	$qstart = $row->{qstart};
      	$qstop = $row->{qstop};
	$hsp_num = $row->{hsp_num};
      	$query_seq =~ s/-//g;
      	$fasta .= ">HSP".$hsp_num."_".$name.", Subsequence: ".$qstart."-".$qstop."\n".$query_seq."\n";
      }
    open(NEW,"> $TEMPDIR/query_fasta_$filename.txt");
    print  NEW $fasta;
    close NEW;
    return "$TEMPURL/query_fasta_$filename.txt";
      	
}

sub export_hsp_subject_fasta
  {
    my $filename = shift;
    my $dna = 1 if shift;
    my $cogeweb = initialize_basefile(basename=>$filename, prog=>"CoGeBlast");
    my $dbh = DBI->connect("dbi:SQLite:dbname=".$cogeweb->sqlitefile,"","");

    my $sth = $dbh->prepare(qq{SELECT * FROM hsp_data});
    $sth->execute();
    my ($fasta,$subject_seq,$name,$sstart,$sstop, $hsp_num);
    while (my $row = $sth->fetchrow_hashref())
      {
      	$name = $row->{sname};
	my ($chr, $dsid) = $name =~ /chromosome:\s+(\S+?),.*?id:\s+(\d+)/;
	my $strand = $row->{strand} =~ /\+/ ? 1 : -1;
      	$sstart = $row->{sstart};
      	$sstop = $row->{sstop};
	$hsp_num = $row->{hsp_num};
      	$name =~ s/^\s+//;
      	$name =~ s/\s+$//;
	$name =~ s/\s+/_/g;
	if ($dna)
	  {
	    my ($ds) = $coge->resultset('Dataset')->find($dsid);
	    $subject_seq = $ds->get_genomic_sequence(chr=>$chr, strand=>$strand, start=>$sstart, stop=>$sstop);
	  }
	else
	  {
	    $subject_seq = $row->{salign};
	    $subject_seq =~ s/-//g;
	  }
	$fasta .= ">HSP".$hsp_num."_".$name.", Location: ".$sstart."-".$sstop."($strand)"."\n".$subject_seq."\n";
      }
    open(NEW,"> $TEMPDIR/subject_fasta_$filename.txt");
    print  NEW $fasta;
    close NEW;
    return "$TEMPURL/subject_fasta_$filename.txt";
      	
}

sub export_alignment_file
{
	my $filename = shift;
	my $cogeweb = initialize_basefile(basename=>$filename, prog=>"CoGeBlast");
    my $dbh = DBI->connect("dbi:SQLite:dbname=".$cogeweb->sqlitefile,"","");

    my $sth = $dbh->prepare(qq{SELECT * FROM hsp_data ORDER BY org,hsp_num});
    $sth->execute();
    my ($sname,$qname,$qstop,$qstart,$sstart,$sstop,$align,$qseq,$sseq,$str,$hsp_num);
    while (my $row = $sth->fetchrow_hashref())
    {
    	$hsp_num = $row->{hsp_num};
    	$sseq = $row->{salign};
      	$sname = $row->{sname};
      	$sstart = $row->{sstart};
      	$sstop = $row->{sstop};
      	$qseq = $row->{qalign};
      	$qname = $row->{qname};
      	$qname =~ s/^\s+//;
      	$qname =~ s/\s+$//;
        $sname =~ s/^\s+//;
      	$sname =~ s/\s+$//;
      	$qstart = $row->{qstart};
      	$qstop = $row->{qstop};
      	$align = $row->{align};
      	$str .= "HSP: ".$hsp_num."\n>Query: ".$qname.", Subsequence: ".$qstart."-".$qstop."\n>Subject: ".$sname.", Location: ".$sstart."-".$sstop."\n".$qseq."\n".$align."\n".$sseq."\n\n";
    }
    $str =~ s/\n+$//;
    open(NEW,"> $TEMPDIR/alignment_file_$filename.txt");
    print  NEW $str;
    close NEW;
    return "$TEMPURL/alignment_file_$filename.txt";
}

sub save_settings_cogeblast
  {
    my %opts = @_;
    my $opts = Dumper \%opts;
    my $item = save_settings(opts=>$opts, user=>$USER, page=>$PAGE_NAME);
  }

sub color_pallet
  {
    my %opts = @_;
    my $start = $opts{start} || [0,200,0];
    my $offset = $opts{offset} || 75;
    my $num_seqs = $opts{num_seqs};

    my @colors;
    my %set = (HSP_NUM=>1,
	       RED=>$start->[0],
	       GREEN=>$start->[1],
	       BLUE=>$start->[2]);
#    push @colors, \%set;
    my $temp = [@$start];
    for (my $i = 1; $i <= $num_seqs; $i++)
      {
	my @color;
	@color = @$temp;
	push @colors, [$color[0],
		       $color[1],
		       $color[2],
		      ];
	foreach (@$temp)
	  {
	    $_ = 200 if $_ < 0;
	  }
	unless ($i%3)
	  {
	    $temp = [map {int($_/2)} @color];
	  }
	else
	  {
	    $temp =[$temp->[2], $temp->[0], $temp->[1]];
	  }
	$temp = [map {$_-1} @$temp] unless ($i%6);

      }
    return wantarray ? @colors : \@colors;
  }

