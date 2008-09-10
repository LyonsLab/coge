#! /usr/bin/perl -w
use strict;
use CGI;
use CGI::Carp 'fatalsToBrowser';
use CGI::Ajax;
use CoGe::Accessory::LogUser;
use CoGe::Accessory::Web;
use CoGe::Accessory::Restricted_orgs;
use CoGeX;
use Data::Dumper;
use HTML::Template;
use Digest::MD5 qw(md5_hex);
use Parallel::ForkManager;
use GD;
use File::Path;
use Mail::Mailer;
use Benchmark;

$ENV{PATH} = "/opt/apache2/CoGe/";
umask(0);
use vars qw( $DATE $DEBUG $DIR $URL $USER $FORM $coge $cogeweb $FORMATDB $BLAST $DATADIR $FASTADIR $BLASTDBDIR $DIAGSDIR $MAX_PROC $DAG_TOOL $PYTHON $TANDEM_FINDER $FILTER_REPETITIVE_MATCHES $RUN_DAGCHAINER $FIND_NEARBY $PLOT_DAG);
$DEBUG = 0;
$DIR = "/opt/apache/CoGe/";
$URL = "/CoGe/";
$FORMATDB = "/usr/bin/formatdb";
$BLAST = "/usr/bin/blast -a 8 -K 80 -m 8 -e 0.05";
$DATADIR = "/home/apache/data/";
$DIAGSDIR = "$DATADIR/diags";
$FASTADIR = $DATADIR.'/fasta/';
$BLASTDBDIR = $DATADIR.'/blast/db/';
$MAX_PROC=8;
$PYTHON = "/usr/bin/python";
$DAG_TOOL = $DIR."/bin/parepair/dag_tools.py";
$TANDEM_FINDER = $DIR."/bin/parepair/tandems.py -d 5 -s -r"; #-d option is the distance (in genes) between dups -- not sure if the -s and -r options are needed -- they create dups files based on the input file name
$FILTER_REPETITIVE_MATCHES = $DIR."/bin/parepair/filter_repetitive_matches.pl 200000"; #filters multiple matches to the same accession within the "window size", here set to 80000 nt
$RUN_DAGCHAINER = $DIR."/bin/parepair/run_dagchainer.pl -E 0.05";
$FIND_NEARBY = $DIR."/bin/parepair/find_nearby.py -d 200000";
$PLOT_DAG = $PYTHON ." ".$DIR."/bin/parepair/plot_dag.py";

$| = 1; # turn off buffering
$DATE = sprintf( "%04d-%02d-%02d %02d:%02d:%02d",
                 sub { ($_[5]+1900, $_[4]+1, $_[3]),$_[2],$_[1],$_[0] }->(localtime));
$FORM = new CGI;
($USER) = CoGe::Accessory::LogUser->get_user();
my %ajax = CoGe::Accessory::Web::ajax_func();

$coge = CoGeX->dbconnect();
#$coge->storage->debugobj(new DBIxProfiler());
#$coge->storage->debug(1);
my $pj = new CGI::Ajax(
		       get_orgs => \&get_orgs,
		       get_datasets => \&get_datasets,
		       get_previous_analyses=>\&get_previous_analyses,
		       get_pair_info=> \&get_pair_info,
		       go=>\&go,
		       check_address_validity=>\&check_address_validity,
		       generate_basefile=>\&generate_basefile,
		       get_iframe=>\&get_iframe,
		       %ajax,
		      );
print $pj->build_html($FORM, \&gen_html);
#print "Content-Type: text/html\n\n";print gen_html($FORM);


sub gen_html
  {
    my $html;
    my ($body) = gen_body();
    my $template = HTML::Template->new(filename=>'/opt/apache/CoGe/tmpl/generic_page.tmpl');
    $template->param(TITLE=>'SynMap: Powered by <a href=http://dagchainer.sourceforge.net/ target=_new>DAGChainer</a>');
    $template->param(HEAD=>qq{});
    my $name = $USER->user_name;
    $name = $USER->first_name if $USER->first_name;
    $name .= " ".$USER->last_name if $USER->first_name && $USER->last_name;
    $template->param(USER=>$name);
    
    $template->param(LOGON=>1) unless $USER->user_name eq "public";
    $template->param(DATE=>$DATE);
    #$template->param(ADJUST_BOX=>1);
    $template->param(LOGO_PNG=>"SynMap-logo.png");
    $template->param(BODY=>$body);
    $html .= $template->output;
    return $html;
  }

sub gen_body
  {
    my $form = shift || $FORM;
    my $template = HTML::Template->new(filename=>'/opt/apache/CoGe/tmpl/SynMap.tmpl');
    for (my $i=1; $i<=2;$i++)
      {
	my $name = $form->param('org_name'.$i);
	my $desc = $form->param('org_desc'.$i);
	my $oid = $form->param('oid'.$i);
	my ($org) = $coge->resultset('Organism')->resolve($oid) if $oid;
	$name = $org->name if $org;
	$name = "Search" unless $name;
	$template->param('ORG_NAME'.$i=>$name);
	$desc = "Search" unless $desc;
	$template->param('ORG_DESC'.$i=>$desc);
	$name = "" if $name =~ /Search/;
	$template->param('ORG_LIST'.$i=>get_orgs(name=>$name,i=>$i));
      }
    my $file = $form->param('file');
    if($file)
    {
    	my $results = read_file($file);
    	$template->param(RESULTS=>$results);
    }
    return $template->output;
  }
  
  sub read_file
  {
  	my $file = shift;

	my $html;
	open (IN, "/opt/apache/CoGe/tmp/SynMap/".$file) || die "can't open $file for reading: $!";
    while (<IN>)
      {
		$html .= $_;
      }
    close IN;
    return $html;
  }

sub get_orgs
  {
    my %opts = @_;
    my $name = $opts{name};
    my $desc = $opts{desc};
    my $i = $opts{i};
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
    unless (@opts && ($name || $desc)) 
      {
	$html .=  qq{<input type = hidden name="org_id$i" id="org_id$i">};
	return $html;
      }

    $html .= qq{<SELECT id="org_id$i" SIZE="5" MULTIPLE onChange="\$('#ds_info'+$i).html('<div class=dna_small class=loading class=small>loading. . .</div>'); get_datasets(['args__oid','org_id$i', 'args__masked','masked$i', 'args__seq_type','seq_type$i'],['ds_info$i']);check_previous_analyses();" >\n};
    $html .= join ("\n", @opts);
    $html .= "\n</SELECT>\n";
    $html =~ s/OPTION/OPTION SELECTED/;
    return $html;
  }

sub get_datasets
  {
    my %opts = @_;
    my $oid = $opts{oid};
    my $masked = $opts{masked};
    my $seq_type = $opts{seq_type};
    my $html; 
    my ($org) = $coge->resultset("Organism")->resolve($oid);
    return unless $org;
    my $orgname = $org->name;
    $orgname .= ": ".$org->description if $org->description;
    $html .= qq{<div><span class="oblique">Organism:</span> $orgname\n</div>};
    $html .= qq{<div><span class="oblique">Current Datasets: </span>};
    my $i =0;
    my @ds = $org->current_datasets(type=>$masked);
    unless (@ds)
      {
	$masked=1; #unmasked sequence
	@ds = $org->current_datasets(type=>$masked);
      }
    my ($type) = $coge->resultset('GenomicSequenceType')->find($masked);
    $html.= "<span class=alert> ".$type->name."</span>" if $type;
    $html.= ", ".$type->description if $type && $type->description;
    $html .= "</DIV>";
    my $has_cds=0;
    foreach my $ds (@ds)
      {
	my $name = $ds->name;
        $name .= ": ".$ds->description if $ds->description;
        $name = "<a href=GenomeView.pl?dsid=".$ds->id." target=_new>".$name."</a>";
        my $source = $ds->datasource->name;
        $source .= ": ".$ds->datasource->description if $ds->datasource->description;
        $source = "<a href=".$ds->datasource->link." target=_new>".$source."</a>" if $ds->datasource->link;
        $source =~ s/href=/href=http:\/\// unless $source =~ /http/;
        $html .= "<div";
        $html .= " class='even'" if $i % 2 == 0;
        $html .= ">".join (": ", $name, $source)."</div>\n";
	if ($seq_type == 1) #want to use CDS.  Let's see if any exist for this dataset
	  {
	    foreach my $ft ($coge->resultset('Feature')->search(
								{
								 feature_type_id=>3,
								 dataset_id=>$ds->id,
								},
								{
								 limit=>1
								}
							       ))
	      {
		$has_cds=1;
	      }
	  }

        $i++;
      }
    $html .= "<span class=alert>WARNING: datasets contain no CDS features</span>" if ($seq_type == 1 &&  !$has_cds);
    return $html;
  }

sub gen_fasta
  {
    my %opts = @_;
    my $oid = $opts{oid};
    my $masked = $opts{masked};
    my $seq_type = $opts{seq_type};
    my ($org_name, $md5,$ds_list);
    ($org_name, $md5,$ds_list,$masked) = gen_org_name(oid=>$oid, masked=>$masked, seq_type=>$seq_type);
    my $file = $FASTADIR."/$md5.fasta";
    my $res;
    if (-r $file)
      {
	write_log("fasta file for *".$org_name."* ($md5) exists", $cogeweb->logfile);
	$res = 1;
      }
    else
      {
	$res = generate_fasta(dslist=>$ds_list, file=>$file, type=>$seq_type) unless -r $file;
      }
    return $file, $md5, $org_name, $masked if $res;
    return 0;
  }
    
sub gen_org_name
  {
    my %opts = @_;
    my $oid = $opts{oid};
    my $masked = $opts{masked} || 1;
    my $seq_type = $opts{seq_type} || 1;
    my ($org) = $coge->resultset('Organism')->resolve($oid);
    my @ds = $org->current_datasets(type=>$masked);
    unless (@ds)
      {
	$masked=1; #unmasked sequence
	@ds = $org->current_datasets(type=>$masked);
      }
    @ds = sort {$a->id <=> $b->id }@ds;
    return unless @ds;
    my $org_name = $ds[0]->organism->name;
    my $title = $org_name ." $seq_type (";
#    my $title = $org_name ." CDS (";
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
    return ($org_name, $md5, \@ds, $masked);
  }

sub generate_fasta
  {
    my %opts = @_;
    my $dslist = $opts{dslist};
    my $file = $opts{file};
    my $type = $opts{type};
    $file = $FASTADIR."/$file" unless $file =~ /$FASTADIR/;
    write_log("creating fasta file.", $cogeweb->logfile);
    open (OUT, ">$file") || die "Can't open $file for writing: $!";;
    foreach my $ds (@$dslist)
      {
	if ($type eq "CDS")
	  {
	    foreach my $feat (sort {$a->chromosome cmp $b->chromosome || $a->start <=> $b->start} $ds->features({feature_type_id=>3}))
	      {
		my ($chr) = $feat->chromosome;#=~/(\d+)/;
		my $name;
		foreach my $n ($feat->names)
		  {
		    $name = $n;
		    last unless $name =~ /\s/;
		  }
		$name =~ s/\s+/_/g;
		my $title = join ("||",$chr, $feat->start, $feat->stop, $name, $feat->strand, $feat->type->name, $feat->id);
		my $seq = $feat->genomic_sequence;
		next unless $seq;
		print OUT ">".$title."\n";
		print OUT $seq,"\n";
	      }
	  }
	else
	  {
	    foreach my $chr (sort $ds->get_chromosomes)
	      {
#		my $title = join ("||",$chr, 1, $ds->last_chromosome_position($chr), "Chr_$chr",1, "genomic", "N/A");
		my $seq = $ds->get_genomic_sequence(chr=>$chr);
		next unless $seq;
		print OUT ">".$chr."\n";
		print OUT $seq,"\n";
	      }
	  }
      }
    close OUT;
    return 1 if -r $file;
    write_log("Error with fasta file creation", $cogeweb->logfile);
    return 0;
  }

sub gen_blastdb
  {
    my %opts = @_;
    my $md5 = $opts{md5};
    my $fasta = $opts{fasta};
    my $org_name = $opts{org_name};
    my $blastdb = "$BLASTDBDIR/$md5";
    my $res = 0;
    if (-r $blastdb.".nsq")
      {
	write_log("blastdb file for *".$org_name."* ($md5) exists", $cogeweb->logfile);
	$res = 1;
      }
    else
      {
	$res = generate_blast_db(fasta=>$fasta, blastdb=>$blastdb, org=>$org_name);
      }
    return $blastdb if $res;
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
    write_log("creating blastdb for *".$org."* ($blastdb)",$cogeweb->logfile);
    `$command`;
    return 1 if -r "$blastdb.nsq";
    write_log("error creating blastdb for $org ($blastdb)",$cogeweb->logfile);
    return 0;
  }
  
  sub generate_basefile
{
	$cogeweb = initialize_basefile(prog=>"SynMap");
	#print STDERR $cogeweb->basefilename,"\n";
	return $cogeweb->basefilename;
}

sub run_blast
  {
    my %opts = @_;
    my $fasta = $opts{fasta};
    my $blastdb = $opts{blastdb};
    my $outfile = $opts{outfile};
    my $prog = $opts{prog};
    $prog = "blastn" unless $prog;
    if (-r $outfile)
      {
	write_log("blastfile $outfile already exists",$cogeweb->logfile);
	return 1;
      }
    my $pre_command = "$BLAST -p $prog -o $outfile -i $fasta -d $blastdb";
    my $x;
    ($x, $pre_command) = check_taint($pre_command);
    write_log("running $pre_command" ,$cogeweb->logfile);
    `$pre_command`;
    return 1 if -r $outfile;
  }

sub run_dag_tools
    {
      my %opts = @_;
      my $query = $opts{query};
      my $subject = $opts{subject};
      my $blast = $opts{blast};
      my $outfile = $opts{outfile};
      my $seq_type1 = $opts{seq_type1};
      my $seq_type2 = $opts{seq_type2};
      if (-r $outfile)
      {
	write_log("run dag_tools: file $outfile already exists",$cogeweb->logfile);
	return 1;
      }
      my $query_dup_file= $opts{query_dup_files};
      my $subject_dup_file= $opts{subject_dup_files};
      my $cmd = "$PYTHON $DAG_TOOL -q \"$query\" -s \"$subject\" -b $blast";
      $cmd .= " -c" if $seq_type1 eq "genomic" && $seq_type2 eq "genomic";
      $cmd .= " --query_dups $query_dup_file" if $query_dup_file;
      $cmd .= " --subject_dups $subject_dup_file" if $subject_dup_file;
      $cmd .=  " > $outfile";
      write_log("run dag_tools: running $cmd",$cogeweb->logfile);
      `$cmd`;
      return 1 if -r $outfile;
    }

sub run_tandem_finder
  {
    my %opts = @_;
    my $infile = $opts{infile}; #dag file produced by dat_tools.py
    my $outfile = $opts{outfile};
    if (-r $outfile)
      {
	write_log("run_tandem_filter: file $outfile already exists",$cogeweb->logfile);
	return 1;
      }
    my $cmd = "$PYTHON $TANDEM_FINDER -i $infile > $outfile";
    write_log("run_tandem_filter: running $cmd", $cogeweb->logfile);
    `$cmd`;
    return 1 if -r $outfile;
  }

sub run_filter_repetitive_matches
  {
    my %opts = @_;
    my $infile = $opts{infile};
    my $outfile = $opts{outfile};
    if (-r $outfile)
      {
	write_log("run_filter_repetitive+_matches: file $outfile already exists",$cogeweb->logfile);
	return 1;
      }
    my $cmd = "$FILTER_REPETITIVE_MATCHES < $infile > $outfile";
    write_log("run_filter_repetitive+_matches: running $cmd", $cogeweb->logfile);
    `$cmd`;
    return 1 if -r $outfile;
  }

sub run_dagchainer
  {
    my %opts = @_;
    my $infile = $opts{infile};
    my $D = $opts{D}; #distance allowed between two matches in basepairs
    my $g = $opts{g}; #length of a gap in bp (ave distance expected between two syntenic genes)
    my $A = $opts{A}; #Minium number of Aligned Pairs 
    my $outfile = $infile;
    $outfile .= "_D$D" if $D;
    $outfile .= "_g$g" if $g;
    $outfile .= "_A$A" if $A;
    $outfile .= ".aligncoords";
    if (-r $outfile)
      {
	write_log("run dagchainer: file $outfile already exists",$cogeweb->logfile);
	return $outfile;
      }
    my $cmd = "$RUN_DAGCHAINER -i $infile";
    $cmd .= " -D $D" if $D;
    $cmd .= " -g $g" if $g;
    $cmd .= " -A $A" if $A;
    write_log("run dagchainer: running $cmd", $cogeweb->logfile);
    `$cmd`;
    `mv $infile.aligncoords $outfile`;
    return $outfile
  }

sub run_find_nearby
  {
    my %opts = @_;
    my $infile = $opts{infile};
    my $dag_all_file = $opts{dag_all_file};
    my $outfile = $opts{outfile};
    if (-r $outfile)
      {
	write_log("run find_nearby: file $outfile already exists",$cogeweb->logfile);
	return 1;
      }
    my $cmd = "$PYTHON $FIND_NEARBY --diags=$infile --all=$dag_all_file > $outfile";
    write_log("run find_nearby: running $cmd", $cogeweb->logfile);
    `$cmd`;
    return 1 if -r $outfile;
  }

sub add_GEvo_links
  {
    my %opts = @_;
    my $infile = $opts{infile};
    open (IN, $infile);
    open (OUT,">$infile.tmp");
    while (<IN>)
      {
	chomp;
	if (/^#/ || /GEvo/)
	  {
	    s/toxic/synteny.cnr/;
	    print OUT $_,"\n";
	    next;
	  }
	s/^\s+//;
	next unless $_;
	my @line = split/\t/;
	my @feat1 = split/\|\|/,$line[1];
	my @feat2 = split/\|\|/,$line[5];
	my $link = "http://synteny.cnr.berkeley.edu/CoGe/GEvo.pl?accn1=".$feat1[3]."&fid1=".$feat1[6]."&accn2=".$feat2[3]."&fid2=".$feat2[6] if $feat1[3] && $feat1[6] && $feat2[3] && $feat2[6];
	print OUT $_;
	print OUT "\t",$link if $link;
	print OUT "\n";
      }
    close IN;
    close OUT;
    my $cmd = "/bin/mv $infile.tmp $infile";
    `$cmd`;
  }

sub add_reverse_match#this is for when there is a self-self comparison.  DAGchainer, for some reason, is only adding one diag.  For example, if chr_14 vs chr_10 have a diag, chr_10 vs chr_14 does not.
  {
    my %opts = @_;
    my $infile = $opts{infile};
    open (IN, $infile);
    my $stuff;
    my $skip =0;
    while (<IN>)
      {
	chomp;
	s/^\s+//;
	$skip = 1 if /GEvo\.pl/; #GEvo links have been added, this file was generated on a previous run.  Skip!
	last if ($skip);
	next unless $_;
	my @line = split/\s+/;
	if (/^#/)
	  {
	    my $chr1 = $line[2];
	    my $chr2 = $line[4];
	    $chr1 =~ s/^a//;
	    $chr2 =~ s/^b//;
	    next if $chr1 eq $chr2;
	    $line[2] = "b".$chr2;
	    $line[4] = "a".$chr1;
	    $stuff .= join (" ",@line)."\n";
	    next;
	  }
	my $chr1 = $line[0];
	my $chr2 = $line[4];
	$chr1 =~ s/^a//;
	$chr2 =~ s/^b//;
	next if $chr1 eq $chr2;
	my @tmp1 = @line[1..3];
	my @tmp2 = @line[5..7];
	@line[1..3] = @tmp2;
	@line[5..7] = @tmp1;
	$line[0] = "a".$chr2;
	$line[4] = "b".$chr1;
	$stuff .= join ("\t",@line)."\n";
      }
    return if $skip;
    close IN;
    open (OUT, ">>$infile");
    print OUT $stuff;
    close OUT;

  }

sub retrieve_chrs_with_diags
  {
    my %opts = @_;
    my $infile = $opts{infile};
    my %data;
    open (IN, $infile);
    while (<IN>)
      {
	next unless /^#/;
	my @line = split/\s+/;
	my $chr1 = $line[2];
	my $chr2 = $line[4];
	$chr1 =~ s/^a//;
	$chr2 =~ s/^a//;
	$chr1 =~ s/^b//;
	$chr2 =~ s/^b//;
	$data{$chr1}{$chr2}=1;
	$data{$chr2}{$chr1}=1;
      }
    close IN;
    return \%data;
  }

sub generate_dotplot
  {
    my %opts = @_;
    my $dag = $opts{dag};
    my $coords = $opts{coords};
    my $outfile = $opts{outfile};
    my $qchr = $opts{qchr};
    my $schr = $opts{schr};
    my $q_dsid = $opts{qdsid};
    my $s_dsid = $opts{sdsid};
    my $q_label = $opts{qlabel};
    my $s_label = $opts{slabel};
    my $q_max = $opts{"q_max"};
    my $s_max = $opts{"s_max"};
    my $regen_images = $opts{regen_images}=~/true/i ? 1 : 0;
    if (-r $outfile && !$regen_images)
      {
	write_log("generate dotplot: file $outfile already exists",$cogeweb->logfile);
	return 1;
      }
    my $cmd = "$PLOT_DAG -d $dag -a $coords --html $outfile --q_dsid $q_dsid --s_dsid $s_dsid --qchr $qchr --schr $schr --q_max $q_max --s_max $s_max";
    $cmd .= " --q_label=$q_label" if $q_label;
    $cmd .= " --s_label=$s_label" if $s_label;
    write_log("generate dotplot: running $cmd", $cogeweb->logfile);
    `$cmd`;
    return 1 if -r $outfile;
  }

sub go
  {
    my %opts = @_;
    my $oid1 = $opts{oid1};
    my $oid2 = $opts{oid2};
    my $masked1 = $opts{masked1};
    my $masked2 = $opts{masked2};
    my $dagchainer_D = $opts{D};
    my $dagchainer_g = $opts{g};
    my $dagchainer_A = $opts{A};
    my $chr_length_limit = $opts{chr_length_limit};
    my $show_diags_only = $opts{show_diags_only};
    my $regen_images = $opts{regen_images};
    my $email = $opts{email};
    my $job_title = $opts{jobtitle};
    
    my $basename = $opts{basename};
    $cogeweb = initialize_basefile(basename=>$basename, prog=>"SynMap");
    
    $email = 0 if check_address_validity($email) eq 'invalid';
    $show_diags_only = $show_diags_only eq "true" ? 1 : 0;
    my $flip_output = $opts{flip_output};
    $flip_output = $flip_output eq "true" ? 1 : 0;
    my $blast = $opts{blast};
    $blast = $blast == 2 ? "tblastx" : "blastn";
    my $seq_type1 = $opts{seq_type1};
    $seq_type1 = $seq_type1 == 2 ? "genomic" : "CDS";
    my $seq_type2 = $opts{seq_type2};
    $seq_type2 = $seq_type2 == 2 ? "genomic" : "CDS";
    unless ($oid1 && $oid2)
      {
	return "<span class=alert>You must select two organisms.</span>"
      }
    #$cogeweb = initialize_basefile(prog=>"SynMap");
    ##generate fasta files and blastdbs
    my $pm = new Parallel::ForkManager($MAX_PROC);
    my @oids = ([$oid1,$masked1, $seq_type1]);
    push @oids, [$oid2,$masked2,$seq_type2] unless $oid1 == $oid2 && $masked1 == $masked2 && $seq_type1 eq $seq_type2;
    foreach my $item (@oids)
      {
	$pm->start and next;
	my $oid = $item->[0];
	my $masked = $item->[1];
	my $seq_type = $item->[2];
	my ($fasta,$md5,$org_name) = gen_fasta(oid=>$oid, masked=>$masked, seq_type=>$seq_type);
	gen_blastdb(md5=>$md5,fasta=>$fasta,org_name=>$org_name);
	$pm->finish;
      }
    $pm->wait_all_children();
    my ($fasta1,$md51,$org_name1);
    my ($fasta2,$md52,$org_name2);
    ($fasta1,$md51,$org_name1,$masked1) = gen_fasta(oid=>$oid1, masked=>$masked1, seq_type=>$seq_type1);
    ($fasta2,$md52,$org_name2,$masked2) = gen_fasta(oid=>$oid2, masked=>$masked2, seq_type=>$seq_type2);
    ($oid1, $masked1, $md51,$org_name1,$fasta1,$seq_type1, $oid2,
    $masked2, $md52,$org_name2,$fasta2, $seq_type2) = ($oid2,
    $masked2, $md52,$org_name2,$fasta2, $seq_type2,$oid1, $masked1,
    $md51,$org_name1,$fasta1, $seq_type1) if ($md52 lt $md51);
    my ($org1) = $coge->resultset('Organism')->resolve($oid1);
    my ($org2) = $coge->resultset('Organism')->resolve($oid2);
    unless ($fasta1 && $fasta2)
      {
	return "<span class=alert>Something went wrong generating the fasta files: ".$cogeweb->logfile."</span>";
      }
    else{
    	write_log("Completed fasta creation", $cogeweb->logfile);
    }
    
    my ($blastdb1) = gen_blastdb(md5=>$md51,fasta=>$fasta1,org_name=>$org_name1);
    my ($blastdb2) = gen_blastdb(md5=>$md52,fasta=>$fasta2,org_name=>$org_name2);
    unless ($blastdb1 && $blastdb2)
      {
	return "<span class=alert>Something went wrong generating the blastdb files: ".$cogeweb->logfile."</span>";
      }
    else{
    	write_log("Completed blastdb creation", $cogeweb->logfile);
    }
    my $html;

    #need to blast each org against itself for finding local dups, then to one another
    my %org_dirs = ( 
		    $org_name1."_".$org_name2=>{fasta=>$fasta1,
						db=>$blastdb2,
						basename=>$md51."_".$md52.".$masked1-$masked2.$seq_type1-$seq_type2.$blast",
						},
		    $org_name1."_".$org_name1=>{fasta=>$fasta1,
						db=>$blastdb1,
						basename=>$md51."_".$md51.".$masked1-$masked1.$seq_type1-$seq_type1.$blast",
						},
		    $org_name2."_".$org_name2=>{fasta=>$fasta2,
						db=>$blastdb2,
						basename=>$md52."_".$md52.".$masked2-$masked2.$seq_type2-$seq_type2.$blast",
						},
		    );
    foreach my $org_dir (keys %org_dirs)
      {
	my $tmp = $org_dir;
	$tmp =~ s/\///g;
	$tmp =~ s/\s+/_/g;
	$tmp =~ s/\(//g;
	$tmp =~ s/\)//g;
	my $outfile = $DIAGSDIR."/".$tmp;
	mkpath ($outfile,0,0777) unless -d $outfile;
	$org_dirs{$org_dir}{dir}=$outfile;
	warn "didn't make $outfile: $!" unless -d $outfile;
	$outfile .= "/".$org_dirs{$org_dir}{basename};
	$org_dirs{$org_dir}{blastfile}=$outfile.".blast";
      }
    #blast!
    foreach my $key (keys %org_dirs)
      {
	$pm->start and next;
	my $fasta = $org_dirs{$key}{fasta};
	my $db = $org_dirs{$key}{db};
	my $outfile = $org_dirs{$key}{blastfile};
	run_blast(fasta=>$fasta, blastdb=>$db, outfile=>$outfile, prog=>$blast) unless -r $outfile;
	$pm->finish;
      }
    $pm->wait_all_children();
	write_log("Completed blast run", $cogeweb->logfile);
    #Find local dups
    my $dag_file11 = $org_dirs{$org_name1."_".$org_name1}{dir}."/".$org_dirs{$org_name1."_".$org_name1}{basename}.".dag";
    run_dag_tools(query=>"a".$md51, subject=>"b".$md51, blast=>$org_dirs{$org_name1."_".$org_name1}{blastfile}, outfile=>$dag_file11, seq_type1=>$seq_type1, seq_type1=>$seq_type1);
    my $dag_file22 = $org_dirs{$org_name2."_".$org_name2}{dir}."/".$org_dirs{$org_name2."_".$org_name2}{basename}.".dag";
    run_dag_tools(query=>"a".$md52, subject=>"b".$md52, blast=>$org_dirs{$org_name2."_".$org_name2}{blastfile}, outfile=>$dag_file22, seq_type1=>$seq_type2, seq_type1=>$seq_type2);
    my $dup_file1  = $org_dirs{$org_name1."_".$org_name1}{dir}."/".$org_dirs{$org_name1."_".$org_name1}{basename}.".dups";
    run_tandem_finder(infile=>$dag_file11,outfile=>$dup_file1);
    my $dup_file2  = $org_dirs{$org_name2."_".$org_name2}{dir}."/".$org_dirs{$org_name2."_".$org_name2}{basename}.".dups";
    run_tandem_finder(infile=>$dag_file22,outfile=>$dup_file2);

    #prepare dag for synteny analysis
    my $dag_file12 = $org_dirs{$org_name1."_".$org_name2}{dir}."/".$org_dirs{$org_name1."_".$org_name2}{basename}.".dag";
    run_dag_tools(query=>"a".$md51, subject=>"b".$md52, blast=>$org_dirs{$org_name1."_".$org_name2}{blastfile}, outfile=>$dag_file12.".all", query_dup_file=>$dup_file1,subject_dup_file=>$dup_file2, seq_type1=>$seq_type1, seq_type2=>$seq_type2);
    #remove repetitive matches
    run_filter_repetitive_matches(infile=>$dag_file12.".all",outfile=>$dag_file12);
    #run dagchainer
    my $dagchainer_file = run_dagchainer(infile=>$dag_file12, D=>$dagchainer_D, g=>$dagchainer_g,A=>$dagchainer_A);
    write_log("Completed dagchainer run", $cogeweb->logfile);
    if (-r $dagchainer_file)
      {
	#add pairs that were skipped by dagchainer
	my $tmp = $dagchainer_file;
	$tmp =~ s/aligncoords/all\.aligncoords/;
	run_find_nearby(infile=>$dagchainer_file, dag_all_file=>$dag_file12.".all", outfile=>$tmp);
	add_reverse_match(infile=>$tmp) if $md51 eq $md52;
	my $chrs_w_diags = retrieve_chrs_with_diags(infile=>$tmp) if $show_diags_only;
	#generate dotplot images
	my %chr1;
	my %chr2;
	my @ds1 = $org1->current_datasets(type=>$masked1);
	unless (@ds1)
	  {
	    $masked1=1; #unmasked sequence
	    @ds1 = $org1->current_datasets(type=>$masked1);
	  }
	my @ds2 = $org2->current_datasets(type=>$masked2);
	unless (@ds2)
	  {
	    $masked2=1; #unmasked sequence
	    @ds2 = $org2->current_datasets(type=>$masked2);
	  }

	foreach my $ds (@ds1)
	  {
	    foreach my $chr ($ds->get_chromosomes)
	      {
		next if $chr =~ /random/;
		my $last = $ds->last_chromosome_position($chr);
		next if $last < $chr_length_limit;
		$chr1{$md51."_".$chr}= {dsid => $ds->id,
					chr_end=>$last};
	      }
	  }
	foreach my $ds (@ds2)
	  {
	    foreach my $chr ($ds->get_chromosomes)
	      {
		next if $chr =~ /random/;
		my $last = $ds->last_chromosome_position($chr);
		next if $last < $chr_length_limit;
		$chr2{$md52."_".$chr}= {dsid => $ds->id,
					chr_end=>$last};
	      }
	  }

	my $name1 = $org_name1;
	$name1 =~ s/\s+/_/g;
	$name1 =~ s/\(//g;
	$name1 =~ s/\)//g;

	my $name2 = $org_name2;
	$name2 =~ s/\s+/_/g;
	$name2 =~ s/\(//g;
	$name2 =~ s/\)//g;
	my $qlead = "a";
	my $slead = "b";
	my $qdsid = "qdsid";
	my $sdsid = "sdsid";
	if ($flip_output)
	  {
	    my %tmp = %chr2;
	    %chr2 = %chr1;
	    %chr1 = %tmp;
	    ($name1, $name2) = ($name2, $name1);
	    $qlead = "b";
	    $slead = "a";
	    $qdsid = "sdsid";
	    $sdsid = "qdsid";
	  }

	foreach my $chr2 (sort keys %chr2)
	  {
	    my @chr2 = split/_/,$chr2;
	    foreach my $chr1 (sort keys %chr1)
	      {
		my @chr1 = split/_/,$chr1;
		$pm->start and next;
		my $out = $org_dirs{$org_name1."_".$org_name2}{dir}."/html/";
		mkpath ($out,0,0777) unless -d $out;
		$out .= "$chr1-$masked1-$seq_type1"."_".$chr2."-$masked2-$seq_type2"."_D$dagchainer_D"."_g$dagchainer_g"."_A$dagchainer_A".".$blast.html";
		generate_dotplot(dag=>$dag_file12.".all", coords=>$tmp, outfile=>$out, qchr=>$qlead.$chr1, schr=>$slead.$chr2, $qdsid=>$chr1{$chr1}{dsid}, $sdsid=>$chr2{$chr2}{dsid},qlabel=>$name1.":".$chr1[-1],slabel=>$name2.":".$chr2[-1], q_max=>$chr1{$chr1}{chr_end}, s_max=>$chr2{$chr2}{chr_end}, regen_images=>$regen_images);
		$pm->finish;
	      }
	  }
	$pm->wait_all_children();
	my $count = scalar (keys %chr1) * scalar (keys%chr2);
#	my $bm1 = new Benchmark;
	$html .= "<table cellspacing=0 cellpadding=0>";
	my $id_num =1;
	foreach my $tchr2 (sort keys %chr2)
	  {
	    $html .= "<tr align=center>";
	    foreach my $tchr1 (sort keys %chr1)
	      {
		my ($chr1, $chr2) = ($tchr1, $tchr2);
		my $out = $org_dirs{$org_name1."_".$org_name2}{dir}."/html/$chr1-$masked1-$seq_type1"."_".$chr2."-$masked2-$seq_type2"."_D$dagchainer_D"."_g$dagchainer_g"."_A$dagchainer_A".".$blast.html";
		if ($show_diags_only)
		  {
		    unless ($chrs_w_diags->{$chr1}{$chr2})
		      {
			next;
		      }
		  }

		my $png = $out;
		$png =~ s/html$/png/;
		#image dimensions hard-coded because of the time it takes to get w and h from image directly
		my ($w, $h) = (576,504);
#		if (-r $png)
#		  {
#		    my $img = GD::Image->new($png);
#		    ($w,$h) = $img->getBounds();
#		  }
		my $stuff;
		if (-r $out)
		  {
		    $png =~ s/$DATADIR/$URL\/data/;
		    if ($count > 4)
		      {
			$stuff = qq{

<img id=close$id_num src="picts/delete.png" style="display: none; float:right; position: relative; top: 60px; right: 40px;" onclick="\$('#img$id_num').toggle();\$('#close$id_num').toggle();\$('#iframe$id_num').html('');" valign=top \>
<div id=iframe$id_num></div>
};
			my ($tmpw, $tmph) = (sprintf("%0f",$w/4),sprintf("%0f",$h/4));
			$stuff .= qq{
<img id= img$id_num src=$png width=$tmpw height=$tmph onclick="\$('#img$id_num').toggle();get_iframe(['args__src','args__$out'],['iframe$id_num']);\$('#close$id_num').toggle();" \>
};
		      }
		    else
		      {
			$out =~ s/$DATADIR/$URL\/data/;
			$html = qq{<iframe src=$out frameborder=0 width=$w height=$h scrolling=no></iframe>};
		      }
		  }
		else {$stuff = " ";}
		$html .= "<td id=cell$id_num>$stuff";
		$id_num++;
	      }
	  }
	$html .= "</table>";
#	my $bm2 = new Benchmark;
#	my $diff = timediff($bm2, $bm1);
#	print STDERR "HTML table generation took: ",timestr($diff),"\n";
	add_GEvo_links (infile=>$tmp);
	$tmp =~ s/$DATADIR/$URL\/data/;
	$html .= "<br><a href=$tmp target=_new>Syntolog file with GEvo links</a><br>";
      }
    

    foreach my $org_dir (keys %org_dirs)
      {
	my $output = $org_dirs{$org_dir}{blastfile};
	$output =~ s/$DIR/$URL/;
	$html .= "<a href=$output target=_new>Blast results for $org_dir</a><br>";;	
      }
    my $log = $cogeweb->logfile;
    $log =~ s/$DIR/$URL/;
    $html .= "<a href=$log target=_new>log</a><br>";

    email_results(email=>$email,html=>$html,org1=>$org_name1,org2=>$org_name2, jobtitle=>$job_title) if $email;

    return $html;
  }

sub get_previous_analyses
  {
    my %opts = @_;
    my $oid1 = $opts{oid1};
    my $oid2 = $opts{oid2};
    return unless $oid1 && $oid2;
    my $masked1 = $opts{masked1};
    my $masked2 = $opts{masked2};
    my $seq_type1 = $opts{seq_type1};
    $seq_type1 = $seq_type1 == 2 ? "genomic" : "CDS";
    my $seq_type2 = $opts{seq_type2};
    $seq_type2 = $seq_type2 == 2 ? "genomic" : "CDS";
    my ($org_name1, $md51) = gen_org_name(oid=>$oid1, masked=>$masked1, seq_type=>$seq_type1);
    my ($org_name2, $md52) = gen_org_name(oid=>$oid2, masked=>$masked2, seq_type=>$seq_type2);
    ($oid1, $masked1, $md51,$org_name1,$seq_type1, $oid2, $masked2, $md52,$org_name2,$seq_type2) = ($oid2, $masked2, $md52,$org_name2,$seq_type2,$oid1, $masked1, $md51,$org_name1,$seq_type1) if ($md52 lt $md51);
    my $dir = "$org_name1"."_".$org_name2;
    $dir =~ s/\///g;
    $dir =~ s/\s+/_/g;
    $dir =~ s/\(//g;
    $dir =~ s/\)//g;
    $dir = "$DIAGSDIR/".$dir;
    my @items;
    if (-d $dir)
      {
	opendir (DIR, $dir);
	while (my $file = readdir(DIR))
	  {
	    next unless $file =~ /all\.aligncoords/;
	    my ($D, $g, $A) = $file =~ /D(\d+)_g(\d+)_A(\d+)/;
	    my $blast = $file =~ /blastn/ ? "BlastN" : "TBlastX";
	    my ($mask1, $mask2, $type1, $type2) = $file =~ /\.(\d)-(\d)\.(\w+)-(\w+)/;
	    next unless ($mask1 && $mask2 && $type1 && $type2);
	    my %data = (D=>$D,
			g=>$g,
			A=>$A,
			blast=>$blast);
	    ($mask1, $type1, $mask2, $type2) = ($mask2, $type2, $mask1, $type1)  if ($md52 lt $md51);
	    $data{mask1} = $mask1-1;
	    $data{mask2} = $mask2-1;
	    $mask1 = $mask1 == "1" ? "unmasked" : "masked";
	    $mask2 = $mask2 == "1" ? "unmasked" : "masked";
	    $data{name} = "$org_name1: $mask1-$type1, $org_name2: $mask2-$type2";
	    $type1 = $type1 eq "CDS" ? 0 : 1; 
	    $type2 = $type2 eq "CDS" ? 0 : 1; 
	    $data{type1} = $type1;
	    $data{type2} = $type2;
	    push @items, \%data;
	  }
      }
    return unless @items;
    my $html = qq{
<select id="prev_params" size=4 multiple onChange="update_params();">
};
    
    foreach (sort {$a->{g}<=>$b->{g} } @items)
      {
	my $blast = $_->{blast} =~ /^blastn$/i ? 0 : 1;
	my $val = join ("_",$_->{g},$_->{D},$_->{A}, $oid1, $_->{mask1},$_->{type1},$oid2, $_->{mask2},$_->{type2}, $blast);
	my $name = $_->{blast}.": g:".$_->{g}." D:".$_->{D}." A:".$_->{A}." ".$_->{name};
	$html .= qq{
 <option value="$val">$name
};
      }
    $html .= "</select>";
    return "  Previously run parameters:<br>$html";
  }

sub get_pair_info
  {
    my @anno;
    foreach my $fid (@_)
      {
	unless ($fid =~ /^\d+$/)
	  {
	    push @anno, $fid;
	    next;
	  }
	my $feat = $coge->resultset('Feature')->find($fid);
	my $anno = "Name: ".join (", ",map {"<a class=\"data link\" href=\"/CoGe/FeatView.pl?accn=".$_."\" target=_new>".$_."</a>"} $feat->names);
	my $location = "Chr ".$feat->chromosome." ";
	$location .= $feat->start."-".$feat->stop;
	$location .="(".$feat->strand.")";
	push @anno, $anno."<br>".$location;
      }
    return unless @anno;
    return "<table class=small><tr>".join ("<td>",@anno)."</table>";
  }

  
  sub check_address_validity {
	my $address = shift;
	return 'valid' unless $address;
	my $validity = $address =~/^[_a-zA-Z0-9-]+(\.[_a-zA-Z0-9-]+)*@[a-zA-Z0-9-]+(\.[a-zA-Z0-9-]+)*\.(([0-9]{1,3})|([a-zA-Z]{2,3})|(aero|coop|info|museum|name))$/ ? 'valid' : 'invalid';
	return $validity;
}

sub email_results {
	my %opts = @_;
	my $email_address = $opts{email};
	my $html = $opts{html};
	my $org_name1 = $opts{org1};
	my $org_name2 = $opts{org2};
	my $job_title = $opts{jobtitle};

	my $file = $cogeweb->basefile."_results.data";
    open(NEW,"> $file") || die "Cannot Save $!\n";
    print NEW $html;
    close NEW;
    
    my $subject = "SynMap Results";
    $subject .= ": $job_title" if $job_title;
    
    ($file) = $file =~/SynMap\/(.+\.data)/;
    
	my $server = $ENV{HTTP_HOST};
	
	my $url = "http://".$server."/CoGe/SynMap.pl?file=".$file;
	
	my $mailer = Mail::Mailer->new("sendmail");
	$mailer->open({From	=> 'CoGE <coge_results@synteny.cnr.berkeley.edu>',
		       To	=> $email_address,
		       Subject	=> $subject,
		      })
	  or die "Can't open: $!\n";
	my $username = $USER->user_name;
    $username = $USER->first_name if $USER->first_name;
    $username .= " ".$USER->last_name if $USER->first_name && $USER->last_name;
	my $body = qq{Dear $username,
		
Thank you for using SynMap! The results from your latest analysis between $org_name1 and $org_name2 are ready, and can be viewed here:
	
$url

These results will remain on our servers for approximately 24 hours; please save them to your own computer before they are removed.
			
Thank you for using the CoGe Software Package.
	
- The CoGe Team
};
	
	print $mailer $body;
	$mailer->close();
}

sub get_iframe
  {
    my %args = @_;
    my $src = $args{src};
#    my $png = $src;
#    $png =~ s/html$/png/;
    $src =~ s/$DATADIR/$URL\/data/;
#    my $img = GD::Image->new($png);
#    my ($w,$h) = $img->getBounds();
    my ($w, $h) = (576,504);
    my $html = qq{<iframe src=$src frameborder=0 width=$w height=$h scrolling=no></iframe>};
    return $html;
  }
