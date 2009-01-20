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
use LWP::Simple;
$ENV{PATH} = "/opt/apache2/CoGe/";
umask(0);
use vars qw( $DATE $DEBUG $DIR $URL $USER $FORM $coge $cogeweb $FORMATDB $BLAST $DATADIR $FASTADIR $BLASTDBDIR $DIAGSDIR $MAX_PROC $DAG_TOOL $PYTHON $TANDEM_FINDER $FILTER_REPETITIVE_MATCHES $RUN_DAGCHAINER $FIND_NEARBY $CONVERT_TO_GENE_ORDER $DOTPLOT);
$DEBUG = 0;
$DIR = "/opt/apache/CoGe/";
$URL = "/CoGe/";
$FORMATDB = "/usr/bin/formatdb";
$BLAST = "/usr/bin/blast -a 8 -K 80 -m 8 -e 0.05";
$DATADIR = "$DIR/data/";
$DIAGSDIR = "$DATADIR/diags";
$FASTADIR = $DATADIR.'/fasta/';
$BLASTDBDIR = $DATADIR.'/blast/db/';
$MAX_PROC=8;
$PYTHON = "/usr/bin/python";
$DAG_TOOL = $DIR."/bin/dagchainer/dag_tools.py";
$TANDEM_FINDER = $DIR."/bin/dagchainer/tandems.py -d 5 -s -r"; #-d option is the distance (in genes) between dups -- not sure if the -s and -r options are needed -- they create dups files based on the input file name
$FILTER_REPETITIVE_MATCHES = $DIR."/bin/dagchainer/filter_repetitive_matches.pl 200000"; #filters multiple matches to the same accession within the "window size", here set to 80000 nt
$RUN_DAGCHAINER = $DIR."/bin/dagchainer/DAGCHAINER/run_DAG_chainer.pl -E 0.05";
$FIND_NEARBY = $DIR."/bin/dagchainer/find_nearby.py -d 200000";
$DOTPLOT = $DIR."/bin/dotplot.pl";#Eric gives up waiting for new and improved to really work, and writes his own.
$CONVERT_TO_GENE_ORDER = $DIR."/bin/dagchainer/convert_to_gene_order.pl";  #this needs to be implemented
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
		       get_dotplot=>\&get_dotplot,
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

    $html .= qq{<SELECT id="org_id$i" SIZE="5" MULTIPLE onChange="\$('#ds_info'+$i).html('<div class=dna_small class=loading class=small>loading. . .</div>'); check_previous_analyses();get_datasets(['args__oid','org_id$i', 'args__masked','masked$i', 'args__seq_type','seq_type$i'],['ds_info$i']);" >\n};
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
    return unless $oid;
    my $html; 
    my ($org) = $coge->resultset("Organism")->resolve($oid);
    return unless $org;
    my $orgname = $org->name;
    $orgname = "<a href=\"GenomeView.pl?oid=".$org->id."\" target=_new>$orgname</a>: ".$org->description if $org->description;
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
    my $total_length=0;
    foreach my $ds (sort {$a->name cmp $b->name} @ds)
      {
	my $name = $ds->name;
        #$name .= ": ".$ds->description if $ds->description;
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
	$length = reverse $length;
	$length =~ s/(\d\d\d)(?=\d)(?!\d*\.)/$1,/g;
	$length = reverse $length;
        $html .= "<div";
        $html .= " class='small";# if $i % 2 == 0;
        $html .= " even'" if $i % 2 == 0;
        $html .= "'>";
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
	$html .= "</div>\n";
#	if ($seq_type == 1) #want to use CDS.  Let's see if any exist for this dataset
#	  {
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
#	  }

        $i++;
      }
    $total_length = reverse $total_length;
    $total_length =~ s/(\d\d\d)(?=\d)(?!\d*\.)/$1,/g;
    $total_length = reverse $total_length;
    $html .= "Total Length: $total_length";
    $html .= "<br><span class=alert>WARNING: datasets contain no CDS features</span>" unless $has_cds;# if ($seq_type == 1 &&  !$has_cds);
    return $html, $has_cds;
  }

sub gen_fasta
  {
    my %opts = @_;
    my $oid = $opts{oid};
    my $masked = $opts{masked};
    my $seq_type = $opts{seq_type};
    my $write_log = $opts{write_log} || 0;
    my ($org_name, $md5,$ds_list);
    ($org_name, $md5,$ds_list,$masked) = gen_org_name(oid=>$oid, masked=>$masked, seq_type=>$seq_type, write_log=>$write_log);
    my $file = $FASTADIR."/$md5.fasta";
    my $res;
    while (-e "$file.running")
      {
	sleep 60;
      }
    if (-r $file)
      {
	write_log("fasta file for *".$org_name."* ($md5) exists", $cogeweb->logfile);
	$res = 1;
      }
    else
      {
	system "touch $file.running"; #track that a blast anlaysis is running for this
	$res = generate_fasta(dslist=>$ds_list, file=>$file, type=>$seq_type) unless -r $file;
	system "rm $file.running" if -r "$file.running"; #remove track file
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
    my $write_log = $opts{write_log} || 0;
    my ($org) = $coge->resultset('Organism')->resolve($oid);
    my @ds = $org->current_datasets(type=>$masked);
    unless (@ds)
      {
	$masked=1; #unmasked sequence
	@ds = $org->current_datasets(type=>$masked);
      }
    @ds = sort {$a->id <=> $b->id }@ds;
    return $org->name unless @ds;
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
    write_log("ORGANISM: ".$title, $cogeweb->logfile) if $write_log;
    my $md5 = md5_hex($title);
    return ($org_name, $md5, \@ds, $masked, $title);
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
    while (-e "$blastdb.running")
      {
	sleep 60;
      }
    if (-r $blastdb.".nsq")
      {
	write_log("blastdb file for *".$org_name."* ($md5) exists", $cogeweb->logfile);
	$res = 1;
      }
    else
      {
	system "touch $blastdb.running"; #track that a blast anlaysis is running for this
	$res = generate_blast_db(fasta=>$fasta, blastdb=>$blastdb, org=>$org_name);
	system "rm $blastdb.running" if -r "$blastdb.running"; #remove track file
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
    while (-e "$outfile.running")
      {
	sleep 60;
      }
    if (-r $outfile)
      {
	unless (-s $outfile)
	  {
	    write_log("WARNING: Blast output file ($outfile) contains no data!" ,$cogeweb->logfile);
	    return 0;
	  }
	write_log("blastfile $outfile already exists",$cogeweb->logfile);
	return 1;
      }
    my $pre_command = "$BLAST -p $prog -o $outfile -i $fasta -d $blastdb";
    my $x;
    system "touch $outfile.running"; #track that a blast anlaysis is running for this
    ($x, $pre_command) = check_taint($pre_command);
    write_log("running $pre_command" ,$cogeweb->logfile);
    `$pre_command`;
    system "rm $outfile.running" if -r "$outfile.running"; #remove track file
    unless (-s $outfile)
      {
	    write_log("WARNING: Problem running $pre_command command.  Blast output file contains no data!" ,$cogeweb->logfile);
	    return 0;
      }
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
      while (-e "$outfile.running")
      {
	sleep 60;
      }
      unless (-r $blast && -s $blast)
	{
	  write_log("WARNING:   Cannot create input file for DAGChainer! Blast output file ($blast) contains no data!" ,$cogeweb->logfile);
	  return 0;
	}
      if (-r $outfile)
      {
	write_log("run dag_tools: file $outfile already exists",$cogeweb->logfile);
	return 1;
      }
      my $query_dup_file= $opts{query_dup_files};
      my $subject_dup_file= $opts{subject_dup_files};
      my $cmd = "$PYTHON $DAG_TOOL -q \"$query\" -s \"$subject\" -b $blast";
      $cmd .= " -c";# if $seq_type1 eq "genomic" && $seq_type2 eq "genomic";
      $cmd .= " --query_dups $query_dup_file" if $query_dup_file;
      $cmd .= " --subject_dups $subject_dup_file" if $subject_dup_file;
      $cmd .=  " > $outfile";
      system "touch $outfile.running"; #track that a blast anlaysis is running for this
      write_log("run dag_tools: running $cmd",$cogeweb->logfile);
      `$cmd`;
      system "rm $outfile.running" if -r "$outfile.running"; #remove track file
      unless (-s $outfile)
	{
	  write_log("WARNING: DAGChainer input file ($outfile) contains no data!" ,$cogeweb->logfile);
	  return 0;
	}
      return 1 if -r $outfile;
    }

sub run_tandem_finder
  {
    my %opts = @_;
    my $infile = $opts{infile}; #dag file produced by dat_tools.py
    my $outfile = $opts{outfile};
    while (-e "$outfile.running")
      {
	sleep 60;
      }
    unless (-r $infile && -s $infile)
	{
	  write_log("WARNING:   Cannot run tandem finder! DAGChainer input file ($infile) contains no data!" ,$cogeweb->logfile);
	  return 0;
	}
    if (-r $outfile)
      {
	write_log("run_tandem_filter: file $outfile already exists",$cogeweb->logfile);
	return 1;
      }
    my $cmd = "$PYTHON $TANDEM_FINDER -i $infile > $outfile";
    system "touch $outfile.running"; #track that a blast anlaysis is running for this
    write_log("run_tandem_filter: running $cmd", $cogeweb->logfile);
    `$cmd`;
    system "rm $outfile.running" if -r "$outfile.running"; #remove track file
    return 1 if -r $outfile;
  }

sub run_filter_repetitive_matches
  {
    my %opts = @_;
    my $infile = $opts{infile};
    my $outfile = $opts{outfile};
    while (-e "$outfile.running")
      {
	sleep 60;
      }
    unless (-r $infile && -s $infile)
	{
	  write_log("WARNING:   Cannot filter repetitive matches! DAGChainer input file ($infile) contains no data!" ,$cogeweb->logfile);
	  return 0;
	}
    if (-r $outfile)
      {
	write_log("run_filter_repetitive+_matches: file $outfile already exists",$cogeweb->logfile);
	return 1;
      }
    my $cmd = "$FILTER_REPETITIVE_MATCHES < $infile > $outfile";
    system "touch $outfile.running"; #track that a blast anlaysis is running for this
    write_log("run_filter_repetitive_matches: running $cmd", $cogeweb->logfile);
    `$cmd`;
    system "rm $outfile.running" if -r "$outfile.running";; #remove track file
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
    while (-e "$outfile.running")
      {
	sleep 60;
      }
    unless (-r $infile && -s $infile)
      {
	  write_log("WARNING:   Cannot run DAGChainer! DAGChainer input file ($infile) contains no data!" ,$cogeweb->logfile);
	  return 0;
	}
    if (-r $outfile)
      {
	write_log("run dagchainer: file $outfile already exists",$cogeweb->logfile);
	return $outfile;
      }
    my $cmd = "$RUN_DAGCHAINER -i $infile";
    $cmd .= " -D $D" if $D;
    $cmd .= " -g $g" if $g;
    $cmd .= " -A $A" if $A;
    system "touch $outfile.running"; #track that a blast anlaysis is running for this
    write_log("run dagchainer: running $cmd", $cogeweb->logfile);
    `$cmd`;
    `mv $infile.aligncoords $outfile`;
    system "rm $outfile.running" if -r "$outfile.running";; #remove track file
    return $outfile
  }

sub run_find_nearby
  {
    my %opts = @_;
    my $infile = $opts{infile};
    my $dag_all_file = $opts{dag_all_file};
    my $outfile = $opts{outfile};
    while (-e "$outfile.running")
      {
	sleep 60;
      }
    if (-r $outfile)
      {
	write_log("run find_nearby: file $outfile already exists",$cogeweb->logfile);
	return 1;
      }
    my $cmd = "$PYTHON $FIND_NEARBY --diags=$infile --all=$dag_all_file > $outfile";
    system "touch $outfile.running"; #track that a blast anlaysis is running for this
    write_log("run find_nearby: running $cmd", $cogeweb->logfile);
    `$cmd`;
    system "rm $outfile.running" if -r "$outfile.running";; #remove track file
    return 1 if -r $outfile;
  }

sub add_GEvo_links
  {
    my %opts = @_;
    my $infile = $opts{infile};
    my $chr1 = $opts{chr1};
    my $chr2 = $opts{chr2};
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
	my $link = "http://synteny.cnr.berkeley.edu/CoGe/GEvo.pl?";
	if ($feat1[6])
	  {
	    $link .= "fid1=".$feat1[6];
	  }
	else
	  {
	    my ($xmin) = sort ($feat1[1], $feat1[2]);
	    my $x = sprintf("%.0f", $xmin+abs($feat1[1]-$feat1[2])/2);
	    $link .= "chr1=".$feat1[0].";dsid1=".$chr1->{$feat1[0]}->id.";x1=".$x;
	  }
	if ($feat2[6])
	  {
	    $link .= ";fid2=".$feat2[6];
	  }
	else
	  {
	    my ($xmin) = sort ($feat2[1], $feat2[2]);
	    my $x = sprintf("%.0f", $xmin+abs($feat2[1]-$feat2[2])/2);
	    $link .= ";chr2=".$feat2[0].";dsid2=".$chr2->{$feat2[0]}->id.";x2=".$x;
	  }
#	accn1=".$feat1[3]."&fid1=".$feat1[6]."&accn2=".$feat2[3]."&fid2=".$feat2[6] if $feat1[3] && $feat1[6] && $feat2[3] && $feat2[6];
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
    my $oid1 = $opts{oid1};
    my $oid2 = $opts{oid2};
    my ($basename) = $coords =~ /([^\/]*).all.aligncoords/;
    my $regen_images = $opts{regen_images}=~/true/i ? 1 : 0;
    my $width = $opts{width} || 1000;
    while (-e "$outfile.running")
      {
	sleep 60;
      }
    if (-r "$outfile.png" && !$regen_images)
      {
	write_log("generate dotplot: file $outfile already exists",$cogeweb->logfile);
	return 1;
      }
    my $cmd = qq{$DOTPLOT -d $dag -a $coords -b $outfile -l 'javascript:synteny_zoom("$oid1","$oid2","$basename","XCHR","YCHR")' -o1 $oid1 -o2 $oid2 -w $width -lt 2};
    system "touch $outfile.running"; #track that a blast anlaysis is running for this
    write_log("generate dotplot: running $cmd", $cogeweb->logfile);
    system "rm $outfile.running" if -r "$outfile.running";; #remove track file
    `$cmd`;
    return 1 if -r "$outfile.html";
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
    my $regen_images = $opts{regen_images};
    my $email = $opts{email};
    my $job_title = $opts{jobtitle};
    my $width = $opts{width};
    my $basename = $opts{basename};
    $cogeweb = initialize_basefile(basename=>$basename, prog=>"SynMap");
    
    $email = 0 if check_address_validity($email) eq 'invalid';
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
    my @oids = ([$oid1,$masked1,$seq_type1]);
    push @oids, [$oid2,$masked2,$seq_type2] unless $oid1 == $oid2 && $masked1 == $masked2 && $seq_type1 eq $seq_type2;
    foreach my $item (@oids)
      {
	$pm->start and next;
	my $oid = $item->[0];
	my $masked = $item->[1];
	my $seq_type = $item->[2];
	my ($fasta,$md5,$org_name) = gen_fasta(oid=>$oid, masked=>$masked, seq_type=>$seq_type, write_log=>1);
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
    $md51,$org_name1,$fasta1, $seq_type1) if ($org_name2 lt $org_name1);
    my ($org1) = $coge->resultset('Organism')->resolve($oid1);
    my ($org2) = $coge->resultset('Organism')->resolve($oid2);
    unless ($fasta1 && $fasta2)
      {
	my $log = $cogeweb->logfile;
	$log =~ s/$DIR/$URL/;
	return "<span class=alert>Something went wrong generating the fasta files: <a href=$log>log file</a></span>";
      }
    else{
    	write_log("Completed fasta creation", $cogeweb->logfile);
    	write_log("", $cogeweb->logfile);
    }
    
    my ($blastdb1) = gen_blastdb(md5=>$md51,fasta=>$fasta1,org_name=>$org_name1);
    my ($blastdb2) = gen_blastdb(md5=>$md52,fasta=>$fasta2,org_name=>$org_name2);
    unless ($blastdb1 && $blastdb2)
      {
	my $log = $cogeweb->logfile;
	$log =~ s/$DIR/$URL/;
	return "<span class=alert>Something went wrong generating the blastdb files: <a href=$log>log file</a></span>";
      }
    else{
    	write_log("Completed blastdb creation", $cogeweb->logfile);
    	write_log("", $cogeweb->logfile);
    }
    my $html;

    #need to blast each org against itself for finding local dups, then to one another
    my $tmp1 = $org_name1;
    my $tmp2 = $org_name2;
    foreach my $tmp ($tmp1, $tmp2)
      {
	$tmp =~ s/\///g;
	$tmp =~ s/\s+/_/g;
	$tmp =~ s/\(//g;
	$tmp =~ s/\)//g;
	$tmp =~ s/://g;
      }
    my $orgkey1 = $org_name1.$masked1.$seq_type1;
    my $orgkey2 = $org_name2.$masked2.$seq_type2;
    my %org_dirs = (
		    $orgkey1."_".$orgkey2=>{fasta=>$fasta1,
					    db=>$blastdb2,
					    basename=>$md51."_".$md52.".$masked1-$masked2.$seq_type1-$seq_type2.$blast",
					    dir=>$DIAGSDIR."/".$tmp1."/".$tmp2,
					   },
		    $orgkey1."_".$orgkey1=>{fasta=>$fasta1,
					    db=>$blastdb1,
					    basename=>$md51."_".$md51.".$masked1-$masked1.$seq_type1-$seq_type1.$blast",
					    dir=>$DIAGSDIR."/".$tmp1."/".$tmp1,
					   },
		    $orgkey2."_".$orgkey2=>{fasta=>$fasta2,
					    db=>$blastdb2,
					    basename=>$md52."_".$md52.".$masked2-$masked2.$seq_type2-$seq_type2.$blast",
					    dir=>$DIAGSDIR."/".$tmp2."/".$tmp2,
					   },
		   );
    foreach my $org_dir (keys %org_dirs)
      {
	my $outfile = $org_dirs{$org_dir}{dir};
	mkpath ($outfile,0,0777) unless -d $outfile;
	warn "didn't create path $outfile: $!" unless -d $outfile;
	$outfile .= "/".$org_dirs{$org_dir}{basename};
	$org_dirs{$org_dir}{blastfile}=$outfile.".blast";
      }
    #blast! use Parallel::ForkManager
    foreach my $key (keys %org_dirs)
      {
	$pm->start and next;
	my $fasta = $org_dirs{$key}{fasta};
	my $db = $org_dirs{$key}{db};
	my $outfile = $org_dirs{$key}{blastfile};
	run_blast(fasta=>$fasta, blastdb=>$db, outfile=>$outfile, prog=>$blast);# unless -r $outfile;
	$pm->finish;
      }
    $pm->wait_all_children();
    #check blast runs for problems;  Not forked in order to keep variables
    my $problem =0;
    foreach my $key (keys %org_dirs)
      {
	my $fasta = $org_dirs{$key}{fasta};
	my $db = $org_dirs{$key}{db};
	my $outfile = $org_dirs{$key}{blastfile};
	my $blast_run = run_blast(fasta=>$fasta, blastdb=>$db, outfile=>$outfile, prog=>$blast);
	$problem=1 unless $blast_run;
      }
    write_log("Completed blast run(s)", $cogeweb->logfile);
    write_log("", $cogeweb->logfile);
    #Find local dups
    my $dag_file11 = $org_dirs{$orgkey1."_".$orgkey1}{dir}."/".$org_dirs{$orgkey1."_".$orgkey1}{basename}.".dag";
    $problem=1 unless (run_dag_tools(query=>"a".$md51, subject=>"b".$md51, blast=>$org_dirs{$orgkey1."_".$orgkey1}{blastfile}, outfile=>$dag_file11, seq_type1=>$seq_type1, seq_type1=>$seq_type1));
    my $dag_file22 = $org_dirs{$orgkey2."_".$orgkey2}{dir}."/".$org_dirs{$orgkey2."_".$orgkey2}{basename}.".dag";
    $problem=1 unless run_dag_tools(query=>"a".$md52, subject=>"b".$md52, blast=>$org_dirs{$orgkey2."_".$orgkey2}{blastfile}, outfile=>$dag_file22, seq_type1=>$seq_type2, seq_type1=>$seq_type2);
    my $dup_file1  = $org_dirs{$orgkey1."_".$orgkey1}{dir}."/".$org_dirs{$orgkey1."_".$orgkey1}{basename}.".dups";
    run_tandem_finder(infile=>$dag_file11,outfile=>$dup_file1);
    my $dup_file2  = $org_dirs{$orgkey2."_".$orgkey2}{dir}."/".$org_dirs{$orgkey2."_".$orgkey2}{basename}.".dups";
    run_tandem_finder(infile=>$dag_file22,outfile=>$dup_file2);

    #prepare dag for synteny analysis
    my $dag_file12 = $org_dirs{$orgkey1."_".$orgkey2}{dir}."/".$org_dirs{$orgkey1."_".$orgkey2}{basename}.".dag";
    $problem=1 unless run_dag_tools(query=>"a".$md51, subject=>"b".$md52, blast=>$org_dirs{$orgkey1."_".$orgkey2}{blastfile}, outfile=>$dag_file12.".all", query_dup_file=>$dup_file1,subject_dup_file=>$dup_file2, seq_type1=>$seq_type1, seq_type2=>$seq_type2);
    #remove repetitive matches
    run_filter_repetitive_matches(infile=>$dag_file12.".all",outfile=>$dag_file12);
    #run dagchainer
    my $dagchainer_file = run_dagchainer(infile=>$dag_file12, D=>$dagchainer_D, g=>$dagchainer_g,A=>$dagchainer_A);
    write_log("Completed dagchainer run", $cogeweb->logfile);
    write_log("", $cogeweb->logfile);
    if (-r $dagchainer_file)
      {
	#add pairs that were skipped by dagchainer
	my $tmp = $dagchainer_file;
	$tmp =~ s/aligncoords/all\.aligncoords/;
	run_find_nearby(infile=>$dagchainer_file, dag_all_file=>$dag_file12.".all", outfile=>$tmp);
	#generate dotplot images
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
	my $org1_length =0;
	my $org2_length =0;
	my $chr1_count = 0;
	my $chr2_count = 0;
	my %chr1;
	my %chr2;
	foreach my $ds (@ds1)
	  {
	    foreach my $chr ($ds->get_chromosomes)
	      {
		$chr1{$chr}=$ds;
		$chr1_count++;
		my $last = $ds->last_chromosome_position($chr);
		$org1_length += $last;
	      }
	  }
	foreach my $ds (@ds2)
	  {
	    foreach my $chr ($ds->get_chromosomes)
	      {
		$chr2{$chr}=$ds;
		$chr2_count++;
		my $last = $ds->last_chromosome_position($chr);
		$org2_length += $last;
	      }
	  }
	my $test = $org1_length > $org2_length ? $org1_length : $org2_length;
	$width = int($test/100000) unless $width;
	$width = 1200 if $width > 1200;
	$width = 500 if $width < 500;
	$width = 1200 if $chr1_count > 9 || $chr2_count > 9;
	$width = 2000 if $chr1_count > 100 || $chr2_count > 100;
	my $qlead = "a";
	my $slead = "b";
	my $qdsid = "qdsid";
	my $sdsid = "sdsid";
	my $out = $org_dirs{$orgkey1."_".$orgkey2}{dir}."/html/";
	mkpath ($out,0,0777) unless -d $out;
	$out .="master_".$org_dirs{$orgkey1."_".$orgkey2}{basename};
	$out .= "_D$dagchainer_D" if $dagchainer_D;
	$out .= "_g$dagchainer_g" if $dagchainer_g;
	$out .= "_A$dagchainer_A" if $dagchainer_A;
	$out .= ".w$width";
	generate_dotplot(dag=>$dag_file12.".all", coords=>$tmp, outfile=>"$out", regen_images=>$regen_images, oid1=>$oid1, oid2=>$oid2, width=>$width);
	add_GEvo_links (infile=>$tmp, chr1=>\%chr1, chr2=>\%chr2);
	$tmp =~ s/$DATADIR/$URL\/data/;
	if (-r "$out.html")
	  {
	    open (IN, "$out.html") || warn "problem opening $out.html for reading\n";
	    $html = "<span class='species small'>$org_name2</span><table><tr valign=top><td valign=top>";
	    $/ = "\n";
	    while (<IN>)
	      {
		next if /<\/?html>/;
		$html .= $_;
	      }
	    close IN;
	    $out =~ s/$DATADIR//;
	    $html =~ s/master.*\.png/data\/$out.png/;
	    warn "$out.html did not parse correctly\n" unless $html =~ /map/i;
	    $html .= qq{
<br><span class="species small">$org_name1</span><br>
Zoomed SynMap:
<table class=species>
<tr>
<td> Display Location:
<td><select name=map_loc id=map_loc>
 <option value="window1">New Window 1
 <option value="window2">New Window 2
 <option value="window3">New Window 3
 <option value="1" selected>Area 1
 <option value="2">Area 2
 <option value="3">Area 3
</select>

<tr>
<td>Flip axes?
<td><input type=checkbox id=flip>
<tr>
<td>Image Width
<td><input type=text name=zoom_width id=zoom_width size=6 value="600">
</table>
};
	    $html .= "<br><span class=small><a href=data/$out.png target=_new>Image File</a>";
	    $html .= "<br><span class=small><a href=$tmp target=_new>Syntolog file with GEvo links</a><br>";
	  }
      }
    else
      {
	$problem=1
      }
    foreach my $org_dir (keys %org_dirs)
      {
	my $output = $org_dirs{$org_dir}{blastfile};
	next unless -s $output;
	$output =~ s/$DATADIR/$URL\/data/;
	$html .= "<a href=$output target=_new>Blast results for $org_dir</a><br>";;	
      }
    my $log = $cogeweb->logfile;
    $log =~ s/$DIR/$URL/;
    $html .= "<a href=$log target=_new>log</a></span><br>";
    $html .= "<td valign=top>";
    $html .= "<div id=syn_loc1></div>";
    $html .= "<div id=syn_loc2></div>";
    $html .= "<div id=syn_loc3></div>";
    $html .= "</table>";
    if ($problem)
      {
	$html .= qq{
<span class=alert>There was a problem running your analysis.  Please check the log file for details.</span><br>
	  };
      }
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
    ($oid1, $masked1, $md51,$org_name1,$seq_type1, $oid2, $masked2, $md52,$org_name2,$seq_type2) = ($oid2, $masked2, $md52,$org_name2,$seq_type2,$oid1, $masked1, $md51,$org_name1,$seq_type1) if ($org_name2 lt $org_name1);

    my $tmp1 = $org_name1;
    my $tmp2 = $org_name2;
    foreach my $tmp ($tmp1, $tmp2)
      {
	$tmp =~ s/\///g;
	$tmp =~ s/\s+/_/g;
	$tmp =~ s/\(//g;
	$tmp =~ s/\)//g;
	$tmp =~ s/://g;
      }

    my $dir = $tmp1."/".$tmp2;
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
#	    ($mask1, $type1, $mask2, $type2) = ($mask2, $type2, $mask1, $type1)  if ($md52 lt $md51); #not sure about this comment for now
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
	    push @anno, $fid."<br>genomic";
	    next;
	  }
	my $feat = $coge->resultset('Feature')->find($fid);
	my $anno = "Name: ".join (", ",map {"<a class=\"data link\" href=\"/CoGe/FeatView.pl?accn=".$_."\" target=_new>".$_."</a>"} $feat->names);
	my $location = "Chr ".$feat->chromosome." ";
	$location .= commify($feat->start)." - ".commify($feat->stop);
#	$location .=" (".$feat->strand.")";
	push @anno, $location."<br>".$anno;
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

Thank you for using the CoGe Software Package.
	
- The CoGe Team
};
	
	print $mailer $body;
	$mailer->close();
}

sub get_dotplot
  {
    my %args = @_;
    my $src = $args{src};
    my $loc = $args{loc};
    my $flip = $args{flip} eq "true" ? 1 : 0;
    my $regen = $args{regen_images} eq "true" ? 1 : 0;
    my $width = $args{width};
    $src .= ";flip=$flip" if $flip;
    $src .= ";regen=$regen" if $regen;
    $src .= ";width=$width" if $width;
    my $content = get("http://".$ENV{SERVER_NAME}."/".$src);
    my ($url) = $content =~ /url=(.*?)"/is;
    my $png = $url;
    $png =~ s/html$/png/;
    $png =~ s/$URL\/data/$DATADIR/;
    my $img = GD::Image->new($png);
    my ($w,$h) = $img->getBounds();
    $w+=20;
    $h+=60;
    if ($loc)
      {
	return ($url,$loc, $w, $h);
      }
    my $html = qq{<iframe src=$url frameborder=0 width=$w height=$h scrolling=no></iframe>};
    return $html;
  }

sub commify
  {
    my $text = reverse $_[0];
    $text =~ s/(\d\d\d)(?=\d)(?!\d*\.)/$1,/g;
    return scalar reverse $text;
  }
