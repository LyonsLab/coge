#! /usr/bin/perl -w
use strict;
use CGI;
use CGI::Carp 'fatalsToBrowser';
use CGI::Ajax;
use CoGe::Accessory::LogUser;
use CoGe::Accessory::Web;
use CoGe::Algos::KsCalc;
use CoGeX;
use DBIxProfiler;
use Data::Dumper;
use HTML::Template;
use Parallel::ForkManager;
use GD;
use File::Path;
use Mail::Mailer;
use Benchmark;
use LWP::Simple;
use DBI;

$ENV{PATH} = "/opt/apache2/CoGe/";
umask(0);
use vars qw( $DATE $DEBUG $DIR $URL $USER $FORM $coge $cogeweb $FORMATDB $BLAST $DATADIR $FASTADIR $BLASTDBDIR $DIAGSDIR $MAX_PROC $DAG_TOOL $PYTHON $TANDEM_FINDER $FILTER_REPETITIVE_MATCHES $RUN_DAGCHAINER $FIND_NEARBY $CONVERT_TO_GENE_ORDER $DOTPLOT);
$DEBUG = 0;
$DIR = "/opt/apache/CoGe/";
$URL = "/CoGe/";
$FORMATDB = "/usr/bin/formatdb";
$BLAST = "nice -20 /usr/bin/blast -a 8 -K 80 -m 8 -e 0.05";
$DATADIR = "$DIR/data/";
$DIAGSDIR = "$DIR/diags";
$FASTADIR = $DATADIR.'/fasta/';
$BLASTDBDIR = $DATADIR.'/blast/db/';
$MAX_PROC=8;
$PYTHON = "/usr/bin/python";
$DAG_TOOL = $DIR."/bin/dagchainer/dag_tools.py";
$TANDEM_FINDER = $DIR."/bin/dagchainer/tandems.py -d 5 -s -r"; #-d option is the distance (in genes) between dups -- not sure if the -s and -r options are needed -- they create dups files based on the input file name
$FILTER_REPETITIVE_MATCHES = $DIR."/bin/dagchainer/filter_repetitive_matches.pl 200000"; #filters multiple matches to the same accession within the "window size", here set to 80000 nt
$RUN_DAGCHAINER = $DIR."/bin/dagchainer/DAGCHAINER/run_DAG_chainer.pl -E 0.05 -s";
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
		       get_dataset_group_info => \&get_dataset_group_info,
		       get_previous_analyses=>\&get_previous_analyses,
		       get_pair_info=> \&get_pair_info,
		       go=>\&go,
		       check_address_validity=>\&check_address_validity,
		       generate_basefile=>\&generate_basefile,
		       get_dotplot=>\&get_dotplot,
		       gen_dsg_menu=>\&gen_dsg_menu,
		       %ajax,
		      );
print $pj->build_html($FORM, \&gen_html);
#print "Content-Type: text/html\n\n";print gen_html($FORM);


sub gen_html
  {
    my $html;
    my ($body) = gen_body();
    my $template = HTML::Template->new(filename=>'/opt/apache/CoGe/tmpl/generic_page.tmpl');
    $template->param(PAGE_TITLE=>'SynMap');
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
    $template->param(MAIN=>1);
    
    my $master_width = $FORM->param('w') || 0;
    $template->param(MWIDTH=>$master_width);
    if ($FORM->param('b') && $FORM->param('b') == 2)
      {
	$template->param(TBLASTX_SELECT=>"selected");
      }
    else
      {
	$template->param(BLASTN_SELECT=>"selected");
      }
    my ($D, $g, $A, $dt) = ($FORM->param('D'),$FORM->param('g'),$FORM->param('A'), $FORM->param('dt'));
    my $display_dagchainer_settings;
    if ($D && $g && $A && $dt) 
      {
	my $type;
	if ($dt =~ /gene/i)
	  {
	    $type = " genes";
	    $template->param('DAG_GENE_SELECT'=>'checked');
	  }
	else
	  {
	    $type = " bp";
	    $template->param('DAG_DISTANCE_SELECT'=>'checked');
	  }
	$display_dagchainer_settings = qq{display_dagchainer_settings([$g,$D,$A],'$type');};
      }
    else
      {
	$template->param('DAG_GENE_SELECT'=>'checked');
	$display_dagchainer_settings = qq{display_dagchainer_settings();};
      }
    $template->param('DISPLAY_DAGCHAINER_SETTINGS'=>$display_dagchainer_settings);
    
    #populate organism menus
    for (my $i=1; $i<=2;$i++)
      {
	my $dsgid = $form->param('dsgid'.$i) || 0;
	my $feattype_param = $FORM->param('ft'.$i) if $FORM->param('ft'.$i);
	my $org_menu = gen_org_menu(dsgid=>$dsgid, num=>$i, feattype_param=>$feattype_param);
	$template->param("ORG_MENU".$i=>$org_menu);
      }


    my $file = $form->param('file');
    if($file)
    {
    	my $results = read_file($file);
    	$template->param(RESULTS=>$results);
    }
    return $template->output;
  }
  

sub gen_org_menu
  {
    my %opts = @_;
    my $oid = $opts{oid};
    my $num = $opts{num};
    my $name = $opts{name};
    my $desc = $opts{desc};
    my $dsgid = $opts{dsgid};
    my $feattype_param = $opts{feattype_param};
    $feattype_param = 1 unless $feattype_param;
    my $org;
    if ($dsgid)
      {
	$org = $coge->resultset('DatasetGroup')->find($dsgid)->organism;
	$oid = $org->id;
      }
    if ($USER->user_name =~ /public/i && $org && $org->restricted)
      {
	 $oid = undef;
	 $dsgid = undef;
      }
    $name = "Search" unless $name;
    $desc = "Search" unless $desc;
    my $menu_template = HTML::Template->new(filename=>'/opt/apache/CoGe/tmpl/SynMap.tmpl');
    $menu_template->param(ORG_MENU=>1);
    $menu_template->param(NUM=>$num);
    $menu_template->param('ORG_NAME'=>$name);
    $menu_template->param('ORG_DESC'=>$desc);
    $menu_template->param('ORG_LIST'=>get_orgs(name=>$name,i=>$num, oid=>$oid));
    my ($dsg_menu) = gen_dsg_menu(oid=>$oid, dsgid=>$dsgid, num=>$num);
    $menu_template->param(DSG_MENU=>$dsg_menu);
    if ($dsgid)
      {
	my ($dsg_info, $feattype_menu, $message) = get_dataset_group_info(dsgid=> $dsgid, org_num=>$num, feattype=>$feattype_param);
	$menu_template->param(DSG_INFO=>$dsg_info);
	$menu_template->param(FEATTYPE_MENU=>$feattype_menu);
	$menu_template->param(GENOME_MESSAGE=>$message);
      }
    return $menu_template->output;
  }


sub gen_dsg_menu
  {
    my $t1 = new Benchmark;
    my %opts = @_;
    my $oid = $opts{oid};
    my $num = $opts{num};
    my $dsgid = $opts{dsgid};
    my @dsg_menu;
    my $message;
    foreach my $dsg (sort {$b->version <=> $a->version || $a->type->id <=> $b->type->id} $coge->resultset('DatasetGroup')->search({organism_id=>$oid},{prefetch=>['genomic_sequence_type']}))
      {
#	next if $USER->user_name =~ /public/i && $dsg->organism->restricted;
	push @dsg_menu, [$dsg->id, $dsg->type->name." (v".$dsg->version.")"];
      }

    my $dsg_menu = qq{
   <select id=dsgid$num onChange="\$('#dsg_info$num').html('<div class=dna_small class=loading class=small>loading. . .</div>'); get_dataset_group_info(['args__dsgid','dsgid$num','args__org_num','args__$num'],['dsg_info$num', 'feattype_menu$num','genome_message$num'])">
};
    foreach (@dsg_menu)
      {
	my ($numt, $name) = @$_;
	my $selected = " selected" if $dsgid && $numt == $dsgid;
	$selected = " " unless $selected;
	$dsg_menu .= qq{
   <OPTION VALUE=$numt $selected>$name</option>
};
      }
    $dsg_menu .= "</select>";
    my $t2 = new Benchmark;
    my $time = timestr(timediff($t2,$t1));
#    print STDERR qq{
#-----------------
#sub gen_dsg_menu runtime:  $time
#-----------------
#};
    return ($dsg_menu, $message);
    
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
    my $oid = $opts{oid};
    my $i = $opts{i};
    my @db;
    #get rid of trailing white-space
    $name =~ s/^\s+//g if $name;
    $name =~ s/\s+$//g if $name;
    $desc =~ s/^\s+//g if $desc;
    $desc =~ s/\s+$//g if $desc;

    $name = "" if $name =~ /Search/; #need to clear to get full org count
    if ($oid)
      {
	my $org = $coge->resultset("Organism")->find($oid);
	$name = $org->name if $org;
	push @db, $org if $name;
      }
    elsif ($name) 
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
	$option .= " selected" if $oid && $oid == $item->id;
	$option .= ">".$item->name." (id".$item->id.")</OPTION>";
	push @opts, $option;

      }
    my $html;
    $html .= qq{<FONT CLASS ="small">Organism count: }.scalar @opts.qq{</FONT>\n<BR>\n};
    unless (@opts && ($name || $desc)) 
      {
	$html .=  qq{<input type = hidden name="org_id$i" id="org_id$i">};
	return $html;
      }

    $html .= qq{<SELECT id="org_id$i" SIZE="5" MULTIPLE onChange="get_dataset_group_info_chain($i)" >\n};
    $html .= join ("\n", @opts);
    $html .= "\n</SELECT>\n";
    $html =~ s/OPTION/OPTION SELECTED/ unless $oid;
    return $html;
  }

sub get_dataset_group_info
  {
    my $t1 = new Benchmark;
    my %opts = @_;
    my $dsgid = $opts{dsgid};
    my $org_num = $opts{org_num};
    my $feattype = $opts{feattype};
    $feattype = 1 unless defined $feattype;
    return " "," "," " unless $dsgid;
    my $html_dsg_info; 
    my ($dsg) = $coge->resultset("DatasetGroup")->find({dataset_group_id=>$dsgid},{join=>['organism','genomic_sequences'],prefetch=>['organism','genomic_sequences']});
    return " "," "," " unless $dsg;
    my $org = $dsg->organism;
#    next if $USER->user_name =~ /public/i && $org->restricted;
    my $orgname = $org->name;
    $orgname = "<a href=\"OrganismView.pl?oid=".$org->id."\" target=_new>$orgname</a>: ".$org->description if $org->description;
    $html_dsg_info .= qq{<div><span class="oblique">Organism:</span><span class="small">$orgname</small></div>};
    $html_dsg_info .= qq{<div><span class="oblique link" onclick=window.open('OrganismView.pl?dsgid=$dsgid')>Genome Information: </span><br>};
    my $i =0;
    my $chr_length=0;
    my $chr_count =0;
    my $plasmid =0;
    my $contig = 0;
    my $scaffold = 0;
    my @gs = $dsg->genomic_sequences;

    foreach my $gs (@gs)
      {
	$chr_length += $gs->sequence_length;
	$chr_count++;
	$plasmid =1 if $gs->chromosome =~ /plasmid/i;
	$contig =1 if $gs->chromosome =~ /contig/i;
	$scaffold =1 if $gs->chromosome =~ /scaffold/i;
      }
    $html_dsg_info .= "<span class=small>";
    my ($ds) = $dsg->datasets;
    my $link = $ds->data_source->link;
    $link = "http://".$link unless $link =~ /^http/;
    $html_dsg_info .= "Source:  <a href=".$link." target=_new>".$ds->data_source->name."</a><br>";
    #$html_dsg_info .= $dsg->chr_info(summary=>1);
    $html_dsg_info .= "Chromosome count: $chr_count<br>";
    $html_dsg_info .= "Total length: ".commify($chr_length);
    $html_dsg_info .= "<br>Contains plasmid" if $plasmid;
    $html_dsg_info .= "<br>Contains contigs" if $contig;
    $html_dsg_info .= "<br>Contains scaffolds" if $scaffold;
    $html_dsg_info .= "</span>";
    my $t2 = new Benchmark;
    my $time = timestr(timediff($t2,$t1));
#    print STDERR qq{
#-----------------
#sub get_dataset_group_info runtime:  $time
#-----------------
#};

    my $message;
    
    #create feature type menu
    my $has_cds;
    foreach my $ft ($coge->resultset('FeatureType')->search(
							    {
							     dataset_group_id=>$dsg->id,
							     'me.feature_type_id'=>3},
							    {
							     join =>{features=>{dataset=>'dataset_connectors'}},
							     rows=>1,
							    }
							   )
		   )
      {
	$has_cds = 1;
      }
    my ($cds_selected, $genomic_selected) = (" ", " ");
    $cds_selected = "selected" if $feattype eq 1 || $feattype eq "CDS";
    $genomic_selected = "selected" if $feattype eq 2 || $feattype eq "genomic";

    my $feattype_menu = qq{
  <select id="feat_type$org_num" name ="feat_type$org_num">
#};
    $feattype_menu .= qq{
   <OPTION VALUE=1 $cds_selected>CDS</option>
} if $has_cds;
    $feattype_menu .= qq{
   <OPTION VALUE=2 $genomic_selected>genomic</option>
};
    $feattype_menu .= "</select>";
    $message = "<span class='small alert'>No Coding Sequence in Genome</span>" unless $has_cds;

    return $html_dsg_info, $feattype_menu, $message;
  }

sub gen_fasta
  {
    my %opts = @_;
    my $dsgid = $opts{dsgid};
    my $feat_type = $opts{feat_type};
    my $write_log = $opts{write_log} || 0;
    my ($org_name, $title);
    ($org_name, $title) = gen_org_name(dsgid=>$dsgid, feat_type=>$feat_type, write_log=>$write_log);
    my $file = $FASTADIR."/$dsgid-$feat_type.fasta";
    my $res;
    while (-e "$file.running")
      {
	print STDERR "detecting $file.running.  Waiting. . .\n";
	sleep 60;
      }
    if (-r $file)
      {
	write_log("fasta file for *".$org_name."* ($file) exists", $cogeweb->logfile);
	$res = 1;
      }
    else
      {
	system "touch $file.running"; #track that a blast anlaysis is running for this
	$res = generate_fasta(dsgid=>$dsgid, file=>$file, type=>$feat_type) unless -r $file;
	system "rm $file.running" if -r "$file.running"; #remove track file
      }
    return $file, $org_name, $title if $res;
    return 0;
  }
    
sub gen_org_name
  {
    my %opts = @_;
    my $dsgid = $opts{dsgid};
    my $feat_type = $opts{feat_type} || 1;
    my $write_log = $opts{write_log} || 0;
    my ($dsg) = $coge->resultset('DatasetGroup')->search({dataset_group_id=>$dsgid}, {join=>'organism',prefetch=>'organism'});
    my $org_name = $dsg->organism->name;
    my $title = $org_name ." (v".$dsg->version.", dsgid".$dsgid.")".$feat_type;
    $title =~ s/(`|')//g;
    write_log("ORGANISM: ".$title, $cogeweb->logfile) if $write_log;
    return ($org_name, $title);
  }

sub generate_fasta
  {
    my %opts = @_;
    my $dsgid = $opts{dsgid};
    my $file = $opts{file};
    my $type = $opts{type};
    my ($dsg) = $coge->resultset('DatasetGroup')->search({"me.dataset_group_id"=>$dsgid},{join=>'genomic_sequences',prefetch=>'genomic_sequences'});

    $file = $FASTADIR."/$file" unless $file =~ /$FASTADIR/;
    write_log("creating fasta file.", $cogeweb->logfile);
    open (OUT, ">$file") || die "Can't open $file for writing: $!";;
    if ($type eq "CDS")
      {
	foreach my $feat (sort {$a->chromosome cmp $b->chromosome || $a->start <=> $b->start} 
			  $coge->resultset('Feature')->search(
							      {
							       feature_type_id=>3, 
							       dataset_group_id=>$dsgid
							      },{
								 join=>[{dataset=>'dataset_connectors'}], 
								 prefetch=>['feature_names']}
							     ))
	  {
	    my ($chr) = $feat->chromosome;#=~/(\d+)/;
	    my $name;
	    foreach my $n ($feat->names)
	      {
		$name = $n;
		last unless $name =~ /\s/;
	      }
	    $name =~ s/\s+/_/g;
	    my $title = join ("||",$chr, $feat->start, $feat->stop, $name, $feat->strand, $type, $feat->id);
	    my $seq = $feat->genomic_sequence(dsgid=>$dsg);
	    next unless $seq;
	    print OUT ">".$title."\n";
	    print OUT $seq,"\n";
	  }
      }
    else
      {
	foreach my $chr (sort $dsg->get_chromosomes)
	  {
	    #		my $title = join ("||",$chr, 1, $ds->last_chromosome_position($chr), "Chr_$chr",1, "genomic", "N/A");
	    my $seq = $dsg->get_genomic_sequence(chr=>$chr);
	    next unless $seq;
	    print OUT ">".$chr."\n";
	    print OUT $seq,"\n";
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
    my $dbname = $opts{dbname};
    my $fasta = $opts{fasta};
    my $org_name = $opts{org_name};
    my $blastdb = "$BLASTDBDIR/$dbname";
    my $res = 0;
    while (-e "$blastdb.running")
      {
	
	print STDERR "detecting $blastdb.running.  Waiting. . .\n";
	sleep 60;
      }
    if (-r $blastdb.".nsq")
      {
	write_log("blastdb file for *".$org_name."* ($dbname) exists", $cogeweb->logfile);
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
	print STDERR "detecting $outfile.running.  Waiting. . .\n";
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
      my $feat_type1 = $opts{feat_type1};
      my $feat_type2 = $opts{feat_type2};
      while (-e "$outfile.running")
      {
	print STDERR "detecting $outfile.running.  Waiting. . .\n";
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
      $cmd .= " -c";# if $feat_type1 eq "genomic" && $feat_type2 eq "genomic";
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
	print STDERR "detecting $outfile.running.  Waiting. . .\n";
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
	print STDERR "detecting $outfile.running.  Waiting. . .\n";
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

sub run_convert_to_gene_order
  {
    my %opts = @_;
    my $infile = $opts{infile};
    my $dsgid1 = $opts{dsgid1};
    my $dsgid2 = $opts{dsgid2};
    my $ft1 = $opts{ft1};
    my $ft2 = $opts{ft2};
    my $ftid1 = $ft1 eq "CDS" ? 3 : 0; #set feat type id to 3 if CDS
    my $ftid2 = $ft2 eq "CDS" ? 3 : 0; #set feat type id to 3 if CDS
    my $outfile = $infile."_geneorder";
    while (-e "$outfile.running")
      {
	print STDERR "detecting $outfile.running.  Waiting. . .\n";
	sleep 60;
      }
    unless (-r $infile && -s $infile)
	{
	  write_log("WARNING:   Cannot convert to gene order! DAGChainer input file ($infile) contains no data!" ,$cogeweb->logfile);
	  return 0;
	}
    if (-r $outfile)
      {
	write_log("run_filter_repetitive+_matches: file $outfile already exists",$cogeweb->logfile);
	return $outfile;
      }
    my $cmd = $CONVERT_TO_GENE_ORDER ." -i $infile -dsg1 $dsgid1 -dsg2 $dsgid2 -ft1 $ftid1 -ft2 $ftid2  > $outfile";
    system "touch $outfile.running"; #track that a blast anlaysis is running for this
    write_log("run_convert_to_gene_order: running $cmd", $cogeweb->logfile);
    `$cmd`;
    write_log("Completed conversion!", $cogeweb->logfile);
    system "rm $outfile.running" if -r "$outfile.running";; #remove track filereturn $outfile;
    return $outfile;
  }

sub replace_gene_order_with_genomic_positions
  {
    my %opts = @_;
    my $file = $opts{file};
    #must convert file's coordinates back to genomic
    while (-e "$file.orig.running")
      {
	print STDERR "detecting $file.orig.running.  Waiting. . .\n";
	sleep 60;
      }
    if (-r "$file.orig" && -s "$file.orig")
      {
	write_log("  no cnoversion for $file back to genomic coordinates needed, convered file exists", $cogeweb->logfile);
	return;
      }
    system "touch $file.orig.running"; #track that a blast anlaysis is running for this
    write_log("  converting $file back to genomic coordinates", $cogeweb->logfile);
    `mv $file $file.orig`;
    $/="\n"; #just in case
    open (IN,  "$file.orig");
    open (OUT, ">$file");
    while (<IN>)
      {
	if (/^#/){print OUT $_; next;}
	my @line = split /\t/;
	my @item1 = split /\|\|/, $line[1];
	my @item2 = split /\|\|/, $line[5];
	my ($start, $stop) =($item1[1], $item1[2]);
	($start, $stop) = ($stop, $start) if $item1[4] && $item1[4]=~/-/;
	$line[2] = $start;
	$line[3] = $stop;
	($start, $stop) =($item2[1], $item2[2]);
	($start, $stop) = ($stop, $start) if $item2[4] && $item2[4]=~/-/;
	$line[6] = $start;
	$line[7] = $stop;
	print OUT join ("\t", @line);
      }
    close IN;
    close OUT;
    system "rm $file.orig.running" if -r "$file.orig.running"; 
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
	print STDERR "detecting $outfile.running.  Waiting. . .\n";
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
	print STDERR "detecting $outfile.running.  Waiting. . .\n";
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

sub gen_ks_db
  {
    my %opts = @_;
    my $infile = $opts{infile};
    my ($outfile) = $infile =~ /^(.*?CDS-CDS)/;
    return unless $outfile;
    $outfile .= ".sqlite";
    unless (-r $outfile)
      {
	my $create = qq{
CREATE TABLE ks_data
(
id INTEGER PRIMARY KEY,
fid1 integer,
fid2 integer,
dS varchar,
dN varchar,
dN_dS varchar
)
};

	my $dbh = DBI->connect("dbi:SQLite:dbname=$outfile","","");
	$dbh->do($create) if $create;
	$dbh->do('create INDEX fid1 ON ks_data (fid1)');
	$dbh->do('create INDEX fid2 ON ks_data (fid2)');
	$dbh->disconnect;
      }
    write_log("generating ks database", $cogeweb->logfile);
    write_log("\tconnecting to ks database $outfile", $cogeweb->logfile);
    my %ksdata;
    my $select = "select * from ks_data";
    my $dbh = DBI->connect("dbi:SQLite:dbname=$outfile","","");
    my $sth = $dbh->prepare($select);
    $sth->execute();
    write_log("\texecuting select all from ks database $outfile", $cogeweb->logfile);
    while (my $data = $sth->fetchrow_arrayref)
      {
	$ksdata{$data->[1]}{$data->[2]}=1;# unless $data->[3] eq "";
	print STDERR $data->[1],"\t", $data->[2]."\n" if $data->[3] eq "";
      }
    write_log("\tgathered data from ks database $outfile", $cogeweb->logfile);
    $dbh->disconnect();
    write_log("\tdisconnecting from ks database $outfile", $cogeweb->logfile);
    open (IN, $infile);
    my @data;
    while (<IN>)
      {
	
	next if /^#/;
	chomp;
	my @line = split /\t/;
	my @item1 = split /\|\|/,$line[1];
	my @item2 = split /\|\|/,$line[5];
	unless ($item1[6] && $item2[6])
	  {
	    warn "Line does not appear to contain coge feature ids:  $_\n";
	    next;
	  }
	next if $ksdata{$item1[6]}{$item2[6]};
	push @data,[$line[1],$line[5],$item1[6],$item2[6]];
      }
    close IN;
    print STDERR "generating synonymous substitution values for ".scalar @data." pairs of genes\n";
    my $pm = new Parallel::ForkManager($MAX_PROC);
    foreach my $item (@data)
      {	
	$pm->start and next;

	my ($fid1) = $item->[2] =~ /(^\d+$)/;
	my ($fid2) = $item->[3] =~ /(^\d+$)/;
	my ($feat1) = $coge->resultset('Feature')->find($fid1);
	my ($feat2) = $coge->resultset('Feature')->find($fid2);
	my $ks = new CoGe::Algos::KsCalc();
	$ks->feat1($feat1);
	$ks->feat2($feat2);
	my $res = $ks->KsCalc();
	unless ($res)
	  {
	    print STDERR "Failed KS calculation: $fid1\t$fid2\n";
	    $res = {};
	  }
	my ($dS, $dN, $dNS);
	if (keys %$res)
	  {
	    $dS = $res->{dS};
	    $dN = $res->{dN};
	    $dNS = $res->{'dN/dS'};
	  }
	my $insert = qq{
INSERT INTO ks_data (fid1, fid2, dS, dN, dN_dS) values ($fid1, $fid2, "$dS", "$dN", "$dNS")
};
	my $dbh = DBI->connect("dbi:SQLite:dbname=$outfile","","");
	my $insert_success = 0;
	while (!$insert_success)
	  {
	    $insert_success = $dbh->do($insert);
	    unless ($insert_success)
	      {
		print STDERR $insert;
		sleep .1;
	      }
	  }

	$dbh->disconnect();

	$pm->finish;
      }
    $pm->wait_all_children();

    system "rm $outfile.running" if -r "$outfile.running";; #remove track file
    return $outfile;
  }

sub add_GEvo_links
  {
    my %opts = @_;
    my $infile = $opts{infile};
    my $dsgid1 = $opts{dsgid1};
    my $dsgid2 = $opts{dsgid2};
    open (IN, $infile);
    open (OUT,">$infile.tmp");
    my %condensed;
    my %names;
    while (<IN>)
      {
	my $prev_link =0;
	chomp;
	if (/^#/)
	  {
	    print OUT $_,"\n";
	    next;
	  }
	if (/GEvo/)
	  {
	    s/toxic/synteny.cnr/;
#	    print OUT $_,"\n";
	    $prev_link =1;
	  }
	s/^\s+//;
	next unless $_;
	my @line = split/\t/;
	my @feat1 = split/\|\|/,$line[1];
	my @feat2 = split/\|\|/,$line[5];
	my $link = "http://synteny.cnr.berkeley.edu/CoGe/GEvo.pl?";
	my ($fid1, $fid2);
	if ($feat1[6])
	  {
	    $fid1 = $feat1[6];
	    $link .= "fid1=".$fid1;
	  }
	else
	  {
	    my ($xmin) = sort ($feat1[1], $feat1[2]);
	    my $x = sprintf("%.0f", $xmin+abs($feat1[1]-$feat1[2])/2);
	    $link .= "chr1=".$feat1[0].";x1=".$x;
	  }
	if ($feat2[6])
	  {
	    $fid2 = $feat2[6];
	    $link .= ";fid2=".$fid2;
	  }
	else
	  {
	    my ($xmin) = sort ($feat2[1], $feat2[2]);
	    my $x = sprintf("%.0f", $xmin+abs($feat2[1]-$feat2[2])/2);
	    $link .= ";chr2=".$feat2[0].";x2=".$x;
	  }
	$link .= ";dsgid1=".$dsgid1;
	$link .= ";dsgid2=".$dsgid2;
	
	if ($fid1 && $fid2)
	  {
	    $condensed{$fid1."_".$dsgid1}{$fid2."_".$dsgid2} = 1;
	    $condensed{$fid2."_".$dsgid2}{$fid1."_".$dsgid1} = 1;
	    $names{$fid1} = $feat1[3];
	    $names{$fid2} = $feat2[3];
	  }
#	accn1=".$feat1[3]."&fid1=".$feat1[6]."&accn2=".$feat2[3]."&fid2=".$feat2[6] if $feat1[3] && $feat1[6] && $feat2[3] && $feat2[6];
	print OUT $_;
	print OUT "\t",$link if $link && !$prev_link;
	print OUT "\n";
      }
    close IN;
    close OUT;
    open (OUT,">$infile.condensed");
    foreach my $id1 (sort keys %condensed)
      {
	my ($fid1, $dsgid1) = split /_/, $id1;
	my @names = $names{$fid1};
	my $link = "http://synteny.cnr.berkeley.edu/CoGe/GEvo.pl?pad_gs=10000;fid1=$fid1;dsgid1=$dsgid1";
	my $count =2;
	foreach my $id2 (sort keys %{$condensed{$id1}})
	  {
	    my ($fid2, $dsgid2) = split/_/, $id2;
	    $link .= ";fid$count=$fid2;dsgid$count=$dsgid2";
	    push @names, $names{$fid2};
	    $count++;
	  }
	$count--;
	$link .= ";num_seqs=$count";
	print OUT join ("\t", $link, "<a href=$link;autogo=1>AutoGo</a>",@names), "\n";
      }
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
    my $dsgid1 = $opts{dsgid1};
    my $dsgid2 = $opts{dsgid2};
    my $dagtype = $opts{dagtype};
    my $ks_db = $opts{ks_db};
    my $ks_type = $opts{ks_type};
    my ($basename) = $coords =~ /([^\/]*).all.aligncoords/;
    my $regen_images = $opts{regen_images}=~/true/i ? 1 : 0;
    my $width = $opts{width} || 1000;
    my $cmd = $DOTPLOT;
    #add ks_db to dotplot command if requested
    if ($ks_db && -r $ks_db)
      {
	$cmd .= qq{ -ksdb $ks_db -kst $ks_type -log 1};
	$outfile .= ".$ks_type";
      }
    $cmd .= qq{ -d $dag -a $coords -b $outfile -l 'javascript:synteny_zoom("$dsgid1","$dsgid2","$basename","XCHR","YCHR","$ks_db")' -dsg1 $dsgid1 -dsg2 $dsgid2 -w $width -lt 2};


    while (-e "$outfile.running")
      {
	print STDERR "detecting $outfile.running.  Waiting. . .\n";
	sleep 60;
      }
    if (-r "$outfile.png" && !$regen_images)
      {
	write_log("generate dotplot: file $outfile already exists",$cogeweb->logfile);
	return $outfile;
      }
    system "touch $outfile.running"; #track that a blast anlaysis is running for this
    write_log("generate dotplot: running $cmd", $cogeweb->logfile);
    system "rm $outfile.running" if -r "$outfile.running";; #remove track file
    `$cmd`;
    return $outfile if -r "$outfile.html";
  }

sub go
  {
    my %opts = @_;
    my $dagchainer_D= $opts{D};
    my $dagchainer_g = $opts{g};
    my $dagchainer_A = $opts{A};
    my $regen_images = $opts{regen_images};
    my $email = $opts{email};
    my $job_title = $opts{jobtitle};
    my $width = $opts{width};
    my $basename = $opts{basename};
    my $blast = $opts{blast};
    my $feat_type1 = $opts{feat_type1};
    my $feat_type2 = $opts{feat_type2};
    my $dsgid1 = $opts{dsgid1};
    my $dsgid2 = $opts{dsgid2};
    my $ks_type = $opts{ks_type};

    my $dagchainer_type = $opts{dagchainer_type};
    $dagchainer_type = $dagchainer_type eq "true" ? "geneorder" : "distance";

    unless ($dsgid1 && $dsgid2)
      {
	return "<span class=alert>You must select two genomes.</span>"
      }
    my ($dsg1) = $coge->resultset('DatasetGroup')->find($dsgid1);
    my ($dsg2) = $coge->resultset('DatasetGroup')->find($dsgid2);
    unless ($dsg1 && $dsg2)
      {
	return "<span class=alert>Problem generating dataset group objects for ids:  $dsgid1, $dsgid2.</span>";
      }
    $cogeweb = initialize_basefile(basename=>$basename, prog=>"SynMap");
    my $synmap_link = "SynMap.pl?dsgid1=$dsgid1;dsgid2=$dsgid2;D=$dagchainer_D;g=$dagchainer_g;A=$dagchainer_A;w=$width;b=$blast;ft1=$feat_type1;ft2=$feat_type2";
    $email = 0 if check_address_validity($email) eq 'invalid';
    $blast = $blast == 2 ? "tblastx" : "blastn";
    $feat_type1 = $feat_type1 == 2 ? "genomic" : "CDS";
    $feat_type2 = $feat_type2 == 2 ? "genomic" : "CDS";

    $synmap_link .=";dt=$dagchainer_type"; 

    ##generate fasta files and blastdbs
    my $t0 = new Benchmark;
    my $pm = new Parallel::ForkManager($MAX_PROC);
    my @dsgs = ([$dsgid1, $feat_type1]);
    push @dsgs, [$dsgid2, $feat_type2] unless $dsgid1 == $dsgid2 && $feat_type1 eq $feat_type2;
    foreach my $item (@dsgs)
      {
	$pm->start and next;
	my $dsgid = $item->[0];

	my $feat_type = $item->[1];

	my ($fasta,$org_name) = gen_fasta(dsgid=>$dsgid, feat_type=>$feat_type, write_log=>1);

	gen_blastdb(dbname=>"$dsgid-$feat_type",fasta=>$fasta,org_name=>$org_name);
	$pm->finish;
      }
    $pm->wait_all_children();
    my ($fasta1,$org_name1, $title1);
    my ($fasta2,$org_name2, $title2);
    ($fasta1,$org_name1, $title1) = gen_fasta(dsgid=>$dsgid1, feat_type=>$feat_type1);
    ($fasta2,$org_name2, $title2) = gen_fasta(dsgid=>$dsgid2, feat_type=>$feat_type2);
    ($dsgid1, $org_name1,$fasta1,$feat_type1, $dsgid2,
    $org_name2,$fasta2, $feat_type2) = ($dsgid2,
    $org_name2,$fasta2, $feat_type2, $dsgid1,
    $org_name1,$fasta1, $feat_type1) if ($org_name2 lt $org_name1);
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
    
    my ($blastdb1) = gen_blastdb(dbname=>"$dsgid1-$feat_type1", fasta=>$fasta1,org_name=>$org_name1);
    my ($blastdb2) = gen_blastdb(dbname=>"$dsgid2-$feat_type2", fasta=>$fasta2,org_name=>$org_name2);
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
 	$tmp =~ s/;//g;
       }
     my $orgkey1 = $title1;
     my $orgkey2 = $title2;
     my %org_dirs = (
 		    $orgkey1."_".$orgkey2=>{fasta=>$fasta1,
 					    db=>$blastdb2,
 					    basename=>$dsgid1."_".$dsgid2.".$feat_type1-$feat_type2.$blast",
 					    dir=>$DIAGSDIR."/".$tmp1."/".$tmp2,
 					   },
#  		    $orgkey1."_".$orgkey1=>{fasta=>$fasta1,
#  					    db=>$blastdb1,
#  					    basename=>$dsgid1."_".$dsgid1.".$feat_type1-$feat_type1.$blast",
#  					    dir=>$DIAGSDIR."/".$tmp1."/".$tmp1,
#  					   },
#  		    $orgkey2."_".$orgkey2=>{fasta=>$fasta2,
#  					    db=>$blastdb2,
#  					    basename=>$dsgid2."_".$dsgid2.".$feat_type2-$feat_type2.$blast",
#  					    dir=>$DIAGSDIR."/".$tmp2."/".$tmp2,
#  					   },
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
    my $t1 = new Benchmark;
    my $blast_time = timestr(timediff($t1,$t0));


     #Find local dups
#     my $dag_file11 = $org_dirs{$orgkey1."_".$orgkey1}{dir}."/".$org_dirs{$orgkey1."_".$orgkey1}{basename}.".dag";
#     $problem=1 unless (run_dag_tools(query=>"a".$dsgid1, subject=>"b".$dsgid1, blast=>$org_dirs{$orgkey1."_".$orgkey1}{blastfile}, outfile=>$dag_file11, feat_type1=>$feat_type1, feat_type1=>$feat_type1));
#     my $dag_file22 = $org_dirs{$orgkey2."_".$orgkey2}{dir}."/".$org_dirs{$orgkey2."_".$orgkey2}{basename}.".dag";
#     $problem=1 unless run_dag_tools(query=>"a".$dsgid2, subject=>"b".$dsgid2, blast=>$org_dirs{$orgkey2."_".$orgkey2}{blastfile}, outfile=>$dag_file22, feat_type1=>$feat_type2, feat_type1=>$feat_type2);
#     my $dup_file1  = $org_dirs{$orgkey1."_".$orgkey1}{dir}."/".$org_dirs{$orgkey1."_".$orgkey1}{basename}.".dups";
#     run_tandem_finder(infile=>$dag_file11,outfile=>$dup_file1);
#     my $dup_file2  = $org_dirs{$orgkey2."_".$orgkey2}{dir}."/".$org_dirs{$orgkey2."_".$orgkey2}{basename}.".dups";
#     run_tandem_finder(infile=>$dag_file22,outfile=>$dup_file2);
    my $t2 = new Benchmark;
    my $local_dup_time = timestr(timediff($t2,$t1));

     #prepare dag for synteny analysis
     my $dag_file12 = $org_dirs{$orgkey1."_".$orgkey2}{dir}."/".$org_dirs{$orgkey1."_".$orgkey2}{basename}.".dag";
#     $problem=1 unless run_dag_tools(query=>"a".$dsgid1, subject=>"b".$dsgid2, blast=>$org_dirs{$orgkey1."_".$orgkey2}{blastfile}, outfile=>$dag_file12.".all", query_dup_file=>$dup_file1,subject_dup_file=>$dup_file2, feat_type1=>$feat_type1, feat_type2=>$feat_type2);
     $problem=1 unless run_dag_tools(query=>"a".$dsgid1, subject=>"b".$dsgid2, blast=>$org_dirs{$orgkey1."_".$orgkey2}{blastfile}, outfile=>$dag_file12.".all", feat_type1=>$feat_type1, feat_type2=>$feat_type2);
     my $dag_file12_all = $dag_file12.".all";
     #remove repetitive matches
     run_filter_repetitive_matches(infile=>$dag_file12_all,outfile=>$dag_file12);
     #is this an ordered gene run?
     $dag_file12 = run_convert_to_gene_order(infile=>$dag_file12, dsgid1=>$dsgid1, dsgid2=>$dsgid2, ft1=>$feat_type1, ft2=>$feat_type2) if $dagchainer_type eq "geneorder";
    my $t3 = new Benchmark;
    my $convert_to_gene_order_time = timestr(timediff($t3,$t2));

     #run dagchainer
     my $dagchainer_file = run_dagchainer(infile=>$dag_file12, D=>$dagchainer_D, g=>$dagchainer_g,A=>$dagchainer_A, type=>$dagchainer_type);
     write_log("Completed dagchainer run", $cogeweb->logfile);
     write_log("", $cogeweb->logfile);
    my $t4 = new Benchmark;
    my $run_dagchainer_time = timestr(timediff($t4,$t3));
    my ($find_nearby_time, $gen_ks_db_time, $dotplot_time, $add_gevo_links_time);
     if (-r $dagchainer_file)
       {
 	my $tmp = $dagchainer_file; #temp file name for the final post-processed data
 	$tmp =~ s/aligncoords/all\.aligncoords/;
 	#convert to genomic coordinates if gene order was used
 	if ($dagchainer_type eq "geneorder")
 	  {
 	    replace_gene_order_with_genomic_positions(file=>$dagchainer_file);
 	  }
 	#add pairs that were skipped by dagchainer
 	run_find_nearby(infile=>$dagchainer_file, dag_all_file=>$dag_file12_all, outfile=>$tmp);
	my $t5 = new Benchmark;
	$find_nearby_time = timestr(timediff($t5,$t4));

 	#generate dotplot images
 	my $org1_length =0;
 	my $org2_length =0;
 	my $chr1_count = 0;
 	my $chr2_count = 0;
	foreach my $gs ($dsg1->genomic_sequences)
	  {
	    $chr1_count++;
	    $org1_length+=$gs->sequence_length;
	  }
	foreach my $gs ($dsg2->genomic_sequences)
	  {
	    $chr2_count++;
	    $org2_length+=$gs->sequence_length;
	  }
 	my $test = $org1_length > $org2_length ? $org1_length : $org2_length;
 	unless ($width)
 	  {
 	    $width = int($test/100000);
 	    $width = 1200 if $width > 1200;
 	    $width = 500 if $width < 500;
 	    $width = 1200 if $chr1_count > 9 || $chr2_count > 9;
 	    $width = 2000 if ($chr1_count > 100 || $chr2_count > 100);
 	  }
 	my $qlead = "a";
 	my $slead = "b";
 	my $out = $org_dirs{$orgkey1."_".$orgkey2}{dir}."/html/";
 	mkpath ($out,0,0777) unless -d $out;
 	$out .="master_".$org_dirs{$orgkey1."_".$orgkey2}{basename};
 	$out .= "_$dagchainer_type";
 	$out .= "_D$dagchainer_D" if $dagchainer_D;
 	$out .= "_g$dagchainer_g" if $dagchainer_g;
 	$out .= "_A$dagchainer_A" if $dagchainer_A;
 	$out .= ".w$width";
 	#deactivation ks calculations due to how slow it is
 	my $ks_db = gen_ks_db(infile=>$tmp) if $ks_type;
	my $t6 = new Benchmark;
	$gen_ks_db_time = timestr(timediff($t6,$t5));

 	$out = generate_dotplot(dag=>$dag_file12_all, coords=>$tmp, outfile=>"$out", regen_images=>$regen_images, dsgid1=>$dsgid1, dsgid2=>$dsgid2, width=>$width, dagtype=>$dagchainer_type, ks_db=>$ks_db, ks_type=>$ks_type);
	my $hist = $out.".hist.png";
	my $t7 = new Benchmark;
	$dotplot_time = timestr(timediff($t7,$t6));

 	add_GEvo_links (infile=>$tmp, dsgid1=>$dsgid1, dsgid2=>$dsgid2);
	my $t8 = new Benchmark;
	$add_gevo_links_time = timestr(timediff($t8,$t7));

 	$tmp =~ s/$DIR/$URL/;
 	if (-r "$out.html")
 	  {
	    $html .= qq{
<div class="ui-widget-content" id="synmap_zoom_box">
 Zoomed SynMap:
 <table class=small>
 <tr>
 <td> Display Location:
 <td><select name=map_loc id=map_loc>
  <option value="window1" selected>New Window
  <option value="1">Area 1
  <option value="2">Area 2
  <option value="3">Area 3
 </select>
 <tr>
 <td>Flip axes?
 <td><input type=checkbox id=flip>
 <tr>
 <td>Image Width
 <td><input class="backbox" type=text name=zoom_width id=zoom_width size=6 value="600">
 <tr>
 <td>kS, kN, kN/kS cutoffs: 
 <td>Min: <input class="backbox" type=text name=zoom_min id=zoom_min size=6 value="">
 <td>Max: <input class="backbox" type=text name=zoom_max id=zoom_max size=6 value="">
 </table>
</div>
 };


 	    open (IN, "$out.html") || warn "problem opening $out.html for reading\n";
#	    print STDERR "$out.html\n";
 	    $html .= "<span class='species small'>y-axis: $org_name2</span><table><tr valign=top><td valign=top>";
 	    $/ = "\n";
 	    while (<IN>)
 	      {
 		next if /<\/?html>/;
 		$html .= $_;
 	      }
 	    close IN;
 	    $out =~ s/$DIR/$URL/;
#	    print STDERR $out,"!!\n";
 	    $html =~ s/master.*\.png/$out.png/;
 	    warn "$out.html did not parse correctly\n" unless $html =~ /map/i;
 	    $html .= qq{
 <br><span class="species small">x-axis: $org_name1</span><br>
};

	    $html .= "<div><img src='$out.hist.png'></div>" if -r $hist;

 	    $html .= "<br><span class=small><a href=$out.png target=_new>Image File</a>";
 	    $html .= "<br><span class=small><a href=$out.hist.png target=_new>Histogram of synonymous substitutions</a>" if -r $hist;
 	    $html .= "<br><span class=small><a href=$tmp target=_new>DAGChainer syntelog file with GEvo links</a><br>";
 	    $html .= "<span class=small><a href=$tmp.condensed target=_new>Condensed syntelog file with GEvo links</a><br>";
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
 	$output =~ s/$DIR/$URL/;
 	$html .= "<a href=$output target=_new>Blast results</a><br>";;	
       }
     my $log = $cogeweb->logfile;
     $log =~ s/$DIR/$URL/;
#     my $tiny = get("http://tinyurl.com/create.php?url=http://".$ENV{SERVER_NAME}."/CoGe/$synmap_link");
#     ($tiny) = $tiny =~ /<b>(http:\/\/tinyurl.com\/\w+)<\/b>/;
     write_log("\nLink: $synmap_link", $cogeweb->logfile);
#     write_log("tinyurl: $tiny", $cogeweb->logfile);
     $html .= "<a href=$log target=_new>log</a><br>";
    $html .= "<a href='$synmap_link' target=_new>SynMap Link</a></span>";

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
    email_results(email=>$email,html=>$html,org1=>$org_name1,org2=>$org_name2, jobtitle=>$job_title, link=>$synmap_link) if $email;
    my $benchmarks = qq{
Benchmarks:
Blast:                    $blast_time
Find Local Dups:          $local_dup_time
Convert Gene Order:       $convert_to_gene_order_time
DAGChainer:               $run_dagchainer_time
find nearby:              $find_nearby_time
kS calculations:          $gen_ks_db_time
Dotplot:                  $dotplot_time
GEvo links:               $add_gevo_links_time
};
    print STDERR $benchmarks;
    write_log($benchmarks, $cogeweb->logfile);
    return $html;
   }

sub get_previous_analyses
  {
    my %opts = @_;
    my $oid1 = $opts{oid1};
    my $oid2 = $opts{oid2};
    return unless $oid1 && $oid2;
    my ($org1) = $coge->resultset('Organism')->find($oid1);
    my ($org2) = $coge->resultset('Organism')->find($oid2);
    return if ($USER->user_name =~ /public/i && ($org1->restricted || $org2->restricted));
    my ($org_name1) = $org1->name;
    my ($org_name2) = $org2->name;
    ($oid1, $org_name1, $oid2, $org_name2) = ($oid2, $org_name2, $oid1, $org_name1) if ($org_name2 lt $org_name1);

    my $tmp1 = $org_name1;
    my $tmp2 = $org_name2;
    foreach my $tmp ($tmp1, $tmp2)
      {
	$tmp =~ s/\///g;
	$tmp =~ s/\s+/_/g;
	$tmp =~ s/\(//g;
	$tmp =~ s/\)//g;
	$tmp =~ s/://g;
	$tmp =~ s/;//g;
      }

    my $dir = $tmp1."/".$tmp2;
    $dir = "$DIAGSDIR/".$dir;
    my $sqlite =0;
    my @items;
    if (-d $dir)
      {
	opendir (DIR, $dir);
	while (my $file = readdir(DIR))
	  {
	    $sqlite = 1 if $file =~ /sqlite$/;
	    next unless $file =~ /all\.aligncoords$/;
	    my ($D, $g, $A) = $file =~ /D(\d+)_g(\d+)_A(\d+)/;
	    next unless ($D && $g && $A);
	    my $blast = $file =~ /blastn/ ? "BlastN" : "TBlastX";
	    my ($dsgid1, $dsgid2, $type1, $type2) = $file =~ /^(\d+)_(\d+)\.(\w+)-(\w+)/;
	    next unless ($dsgid1 && $dsgid2 && $type1 && $type2);
	    my %data = (D=>$D,
			g=>$g,
			A=>$A,
			blast=>$blast,
			dsgid1=>$dsgid1,
			dsgid2=>$dsgid2);
	    my $geneorder = $file =~ /geneorder/;
	    my $dsg1 = $coge->resultset('DatasetGroup')->find($dsgid1);
	    next unless $dsg1;
	    my ($ds1) = $dsg1->datasets;
	    my $dsg2 = $coge->resultset('DatasetGroup')->find($dsgid2);
	    next unless $dsg2;
	    my ($ds2) = $dsg2->datasets;
	    my $name .= $dsg1->type->name." ".$type1." (".$ds1->data_source->name." v".$dsg1->version.") vs ";
	    $name .= $dsg2->type->name." ".$type2." (".$ds2->data_source->name." v".$dsg2->version.")";
	    $data{name} = $name;
	    $type1 = $type1 eq "CDS" ? 1 : 2; 
	    $type2 = $type2 eq "CDS" ? 1 : 2; 
	    $data{type1} = $type1;
	    $data{type2} = $type2;
	    $data{dagtype} = $geneorder ? "Ordered genes" : "Distance";
	    push @items, \%data;
	  }
      }
    return unless @items;
    my $size = scalar @items;
    $size = 8 if $size > 8;
    my $html = qq{
<select id="prev_params" size=4 multiple onChange="update_params();">
};
    foreach (sort {$a->{g}<=>$b->{g} } @items)
      {
	my $blast = $_->{blast} =~ /^blastn$/i ? 0 : 1;
	my $val = join ("_",$_->{g},$_->{D},$_->{A}, $oid1, $_->{dsgid1}, $_->{type1},$oid2, $_->{dsgid2}, $_->{type2}, $blast, $_->{dagtype});
	my $name =$_->{name}.": ".$_->{blast}.", ".$_->{dagtype}.", g".$_->{g}." D".$_->{D}." A".$_->{A};
	$html .= qq{
 <option value="$val">$name
};
      }
    $html .= "</select>";
    $html .= "<br><span class=small>Synonymous substitution rates previously calculated</span>" if $sqlite;
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
    my $ksdb = $args{ksdb};    
    my $kstype = $args{kstype};
    my $max = $args{max};
    my $min = $args{min};
# base=8_8.CDS-CDS.blastn.dag_geneorder_D60_g30_A5;
    
    $src .= ";flip=$flip" if $flip;
    $src .= ";regen=$regen" if $regen;
    $src .= ";width=$width" if $width;
    $src .= ";ksdb=$ksdb" if $ksdb;
    $src .= ";kstype=$kstype" if $kstype;
    $src .=  ";log=1" if $kstype;
    $src .=  ";min=$min" if defined $min;
    $src .=  ";max=$max" if defined $max;
    my $content = get("http://".$ENV{SERVER_NAME}."/".$src);
    my ($url) = $content =~ /url=(.*?)"/is;
    my $png = $url;
    $png =~ s/html$/png/;
    $png =~ s/$URL/$DIR/;
    my $img = GD::Image->new($png);
    my ($w,$h) = $img->getBounds();
    $w+=500;
    $h+=150;
    if ($loc)
      {
	return ($url,$w, $h);
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
