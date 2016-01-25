package CoGe::Builder::Tools::SynMap;

use Moose;

use CoGe::Accessory::Jex;
use CoGe::Accessory::Web qw( get_defaults );
use CoGe::Accessory::Workflow;
use CoGe::Builder::CommonTasks qw( create_gff_generation_job );
use CoGe::Core::Storage qw( get_workflow_paths );
use Data::Dumper;
use File::Path;
use File::Spec::Functions;
use JSON qw( encode_json );
use POSIX;

BEGIN {
	use Exporter 'import';
	our @EXPORT_OK =
	  qw( algo_lookup check_address_validity gen_org_name generate_pseudo_assembly get_logfile get_query_link go );
}

sub add_jobs {
	my %opts     = @_;
	my $workflow = $opts{workflow};
	my $db       = $opts{db};
	my $config   = $opts{config};
	$workflow->logfile(get_logfile($db, $config, %opts));

	my $dsgid1 = $opts{dsgid1} || $opts{genome_id1};
	my $dsgid2 = $opts{dsgid2} || $opts{genome_id2};
	my ( $dir1, $dir2 ) = sort ( $dsgid1, $dsgid2 );

	my $SEQUENCE_SIZE_LIMIT =
	  50_000_000;    # Limit the maximum genome size for genomic-genomic
	my $DIAGSDIR      = $config->{DIAGSDIR};
	my $SCRIPTDIR     = catdir( $config->{SCRIPTDIR}, 'synmap' );
	my $GEVO_LINKS    = catfile( $SCRIPTDIR, 'gevo_links.pl' );
	my $DAG_TOOL      = 'nice ' . catfile( $SCRIPTDIR, 'dag_tools.py' );
	my $BLAST2BED     = catfile( $SCRIPTDIR, 'blast2bed.pl' );
	my $GENE_ORDER    = catfile( $SCRIPTDIR, 'gene_order.py' );
	my $KSCALC        = catfile( $SCRIPTDIR, 'kscalc.pl' );
	my $GEN_FASTA     = catfile( $SCRIPTDIR, 'generate_fasta.pl' );
	my $RUN_ALIGNMENT = catfile( $SCRIPTDIR, 'quota_align_merge.pl' );
	my $RUN_COVERAGE  = catfile( $SCRIPTDIR, 'quota_align_coverage.pl' );
	my $PROCESS_DUPS  = catfile( $SCRIPTDIR, 'process_dups.pl' );
	my $DOTPLOT       = $config->{DOTPLOT} . " -cf " . $config->{_CONFIG_PATH};
	my $SVG_DOTPLOT   = catfile( $SCRIPTDIR, 'dotplot.py' );

#$RUN_DAGHAINER = $DIR."/bin/dagchainer/DAGCHAINER/run_DAG_chainer.pl -E 0.05 -s";
	my $RUN_DAGCHAINER =
	  'nice ' . $config->{PYTHON} . ' ' . $config->{DAGCHAINER};
	my $BLAST2RAW  = $config->{BLAST2RAW};    #find local duplicates
	my $FORMATDB   = $config->{FORMATDB};
	my $BLASTDBDIR = $config->{BLASTDB};
	my $FASTADIR   = $config->{FASTADIR};

	############################################################################
	# Fetch organism name and title
	############################################################################
	my $feat_type1 = $opts{feat_type1};
	my $feat_type2 = $opts{feat_type2};

	my ($genome1) = $db->resultset('Genome')->find($dsgid1);
	my ($genome2) = $db->resultset('Genome')->find($dsgid2);

	# Block large genomic-genomic jobs from running
	if (
		(
			   $feat_type1 == 2
			&& $genome1->length > $SEQUENCE_SIZE_LIMIT
			&& !$genome1->type->name =~ /hard/i
		)
		&& (   $feat_type2 == 2
			&& $genome2->length > $SEQUENCE_SIZE_LIMIT
			&& !$genome2->type->name =~ /hard/i )
	  )
	{

		return encode_json(
			{
				success => JSON::false,
				error   => "The analysis was blocked: "
				  . "a comparison of two unmasked and unannotated genomes larger than 50Mb requires many days to weeks to finish. "
				  . "Please use at least one annotated genome in the analysis."
			}
		);
	}

	#	my $basename = $opts{basename};
	my ( $org_name1, $title1 ) = gen_org_name(
		db        => $db,
		dsgid     => $dsgid1,
		feat_type => $feat_type1
	);

	my ( $org_name2, $title2 ) = gen_org_name(
		db        => $db,
		dsgid     => $dsgid2,
		feat_type => $feat_type2
	);

	my $ks_type = $opts{ks_type};
	############################################################################
	# Parameters
	############################################################################

	# blast options
	my $blast = $opts{blast};

	# blast2bed options

	# dagchainer options
	#my $dagchainer_g = $opts{g}; #depreciated -- will be a factor of -D
	my $dagchainer_D = $opts{D};
	my $dagchainer_A = $opts{A};
	my $Dm           = $opts{Dm};
	my $gm           = $opts{gm};
	($Dm) = $Dm =~ /(\d+)/;
	($gm) = $gm =~ /(\d+)/;

	#$dagchainer_type = $dagchainer_type eq "true" ? "geneorder" : "distance";

	#my $repeat_filter_cvalue = $opts{c};
	##parameter to be passed to run_adjust_dagchainer_evals

	#c-score for filtering low quality blast hits, fed to blast to raw
	my $cscore = $opts{csco};

	#tandem duplication distance, fed to blast to raw
	my $dupdist = defined( $opts{tdd} ) ? $opts{tdd} : 10;

	# dotplot options
	my $regen             = $opts{regen_images};
	my $regen_images      = ( $regen && $regen eq "true" ) ? 1 : 0;
	my $job_title         = $opts{jobtitle};
	my $width             = $opts{width};
	my $axis_metric       = $opts{axis_metric};
	my $axis_relationship = $opts{axis_relationship};
	my $min_chr_size      = $opts{min_chr_size};
	my $dagchainer_type   = $opts{dagchainer_type};
	my $color_type        = $opts{color_type};
	my $merge_algo        = $opts{merge_algo}; #is there a merging function? will non-syntenic dots be shown?
	my $snsd              = $opts{show_non_syn_dots} =~ /true/i ? 1 : 0;

	#	my $algo_name = $ALGO_LOOKUP->{$blast}{displayname};

	#will the axis be flipped?
	my $flip = $opts{flip} =~ /true/i ? 1 : 0;

	#are axes labeled?
	my $clabel = $opts{clabel} =~ /true/i ? 1 : 0;

	#are random chr skipped
	my $skip_rand = $opts{skip_rand} =~ /true/i ? 1 : 0;

	#which color scheme for ks/kn dots?
	my $color_scheme = $opts{color_scheme};

	#fids that are passed in for highlighting the pair in the dotplot
	my $fid1 = $opts{fid1};
	my $fid2 = $opts{fid2};

	#draw a box around identified diagonals?
	my $box_diags = $opts{box_diags};
	$box_diags = $box_diags eq "true" ? 1 : 0;

	#how are the chromosomes to be sorted?
	my $chr_sort_order = $opts{chr_sort_order};

	#codeml min and max calues
	my $codeml_min = $opts{codeml_min};
	$codeml_min = undef
	  unless $codeml_min =~ /\d/ && $codeml_min =~ /^-?\d*.?\d*$/;
	my $codeml_max = $opts{codeml_max};
	$codeml_max = undef
	  unless $codeml_max =~ /\d/ && $codeml_max =~ /^-?\d*.?\d*$/;
	my $logks = $opts{logks};
	$logks = $logks eq "true" ? 1 : 0;

	my $assemble = $opts{assemble} =~ /true/i ? 1 : 0;
	$assemble = 2 if $assemble && $opts{show_non_syn} =~ /true/i;
	$assemble *= -1 if $assemble && $opts{spa_ref_genome} < 0;

	#options for finding syntenic depth coverage by quota align (Bao's algo)
	my $depth_algo        = $opts{depth_algo};
	my $depth_org_1_ratio = $opts{depth_org_1_ratio};
	my $depth_org_2_ratio = $opts{depth_org_2_ratio};
	my $depth_overlap     = $opts{depth_overlap};

	$feat_type1 = $feat_type1 == 2 ? "genomic" : "CDS";
	$feat_type2 = $feat_type2 == 2 ? "genomic" : "CDS";
	$feat_type1 = "protein" if $blast == 5 && $feat_type1 eq "CDS"; #blastp time
	$feat_type2 = "protein" if $blast == 5 && $feat_type2 eq "CDS"; #blastp time

	# Sort by genome id
	(
		$dsgid1,     $genome1,           $org_name1,
		$feat_type1, $depth_org_1_ratio, $dsgid2,
		$genome2,    $org_name2,         $feat_type2,
		$depth_org_2_ratio
	  )
	  = (
		$dsgid2,     $genome2,           $org_name2,
		$feat_type2, $depth_org_2_ratio, $dsgid1,
		$genome1,    $org_name1,         $feat_type1,
		$depth_org_1_ratio
	  ) if ( $dsgid2 lt $dsgid1 );

	############################################################################
	# Generate Fasta files
	############################################################################
	my ( $fasta1, $fasta2 );

	#my @dsgs = ([$dsgid1, $feat_type1]);
	#push @dsgs, [$dsgid2, $feat_type2]
	#  unless $dsgid1 == $dsgid2 && $feat_type1 eq $feat_type2;
	#
	#    foreach my $item (@dsgs) {
	#        my $dsgid = $item->[0];
	#        my $feat_type = $item->[1];
	#
	#        #TODO: Schedule fasta generation only if feat_type not genomic
	#        #if ($feat_type eq "genomic") {
	#        #    my $genome = $coge->resultset('Genome')->find($gid);
	#        #    $file = $genome->file_path;
	#        #} else {
	#        #    $file = $FASTADIR . "/$gid-$feat_type.fasta";
	#        #}
	#    }
	if ( $feat_type1 eq "genomic" ) {
		$fasta1 = $genome1->file_path;

		$workflow->log( "Fetched fasta file for:" );
		$workflow->log( " " x (2) . $org_name1 );
	}
	else {
		my @fasta1args = ();
		$fasta1 = $FASTADIR . "/$dsgid1-$feat_type1.fasta";
		push @fasta1args, [ "--config", $config->{_CONFIG_PATH}, 0 ];
		push @fasta1args, [ "--genome_id",    $dsgid1,     1 ];
		push @fasta1args, [ "--feature_type", $feat_type1, 1 ];
		push @fasta1args, [ "--fasta",        $fasta1,     1 ];

		$workflow->add_job(
			{
				cmd         => $GEN_FASTA,
				script      => undef,
				args        => \@fasta1args,
				inputs      => undef,
				outputs     => [$fasta1],
				description => "Generating fasta file...",
			}
		);

		$workflow->log( "Added fasta file generation for:" );
		$workflow->log( " " x (2) . $org_name1 );
	}

	if ( $feat_type2 eq "genomic" ) {
		$fasta2 = $genome2->file_path;

		$workflow->log( "" );
		$workflow->log( "Fetched fasta file for:" );
		$workflow->log( " " x (2) . $org_name2 );
	}
	else {
		$fasta2 = $FASTADIR . "/$dsgid2-$feat_type2.fasta";
		my @fasta2args = ();
		push @fasta2args, [ "--config", $config->{_CONFIG_PATH}, 0 ];
		push @fasta2args, [ "--genome_id",    $dsgid2,     1 ];
		push @fasta2args, [ "--feature_type", $feat_type2, 1 ];
		push @fasta2args, [ "--fasta",        $fasta2,     1 ];

		$workflow->add_job(
			{
				cmd         => $GEN_FASTA,
				script      => undef,
				args        => \@fasta2args,
				inputs      => undef,
				outputs     => [$fasta2],
				description => "Generating fasta file...",
			}
		);

		$workflow->log( "" );
		$workflow->log( "Added fasta file generation for:" );
		$workflow->log( " " x (2) . $org_name2 );
	}

	############################################################################
	# Generate blastdb files
	############################################################################
	my ( $blastdb, @blastdb_files );
	my $ALGO_LOOKUP = algo_lookup();
	if ( $ALGO_LOOKUP->{$blast}{formatdb} ) {

		my $basename = "$BLASTDBDIR/$dsgid2-$feat_type2";

		my @blastdbargs = ();
		push @blastdbargs, [ '-p', $feat_type2 eq "protein" ? "T" : "F", 1 ];
		push @blastdbargs, [ '-i', $fasta2, 1 ];
		push @blastdbargs, [ '-t', '"' . $org_name2 . '"', 1 ];
		push @blastdbargs, [ '-n', $basename, 1 ];

		$blastdb = $basename;
		$basename .= $feat_type2 eq "protein" ? ".p" : ".n";

		push @blastdb_files, "$basename" . "sq";
		push @blastdb_files, "$basename" . "in";
		push @blastdb_files, "$basename" . "hr";

		$workflow->add_job(
			{
				cmd         => $FORMATDB,
				script      => undef,
				args        => \@blastdbargs,
				inputs      => [$fasta2],
				outputs     => \@blastdb_files,
				description => "Generating BlastDB...",
			}
		);

		$workflow->log( "" );
		$workflow->log( "Added BlastDB generation" );
		$workflow->log( $blastdb );
	}
	else {
		$blastdb = $fasta2;
		push @blastdb_files, $blastdb;
	}
	my ( $orgkey1, $orgkey2 ) = ( $title1, $title2 );
	my %org_dirs = (
		$orgkey1 . "_"
		  . $orgkey2 => {
			fasta    => $fasta1,
			db       => $blastdb,
			basename => $dsgid1 . "_" 
			  . $dsgid2
			  . ".$feat_type1-$feat_type2."
			  . $ALGO_LOOKUP->{$blast}{filename},
			dir => "$DIAGSDIR/$dir1/$dir2",
		  },
	);

	foreach my $org_dir ( keys %org_dirs ) {
		my $outfile = $org_dirs{$org_dir}{dir};
		mkpath( $outfile, 0, 0777 ) unless -d $outfile;
		warn "didn't create path $outfile: $!" unless -d $outfile;
		$outfile .= "/" . $org_dirs{$org_dir}{basename};
		$org_dirs{$org_dir}{blastfile} = $outfile;    #.".blast";
	}

	############################################################################
	# Run Blast
	############################################################################

	#check blast runs for problems;  Not forked in order to keep variables
	my $raw_blastfile = $org_dirs{ $orgkey1 . "_" . $orgkey2 }{blastfile};

	foreach my $key ( keys %org_dirs ) {
		my $cmd = 'nice '
		  . $ALGO_LOOKUP->{$blast}
		  {algo};    #$prog =~ /tblastx/i ? $TBLASTX : $BLASTN;

		my $fasta   = $org_dirs{$key}{fasta};
		my $db      = $org_dirs{$key}{db};
		my $outfile = $org_dirs{$key}{blastfile};
		my @blastargs;

		if ( $cmd =~ /lastz/i ) {
			push @blastargs, [ "-i", $fasta,   0 ];
			push @blastargs, [ "-d", $db,      0 ];
			push @blastargs, [ "-o", $outfile, 1 ];
		}
		elsif ( $cmd =~ /last_wrapper/i ) {

			# mdb added 9/20/13 issue 213
			my $dbpath = $config->{LASTDB} . '/' . $dsgid2;
			mkpath( $dbpath, 0, 0777 );
			push @blastargs, [ "--dbpath", $dbpath, 0 ];

			push @blastargs, [ "",   $db,      0 ];
			push @blastargs, [ "",   $fasta,   0 ];
			push @blastargs, [ "-o", $outfile, 1 ];
		}
		else {
			push @blastargs, [ "-out",   $outfile, 1 ];
			push @blastargs, [ "-query", $fasta,   0 ];
			push @blastargs, [ "-db",    $db,      0 ];
		}

		( undef, $cmd ) = CoGe::Accessory::Web::check_taint($cmd);
		push @blastdb_files, $fasta;
		$workflow->add_job(
			{
				cmd         => $cmd,
				script      => undef,
				args        => \@blastargs,
				inputs      => \@blastdb_files,
				outputs     => [$outfile],
				description => "Running genome comparison...",
			}
		);
	}

	$workflow->log( "" );
	$workflow->log(
		"Added genome comparison (algorithm: "
		  . $ALGO_LOOKUP->{$blast}{displayname} . ")");

	###########################################################################
	# Converting blast to bed and finding local duplications
	###########################################################################
	#
	# NOTES: The blast2bed program will take the input rawblast file and
	# filter it and creating a new rawblast and moving the old to
	# rawblast.orig
	#
	my $query_bed   = $raw_blastfile . ".q.bed";
	my $subject_bed = $raw_blastfile . ".s.bed";

	my @blastargs = ();
	push @blastargs, [ '-infile',   $raw_blastfile, 1 ];
	push @blastargs, [ '-outfile1', $query_bed,     1 ];
	push @blastargs, [ '-outfile2', $subject_bed,   1 ];

	my @bedoutputs = ();
	push @bedoutputs, $query_bed;
	push @bedoutputs, $subject_bed;
	push @bedoutputs, $raw_blastfile if ( $raw_blastfile =~ /genomic/ );

#    push @bedoutputs, "$raw_blastfile.orig" if ( $raw_blastfile =~ /genomic/ );
	push @bedoutputs, "$raw_blastfile.new" if ( $raw_blastfile =~ /genomic/ );
	$workflow->add_job(
		{
			cmd         => $BLAST2BED,
			script      => undef,
			args        => \@blastargs,
			inputs      => [$raw_blastfile],
			outputs     => \@bedoutputs,
			description => "Creating .bed files...",
		}
	);

	$workflow->log( "" );
	$workflow->log( "Added .bed files creation" );

	###########################################################################
	# Converting blast to raw and finding local duplications
	###########################################################################
	my $qlocaldups = $raw_blastfile . ".q.localdups";
	my $slocaldups = $raw_blastfile . ".s.localdups";
	$raw_blastfile .= ".new" if $raw_blastfile =~ /genomic/;
	my $filtered_blastfile = $raw_blastfile;

	#    $filtered_blastfile .=".new" if $raw_blastfile =~ /genomic/;
	$filtered_blastfile .= ".tdd$dupdist";
	$filtered_blastfile .= ".cs$cscore" if $cscore < 1;
	$filtered_blastfile .= ".filtered";

	my @rawargs = ();
	push @rawargs, [ "", $raw_blastfile, 1 ];
	push @rawargs, [ "--localdups", "", 1 ];
	push @rawargs, [ "--qbed",        $query_bed,          1 ];
	push @rawargs, [ "--sbed",        $subject_bed,        1 ];
	push @rawargs, [ "--tandem_Nmax", $dupdist,            1 ];
	push @rawargs, [ "--cscore",      $cscore,             1 ] if $cscore < 1;
	push @rawargs, [ ">",             $filtered_blastfile, 1 ];

	my @rawoutputs = ();
	push @rawoutputs, $filtered_blastfile;
	push @rawoutputs, $qlocaldups;
	push @rawoutputs, $slocaldups;

	$workflow->add_job(
		{
			cmd         => $BLAST2RAW,
			script      => undef,
			args        => \@rawargs,
			inputs      => [ $raw_blastfile, $query_bed, $subject_bed ],
			outputs     => \@rawoutputs,
			description => "Filtering tandem dups...",
		}
	);

	$workflow->log( "" );
	my $msg = "Added Filtering results of tandem ";
	$msg .= "duplicates (tandem duplication distance: $dupdist";
	$msg .= ", c-score: $cscore" if $cscore;
	$msg .= ")";

	$workflow->log( $msg );

 #TODO: This feature is currently disabled
 #needed to comment out as the bed files and blast files have changed in SynFind
 #   my $synteny_score_db = run_synteny_score(
 #       blastfile => $filtered_blastfile,
 #       bedfile1  => $bedfile1,
 #       bedfile2  => $bedfile2,
 #       outfile   => $org_dirs{$orgkey1 . "_" . $orgkey2}{dir} . "/"
 #         . $dsgid1 . "_"
 #         . $dsgid2
 #         . ".$feat_type1-$feat_type2");

	############################################################################
	# Run dag tools - Prepares DAG for syntenic analysis.
	############################################################################
	my $dag_file12       = $filtered_blastfile . ".dag";
	my $dag_file12_all   = $dag_file12 . ".all";
	my $query_dup_file   = $opts{query_dup_files};
	my $subject_dup_file = $opts{subject_dup_files};
	my $query            = "a" . $dsgid1;
	my $subject          = "b" . $dsgid2;

	my @dagtoolargs = ();
	push @dagtoolargs, [ '-q', $query,              1 ];
	push @dagtoolargs, [ '-s', $subject,            1 ];
	push @dagtoolargs, [ '-b', $filtered_blastfile, 1 ];
	push @dagtoolargs, [ '-c', "", 1 ];

	push @dagtoolargs, [ '--query_dups', $query_dup_file, 1 ]
	  if $query_dup_file;
	push @dagtoolargs, [ '--subject_dups', $subject_dup_file, 1 ]
	  if $subject_dup_file;
	push @dagtoolargs, [ '>', $dag_file12_all, 1 ];

	$workflow->add_job(
		{
			cmd         => $DAG_TOOL,
			script      => undef,
			args        => \@dagtoolargs,
			inputs      => [$filtered_blastfile],
			outputs     => [$dag_file12_all],
			description => "Formatting for DAGChainer...",
		}
	);

	$workflow->log( "" );
	$workflow->log( "Added convertion of blast file to dagchainer input file" );

	############################################################################
	# Convert to gene order
	############################################################################
	my $dag_file12_all_geneorder = "$dag_file12_all.go";
	my $all_file;

	if ( $dagchainer_type eq "geneorder" ) {

		$workflow->log( "" );
		$workflow->log("Added convertion of dagchainer input into gene order coordinates");

		my @geneorderargs = ();
		push @geneorderargs, [ "",           $dag_file12_all,           1 ];
		push @geneorderargs, [ "",           $dag_file12_all_geneorder, 1 ];
		push @geneorderargs, [ "--gid1",     $dsgid1,                   1 ];
		push @geneorderargs, [ "--gid2",     $dsgid2,                   1 ];
		push @geneorderargs, [ "--feature1", $feat_type1,               1 ];
		push @geneorderargs, [ "--feature2", $feat_type2,               1 ];

		$workflow->add_job(
			{
				cmd         => $GENE_ORDER,
				script      => undef,
				args        => \@geneorderargs,
				inputs      => [$dag_file12_all],
				outputs     => [$dag_file12_all_geneorder],
				description => "Converting to genomic order...",
			}
		);

		$all_file = $dag_file12_all_geneorder;
		$dag_file12 .= ".go";
	}
	else {
		$all_file = $dag_file12_all;
	}

  #B Pedersen's program for automatically adjusting the evals in the dag file to
  # remove bias from local gene duplicates and transposons
  #   $dag_file12 .= "_c" . $repeat_filter_cvalue;
  #   $workflow->log( "#" x (20) );
  #   $workflow->log( "Adjusting evalue of blast hits to correct
  #   for repeat sequences" );
  #   run_adjust_dagchainer_evals( infile => $all_file, outfile => $dag_file12,
  #   cvalue => $repeat_filter_cvalue );
  #   $workflow->log( "#" x (20) );
  #   $workflow->log( "" );

	# This step will fail if the dag_file_all is larger than the system memory
	# limit. If this file does not exist, let's send a warning to the log file
	# and continue with the analysis using the dag_all file
	unless ( ( -r $dag_file12 && -s $dag_file12 )
		|| ( -r $dag_file12 . ".gz" && -s $dag_file12 . ".gz" ) )
	{
		$dag_file12 = $all_file;
		$workflow->log( "" );
		$workflow->log(
			"WARNING: sub run_adjust_dagchainer_evals failed. "
			  . "Perhaps due to Out of Memory error. "
			  . "Proceeding without this step!");
	}

	############################################################################
	# Run dagchainer
	############################################################################
	my ( $dagchainer_file, $merged_dagchainer_file );

	#this is for using dagchainer's merge function
	my $dag_merge_enabled = ( $merge_algo == 2 ) ? 1 : 0;
	my $self_comparision = ( $dsgid1 eq $dsgid2 ) ? 1 : 0;

	#length of a gap (average distance expected between two syntenic genes)
	my $gap = defined( $opts{g} ) ? $opts{g} : floor( $dagchainer_D / 2 );
	$gap = 1 if $gap < 1;
	$dagchainer_file = $dag_file12;
	$dagchainer_file .= "_D$dagchainer_D" if defined $dagchainer_D;
	$dagchainer_file .= "_g$gap"          if defined $gap;
	$dagchainer_file .= "_A$dagchainer_A" if defined $dagchainer_A;
	$dagchainer_file .= "_Dm$Dm"          if $dag_merge_enabled;
	$dagchainer_file .= "_gm$gm"          if $dag_merge_enabled;
	$dagchainer_file .= ".aligncoords";
	$dagchainer_file .= ".ma2.dag"        if $dag_merge_enabled;

	my @dagargs = ();
	push @dagargs, [ "-E", "0.05", 1 ];
	push @dagargs, [ "-i",   $dag_file12,   1 ];
	push @dagargs, [ "-D",   $dagchainer_D, 1 ] if defined $dagchainer_D;
	push @dagargs, [ "-g",   $gap,          1 ] if defined $gap;
	push @dagargs, [ "-A",   $dagchainer_A, 1 ] if defined $dagchainer_A;
	push @dagargs, [ "--Dm", $Dm,           1 ] if $dag_merge_enabled;
	push @dagargs, [ "--m",  $gm,           1 ] if $dag_merge_enabled;
	push @dagargs, [ "--new_behavior", "", 1 ] if $self_comparision;

	###MERGING OF DIAGONALS FUNCTION
	# --merge $outfile
	# this will automatically cause the merging of diagonals to happen.
	# $outfile will be created and $outfile.merge (which is the inappropriately
	# named merge of diagonals file.
	# --gm  average distance between diagonals
	# --Dm  max distance between diagonals
	##Both of these parameters' default values is 4x -g and -D respectively.
	my $post_dagchainer_file;

	if ($dag_merge_enabled) {
		$merged_dagchainer_file = "$dagchainer_file.merged";
		push @dagargs, [ "--merge", $merged_dagchainer_file, 1 ];

		$workflow->add_job(
			{
				cmd         => $RUN_DAGCHAINER,
				script      => undef,
				args        => \@dagargs,
				inputs      => [$dag_file12],
				outputs     => [$merged_dagchainer_file],
				description => "Running DAGChainer (with merge)...",
			}
		);
		$post_dagchainer_file = $merged_dagchainer_file;
		$workflow->log( "" );
		$workflow->log(
			"Added DagChainer (merge: enabled,"
			  . " Maximum distance between two blocks: $Dm genes, "
			  . " Average distance expected between syntenic blocks: $gm genes)");
	}
	else {
		push @dagargs, [ ">", $dagchainer_file, 1 ];

		$workflow->add_job(
			{
				cmd         => $RUN_DAGCHAINER,
				script      => undef,
				args        => \@dagargs,
				inputs      => [$dag_file12],
				outputs     => [$dagchainer_file],
				description => "Running DAGChainer...",
			}
		);

		$post_dagchainer_file = $dagchainer_file;
		$workflow->log( "" );
		$workflow->log( "Added DagChainer (merge: disabled)" );
	}

	############################################################################
	# Run quota align merge
	############################################################################
	#id 1 is to specify quota align as a merge algo
	if ( $merge_algo == 1 ) {
		$merged_dagchainer_file = "$dagchainer_file.Dm$Dm.ma1";

		my @mergeargs = ();
		push @mergeargs, [ '--config', $config->{_CONFIG_PATH}, 0 ];
		push @mergeargs, [ '--infile',       $dagchainer_file,        1 ];
		push @mergeargs, [ '--outfile',      $merged_dagchainer_file, 1 ];
		push @mergeargs, [ '--max_distance', $Dm,                     1 ];

		$workflow->add_job(
			{
				cmd         => $RUN_ALIGNMENT,
				script      => undef,
				args        => \@mergeargs,
				inputs      => [$dagchainer_file],
				outputs     => [$merged_dagchainer_file],
				description => "Merging Syntenic Blocks...",
			}
		);
		$workflow->log( "" );
		$workflow->log(
			"Added Merge Syntenic Blocks"
			  . " (Maximum distance between two blocks: $Dm genes)");

		$post_dagchainer_file = $merged_dagchainer_file;
	}
	my $post_dagchainer_file_w_nearby = $post_dagchainer_file;
	$post_dagchainer_file_w_nearby =~ s/aligncoords/all\.aligncoords/;

	#add pairs that were skipped by dagchainer
	$post_dagchainer_file_w_nearby = $post_dagchainer_file; # are the previous two lines necessary if we just overwrite here?

	#TODO: Run find nearby is currently disabled
	#run_find_nearby(
	#    infile       => $post_dagchainer_file,
	#    dag_all_file => $all_file,
	#    outfile      => $post_dagchainer_file_w_nearby
	#);    #program is not working correctly.

	############################################################################
	# Run quota align coverage
	############################################################################
	my ( $quota_align_coverage, $grimm_stuff, $final_dagchainer_file );

	if ( $depth_algo == 1 )    #id 1 is to specify quota align
	{
		$quota_align_coverage = $post_dagchainer_file_w_nearby;
		$quota_align_coverage .= ".qac" . $depth_org_1_ratio . ".";
		$quota_align_coverage .= $depth_org_2_ratio . "." . $depth_overlap;

		my @depthargs = ();
		push @depthargs, [ '--config', $config->{_CONFIG_PATH}, 0 ];
		push @depthargs, [ '--infile',  $post_dagchainer_file_w_nearby, 1 ];
		push @depthargs, [ '--outfile', $quota_align_coverage,          1 ];
		push @depthargs, [ '--depth_ratio_org1', $depth_org_1_ratio, 1 ];
		push @depthargs, [ '--depth_ratio_org2', $depth_org_2_ratio, 1 ];
		push @depthargs, [ '--depth_overlap',    $depth_overlap,     1 ];

		$workflow->add_job(
			{
				cmd         => $RUN_COVERAGE,
				script      => undef,
				args        => \@depthargs,
				inputs      => [$post_dagchainer_file_w_nearby],
				outputs     => [$quota_align_coverage],
				description => "Calculating Syntenic depth...",
			}
		);

		$workflow->log( "" );
		$workflow->log(
			"Added Syntenic Depth "
			  . "(ratio: $depth_org_1_ratio/$depth_org_2_ratio, "
			  . "overlap: $depth_overlap)");
		$workflow->log( "$org_name1 -vs- $org_name2" );
		$final_dagchainer_file = $quota_align_coverage;
	}
	else {
		$final_dagchainer_file = $post_dagchainer_file_w_nearby;
	}

	if ( $dagchainer_type eq "geneorder" ) {
		my @positionargs = ();
		push @positionargs, [ '', $final_dagchainer_file, 1 ];
		push @positionargs, [ '', "$final_dagchainer_file.gcoords", 1 ];
		push @positionargs, [ "--positional", '', 1 ];

		$workflow->add_job(
			{
				cmd         => $GENE_ORDER,
				script      => undef,
				args        => \@positionargs,
				inputs      => [$final_dagchainer_file],
				outputs     => ["$final_dagchainer_file.gcoords"],
				description => "Converting to genomic coordinates...",
			}
		);

		$workflow->log( "" );
		$workflow->log("Added conversion gene order coordinates back to genomic coordinates");

		$final_dagchainer_file = $final_dagchainer_file . ".gcoords";
	}

	#generate dotplot images
	unless ($width) {
		my $chr1_count = $genome1->chromosome_count;
		my $chr2_count = $genome2->chromosome_count;
		my $max_chr    = $chr1_count;
		$max_chr = $chr2_count if $chr2_count > $chr1_count;
		$width   = int( $max_chr * 100 );
		$width   = 1000 if $width > 1000;
		$width   = 500 if $width < 500;
	}

	############################################################################
	# Create html output directory
	############################################################################
#	my ( $qlead, $slead ) = ( "a", "b" );
	my $out = $org_dirs{ $orgkey1 . "_" . $orgkey2 }{dir} . "/html/";
	mkpath( $out, 0, 0777 ) unless -d $out;
	$out .= "master_";
	my ($base) = $final_dagchainer_file =~ /([^\/]*$)/;
	$out .= $base;

	#	my $json_basename = "$out";
	$out .= "_ct$color_type" if defined $color_type;
	$out .= ".w$width";

	############################################################################
	# KS Calculations
	############################################################################
	my ( $ks_db, $ks_blocks_file, $svg_file );
	my $check_ks = $final_dagchainer_file =~ /^(.*?CDS-CDS)/;
	$check_ks = $final_dagchainer_file =~ /^(.*?protein-protein)/
	  unless $check_ks;

	if ( $ks_type and $check_ks ) {
		$ks_db          = "$final_dagchainer_file.sqlite";
		$ks_blocks_file = "$final_dagchainer_file.ks";

		my @ksargs = ();
		push @ksargs, [ '--config', $config->{_CONFIG_PATH}, 0 ];
		push @ksargs, [ '--infile',    $final_dagchainer_file, 1 ];
		push @ksargs, [ '--dbfile',    $ks_db,                 1 ];
		push @ksargs, [ '--blockfile', $ks_blocks_file,        1 ];

		$workflow->add_job(
			{
				cmd         => $KSCALC,
				script      => undef,
				args        => \@ksargs,
				inputs      => [$final_dagchainer_file],
				outputs     => [ $ks_blocks_file, $ks_db ],
				description => "Calculating synonymous changes (slow)...",
			}
		);

		$workflow->log( "" );
		$workflow->log("Added ($ks_type) calculation of syntenic CDS pairs and color dots");

		####################################################################
		# Generate svg dotplot
		####################################################################

		my @svgargs = ();
		push @svgargs, [ '--dag_file', $ks_blocks_file, 1 ];
		push @svgargs, [ '--flip', "", 1 ] if $flip;
		push @svgargs, [ '--xhead', '"' . $org_name1 . '"', 1 ]
		  if $org_name1;
		push @svgargs, [ '--yhead', '"' . $org_name2 . '"', 1 ]
		  if $org_name2;
		push @svgargs, [ '--output', $ks_blocks_file, 1 ];

		$svg_file = $ks_blocks_file . ".svg";
		$workflow->add_job(
			{
				cmd         => $SVG_DOTPLOT,
				script      => undef,
				args        => \@svgargs,
				inputs      => [$ks_blocks_file],
				outputs     => [$svg_file],
				description => "Generating svg image...",
			}
		);

		$workflow->log( "" );
		$workflow->log( "Added generation of svg dotplot" );
	}
	else {
		$ks_type = undef;

		####################################################################
		# Generate svg dotplot
		####################################################################

		my @svgargs = ();
		push @svgargs, [ '--dag_file', $final_dagchainer_file, 1 ];
		push @svgargs, [ '--flip', "", 1 ] if $flip;
		push @svgargs, [ '--xhead', '"' . $org_name1 . '"', 1 ]
		  if $org_name1;
		push @svgargs, [ '--yhead', '"' . $org_name2 . '"', 1 ]
		  if $org_name2;
		push @svgargs, [ '--output', $final_dagchainer_file, 1 ];
		push @svgargs, [ '--no-ks', "", 1 ];

		$svg_file = $final_dagchainer_file . ".svg";

		$workflow->add_job(
			{
				cmd         => $SVG_DOTPLOT,
				script      => undef,
				args        => \@svgargs,
				inputs      => [$final_dagchainer_file],
				outputs     => [$svg_file],
				description => "Generating svg image...",
			}
		);

		$workflow->log( "" );
		$workflow->log( "Added generation of svg dotplot" );

	}

	############################################################################
	# Generate dot plot
	############################################################################
	my @plotargs   = ();
	my @plotinputs = ();

	my ($basename) =
	  $final_dagchainer_file =~ /([^\/]*aligncoords.*)/;    #.all.aligncoords/;
	$width = 1000 unless defined($width);

	my $dotfile = "$out";
	$dotfile .= ".spa$assemble"     if $assemble;
	$dotfile .= ".gene"             if $axis_metric =~ /gene/i;
	$dotfile .= ".s"                if $axis_relationship =~ /s/i;
	$dotfile .= ".mcs$min_chr_size" if $min_chr_size;
	$dotfile .= ".$fid1"            if $fid1;
	$dotfile .= ".$fid2"            if $fid2;

	#add ks_db to dotplot command if requested
	if ($ks_type) {
		$dotfile .= ".$ks_type";

		push @plotargs, [ '-ksdb', $ks_db,   1 ];
		push @plotargs, [ '-kst',  $ks_type, 1 ];
		push @plotargs, [ '-log',  $logks,   1 ];
		push @plotinputs, $ks_db;
	}
	$dotfile .= ".box"                if $box_diags;
	$dotfile .= ".flip"               if $flip;
	$dotfile .= ".c0"                 if $clabel eq 0;
	$dotfile .= ".sr"                 if $skip_rand;
	$dotfile .= ".cs$color_scheme"    if defined $color_scheme;
	$dotfile .= ".cso$chr_sort_order" if defined $chr_sort_order;
	$dotfile .= ".min$codeml_min"     if defined $codeml_min;
	$dotfile .= ".max$codeml_max"     if defined $codeml_max;
	$dotfile .= ".log"                if $logks;

	#are non-syntenic dots being displayed
	if ($snsd) {
		push @plotargs, [ '-d', $dag_file12_all, 0 ];
		push @plotinputs, $dag_file12_all;
	}
	else {
		$dotfile .= ".nsd";    #no syntenic dots, yes, nomicalture is confusing.
	}

	push @plotargs, [ '-a', $final_dagchainer_file, 1 ];
	push @plotargs, [ '-b', $dotfile, 1 ];

	my $jsoption = "";
	$jsoption .= qq{'javascript:synteny_zoom(};
	$jsoption .= qq{"$dsgid1",};
	$jsoption .= qq{"$dsgid2",};
	$jsoption .= qq{"$basename",};
	$jsoption .= $flip ? qq{"YCHR","XCHR"} : qq{"XCHR","YCHR"};
	$jsoption .= qq{,"$ks_db"} if $ks_db;
	$jsoption .= qq{)'};

	push @plotargs, [ '-l',    $jsoption, 0 ];
	push @plotargs, [ '-dsg1', $dsgid1,   1 ];
	push @plotargs, [ '-dsg2', $dsgid2,   1 ];
	push @plotargs, [ '-w',    $width,    1 ];
	push @plotargs, [ '-lt', 2, 1 ];
	push @plotargs, [ '-assemble', $assemble,    1 ] if $assemble;
	push @plotargs, [ '-am',       $axis_metric, 1 ] if $axis_metric;
	push @plotargs, [ '-fb', '', 1 ]
	  if $axis_relationship && $axis_relationship =~ /s/;
	push @plotargs, [ '-mcs', $min_chr_size, 1 ] if $min_chr_size;
	push @plotargs, [ '-cdt', $color_type,   1 ] if $color_type;
	push @plotargs, [ '-bd', 1, 1 ] if $box_diags;
	push @plotargs, [ '-fid1', $fid1, 1 ] if $fid1;
	push @plotargs, [ '-fid2', $fid2, 1 ] if $fid2;
	push @plotargs, [ '-f',      1, 1 ] if $flip;
	push @plotargs, [ '-labels', 0, 1 ] if $clabel eq 0;
	push @plotargs, [ '-sr',     1, 1 ] if $skip_rand;
	push @plotargs, [ '-color_scheme', $color_scheme, 1 ]
	  if defined $color_scheme;
	push @plotargs, [ '-chr_sort_order', $chr_sort_order, 1 ]
	  if defined $chr_sort_order;
	push @plotargs, [ '-min', $codeml_min, 1 ] if defined $codeml_min;
	push @plotargs, [ '-max', $codeml_max, 1 ] if defined $codeml_max;

	push @plotinputs, $final_dagchainer_file;

	my $hist = $dotfile . ".hist.png";

	# this would be generated by the DOTPLOT program is Syntenic path assembly
	# was requested
	my $spa_file = $dotfile . ".spa_info.txt";
	my @plotoutputs = ( "$dotfile.html", "$dotfile.png" );

	push @plotoutputs, ( "$dotfile.x.png", "$dotfile.y.png" ) if $clabel eq "1";
	push @plotoutputs, $hist     if $ks_db;
	push @plotoutputs, $spa_file if $assemble;

	$workflow->add_job(
		{
			cmd         => $DOTPLOT,
			script      => undef,
			args        => \@plotargs,
			inputs      => \@plotinputs,
			outputs     => \@plotoutputs,
			overwrite   => $regen_images,
			description => "Generating images...",
		}
	);

	#    my $dot_args = [
	#        [ '-cf', $config->{_CONFIG_PATH}, 0 ],
	#        [ '-genome1', $dsgid1, 0 ],
	#        [ '-genome2', $dsgid2, 0 ],
	#        [ '-a', $final_dagchainer_file, 1 ],
	#        [ '-b', $json_basename, 1 ],
	#    ];
	#
	#    push @$dot_args, [ '-ksdb', $ks_db,   1 ] if $ks_db;
	#
	#    my $dot_inputs = [
	#        $final_dagchainer_file,
	#    ];
	#
	#    my $dot_outputs = [
	#        "$json_basename.json",
	#    ];
	#
	#    if ($dag_file12_all) {
	#        push @$dot_args, [ '-d', $dag_file12_all, 0 ];
	#        push @$dot_inputs, $dag_file12_all;
	#        push @$dot_outputs, "$json_basename.all.json";
	#    }
	#
	#    if ($ks_db) {
	#        push @$dot_inputs, $ks_db if $ks_db;
	#        push @$dot_outputs, "$json_basename.datasets.json";
	#    }
	#
	#    $workflow->add_job(
	#        cmd         => $DOTPLOT_DOTS,
	#        script      => undef,
	#        args        => $dot_args,
	#        inputs      => $dot_inputs,
	#        outputs     => $dot_outputs,
	#        description => "Generating dotplot dots...",
	#    );

	############################################################################
	# Post Processing
	############################################################################

	$workflow->log( "" );
	$workflow->log( "Final Post Processing" );

	my $subject_dup_args = [
		[ '--config', $config->{_CONFIG_PATH}, 0 ],
		[ '--infile', $slocaldups, 1 ],   #$raw_blastfile . ".s.localdups", 1 ],
		[ '--outfile', $raw_blastfile . ".s.tandems", 1 ],
	];

	$workflow->add_job(
		{
			cmd    => $PROCESS_DUPS,
			script => undef,
			args   => $subject_dup_args,
			inputs => [$slocaldups],       #[$raw_blastfile . ".s.localdups"],
			outputs     => [ $raw_blastfile . ".s.tandems" ],
			description => "Processing Subject Tandem Duplicate File...",
		}
	);

	my $query_dup_args = [
		[ '--config', $config->{_CONFIG_PATH}, 0 ],
		[ '--infile', $qlocaldups, 1 ],   #$raw_blastfile . ".q.localdups", 1 ],
		[ '--outfile', $raw_blastfile . ".q.tandems", 1 ],
	];

	$workflow->add_job(
		{
			cmd    => $PROCESS_DUPS,
			script => undef,
			args   => $query_dup_args,
			inputs => [$qlocaldups],      #[$raw_blastfile . ".q.localdups"],
			outputs     => [ $raw_blastfile . ".q.tandems" ],
			description => "Processing Query Tandem Duplicate File...",
		}
	);
	$workflow->log( "" );
	$workflow->log( "Added Processing of Tandem Duplicate Files" );

	my $condensed = "$final_dagchainer_file.condensed";

	my $link_args = [
		[ '--config',  $config->{_CONFIG_PATH}, 0 ],
		[ '--infile',  $final_dagchainer_file,  1 ],
		[ '--dsgid1',  $dsgid1,                 1 ],
		[ '--dsgid2',  $dsgid2,                 1 ],
		[ '--outfile', $condensed,              1 ],
	];

	$workflow->add_job(
		{
			cmd         => $GEVO_LINKS,
			script      => undef,
			args        => $link_args,
			inputs      => [$final_dagchainer_file],
			outputs     => [$condensed],
			description => "Generating GEvo links...",
		}
	);

	if ( $opts{frac_bias} =~ /true/i ) {
		my $organism_name;
		my $target_id;
		if ( $depth_org_1_ratio < $depth_org_2_ratio ) {
			$organism_name = $genome1->organism->name;
			$target_id     = $dsgid1;
		}
		else {
			$organism_name = $genome2->organism->name;
			$target_id     = $dsgid2;
		}

		my $gff_job = create_gff_generation_job(
			gid           => $target_id,
			organism_name => $organism_name
		);
		$workflow->add_job($gff_job);

		my $output_dir = catdir( $config->{DIAGSDIR}, $dir1, $dir2 );
		$workflow->add_job(
			{
				cmd => 'python '
				  . $config->{SCRIPTDIR}
				  . '/synmap/fractionation_bias_geco.py',
				script => undef,
				args   => [
					[ '--align',      $final_dagchainer_file,   0 ],
					[ '--gff',        $gff_job->{outputs}->[0], 0 ],
					[ '--target',     $target_id,               0 ],
					[ '--windowsize', $opts{fb_window_size},    0 ],
					[
						'--allgenes',
						( $opts{fb_target_genes} eq 'true' ) ? 'False' : 'True',
						0
					],
					[ '--output', $output_dir, 0 ]
				],
				inputs => [
					$final_dagchainer_file, $condensed,
					$gff_job->{outputs}->[0]
				],
				outputs => [
					catfile(
						catdir( $output_dir, 'html' ),
						'fractbias_figure1.png'
					)
				],
				description => "Running Fractination Bias...",
			}
		);
	}

	$workflow->log( "" );
	$workflow->log( "Added GEvo links generation" );
	$workflow->log( "#" x (25) );
	$workflow->log( "" );
}

sub algo_lookup {

#in the web form, each sequence search algorithm has a unique number.  This table identifies those and adds appropriate options
	my $config        = get_defaults();
	my $MAX_PROC      = $config->{MAX_PROC};
	my $blast_options = " -num_threads $MAX_PROC -evalue 0.0001 -outfmt 6";
	my $TBLASTX       = $config->{TBLASTX} . $blast_options;
	my $BLASTN        = $config->{BLASTN} . $blast_options;
	my $BLASTP        = $config->{BLASTP} . $blast_options;
	my $LASTZ =
	    $config->{PYTHON} . " "
	  . $config->{MULTI_LASTZ}
	  . " -A $MAX_PROC --path="
	  . $config->{LASTZ};
	my $LAST =
	  $config->{MULTI_LAST} . " -a $MAX_PROC --path=" . $config->{LAST_PATH};
	return {
		0 => {
			algo => $BLASTN . " -task megablast",    #megablast
			opt         => "MEGA_SELECT",  #select option for html template file
			filename    => "megablast",
			displayname => "MegaBlast",
			html_select_val => 0,
			formatdb        => 1,
		},
		1 => {
			algo => $BLASTN . " -task dc-megablast",   #discontinuous megablast,
			opt  => "DCMEGA_SELECT",
			filename        => "dcmegablast",
			displayname     => "Discontinuous MegaBlast",
			html_select_val => 1,
			formatdb        => 1,
		},
		2 => {
			algo            => $BLASTN . " -task blastn",    #blastn
			opt             => "BLASTN_SELECT",
			filename        => "blastn",
			displayname     => "BlastN",
			html_select_val => 2,
			formatdb        => 1,
		},
		3 => {
			algo            => $TBLASTX,                     #tblastx
			opt             => "TBLASTX_SELECT",
			filename        => "tblastx",
			displayname     => "TBlastX",
			html_select_val => 3,
			formatdb        => 1,
		},
		4 => {
			algo            => $LASTZ,                       #lastz
			opt             => "LASTZ_SELECT",
			filename        => "lastz",
			displayname     => "(B)lastZ",
			html_select_val => 4,
		},
		5 => {
			algo            => $BLASTP . " -task blastp",    #blastn
			opt             => "BLASTP_SELECT",
			filename        => "blastp",
			displayname     => "BlastP",
			html_select_val => 5,
			formatdb        => 1,
		},
		6 => {
			algo            => $LAST,                        #last
			opt             => "LAST_SELECT",
			filename        => "last",
			displayname     => "Last",
			html_select_val => 6,
		}
	};
}

sub build {
	my $self = shift;

	# Validate inputs
	my $genome_id1 = $self->params->{genome_id1};
	return unless $genome_id1;
	my $genome_id2 = $self->params->{genome_id2};
	return unless $genome_id2;

	my %opts = ( %{ defaults() }, %{ $self->params } );
	add_jobs(
		workflow => $self->workflow,
		db       => $self->db,
		config   => $self->conf,
		%opts
	);

	return 1;
}

sub check_address_validity {
	my $address = shift;
	return 'valid' unless $address;
	my $validity =
	  $address =~
/^[_a-zA-Z0-9-]+(\.[_a-zA-Z0-9-]+)*@[a-zA-Z0-9-]+(\.[a-zA-Z0-9-]+)*\.(([0-9]{1,3})|([a-zA-Z]{2,3})|(aero|coop|info|museum|name))$/
	  ? 'valid'
	  : 'invalid';
	return $validity;
}

sub defaults {
	return {
		A                 => 5,
		assemble          => 'false',
		axis_metric       => 'nt',
		axis_relationship => 'r',
		blast             => 6,
		box_diags         => 'false',
		chr_sort_order    => 'S',
		clabel            => 1,
		codeml_max        => '',
		codeml_min        => '',
		color_scheme      => 1,
		csco              => 0,
		D                 => 20,
		dagchainer_type   => 'geneorder',
		depth_algo        => 0,
		Dm                => 0,
		feat_type1        => 1,
		feat_type2        => 1,
		flip              => 'false',
		frac_bias         => 'false',
		gm                => 0,
		ks_type           => 0,
		logks             => 1,
		merge_algo        => 1,
		show_non_syn_dots => 'false',
		skip_rand         => 1
	};
}

sub gen_org_name {
	my %opts      = @_;
	my $db        = $opts{db};
	my $dsgid     = $opts{dsgid};
	my $feat_type = $opts{feat_type} || 1;

	my ($dsg) = $db->resultset('Genome')->search( { genome_id => $dsgid },
		{ join => 'organism', prefetch => 'organism' } );

	my $org_name = $dsg->organism->name;
	my $title =
	    $org_name . " (v"
	  . $dsg->version
	  . ", dsgid"
	  . $dsgid . ") "
	  . $feat_type;
	$title =~ s/(`|')//g;

	return ( $org_name, $title );
}

sub generate_pseudo_assembly {
	my ( $config, $input, $output, $flip ) = @_;
	$flip = 0 unless $flip;

	my $cmd = "synmap/order_contigs_to_chromosome.pl";

	my $JEX = CoGe::Accessory::Jex->new(
		host => $config->{JOBSERVER},
		port => $config->{JOBPORT}
	);
	my $workflow = $JEX->create_workflow( name => "Generate Pseudo Assembly" );

	$workflow->add_job(
		{
			cmd  => catfile( $config->{SCRIPTDIR}, $cmd ),
			args => [
				[ "-cfg",    $config->{_CONFIG_PATH}, 1 ],
				[ "-input",  $input,                  1 ],
				[ "-output", $output,                 1 ],
				[ "-flip",   $flip,                   1 ],
			],
			inputs  => [ $input, $config->{_CONFIG_PATH} ],
			outputs => [$output],
			description => "Generating pseudo assembly"
		}
	);

	my $response = $JEX->submit_workflow($workflow);

	return {
		id      => $response->{id},
		success => $JEX->is_successful($response),
		output  => $output
	};
}

sub get_logfile {
    my $db = shift;
    my $config = shift;
	my %opts     = @_;
	my $dsgid1 = $opts{dsgid1} || $opts{genome_id1};
	my $dsgid2 = $opts{dsgid2} || $opts{genome_id2};
	my ( $dir1, $dir2 ) = sort ( $dsgid1, $dsgid2 );
	my $feat_type1 = $opts{feat_type1};
	my $feat_type2 = $opts{feat_type2};
	my ($genome1) = $db->resultset('Genome')->find($dsgid1);
	my ($genome2) = $db->resultset('Genome')->find($dsgid2);
	my ( $org_name1, $title1 ) = gen_org_name(
		db        => $db,
		dsgid     => $dsgid1,
		feat_type => $feat_type1
	);
	my ( $org_name2, $title2 ) = gen_org_name(
		db        => $db,
		dsgid     => $dsgid2,
		feat_type => $feat_type2
	);
	my $ks_type = $opts{ks_type};
	my $blast = $opts{blast};
	my $dagchainer_D = $opts{D};
	my $dagchainer_A = $opts{A};
	my $Dm           = $opts{Dm};
	my $gm           = $opts{gm};
	($Dm) = $Dm =~ /(\d+)/;
	($gm) = $gm =~ /(\d+)/;
	my $cscore = $opts{csco};
	my $dupdist = defined( $opts{tdd} ) ? $opts{tdd} : 10;
	my $width             = $opts{width};
	my $axis_metric       = $opts{axis_metric};
	my $axis_relationship = $opts{axis_relationship};
	my $min_chr_size      = $opts{min_chr_size};
	my $dagchainer_type   = $opts{dagchainer_type};
	my $color_type        = $opts{color_type};
	my $merge_algo        = $opts{merge_algo}; #is there a merging function? will non-syntenic dots be shown?
	my $snsd              = $opts{show_non_syn_dots} =~ /true/i ? 1 : 0;
	my $flip = $opts{flip} =~ /true/i ? 1 : 0;
	my $clabel = $opts{clabel} =~ /true/i ? 1 : 0;
	my $skip_rand = $opts{skip_rand} =~ /true/i ? 1 : 0;
	my $color_scheme = $opts{color_scheme};
	my $fid1 = $opts{fid1};
	my $fid2 = $opts{fid2};
	my $box_diags = $opts{box_diags};
	$box_diags = $box_diags eq "true" ? 1 : 0;
	my $chr_sort_order = $opts{chr_sort_order};
	my $codeml_min = $opts{codeml_min};
	$codeml_min = undef
	  unless $codeml_min =~ /\d/ && $codeml_min =~ /^-?\d*.?\d*$/;
	my $codeml_max = $opts{codeml_max};
	$codeml_max = undef
	  unless $codeml_max =~ /\d/ && $codeml_max =~ /^-?\d*.?\d*$/;
	my $logks = $opts{logks};
	$logks = $logks eq "true" ? 1 : 0;
	my $assemble = $opts{assemble} =~ /true/i ? 1 : 0;
	$assemble = 2 if $assemble && $opts{show_non_syn} =~ /true/i;
	$assemble *= -1 if $assemble && $opts{spa_ref_genome} < 0;
	my $depth_algo        = $opts{depth_algo};
	my $depth_org_1_ratio = $opts{depth_org_1_ratio};
	my $depth_org_2_ratio = $opts{depth_org_2_ratio};
	my $depth_overlap     = $opts{depth_overlap};
	$feat_type1 = $feat_type1 == 2 ? "genomic" : "CDS";
	$feat_type2 = $feat_type2 == 2 ? "genomic" : "CDS";
	$feat_type1 = "protein" if $blast == 5 && $feat_type1 eq "CDS"; #blastp time
	$feat_type2 = "protein" if $blast == 5 && $feat_type2 eq "CDS"; #blastp time
	(
		$dsgid1,     $genome1,           $org_name1,
		$feat_type1, $depth_org_1_ratio, $dsgid2,
		$genome2,    $org_name2,         $feat_type2,
		$depth_org_2_ratio
	  )
	  = (
		$dsgid2,     $genome2,           $org_name2,
		$feat_type2, $depth_org_2_ratio, $dsgid1,
		$genome1,    $org_name1,         $feat_type1,
		$depth_org_1_ratio
	  ) if ( $dsgid2 lt $dsgid1 );
	my ( $orgkey1, $orgkey2 ) = ( $title1, $title2 );
	my $ALGO_LOOKUP = algo_lookup();
	my $DIAGSDIR      = $config->{DIAGSDIR};
	my %org_dirs = (
		$orgkey1 . "_" . $orgkey2 => {
			basename => $dsgid1 . "_" 
			  . $dsgid2
			  . ".$feat_type1-$feat_type2."
			  . $ALGO_LOOKUP->{$blast}{filename},
			dir => "$DIAGSDIR/$dir1/$dir2",
		  },
	);
	foreach my $org_dir ( keys %org_dirs ) {
		my $outfile = $org_dirs{$org_dir}{dir};
		$outfile .= "/" . $org_dirs{$org_dir}{basename};
		$org_dirs{$org_dir}{blastfile} = $outfile;    #.".blast";
	}
	my $raw_blastfile = $org_dirs{ $orgkey1 . "_" . $orgkey2 }{blastfile};
	$raw_blastfile .= ".new" if $raw_blastfile =~ /genomic/;
	my $filtered_blastfile = $raw_blastfile;
	$filtered_blastfile .= ".tdd$dupdist";
	$filtered_blastfile .= ".cs$cscore" if $cscore < 1;
	$filtered_blastfile .= ".filtered";
	my $dag_file12       = $filtered_blastfile . ".dag";
	my $dag_file12_all   = $dag_file12 . ".all";
	my $dag_file12_all_geneorder = "$dag_file12_all.go";
	my $all_file;
	if ( $dagchainer_type eq "geneorder" ) {
		$all_file = $dag_file12_all_geneorder;
		$dag_file12 .= ".go";
	}
	else {
		$all_file = $dag_file12_all;
	}
	unless ( ( -r $dag_file12 && -s $dag_file12 )
		|| ( -r $dag_file12 . ".gz" && -s $dag_file12 . ".gz" ) )
	{
		$dag_file12 = $all_file;
	}
	my ( $dagchainer_file );
	my $dag_merge_enabled = ( $merge_algo == 2 ) ? 1 : 0;
	my $gap = defined( $opts{g} ) ? $opts{g} : floor( $dagchainer_D / 2 );
	$gap = 1 if $gap < 1;
	$dagchainer_file = $dag_file12;
	$dagchainer_file .= "_D$dagchainer_D" if defined $dagchainer_D;
	$dagchainer_file .= "_g$gap"          if defined $gap;
	$dagchainer_file .= "_A$dagchainer_A" if defined $dagchainer_A;
	$dagchainer_file .= "_Dm$Dm"          if $dag_merge_enabled;
	$dagchainer_file .= "_gm$gm"          if $dag_merge_enabled;
	$dagchainer_file .= ".aligncoords";
	$dagchainer_file .= ".ma2.dag"        if $dag_merge_enabled;
	my $post_dagchainer_file;
	if ($dag_merge_enabled) {
		$post_dagchainer_file = "$dagchainer_file.merged";
	}
	else {
		$post_dagchainer_file = $dagchainer_file;
	}
	if ( $merge_algo == 1 ) {
		$post_dagchainer_file = "$dagchainer_file.Dm$Dm.ma1";
	}
	my $post_dagchainer_file_w_nearby = $post_dagchainer_file;
	$post_dagchainer_file_w_nearby =~ s/aligncoords/all\.aligncoords/;
	$post_dagchainer_file_w_nearby = $post_dagchainer_file; # are the previous two lines necessary if we just overwrite here?
	my ( $quota_align_coverage, $final_dagchainer_file );
	if ( $depth_algo == 1 ) {    #id 1 is to specify quota align
		$quota_align_coverage = $post_dagchainer_file_w_nearby;
		$quota_align_coverage .= ".qac" . $depth_org_1_ratio . ".";
		$quota_align_coverage .= $depth_org_2_ratio . "." . $depth_overlap;
		$final_dagchainer_file = $quota_align_coverage;
	}
	else {
		$final_dagchainer_file = $post_dagchainer_file_w_nearby;
	}
	$final_dagchainer_file .= ".gcoords" if $dagchainer_type eq "geneorder";
	unless ($width) {
		my $chr1_count = $genome1->chromosome_count;
		my $chr2_count = $genome2->chromosome_count;
		my $max_chr    = $chr1_count;
		$max_chr = $chr2_count if $chr2_count > $chr1_count;
		$width   = int( $max_chr * 100 );
		$width   = 1000 if $width > 1000;
		$width   = 500 if $width < 500;
	}
	my $out = $org_dirs{ $orgkey1 . "_" . $orgkey2 }{dir} . "/html/";
	$out .= "master_";
	my ($base) = $final_dagchainer_file =~ /([^\/]*$)/;
	$out .= $base;
	$out .= "_ct$color_type" if defined $color_type;
	$out .= ".w$width";
	my $check_ks = $final_dagchainer_file =~ /^(.*?CDS-CDS)/;
	$check_ks = $final_dagchainer_file =~ /^(.*?protein-protein)/ unless $check_ks;
	$ks_type = undef unless $ks_type and $check_ks;
	my $dotfile = $out;
	$dotfile .= ".spa$assemble"     if $assemble;
	$dotfile .= ".gene"             if $axis_metric =~ /gene/i;
	$dotfile .= ".s"                if $axis_relationship =~ /s/i;
	$dotfile .= ".mcs$min_chr_size" if $min_chr_size;
	$dotfile .= ".$fid1"            if $fid1;
	$dotfile .= ".$fid2"            if $fid2;
	$dotfile .= ".$ks_type"          if $ks_type;
	$dotfile .= ".box"                if $box_diags;
	$dotfile .= ".flip"               if $flip;
	$dotfile .= ".c0"                 if $clabel eq 0;
	$dotfile .= ".sr"                 if $skip_rand;
	$dotfile .= ".cs$color_scheme"    if defined $color_scheme;
	$dotfile .= ".cso$chr_sort_order" if defined $chr_sort_order;
	$dotfile .= ".min$codeml_min"     if defined $codeml_min;
	$dotfile .= ".max$codeml_max"     if defined $codeml_max;
	$dotfile .= ".log"                if $logks;
	$dotfile .= ".nsd" unless $snsd;    #no syntenic dots, yes, nomicalture is confusing.
    return $dotfile . '.log';
}

sub get_name {
	return 'SynMap';
}

sub get_query_link {
	my $config       = shift;
	my $db           = shift;
	my %url_options  = @_;
	my $dagchainer_D = $url_options{D};

   #  my $dagchainer_g = $url_options{g}; #depreciated -- will be a factor of -D
	my $dagchainer_A = $url_options{A};
	my $Dm           = $url_options{Dm};
	my $gm           = $url_options{gm};
	($Dm) = $Dm =~ /(\d+)/;
	($gm) = $gm =~ /(\d+)/;

#   my $repeat_filter_cvalue = $url_options{c}; #parameter to be passed to run_adjust_dagchainer_evals
	my $cscore =
	  $url_options{csco}; #c-score for filtering low quality blast hits, fed to blast to raw
	my $dupdist =
	  $url_options{tdd};    #tandem duplication distance, fed to blast to raw
	my $regen_images = $url_options{regen_images};
	my $job_title    = $url_options{jobtitle};
	my $width        = $url_options{width};

	#	my $basename     = $url_options{basename};
	my $blast = $url_options{blast};

	my $feat_type1 = $url_options{feat_type1};
	my $feat_type2 = $url_options{feat_type2};

	my $dsgid1 = $url_options{dsgid1};
	my $dsgid2 = $url_options{dsgid2};

	unless ( $dsgid1 and $dsgid2 ) {
		return encode_json( { error => "Missing a genome id." } );
	}

	my $assemble = $url_options{assemble} =~ /true/i ? 1 : 0;
	$assemble = 2 if $assemble && $url_options{show_non_syn} =~ /true/i;
	$assemble *= -1 if $assemble && $url_options{spa_ref_genome} < 0;
	my $axis_metric       = $url_options{axis_metric};
	my $axis_relationship = $url_options{axis_relationship};
	my $min_chr_size      = $url_options{min_chr_size};
	my $dagchainer_type   = $url_options{dagchainer_type};
	my $color_type        = $url_options{color_type};
	my $merge_algo = $url_options{merge_algo};    #is there a merging function?

	#options for finding syntenic depth coverage by quota align (Bao's algo)
	my $depth_algo        = $url_options{depth_algo};
	my $depth_org_1_ratio = $url_options{depth_org_1_ratio};
	my $depth_org_2_ratio = $url_options{depth_org_2_ratio};
	my $depth_overlap     = $url_options{depth_overlap};

	#options for fractionation bias
	my $frac_bias       = $url_options{frac_bias} =~ /true/i;
	my $fb_window_size  = $url_options{fb_window_size};
	my $fb_target_genes = $url_options{fb_target_genes};

	#fids that are passed in for highlighting the pair in the dotplot
	#	my $fid1 = $url_options{fid1};
	#	my $fid2 = $url_options{fid2};

	#will non-syntenic dots be shown?
	my $snsd = $url_options{show_non_syn_dots} =~ /true/i ? 1 : 0;
#	my $ALGO_LOOKUP = algo_lookup();
#	my $algo_name = $ALGO_LOOKUP->{$blast}{displayname};

	#will the axis be flipped?
	my $flip = $url_options{flip};
	$flip = $flip =~ /true/i ? 1 : 0;

	#are axes labeled?
	my $clabel = $url_options{clabel};
	$clabel = $clabel =~ /true/i ? 1 : 0;

	#are random chr skipped
	my $skip_rand = $url_options{skip_rand};
	$skip_rand = $skip_rand =~ /true/i ? 1 : 0;

	#which color scheme for ks/kn dots?
	my $color_scheme = $url_options{color_scheme};

	#codeml min and max calues
	my $codeml_min = $url_options{codeml_min};
	$codeml_min = undef
	  unless $codeml_min =~ /\d/ && $codeml_min =~ /^-?\d*.?\d*$/;
	my $codeml_max = $url_options{codeml_max};
	$codeml_max = undef
	  unless $codeml_max =~ /\d/ && $codeml_max =~ /^-?\d*.?\d*$/;
	my $logks = $url_options{logks};
	$logks = $logks eq "true" ? 1 : 0;

	#how are the chromosomes to be sorted?
	my $chr_sort_order = $url_options{chr_sort_order};

	#draw a box around identified diagonals?
	my $box_diags = $url_options{box_diags};
	$box_diags = $box_diags eq "true" ? 1 : 0;

	my ( $org_name1, $titleA ) = gen_org_name(
		db        => $db,
		dsgid     => $dsgid1,
		feat_type => $feat_type1
	);

	my ( $org_name2, $titleB ) = gen_org_name(
		db        => $db,
		dsgid     => $dsgid2,
		feat_type => $feat_type2
	);

	# Sort by genome id
	(
		$dsgid1, $org_name1, $feat_type1, $depth_org_1_ratio, $dsgid2,
		$org_name2, $feat_type2, $depth_org_2_ratio
	  )
	  = (
		$dsgid2, $org_name2, $feat_type2, $depth_org_2_ratio, $dsgid1,
		$org_name1, $feat_type1, $depth_org_1_ratio
	  ) if ( $dsgid2 lt $dsgid1 );

	my $synmap_link =
	    $config->{SERVER}
	  . "SynMap.pl?dsgid1=$dsgid1;dsgid2=$dsgid2"
	  . ";D=$dagchainer_D;A=$dagchainer_A;w=$width;b=$blast;ft1=$feat_type1;"
	  . "ft2=$feat_type2;autogo=1";

	$synmap_link .= ";Dm=$Dm"       if defined $Dm;
	$synmap_link .= ";csco=$cscore" if $cscore;
	$synmap_link .= ";tdd=$dupdist" if defined $dupdist;
	$synmap_link .= ";gm=$gm"       if defined $gm;
	$synmap_link .= ";snsd=$snsd";

	$synmap_link .= ";bd=$box_diags"          if $box_diags;
	$synmap_link .= ";mcs=$min_chr_size"      if $min_chr_size;
	$synmap_link .= ";sp=$assemble"           if $assemble;
	$synmap_link .= ";ma=$merge_algo"         if $merge_algo;
	$synmap_link .= ";da=$depth_algo"         if $depth_algo;
	$synmap_link .= ";do1=$depth_org_1_ratio" if $depth_org_1_ratio;
	$synmap_link .= ";do2=$depth_org_2_ratio" if $depth_org_2_ratio;
	$synmap_link .= ";do=$depth_overlap"      if $depth_overlap;
	$synmap_link .= ";fb=1"                   if $frac_bias;
	$synmap_link .= ";fb_ws=$fb_window_size"  if $fb_window_size;
	$synmap_link .= ";fb_tg=1"                if $fb_target_genes eq 'true';
	$synmap_link .= ";flip=1"                 if $flip;
	$synmap_link .= ";cs=$color_scheme";
	$synmap_link .= ";cmin=$codeml_min"
	  if defined $codeml_min;   #$codeml_min=~/\d/ && $codeml_min=~/^\d*.?\d*$/;
	$synmap_link .= ";cmax=$codeml_max"
	  if defined $codeml_max;   #$codeml_max=~/\d/ && $codeml_max=~/^\d*.?\d*$/;
	$synmap_link .= ";logks=$logks"        if defined $logks;
	$synmap_link .= ";cl=0"                if $clabel eq "0";
	$synmap_link .= ";sr=$skip_rand"       if defined $skip_rand;
	$synmap_link .= ";cso=$chr_sort_order" if $chr_sort_order;

	$feat_type1 = $feat_type1 == 2 ? "genomic" : "CDS";
	$feat_type2 = $feat_type2 == 2 ? "genomic" : "CDS";
	$feat_type1 = "protein" if $blast == 5 && $feat_type1 eq "CDS"; #blastp time
	$feat_type2 = "protein" if $blast == 5 && $feat_type2 eq "CDS"; #blastp time

	$synmap_link .= ";dt=$dagchainer_type";

	my $ks_type = $url_options{ks_type};
	if ($ks_type) {
		my $num;
		if    ( $ks_type eq "ks" )    { $num = 1; }
		elsif ( $ks_type eq "kn" )    { $num = 2; }
		elsif ( $ks_type eq "kn_ks" ) { $num = 3; }
		$synmap_link .= ";ks=$num";
	}
	$synmap_link .= ";am=g" if $axis_metric       && $axis_metric       =~ /g/i;
	$synmap_link .= ";ar=s" if $axis_relationship && $axis_relationship =~ /s/i;
	$synmap_link .= ";ct=$color_type" if $color_type;

	my $tiny_link = CoGe::Accessory::Web::get_tiny_link( url => $synmap_link );

	return $tiny_link;
}

sub go {
	my %opts = @_;
	foreach my $k ( keys %opts ) {
		$opts{$k} =~ s/^\s+//;
		$opts{$k} =~ s/\s+$//;
	}
	my $dsgid1 = $opts{dsgid1};
	my $dsgid2 = $opts{dsgid2};
	return encode_json(
		{
			success => JSON::false,
			error   => "You must select two genomes."
		}
	) unless ( $dsgid1 && $dsgid2 );

	my ( $db, $user, $config ) = CoGe::Accessory::Web->init();

	my ($genome1) = $db->resultset('Genome')->find($dsgid1);
	my ($genome2) = $db->resultset('Genome')->find($dsgid2);

	return encode_json(
		{
			success => JSON::false,
			error   => "The Genome $dsgid1 could not be found."
		}
	) unless $genome1;

	return encode_json(
		{
			success => JSON::false,
			error   => "The Genome $dsgid2 could not be found."
		}
	) unless $genome2;

	return encode_json(
		{
			success => JSON::false,
			error   => "Genome $dsgid1 primary data is missing."
		}
	) unless -r $genome1->file_path;

	return encode_json(
		{
			success => JSON::false,
			error   => "Genome $dsgid2 primary data is missing."
		}
	) unless -r $genome2->file_path;

	my $tiny_link = get_query_link( $config, $db, @_ );
	warn "tiny_link is required for logging." unless defined($tiny_link);
	my ($tiny_id) = $tiny_link =~ /\/(\w+)$/;
	my $workflow_name = "synmap-$tiny_id";

	my $JEX = CoGe::Accessory::Jex->new(
		host => $config->{JOBSERVER},
		port => $config->{JOBPORT}
	);
	my $workflow = $JEX->create_workflow( name => $workflow_name );
	$workflow->log_section( "Creating Workflow" );
	$workflow->log( "Link to Regenerate Analysis" );
	$workflow->log( "$tiny_link" );
	$workflow->log( "" );
	$workflow->log( "Created Workflow: $workflow_name" );
	$workflow->log( "" );

	add_jobs(
		workflow => $workflow,
		db       => $db,
		config   => $config,
		@_
	);

	$workflow->log_section( "Running Workflow" );

	my $response = $JEX->submit_workflow($workflow);

	if ( $response and $response->{id} ) {
		my $feat_type1 = $opts{feat_type1};
		my $feat_type2 = $opts{feat_type2};
		my ( $org_name1, $title1 ) = gen_org_name(
			db        => $db,
			dsgid     => $dsgid1,
			feat_type => $feat_type1
		);
		my ( $org_name2, $title2 ) = gen_org_name(
			db        => $db,
			dsgid     => $dsgid2,
			feat_type => $feat_type2
		);
		my $log_msg = "$org_name1 v. $org_name2"
		  ; #"<a href='OrganismView.pl?dsgid=$dsgid1' target='_blank'>$org_name1</a> v. <a href='OrganismView.pl?dsgid=$dsgid2' target='_blank'>$org_name2</a>";
		$log_msg .= " Ks" if $opts{ks_type};

		my $log = CoGe::Accessory::Web::log_history(
			db          => $db,
			user_id     => $user->id,
			description => $log_msg,
			page        => 'SynMap',
			link        => $tiny_link,
			parent_id   => $response->{id},
			parent_type => 7                  #FIXME magic number
		);

		return encode_json(
			{
				link    => $tiny_link,
				request => "jex/synmap/status/" . $response->{id},
				status  => $response->{status},
				success => $JEX->is_successful($response)
				? JSON::true
				: JSON::false
			}
		);
	}
	else {
		return encode_json( { success => JSON::false } );
	}
}

# FIXME: Currently this feature is disabled.
# @by Evan Briones
# @on 3/1/2013
#sub run_tandem_finder {
#	my %opts    = @_;
#	my $infile  = $opts{infile};    #dag file produced by dat_tools.py
#	my $outfile = $opts{outfile};
#	while ( -e "$outfile.running" ) {
#		print STDERR "detecting $outfile.running.  Waiting. . .\n";
#		sleep 60;
#	}
#	unless ( -r $infile && -s $infile ) {
#		CoGe::Accessory::Web::write_log( "", $cogeweb->logfile );
#		CoGe::Accessory::Web::write_log(
#"WARNING:   Cannot run tandem finder! DAGChainer input file ($infile) contains no data!",
#			$cogeweb->logfile
#		);
#		return 0;
#	}
#	if ( -r $outfile ) {
#		CoGe::Accessory::Web::write_log(
#			"run_tandem_filter: file $outfile already exists",
#			$cogeweb->logfile );
#		return 1;
#	}
#	my $cmd = "$PYTHON $TANDEM_FINDER -i '$infile' > '$outfile'";
#	system "/usr/bin/touch '$outfile.running'"
#	  ;    #track that a blast anlaysis is running for this
#	CoGe::Accessory::Web::write_log( "run_tandem_filter: running\n\t$cmd",
#		$cogeweb->logfile );
#	`$cmd`;
#	system "/bin/rm '$outfile.running'"
#	  if -r "$outfile.running";    #remove track file
#	return 1 if -r $outfile;
#}

#FIXME: Currently this feature is disabled
# @by Evan Briones
# @on 3/1/2013
#sub run_adjust_dagchainer_evals {
#	my %opts    = @_;
#	my $infile  = $opts{infile};
#	my $outfile = $opts{outfile};
#	my $cvalue  = $opts{cvalue};
#	$cvalue = 4 unless defined $cvalue;
#	while ( -e "$outfile.running" ) {
#		print STDERR "detecting $outfile.running.  Waiting. . .\n";
#		sleep 60;
#	}
#	if ( -r $outfile || -r $outfile . ".gz" ) {
#		CoGe::Accessory::Web::write_log(
#			"run_adjust_dagchainer_evals: file $outfile already exists",
#			$cogeweb->logfile );
#		return 1;
#	}
#	CoGe::Accessory::Web::gunzip( $infile . ".gz" ) if -r $infile . ".gz";
#	unless ( -r $infile && -s $infile ) {
#		CoGe::Accessory::Web::write_log(
#"WARNING:   Cannot adjust dagchainer evals! DAGChainer input file ($infile) contains no data!",
#			$cogeweb->logfile
#		);
#		return 0;
#	}
#	my $cmd = "$PYTHON $EVAL_ADJUST -c $cvalue '$infile' > '$outfile'";
#
##There is a parameter that can be passed into this to filter repetitive sequences more or less stringently:
## -c   2 gets rid of more stuff; 10 gets rid of less stuff; default is 4
##consider making this a parameter than can be adjusted from SynMap -- will need to actually play with this value to see how it works
##if implemented, this will require re-naming all the files to account for this parameter
##and updating the auto-SynMap link generator for redoing an analysis
#
#	system "/usr/bin/touch '$outfile.running'"
#	  ;    #track that a blast anlaysis is running for this
#	CoGe::Accessory::Web::write_log(
#		"run_adjust_dagchainer_evals: running\n\t$cmd",
#		$cogeweb->logfile );
#	`$cmd`;
#	system "/bin/rm '$outfile.running'" if -r "$outfile.running";
#	;      #remove track file
#	return 1 if -r $outfile;
#
#}

#FIXME: Currently this feature is disabled
# @by Evan Briones
# @on 6/18/2013
#sub run_find_nearby {
#	my %opts         = @_;
#	my $infile       = $opts{infile};
#	my $dag_all_file = $opts{dag_all_file};
#	my $outfile      = $opts{outfile};
#	while ( -e "$outfile.running" ) {
#		print STDERR "detecting $outfile.running.  Waiting. . .\n";
#		sleep 60;
#	}
#	if ( -r $outfile ) {
#		CoGe::Accessory::Web::write_log(
#			"run find_nearby: file $outfile already exists",
#			$cogeweb->logfile );
#		return 1;
#	}
#	my $cmd =
#"$PYTHON $FIND_NEARBY --diags='$infile' --all='$dag_all_file' > '$outfile'";
#	system "/usr/bin/touch '$outfile.running'"
#	  ;    #track that a blast anlaysis is running for this
#	CoGe::Accessory::Web::write_log( "run find_nearby: running\n\t$cmd",
#		$cogeweb->logfile );
#	`$cmd`;
#	system "/bin/rm '$outfile.running'" if -r "$outfile.running";
#	;      #remove track file
#	return 1 if -r $outfile;
#}

with qw(CoGe::Builder::Buildable);

1;
