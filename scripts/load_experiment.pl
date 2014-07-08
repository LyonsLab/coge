#!/usr/bin/perl -w

use strict;
use Data::Dumper;
use CoGeX;
use Getopt::Long;
use File::Path;
use File::Touch;
use File::Basename qw( basename );
use File::Spec::Functions qw( catdir catfile );
use URI::Escape::JavaScript qw(unescape);
use JSON::XS;
use CoGe::Accessory::Web qw(get_defaults);
use CoGe::Accessory::Utils qw( commify );
use CoGe::Core::Metadata qw( create_annotations );
use CoGe::Accessory::TDS;

use vars qw($staging_dir $result_dir $install_dir $data_file $file_type $log_file
  $name $description $version $restricted $ignore_missing_chr
  $gid $source_name $user_name $config $allow_negative $disable_range_check
  $annotations $types
  $host $port $db $user $pass $P);

#FIXME: use these from Storage.pm instead of redeclaring them
my $DATA_TYPE_QUANT  = 1; # Quantitative data
my $DATA_TYPE_POLY	 = 2; # Polymorphism data
my $DATA_TYPE_ALIGN  = 3; # Alignments
my $DATA_TYPE_MARKER = 4; # Markers

#my $MIN_QUANT_COLUMNS = 5;
#my $MAX_QUANT_COLUMNS = 6;
my $MIN_VCF_COLUMNS = 8;
#my $MAX_VCF_COLUMNS = 10;
my $MIN_GFF_COLUMNS = 9;

GetOptions(
    "staging_dir=s" => \$staging_dir,    # temporary staging path
    "install_dir=s" => \$install_dir,    # final installation path
    "result_dir=s"  => \$result_dir,     # results path
    "data_file=s"   => \$data_file,      # input data file (JS escape)
    "file_type=s"   => \$file_type,		 # input file type
    "name=s"        => \$name,           # experiment name (JS escaped)
    "desc=s"        => \$description,    # experiment description (JS escaped)
    "version=s"     => \$version,        # experiment version (JS escaped)
    "restricted=i"  => \$restricted,     # experiment restricted flag
    "source_name=s" => \$source_name,    # experiment source name (JS escaped)
    "gid=s"         => \$gid,            # genome id
    "user_name=s"   => \$user_name,      # user name
    "annotations=s" => \$annotations,    # optional: semicolon-separated list of locked annotations (link:group:type:text;...)
    "types=s"       => \$types,          # optional: semicolon-separated list of experiment type names
    "config=s"      => \$config,         # configuration file

    # Optional flags for debug and bulk loader
    "ignore-missing-chr=i" => \$ignore_missing_chr,
    "allow_negative=i"     => \$allow_negative,
    "disable_range_check"  => \$disable_range_check, # allow any value in val1 column
    "log_file=s"           => \$log_file
);

# Open log file
$| = 1;
die unless ($staging_dir);
mkpath($staging_dir); # make sure this exists
$log_file = catfile($staging_dir, 'log.txt') unless $log_file;
mkpath($staging_dir, 0, 0777) unless -r $staging_dir;
open( my $log, ">>$log_file" ) or die "Error opening log file $log_file";
$log->autoflush(1);

# Process and verify parameters
$data_file   = unescape($data_file);
$name        = unescape($name);
$description = unescape($description);
$version     = unescape($version);
$source_name = unescape($source_name);
$restricted  = '0' if ( not defined $restricted );

if ($user_name eq 'public') {
	print $log "log: error: not logged in\n";
    exit(-1);
}

# Load config file
unless ($config) {
    print $log "log: error: can't find config file\n";
    print STDERR "can't find config file\n";
    exit(-1);
}
$P    = CoGe::Accessory::Web::get_defaults($config);
$db   = $P->{DBNAME};
$host = $P->{DBHOST};
$port = $P->{DBPORT};
$user = $P->{DBUSER};
$pass = $P->{DBPASS};

my $FASTBIT_LOAD  = $P->{FASTBIT_LOAD};
my $FASTBIT_QUERY = $P->{FASTBIT_QUERY};
my $SAMTOOLS      = $P->{SAMTOOLS};
my $GUNZIP        = $P->{GUNZIP};
if (   not $FASTBIT_LOAD
    or not $FASTBIT_QUERY
    or not $SAMTOOLS
    or not $GUNZIP
    or not -e $FASTBIT_LOAD
    or not -e $FASTBIT_QUERY
    or not -e $SAMTOOLS
    or not -e $GUNZIP )
{
    print $log "log: error: can't find required command(s)\n";
    exit(-1);
}

my $cmd;

# Copy input data file to staging area
# If running via JEX the file will already be there
my ($filename) = basename($data_file);#$data_file =~ /^.+\/([^\/]+)$/;
my $staged_data_file = $staging_dir . '/' . $filename;
unless (-r $staged_data_file) {
    $cmd = "cp -f '$data_file' $staging_dir";
    `$cmd`;
}

# Decompress file if necessary
if ( $staged_data_file =~ /\.gz$/ ) {
    my $cmd = $GUNZIP . ' ' . $staged_data_file;
    #print STDERR "$cmd\n";
    print $log "log: Decompressing '$filename'\n";
    `$cmd`;
    $staged_data_file =~ s/\.gz$//;
}

# Determine file type
my ($file_type, $data_type) = detect_data_type($file_type, $staged_data_file);
if ( !$file_type or !$data_type ) {
    print $log "log: error: unknown or unsupported file type '$file_type'\n";
    exit(-1);
}

# Connect to database
my $connstr = "dbi:mysql:dbname=$db;host=$host;port=$port;";
my $coge = CoGeX->connect( $connstr, $user, $pass );
unless ($coge) {
    print $log "log: couldn't connect to database\n";
    exit(-1);
}

# Retrieve genome
my $genome = $coge->resultset('Genome')->find( { genome_id => $gid } );
unless ($genome) {
    print $log "log: error finding genome id$gid\n";
    exit(-1);
}

# Hash chromosome names
my %genome_chr = map { $_ => 1 } $genome->chromosomes;

# Validate the data file
print $log "log: Validating data file\n";
if (-s $staged_data_file == 0) {
    print $log "log: error: input file '$staged_data_file' is empty\n";
    exit(-1);
}
my $count = 0;
my $pChromosomes;
my $format;
if ( $data_type == $DATA_TYPE_QUANT ) {
    ( $staged_data_file, $format, $count, $pChromosomes ) =
      validate_quant_data_file( file => $staged_data_file, file_type => $file_type, genome_chr => \%genome_chr );
}
elsif ( $data_type == $DATA_TYPE_POLY ) {
    ( $staged_data_file, $format, $count, $pChromosomes ) =
      validate_vcf_data_file( file => $staged_data_file, genome_chr => \%genome_chr );
}
elsif ( $data_type == $DATA_TYPE_ALIGN ) {
	( $staged_data_file, $format, $count, $pChromosomes ) =
      validate_bam_data_file( file => $staged_data_file, genome_chr => \%genome_chr );
}
elsif ( $data_type == $DATA_TYPE_MARKER ) {
    ( $staged_data_file, $format, $count, $pChromosomes ) =
      validate_gff_data_file( file => $staged_data_file, genome_chr => \%genome_chr );
}
if ( not $count ) {
    print $log "log: error: file contains no data\n";
    exit(-1);
}
print $log "log: Successfully read " . commify($count) . " lines\n";

# Verify that chromosome names in input file match those for genome
foreach ( sort keys %genome_chr ) {
    print $log "genome chromosome $_\n";
}
foreach ( sort keys %$pChromosomes ) {
    print $log "input chromosome $_\n";
}
if (not $ignore_missing_chr) {
	my $error = 0;
	foreach ( sort keys %$pChromosomes ) {
	    if ( not defined $genome_chr{$_} ) {
	        print $log "log: chromosome '$_' not found in genome\n";
	        $error++;
	    }
	}
	if ($error) {
	    print $log "log: error: input chromosome names don't match genome\n";
	    exit(-1);
	}
}

# Save data format doc
if ($format) {
    my $format_file = catfile($staging_dir, 'format.json');
    open(my $out, '>', $format_file);
    print $out encode_json($format);
    close($out);
}

# Generate fastbit database/index (doesn't apply to BAM files)
if ( $data_type == $DATA_TYPE_QUANT
     or $data_type == $DATA_TYPE_POLY
     or $data_type == $DATA_TYPE_MARKER )
{
    # Determine data scheme
    my $data_spec = join(',', map { $_->{name} . ':' . $_->{type} } @{$format->{columns}} );

	#TODO redirect fastbit output to log file instead of stderr
	print $log "log: Generating database\n";
	$cmd = "$FASTBIT_LOAD -d $staging_dir -m \"$data_spec\" -t $staged_data_file";
	print $log $cmd, "\n";
	my $rc = system($cmd);
	if ( $rc != 0 ) {
	    print $log "log: error executing ardea command: $rc\n";
	    exit(-1);
	}

	print $log "log: Indexing database (may take a few minutes)\n";
	$cmd =
	"$FASTBIT_QUERY -d $staging_dir -v -b \"<binning precision=2/><encoding equality/>\"";
	print $log $cmd, "\n";
	$rc = system($cmd);
	if ( $rc != 0 ) {
	    print $log "log: error executing ibis command: $rc\n";
	    exit(-1);
	}
}

################################################################################
# If we've made it this far without error then we can feel confident about
# the input data.  Now we can go ahead and create the db entities and
# install the files.
################################################################################

# Create data source
my $data_source =
  $coge->resultset('DataSource')->find_or_create( {
      name => $source_name, description => "" }
  );#, description => "Loaded into CoGe via LoadExperiment" } );
unless ($data_source) {
    print $log "log: error creating data source\n";
    exit(-1);
}

# Create experiment
my $experiment = $coge->resultset('Experiment')->create(
    {
        name        => $name,
        description => $description,
        version     => $version,
        #link				=> $link, #FIXME
        data_source_id => $data_source->id,
        data_type      => $data_type,
        row_count      => $count,
        genome_id      => $gid,
        restricted     => $restricted
    }
);
print $log "experiment id: " . $experiment->id . "\n";
print STDOUT "experiment id: " . $experiment->id . "\n"; # DON'T DELETE: needed by load_batch.pl

# Create types
if ($types) {
    foreach my $type_name ( split(/\s*;\s*/, $types) ) {
        # Try to find a matching type by name, ignoring description
        my $type = $coge->resultset('ExperimentType')->find({ name => $type_name });
        if (!$type) {
            $type = $coge->resultset('ExperimentType')->create({ name => $type_name }); # null description
        }
        unless ($type) {
            print $log "log: error creating experiment type\n";
            exit(-1);
        }
        my $conn = $coge->resultset('ExperimentTypeConnector')->find_or_create({
            experiment_id => $experiment->id,
            experiment_type_id => $type->id
        });
        unless ($conn) {
            print $log "log: error creating experiment type connector\n";
            exit(-1);
        }
    }
}

# Create annotations
if ($annotations) {
    CoGe::Core::Metadata::create_annotations(db => $coge, target => $experiment, annotations => $annotations, locked => 1);
}

# Determine installation path
unless ($install_dir) {
    unless ($P) {
        print $log "log: error: can't determine install directory, set 'install_dir' or 'config' params\n";
        exit(-1);
    }
    $install_dir = $P->{EXPDIR};
}
my $storage_path = catdir($install_dir, CoGe::Core::Storage::get_tiered_path( $experiment->id ));
print $log 'Storage path: ', $storage_path, "\n";

# This is a check for dev server which may be out-of-sync with prod
if ( -e $storage_path ) {
    print $log "log: error: install path already exists\n";
    exit(-1);
}

#TODO create experiment type & connector

# Make user owner of new experiment
my $user = $coge->resultset('User')->find( { user_name => $user_name } );
unless ($user) {
    print $log "log: error finding user '$user_name'\n";
    exit(-1);
}
my $node_types = CoGeX::node_types();
my $conn       = $coge->resultset('UserConnector')->create(
    {
        parent_id   => $user->id,
        parent_type => $node_types->{user},
        child_id    => $experiment->id,
        child_type  => $node_types->{experiment},
        role_id     => 2                            # FIXME hardcoded
    }
);
unless ($conn) {
    print $log "log: error creating user connector\n";
    exit(-1);
}

# Copy files from staging directory to installation directory
mkpath($storage_path);
unless (-r $storage_path) {
	print $log "log: error: could not create installation path\n";
	exit(-1);
}
$cmd = "cp -r $staging_dir/* $storage_path"; #FIXME use perl copy and detect failure
print $log "$cmd\n";
`$cmd`;

# Make sure file permissions are set properly (added for qTeller pipeline)
$cmd = "chmod -R a+r $storage_path";
print $log "$cmd\n";
`$cmd`;

# Save result document
if ($result_dir) {
    mkpath($result_dir);
    CoGe::Accessory::TDS::write(
        catfile($result_dir, '1'),
        {
            experiment_id => int($experiment->id)
        }
    );
}

# Yay!
CoGe::Accessory::Web::log_history(
    db          => $coge,
    user_id     => $user->id,
    page        => "LoadExperiment",
    description => 'load experiment id' . $experiment->id,
    link        => 'ExperimentView.pl?eid=' . $experiment->id
);

# Create "log.done" file to indicate completion to JEX
my $logdonefile = "$staging_dir/log.done";
touch($logdonefile);

print $log "log: All done!\n";
close($log);

exit;

#-------------------------------------------------------------------------------
sub detect_data_type {
    my $filetype = shift;
    my $filepath = shift;
    print $log "detect_data_type: $filepath\n";

    if (!$filetype or $filetype eq 'autodetect') {
        # Try to determine type based on file extension
        print $log "log: Detecting file type\n";
        ($filetype) = lc($filepath) =~ /\.([^\.]+)$/;
    }

    if ( grep { $_ eq $filetype } ('csv', 'tsv', 'bed') ) { #TODO add 'bigbed', 'wig', 'bigwig'
        print $log "log: Detected a quantitative file ($filetype)\n";
        return ($filetype, $DATA_TYPE_QUANT);
    }
    elsif ( $filetype eq 'bam' ) {
        print $log "log: Detected an alignment file ($filetype)\n";
        return ($filetype, $DATA_TYPE_ALIGN);
    }
    elsif ( $filetype eq 'vcf' ) {
        print $log "log: Detected a polymorphism file ($filetype)\n";
        return ($filetype, $DATA_TYPE_POLY);
    }
    elsif ( grep { $_ eq $filetype } ( 'gff', 'gtf' ) ) {
        print $log "log: Detected a marker file ($filetype)\n";
        return ($filetype, $DATA_TYPE_MARKER);
    }
    else {
        print $log "detect_data_type: unknown file ext '$filetype'\n";
        return ($filetype);
    }
}

# Quant file can be .csv or .bed formats
sub validate_quant_data_file {
    my %opts = @_;
    my $filepath = $opts{file};
    my $filetype = $opts{file_type};
    my $genome_chr = $opts{genome_chr};
    my %chromosomes;
    my $line_num = 1;
    my $count    = 0;
    my $hasLabels = 0;
    my $hasVal2   = 0;

    print $log "validate_quant_data_file: $filepath\n";
    open( my $in, $filepath ) || die "can't open $filepath for reading: $!";
    my $outfile = $filepath . ".processed";
    open( my $out, ">$outfile" );
    while ( my $line = <$in> ) {
        $line_num++;
        next if ( $line =~ /^\s*#/ ); # skip comments
        chomp $line;
        next unless $line; # skip blanks

        # Interpret tokens according to file type
        my @tok;
        my ( $chr, $start, $stop, $strand, $val1, $val2, $label );
        if ($filetype eq 'csv') {
        	@tok = split( /,/, $line );
        	( $chr, $start, $stop, $strand, $val1, $val2 ) = @tok;
        }
        elsif ($filetype eq 'tsv') {
        	@tok = split( /\s+/, $line );
        	( $chr, $start, $stop, $strand, $val1, $val2 ) = @tok;
        }
        elsif ($filetype eq 'bed') {
        	next if ( $line =~ /^track/ );
        	@tok = split( /\s+/, $line );
        	( $chr, $start, $stop, $label, $val1, $strand ) = @tok;
        	$val2 = $tok[6] if (@tok >= 7);
        }
        else {
        	die; # sanity check
        }

        # Validate mandatory fields
        if (   not defined $chr
            or not defined $start
            or not defined $stop
            or not defined $strand )
        {
            print $log "log: error at line $line_num: missing value in a column\n";
            return;
        }

        # mdb added 2/19/14 for bulk loading based on user request
        if ($allow_negative and $val1 < 0) {
	       $val1 = abs($val1);
        }
        # mdb added 3/13/14 issue 331
        elsif ($strand eq '.') {
            $strand = ($val1 >= 0 ? 1 : -1);
            $val1 = abs($val1);
        }

        if ( not defined $val1 or (!$disable_range_check and ($val1 < 0 or $val1 > 1)) ) {
            print $log "log: error at line $line_num: value 1 not between 0 and 1\n";
            return;
        }

		$chr = fix_chromosome_id($chr, $genome_chr);
        if (!$chr) {
            print $log "log: error at line $line_num: trouble parsing chromosome\n";
            return;
        }
        $strand = $strand =~ /-/ ? -1 : 1;

        # Build output line
        my @fields  = ( $chr, $start, $stop, $strand, $val1 ); # default fields
        if (defined $val2) {
            $hasVal2 = 1;
            push @fields, $val2;
        }
        if (defined $label) {
            $hasLabels = 1;
            push @fields, $label;
        }
        print $out join( ",", @fields ), "\n";

        $chromosomes{$chr}++;
        $count++;
    }
    close($in);
    close($out);

    #my $format = "chr:key, start:unsigned long, stop:unsigned long, strand:byte, value1:double, value2:double, label:text"; # mdb removed 4/2/14, issue 352
    # mdb added 4/2/14, issue 352
    my $format = {
        columns => [
            { name => 'chr',    type => 'key' },
            { name => 'start',  type => 'unsigned long' },
            { name => 'stop',   type => 'unsigned long' },
            { name => 'strand', type => 'byte' },
            { name => 'value1', type => 'double' }
        ]
    };
    push(@{$format->{columns}}, { name => 'value2', type => 'double' }) if $hasVal2;
    push(@{$format->{columns}}, { name => 'label',  type => 'text' }) if $hasLabels;

    return ( $outfile, $format, $count, \%chromosomes );
}

# For VCF format specification v4.1, see http://www.1000genomes.org/node/101
sub validate_vcf_data_file {
    my %opts       = @_;
    my $filepath   = $opts{file};
    my $genome_chr = $opts{genome_chr};

    my %chromosomes;
    my $line_num = 1;
    my $count    = 0;

    print $log "validate_vcf_data_file: $filepath\n";
    open( my $in, $filepath ) || die "can't open $filepath for reading: $!";
    my $outfile = $filepath . ".processed";
    open( my $out, ">$outfile" );
    while ( my $line = <$in> ) {
        $line_num++;

		#TODO load VCF metadata for storage as experiment annotations in DB (lines that begin with '##')
        next if ( $line =~ /^#/ );
        chomp $line;
        next unless $line;
        my @tok = split( /\s+/, $line );

        # Validate format
        if ( @tok < $MIN_VCF_COLUMNS ) {
            print $log "log: error at line $line_num: more columns expected ("
              . @tok
              . " < $MIN_VCF_COLUMNS)\n";
            return;
        }

        # Validate values and set defaults
        my ( $chr, $pos, $id, $ref, $alt, $qual, undef, $info ) = @tok;
        if (   not defined $chr
            || not defined $pos
            || not defined $ref
            || not defined $alt )
        {
            print $log "log: error at line $line_num: missing required value in a column\n";
            return;
        }
        next if ( $alt eq '.' );    # skip monomorphic sites
        $id   = '.' if ( not defined $id );
        $qual = 0   if ( not defined $qual );
        $info = ''  if ( not defined $info );

        $chr = fix_chromosome_id($chr, $genome_chr);
        if (!$chr) {
            print $log
              "log: error at line $line_num: trouble parsing chromosome\n";
            return;
        }
        $chromosomes{$chr}++;

        # Each line could encode multiple alleles
        my @alleles = split( ',', $alt );
        foreach my $a (@alleles) {

            # Determine site type
            my $type = detect_site_type( $ref, $a );

            # Save to file
            print $out join( ",",
                $chr, $pos, $pos + length($ref) - 1,
                $type, $id, $ref, $a, $qual, $info ),
              "\n";
            $count++;
        }
    }
    close($in);
    close($out);

    #my $format = "chr:key, start:unsigned long, stop:unsigned long, type:key, id:text, ref:key, alt:key, qual:double, info:text"; # mdb removed 4/2/14, issue 352
    # mdb added 4/2/14, issue 352
    my $format = {
        columns => [
            { name => 'chr',   type => 'key' },
            { name => 'start', type => 'unsigned long' },
            { name => 'stop',  type => 'unsigned long' },
            { name => 'type',  type => 'key' },
            { name => 'id',    type => 'text' },
            { name => 'ref',   type => 'key' },
            { name => 'alt',   type => 'key' },
            { name => 'qual',  type => 'double' },
            { name => 'info',  type => 'text' }
        ]
    };

    return ( $outfile, $format, $count, \%chromosomes );
}

sub detect_site_type {
    my $ref = shift;
    my $alt = shift;

    return 'snp' if ( length $ref == 1 and length($alt) == 1 );
    return 'deletion'  if ( length $ref > length $alt );
    return 'insertion' if ( length $ref < length $alt );
    return 'unknown';
}

sub validate_bam_data_file {
    my %opts       = @_;
    my $filepath   = $opts{file};
    my $genome_chr = $opts{genome_chr};

    my %chromosomes;
    my $count = 0;

    print $log "validate_bam_data_file: $filepath\n";

	# Get the number of reads in BAM file
	my $cmd = "$SAMTOOLS view -c $filepath";
	print $log $cmd, "\n";
    my $cmdOut = qx{$cmd};
    if ( $? != 0 ) {
	    print $log "log: error executing samtools view -c command: $?\n";
	    exit(-1);
	}
	if ($cmdOut =~ /\d+/) {
		$count = $cmdOut;
	}

	# Get the BAM file header
	$cmd = "$SAMTOOLS view -H $filepath";
	print $log $cmd, "\n";
    my @header = qx{$cmd};
    print $log "Old header:\n", @header;
    if ( $? != 0 ) {
	    print $log "log: error executing samtools view -H command: $?\n";
	    exit(-1);
	}

	# Parse the chromosome names out of the header
	my %renamed;
	foreach (@header) {
		chomp;
		if ($_ =~ /^\@SQ\s+SN\:(\S+)/) {
			my $chr = $1;
			my $newChr = fix_chromosome_id($chr, $genome_chr);
			$renamed{$chr} = $newChr if ($newChr ne $chr);
			$chromosomes{$newChr}++;
		}
	}

	# Reheader the bam file if chromosome names changed
	my $newfilepath = "$staging_dir/alignment.bam";
	if (keys %renamed) {
		# Replace chromosome names in header
		my @header2;
		foreach my $line (@header) {
			my $match = qr/^(\@SQ\s+SN\:)(\S+)/;
			if ($line =~ $match) {
				my $newChr = $renamed{$2};
				$line =~ s/$match/$1$newChr/;
			}
			push @header2, $line."\n";
		}
		print $log "New header:\n", @header2;

		# Write header to temp file
		my $header_file = "$staging_dir/header.txt";
		open( my $out, ">$header_file" );
		print $out @header2;
		close($out);

		# Run samtools to reformat the bam file header
		$cmd = "$SAMTOOLS reheader $header_file $filepath > $newfilepath";
		print $log $cmd, "\n";
	    qx{$cmd};
	    if ( $? != 0 ) {
		    print $log "log: error executing samtools reheader command: $?\n";
		    exit(-1);
		}

		# Remove the original bam file
		$cmd = "rm -f $filepath";
		qx{$cmd};
	    if ( $? != 0 ) {
		    print $log "log: error removing bam file: $?\n";
		    exit(-1);
		}
	}
	else {
		# Rename original bam file
		$cmd = "mv $filepath $newfilepath";
		qx{$cmd};
	    if ( $? != 0 ) {
		    print $log "log: error renaming bam file: $?\n";
		    exit(-1);
		}
	}

	# Index the bam file
	print $log "log: Indexing file\n";
	$cmd = "$SAMTOOLS index $newfilepath";
	print $log $cmd, "\n";
    qx{$cmd};
    if ( $? != 0 ) {
	    print $log "log: error executing samtools index command: $?\n";
	    exit(-1);
	}

    return ( $newfilepath, undef, $count, \%chromosomes );
}

# http://www.sanger.ac.uk/resources/software/gff/spec.html
sub validate_gff_data_file {
    my %opts       = @_;
    my $filepath   = $opts{file};
    my $genome_chr = $opts{genome_chr};

    my %chromosomes;
    my $line_num = 1;
    my $count    = 0;

    print $log "validate_gff_data_file: $filepath\n";
    open( my $in, $filepath ) || die "can't open $filepath for reading: $!";
    my $outfile = $filepath . ".processed";
    open( my $out, ">$outfile" );
    while ( my $line = <$in> ) {
        $line_num++;
        next if ( $line =~ /^#/ );
        chomp $line;
        next unless $line;
        my @tok = split( /\t/, $line );

        # Validate format
        if ( @tok < $MIN_GFF_COLUMNS ) {
            print $log "log: error at line $line_num: more columns expected (" . @tok . " < $MIN_GFF_COLUMNS)\n";
            return;
        }

        # Validate values and set defaults
        my ( $chr, $source, $type, $start, $stop, $score, $strand, $frame, $attr ) = @tok;
        if (   not defined $chr
            || not defined $type
            || not defined $start
            || not defined $stop )
        {
            print $log "log: error at line $line_num: missing required value in a column\n";
            return;
        }

        $chr = fix_chromosome_id($chr, $genome_chr);
        if (!$chr) {
            print $log "log: error at line $line_num: trouble parsing chromosome\n";
            return;
        }
        $chromosomes{$chr}++;

        $strand = (!defined $strand or $strand ne '-' ? 1 : -1);
        $score = 0 if (!defined $score || $score =~ /\D/);
        $attr = '' if (!defined $attr);

        print $out join(",", $chr, $start, $stop, $strand, $type, $score, $attr), "\n";
        $count++;
    }
    close($in);
    close($out);

    #my $format = "chr:key, start:unsigned long, stop:unsigned long, strand:key, type:key, score:double, attr:text"; # mdb removed 4/2/14, issue 352
    # mdb added 4/2/14, issue 352
    my $format = {
        columns => [
            { name => 'chr',    type => 'key' },
            { name => 'start',  type => 'unsigned long' },
            { name => 'stop',   type => 'unsigned long' },
            { name => 'strand', type => 'key' },
            { name => 'type',   type => 'key' },
            { name => 'score',  type => 'double' },
            { name => 'attr',   type => 'text' }
        ]
    };

    return ( $outfile, $format, $count, \%chromosomes );
}

sub fix_chromosome_id {
	my $chr = shift;
	my $genome_chr = shift;

    # Fix chromosome identifier
    $chr =~ s/^lcl\|//;
    $chr =~ s/chromosome//i;
    $chr =~ s/^chr//i;
    $chr =~ s/^0+//;
    $chr =~ s/^_+//;
    $chr =~ s/\s+/ /;
    $chr =~ s/^\s//;
    $chr =~ s/\s$//;

	# Hack to deal with converting 'chloroplast' and 'mitochondia' to 'C' and 'M' if needed
    if (   $chr =~ /^chloroplast$/i
        && !$genome_chr->{$chr}
        && $genome_chr->{"C"} )
    {
        $chr = "C";
    }
    if (   $chr =~ /^mitochondria$/i
        && !$genome_chr->{$chr}
        && $genome_chr->{"M"} )
    {
        $chr = "M";
    }

	return $chr;
}
