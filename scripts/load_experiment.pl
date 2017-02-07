#!/usr/bin/perl -w

use strict;
use Data::Dumper;
use CoGeX;
use Getopt::Long;
use File::Path;
use File::Touch;
use File::Basename qw(basename);
use File::Spec::Functions qw(catdir catfile);
use URI::Escape::JavaScript qw(unescape);
use JSON::XS;
use CoGe::Accessory::Web qw(get_defaults get_command_path);
use CoGe::Accessory::Utils qw( commify to_pathname to_number );
use CoGe::Accessory::TDS;
use CoGe::Core::Genome qw(fix_chromosome_id);
use CoGe::Core::Storage qw(add_workflow_result $DATA_TYPE_QUANT $DATA_TYPE_ALIGN $DATA_TYPE_POLY $DATA_TYPE_MARKER);
use CoGe::Core::Experiment qw(@SUPPORTED_TYPES detect_data_type);
use CoGe::Core::Metadata qw(create_annotations);

use vars qw($staging_dir $install_dir $data_file $file_type $metadata_file
  $name $description $version $restricted $ignore_missing_chr $creator_id $normalize
  $gid $source_name $user_name $config $allow_negative $disable_range_check
  $exit_without_error_for_empty_input $link
  $user_id $annotations $tags $wid $host $port $dbname $user $pass $P);

#my $MIN_QUANT_COLUMNS = 5;
#my $MAX_QUANT_COLUMNS = 6;
my $MIN_VCF_COLUMNS = 8;
#my $MAX_VCF_COLUMNS = 10;
my $MIN_GFF_COLUMNS = 9;

use constant {
    RC_PARSE_ERROR     => 1,
    RC_PARSE_SKIP_LINE => 2
};

GetOptions(
    "staging_dir=s" => \$staging_dir,    # temporary staging path
    "install_dir=s" => \$install_dir,    # final installation path
#    "result_file=s" => \$result_file,    # results file
    "data_file=s"   => \$data_file,      # input data file (JS escape)
    "metadata_file=s" => \$metadata_file,# input metadata file (output of Data::Dumper)
    "file_type=s"   => \$file_type,		 # input file type
    "name=s"        => \$name,           # experiment name (JS escaped)
    "desc=s"        => \$description,    # experiment description (JS escaped)
    "version=s"     => \$version,        # experiment version (JS escaped)
    "restricted=s"  => \$restricted,     # experiment restricted flag (0|1 or false|true)
    "source_name=s" => \$source_name,    # experiment source name (JS escaped)
    "link=s"        => \$link,           # link (JS escaped)
    "gid=s"         => \$gid,            # genome id
    "wid=s"         => \$wid,            # workflow id
    "user_id=i"     => \$user_id,        # user ID to assign experiment
    "user_name=s"   => \$user_name,      # user name to assign experiment (alternative to user_id)
    "creator_id=i"  => \$creator_id,     # user ID to set as experiment creator
    "annotations=s" => \$annotations,    # optional: semicolon-separated list of locked annotations (link:group:type:text;...)
    "tags=s"        => \$tags,           # optional: semicolon-separated list of experiment tag names
    "normalize=s"   => \$normalize,      # optional: percentage, log10 or loge    
    "config=s"      => \$config,         # configuration file

    # Optional flags (mostly for debug and bulk loader)
    "ignore-missing-chr=i" => \$ignore_missing_chr,
    "allow_negative=i"     => \$allow_negative,
    "disable_range_check"  => \$disable_range_check, # allow any value in val1 column
    "exit_without_error_for_empty_input" => \$exit_without_error_for_empty_input, # added for SNP analysis pipeline where no SNPs found
);

$| = 1;
print STDOUT "Starting $0 (pid $$)\n", qx/ps -o args $$/;

# Setup staging path
unless ($staging_dir) {
    print STDOUT "log: error: staging_dir argument is missing\n";
    exit(-1);
}
mkpath($staging_dir, 0, 0777) unless -r $staging_dir;

# Prevent loading again (issue #417)
my $logdonefile = "$staging_dir/log.done";
if (-e $logdonefile) {
    print STDOUT "log: error: done file already exists: $logdonefile\n";
    exit(-1);
}

# Process and verify parameters
$data_file   = unescape($data_file);
$name        = unescape($name);
$description = unescape($description);
$version     = unescape($version);
$source_name = unescape($source_name);
$link        = unescape($link) if ($link);

unless ($wid) {
    print STDOUT "log: error: required workflow ID not specified\n";
    exit(-1);
}

unless ($data_file && -r $data_file) {
    print STDOUT "log: error: cannot access input data file\n";
    exit(-1);
}

if (not defined $user_id and not defined $user_name) {
    print STDOUT "log: error: user not specified, use user_id or user_name\n";
    exit(-1);
}

if ((defined $user_name and $user_name eq 'public') || (defined $user_id and $user_id eq '0')) {
    print STDOUT "log: error: not logged in\n";
    exit(-1);
}

# Set default parameters
$restricted  = '1' unless (defined $restricted && (lc($restricted) eq 'false' || $restricted eq '0'));
$ignore_missing_chr = '1' unless (defined $ignore_missing_chr); # mdb added 10/6/14 easier just to make this the default

# Load config file
unless ($config) {
    print STDOUT "log: error: can't find config file\n";
    print STDERR "can't find config file\n";
    exit(-1);
}
$P    = CoGe::Accessory::Web::get_defaults($config);
$dbname = $P->{DBNAME};
$host = $P->{DBHOST};
$port = $P->{DBPORT};
$user = $P->{DBUSER};
$pass = $P->{DBPASS};

my $FASTBIT_LOAD  = get_command_path('FASTBIT_LOAD', 'ardea');
my $FASTBIT_QUERY = get_command_path('FASTBIT_QUERY', 'ibis');
my $SAMTOOLS      = get_command_path('SAMTOOLS');
my $GUNZIP        = get_command_path('GUNZIP');
if (   not $FASTBIT_LOAD
    or not $FASTBIT_QUERY
    or not $SAMTOOLS
    or not $GUNZIP )
{
    print STDOUT "log: error: can't find required command(s)\n";
    exit(-1);
}

# Copy input data file to staging area
# If running via JEX the file will already be there
my ($filename) = basename($data_file);
my $staged_data_file = $staging_dir . '/' . $filename;
unless (-r $staged_data_file) {
    my $cmd;
    $cmd = "cp -f '$data_file' $staging_dir";
    `$cmd`;
}

# Decompress file if necessary
if ( $staged_data_file =~ /\.gz$/ ) {
    my $cmd = $GUNZIP . ' ' . $staged_data_file;
    print STDOUT "log: Decompressing '$filename'\n";
    `$cmd`;
    $staged_data_file =~ s/\.gz$//;
}

# Determine file type
my ($file_type, $data_type) = detect_data_type($file_type, $staged_data_file);
if ( !$file_type or !$data_type ) {
    my $types = join ",", sort {$a cmp $b} @SUPPORTED_TYPES;
    print STDOUT "log: error: unknown or unsupported file type '$file_type'\n";
    print STDOUT "log: file must end with one of the following types: $types\n";
    exit(-1);
}

# Connect to database
my $connstr = "dbi:mysql:dbname=$dbname;host=$host;port=$port;";
my $db = CoGeX->connect( $connstr, $user, $pass );
unless ($db) {
    print STDOUT "log: couldn't connect to database\n";
    exit(-1);
}

# Retrieve user
my $user;
if ($user_id) {
    $user = $db->resultset('User')->find($user_id);
}
elsif ($user_name) {
    $user = $db->resultset('User')->find( { user_name => $user_name } );
}
else {
    print STDOUT "log: error user not specified, see user_id or user_name\n";
    exit(-1);
}

unless ($user) {
    print STDOUT "log: error finding user ", ($user_name ? $user_name : $user_id) , "\n";
    exit(-1);
}

# Retrieve creator
my $creator;
if ($creator_id) {
    $creator = $db->resultset('User')->find($creator_id);
    unless ($creator) {
        print STDOUT "log: error finding creator $creator_id\n";
        exit(-1);
    }
}
$creator = $user unless $creator;

# Retrieve genome
my $genome = $db->resultset('Genome')->find( { genome_id => $gid } );
unless ($genome) {
    print STDOUT "log: error finding genome id$gid\n";
    exit(-1);
}

# Hash chromosome names
my %genome_chr = map { $_ => 1 } $genome->chromosomes;

# Validate the data file
print STDOUT "log: Validating data file\n";
if (-s $staged_data_file == 0) {
    print STDOUT "log: error: input file '", basename($staged_data_file), "' is empty\n";
    exit(-1);
}
my ($count, $pChromosomes, $format, $detected_tags, $metadata);
if ( $data_type == $DATA_TYPE_QUANT ) {
    ( $staged_data_file, $format, $count, $pChromosomes, $metadata ) =
      validate_quant_data_file( file => $staged_data_file, file_type => $file_type, genome_chr => \%genome_chr );
}
elsif ( $data_type == $DATA_TYPE_POLY ) {
    ( $staged_data_file, $format, $count, $pChromosomes, $detected_tags ) =
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
if ( !$count ) {
    print STDOUT "log: error: file contains no data\n";
    if ($exit_without_error_for_empty_input) {
        touch($logdonefile);
        print STDOUT "Exited early due to no data and exit_without_error_for_empty_input\n";
        exit(0);
    }
    exit(-1);
}
print STDOUT "log: Successfully read " . commify($count) . " lines\n";

# Verify that chromosome names in input file match those for genome
my $print_limit = 50;
foreach ( sort keys %genome_chr ) {
    print STDOUT "genome chromosome $_\n";
    if ($print_limit-- == 0) {
        print STDOUT "... (stopping here, too many genome chromosomes to show)\n";
        last;
    }
}
$print_limit = 50;
foreach ( sort keys %$pChromosomes ) {
    print STDOUT "input chromosome $_\n";
    if ($print_limit-- == 0) {
        print STDOUT "... stopping here, too many input chromosomes to show\n";
        last;
    }
}
	
my $missing_chr_error = 0;
foreach ( sort keys %$pChromosomes ) {
	if ( not defined $genome_chr{$_} ) { # don't repeat same error message
		if ($missing_chr_error < 5) {
			print STDOUT "log: chromosome '$_' not found in genome, skipping (only showing first 5) ...\n";
		}
	        $missing_chr_error++;
	}
}

if (not $ignore_missing_chr) {
	if ($missing_chr_error) {
	    print STDOUT "log: error: input chromosome names don't match genome\n";
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
	print STDOUT "log: Generating database\n";
	my $cmd = "$FASTBIT_LOAD -d $staging_dir -m \"$data_spec\" -t $staged_data_file";
	print STDOUT $cmd, "\n";
	my $rc = system($cmd);
	if ( $rc != 0 ) {
	    print STDOUT "log: error executing ardea command: $rc\n";
	    exit(-1);
	}

	print STDOUT "log: Indexing database (may take a few minutes)\n";
	$cmd = "$FASTBIT_QUERY -d $staging_dir -v -b \"<binning precision=2/><encoding equality/>\"";
	print STDOUT $cmd, "\n";
	$rc = system($cmd);
	if ( $rc != 0 ) {
	    print STDOUT "log: error executing ibis command: $rc\n";
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
  $db->resultset('DataSource')->find_or_create( {
      name => $source_name, description => "" }
  );#, description => "Loaded into CoGe via LoadExperiment" } );
unless ($data_source) {
    print STDOUT "log: error creating data source\n";
    exit(-1);
}

# Create experiment
my $experiment = $db->resultset('Experiment')->create(
    {
        name        => $name,
        description => $description,
        version     => $version,
        link		=> $link,
        data_source_id => int($data_source->id),
        data_type   => int($data_type),
        row_count   => int($count),
        genome_id   => int($gid),
        creator_id  => int($creator->id),
        restricted  => $restricted
    }
);
print STDOUT "experiment id: " . $experiment->id . "\n";

# Create tags
if ($tags || $detected_tags) {
    foreach my $name ( @$detected_tags, split(/\s*;\s*/, $tags) ) {
        next unless $name;
        
        # Try to find a matching type by name, ignoring description
        my $type = $db->resultset('ExperimentType')->find({ name => $name });
        if (!$type) {
            $type = $db->resultset('ExperimentType')->create({ name => $name }); # null description
        }
        unless ($type) {
            print STDOUT "log: error creating experiment type\n";
            exit(-1);
        }
        my $conn = $db->resultset('ExperimentTypeConnector')->find_or_create({
            experiment_id => $experiment->id,
            experiment_type_id => $type->id
        });
        unless ($conn) {
            print STDOUT "log: error creating experiment type connector\n";
            exit(-1);
        }
    }
}

# Create annotations
if ($annotations || $metadata_file) {
    CoGe::Core::Metadata::create_annotations(
        db => $db,
        user => $user,
        target => $experiment,
        annotations => $annotations,
        anno_file => $metadata_file,
        locked => 1
    );
}

# Determine installation path
unless ($install_dir) {
    unless ($P) {
        print STDOUT "log: error: can't determine install directory, set 'install_dir' or 'config' params\n";
        exit(-1);
    }
    $install_dir = $P->{EXPDIR};
}
my $storage_path = catdir($install_dir, CoGe::Core::Storage::get_tiered_path( $experiment->id ));
print STDOUT 'Storage path: ', $storage_path, "\n";

# This is a check for dev server which may be out-of-sync with prod
if ( -e $storage_path ) {
    print STDOUT "log: error: install path already exists\n";
    exit(-1);
}

#TODO create experiment type & connector

my $node_types = CoGeX::node_types();
my $conn       = $db->resultset('UserConnector')->create(
    {
        parent_id   => $user->id,
        parent_type => $node_types->{user},
        child_id    => $experiment->id,
        child_type  => $node_types->{experiment},
        role_id     => 2                            # FIXME hardcoded
    }
);
unless ($conn) {
    print STDOUT "log: error creating user connector\n";
    exit(-1);
}

# Copy files from staging directory to installation directory
mkpath($storage_path);
unless (-r $storage_path) {
	print STDOUT "log: error: could not create installation path\n";
	exit(-1);
}
my $cmd = "cp -r $staging_dir/* $storage_path"; #FIXME use perl copy and detect failure
print STDOUT "$cmd\n";
`$cmd`;

# Make sure file permissions are set properly (added for Cufflinks pipeline)
$cmd = "chmod -R a+r $storage_path";
print STDOUT "$cmd\n";
`$cmd`;

# Save result
unless (add_workflow_result($user_name, $wid, 
        {
            type        => 'experiment',
            id          => int($experiment->id),
            name        => $name,
            description => $description,
            version     => $version,
            #link       => $link, #FIXME
            source_id   => int($data_source->id),
            data_type   => int($data_type), #FIXME convert from number to string identifier
            row_count   => int($count),
            genome_id   => int($gid),
            restricted  => $restricted
        })
    )
{
    print STDOUT "log: error: could not add workflow result\n";
    exit(-1);
}

# Add experiment ID to log - mdb added 8/19/14, needed after log output was moved to STDOUT for jex
my $logtxtfile = "$staging_dir/log.txt";
open(my $logh, '>', $logtxtfile);
print $logh "experiment id: " . $experiment->id . "\n";
close($logh);

# Save workflow_id in metadata.json file in experiment data path -- #TODO move into own routine in Storage.pm
$metadata = {} unless defined $metadata;
$metadata->{workflow_id} = int($wid);
CoGe::Accessory::TDS::write( catfile($storage_path, 'metadata.json'), $metadata );

# Create "log.done" file to indicate completion to JEX
touch($logdonefile);
print STDOUT "All done!\n";

exit;

#-------------------------------------------------------------------------------

# Parse entire quant data file and determine max and min values
sub max_and_min_of_values {
	my $filepath = shift;
	my $filetype = shift;
	my $line_num = 0;
	my ($max, $min) = (0, 0);
#	my %perChr;
	
    open( my $in, $filepath ) || die "can't open $filepath for reading: $!";
    <$in> if $filetype eq 'seg'; # skip first line of column headings
    while ( my $line = <$in> ) {
        $line_num++;
        next if ( $line =~ /^\s*#/ ); # skip comment lines
        chomp $line;
        next unless $line; # skip blank lines

        my ($status, $chr, undef, undef, $strand, $val1) = parse_quant_line($file_type, $line, $line_num);
        next if ($status == RC_PARSE_SKIP_LINE);
        return if ($status == RC_PARSE_ERROR);
        return unless (defined $chr && defined $strand && defined $val1); # error, let calling function handle it

        $strand = $strand =~ /-/ ? -1 : 1;
        $val1 *= $strand if ($strand == -1);
        
        $max = $val1 if ($val1 > $max);
        $min = $val1 if ($val1 < $min);
        
#        $perChr{$chr}{max} = $max if (!defined($perChr{$chr}{max}) || $val1 > $perChr{$chr}{max});
#        $perChr{$chr}{min} = $min if (!defined($perChr{$chr}{min}) || $val1 < $perChr{$chr}{min});
    }
    close($in);
    print STDOUT "max_and_min_of_values: max=$max min=$min\n";
    return ($max, $min);
}

# Parses multiple line-based file formats for quant data
my $bedType;              # only used for BED formats
my ($stepSpan, $stepChr); # only used for WIG format
sub validate_quant_data_file {
    my %opts = @_;
    my $filepath = $opts{file};
    my $filetype = $opts{file_type};
    my $genome_chr = $opts{genome_chr};
    my %chromosomes;
    my $line_num = 0;
    my $count;
    my $hasLabels = 0;
    my $hasVal2   = 0;
    my $md;

    print STDOUT "validate_quant_data_file: $filepath\n";
    
    # Get min/max values for normalization
    my ($max, $min) = max_and_min_of_values($filepath, $filetype);
    
    # Parse data file line by line
    open( my $in, $filepath ) || die "can't open $filepath for reading: $!";
    my $outfile = $filepath . ".processed";
    open( my $out, ">$outfile" );
    <$in> if $filetype eq 'seg'; # skip first line of column headings
    while ( my $line = <$in> ) {
        $line_num++;
        next if ( $line =~ /^\s*#/ ); # skip comment lines
        chomp $line;
        next unless $line; # skip blank lines
        
        my ($status, $chr, $start, $stop, $strand, $val1, $val2, $label) = parse_quant_line($file_type, $line);
        next if ($status == RC_PARSE_SKIP_LINE);
        return if ($status == RC_PARSE_ERROR);

        # Validate mandatory fields
        if (   not defined $chr
            or not defined $start
            or not defined $stop
            or not defined $strand )
        {
            my $missing;
            $missing = 'chr'    unless $chr;
            $missing = 'start'  unless $start;
            $missing = 'stop'   unless $stop;
            $missing = 'strand' unless $strand;
            log_line("missing value in a column: $missing", $line_num, $line);
            return;
        }

    	# mdb added 2/10/16
    	if ($start =~ /[^\d]/ || $stop =~ /[^\d]/) {
    		log_line('start/stop not integer values', $line_num, $line);
    		return;
    	}

        # mdb added 2/19/14 for bulk loading based on user request
        if ($allow_negative and $val1 < 0) {
	       $val1 = abs($val1);
        }
        # mdb added 3/13/14 issue 331 - set strand based on polarity of value
        elsif ($strand eq '.') {
            $strand = ($val1 >= 0 ? 1 : -1);
            $val1 = abs($val1);
        }
        
        if (!$normalize) {
        	if (not defined $val1 or (!$disable_range_check and ($val1 < 0 or $val1 > 1))) {
	            log_line('value 1 not between 0 and 1', $line_num, $line);
    	        return;
        	}
        }

        # Munge chr name for CoGe
        ($chr) = split(/\s+/, $chr);
		$chr = fix_chromosome_id($chr, $genome_chr);
        unless (defined $chr) {
            log_line("trouble parsing sequence ID", $line_num, $line);
            return;
        }
        $strand = $strand =~ /-/ ? -1 : 1;

        # Build output line
        if ($normalize) {
	        if ($normalize eq 'percentage') {
	        	$val1 /= $max;
	        }
	        elsif ($normalize eq 'log10') {
	        	$val1 = log($val1) / log(10) / $max;
	        }
	        elsif ($normalize eq 'loge') {
	        	$val1 = log($val1) / $max;
	        }
	        else {
	            print STDOUT "log: error: unknown normalization method given: '$normalize'\n";
                return;
	        }
	        
	        # Add min/max values to metadata object -- mdb added 8/31/16 COGE-270
	        $md = {
                max => 1,
                min => -1,
                normalization => $normalize
            };
        }
        else {
            # Add min/max values to metadata object -- mdb added 8/31/16 COGE-270
            $md = {
                max => to_number($max),
                min => to_number($min),
                normalization => undef
            };
        }
        
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

        # Keep track of seen chromosome names for later use
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

    return ( $outfile, $format, $count, \%chromosomes, $md );
}

sub parse_quant_line {
    my ($filetype, $line, $line_num) = @_;
    
    # Interpret tokens according to file type
    my @tok;
    my ( $chr, $start, $stop, $strand, $val1, $val2, $label );
    
    if ($filetype eq 'csv') { # CoGe format, comma-separated
        @tok = split( /,/, $line );
        ( $chr, $start, $stop, $strand, $val1, $val2 ) = @tok;
    }
    elsif ($filetype eq 'tsv') { # CoGe format, tab-separated
        @tok = split( /\s+/, $line );
        ( $chr, $start, $stop, $strand, $val1, $val2 ) = @tok;
    }
    elsif ($filetype eq 'wig') {
        if ( $line =~ /^track/ ) { # ignore "track" line
            return RC_PARSE_SKIP_LINE;
        }
        if ( $line =~ /^variableStep/i ) { # handle step definition line
            if ($line =~ /chrom=(\w+)/i) {
                $stepChr = $1;
            }
            
            $stepSpan = 1;
            if ($line =~ /span=(\d+)/i) {
                $stepSpan = $1;
            }
            return RC_PARSE_SKIP_LINE;
        }
        elsif ( $line =~ /^fixedStep/i ) {
            log_line('fixedStep wiggle format is not currently supported', $line_num, $line);
            return RC_PARSE_ERROR;
        }
        
        if (not defined $stepSpan or not defined $stepChr) {
            log_line('missing or invalid wiggle step definition line', $line_num, $line);
            return RC_PARSE_ERROR;
        }
        
        @tok = split( /\s+/, $line );
        ( $start, $val1 ) = @tok;
        $stop = $start + $stepSpan - 1;
        $chr = $stepChr;
        $strand = '.'; # determine strand by val1 polarity   
    }
    elsif ($filetype eq 'bed') {
        # Check for track type for BED files
        if ( $line =~ /^track/ ) {
            undef $bedType;
            if ($line =~ /type=(\w+)/i) {
                $bedType = lc($1);
            }
            return RC_PARSE_SKIP_LINE;
        }
    
        # Handle different BED formats
        @tok = split( /\s+/, $line );
        if (defined $bedType && $bedType eq 'bedgraph') { # UCSC bedGraph: http://genome.ucsc.edu/goldenPath/help/bedgraph.html
            ( $chr, $start, $stop, $val1 ) = @tok;
            $strand = '.'; # determine strand by val1 polarity
        }
        else { # UCSC standard BED: http://genome.ucsc.edu/FAQ/FAQformat.html#format1
            ( $chr, $start, $stop, $label, $val1, $strand ) = @tok;
            $val2 = $tok[6] if (@tok >= 7); # non-standard CoGe usage
        }
        
        # Adjust coordinates from base-0 to base-1
        if (defined $start and defined $stop) {
            $start += 1;
            $stop += 1;
        }
    }
    elsif ($filetype eq 'seg') {
        @tok = split( /\s+/, $line );
        ( undef, $chr, $start, $stop, undef, $val1 ) = @tok;
        $strand = 1;
    }
    else { # unknown file type (should never happen)
        die "fatal error: unknown file type!";
    }
    
    return (0, $chr, $start, $stop, $strand, $val1, $val2, $label);
}

# For VCF format specification v4.1, see http://www.1000genomes.org/node/101
sub validate_vcf_data_file {
    my %opts       = @_;
    my $filepath   = $opts{file};
    my $genome_chr = $opts{genome_chr};

    my %chromosomes;
    my $isGVCF = 0;
    my $line_num = 1;
    my $count;

    print STDOUT "validate_vcf_data_file: $filepath\n";
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
            log_line('more columns expected ('.@tok.' < '.$MIN_VCF_COLUMNS.')', $line_num, $line);
            return;
        }
        
        # Check for extra genotype columns in GVCF # mdb added 5/4/16
        if ( !$isGVCF && @tok > 10 ) { 
            print STDOUT "Detected multisample GVCF file based on ", scalar(@tok), " columns\n";
            $isGVCF = 1;
        }

        # Validate values and set defaults
        my ( $chr, $pos, $id, $ref, $alt, $qual, undef, $info ) = @tok;
        if (   not defined $chr
            || not defined $pos
            || not defined $ref
            || not defined $alt )
        {
            log_line('missing required value in a column', $line_num, $line);
            return;
        }
        next if ( $alt eq '.' );    # skip monomorphic sites
        $id   = '.' if ( not defined $id );
        $qual = 0   if ( not defined $qual );
        $info = ''  if ( not defined $info );

        $chr = fix_chromosome_id($chr, $genome_chr);
        if (!$chr) {
            log_line('trouble parsing chromosome', $line_num, $line);
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
            { name => 'qual',  type => 'double' }
            #{ name => 'info',  type => 'text' }
        ]
    };
    
    # Add detected tags
    my @tags = ( $isGVCF ? 'Multisample-GVCF' : 'VCF' );

    return ( $outfile, $format, $count, \%chromosomes, \@tags );
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
    my $count;

    print STDOUT "validate_bam_data_file: $filepath\n";

	# Get the number of reads in BAM file
	my $cmd = "$SAMTOOLS view -c $filepath";
	print STDOUT $cmd, "\n";
    my $cmdOut = qx{$cmd};
    if ( $? != 0 ) {
	    print STDOUT "log: error executing samtools view -c command: $?\n";
	    exit(-1);
	}
	if ($cmdOut =~ /\d+/) {
		$count = $cmdOut;
	}

	# Get the BAM file header
	$cmd = "$SAMTOOLS view -H $filepath";
	print STDOUT $cmd, "\n";
    my @header = qx{$cmd};
    #print STDOUT "Old header:\n", @header;
    execute($cmd);

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
				$line =~ s/$match/$1$newChr/ if defined $newChr;
			}
			push @header2, $line."\n";
		}
		print STDOUT "New header:\n", @header2;

		# Write header to temp file
		my $header_file = "$staging_dir/header.txt";
		open( my $out, ">$header_file" );
		print $out @header2;
		close($out);

		# Run samtools to reformat the bam file header
		$cmd = "$SAMTOOLS reheader $header_file $filepath > $newfilepath";
		execute($cmd);

		# Remove the original bam file
		$cmd = "rm -f $filepath";
		execute($cmd);
	}
	elsif ($filepath ne $newfilepath) { # mdb added condition 3/12/15 -- possible that original file is named "alignment.bam"
		# Rename original bam file
		$cmd = "mv $filepath $newfilepath";
		execute($cmd);
	}
	
	# Sort the bam file
	# TODO this can be slow, is it possible to detect if it is sorted already?
	my $sorted_file = "$staging_dir/sorted.bam";
    $cmd = "$SAMTOOLS sort $newfilepath -o $sorted_file"; # mdb changed 1/5/17 -- added -o for SAMtools 1.3.1
    execute($cmd);
    if (-e $sorted_file && -s $sorted_file > 0) {
        # Replace original file with sorted version
        execute("mv $sorted_file $newfilepath");
    }
    else {
        print STDOUT "log: error: samtools sort produced no result\n";
        exit(-1);
    }

	# Index the bam file
	$cmd = "$SAMTOOLS index $newfilepath";
	execute($cmd);

    return ( $newfilepath, undef, $count, \%chromosomes );
}

# http://www.sanger.ac.uk/resources/software/gff/spec.html
sub validate_gff_data_file {
    my %opts       = @_;
    my $filepath   = $opts{file};
    my $genome_chr = $opts{genome_chr};

    my %chromosomes;
    my $line_num = 1;
    my $count;

    print STDOUT "validate_gff_data_file: $filepath\n";
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
            log_line("more columns expected (" . @tok . " < $MIN_GFF_COLUMNS)", $line_num, $line);
            return;
        }

        # Validate values and set defaults
        my ( $chr, $source, $type, $start, $stop, $score, $strand, $frame, $attr ) = @tok;
        if (   not defined $chr
            || not defined $type
            || not defined $start
            || not defined $stop )
        {
            log_line('missing required value in a column', $line_num, $line);
            return;
        }

        $chr = fix_chromosome_id($chr, $genome_chr);
        if (!$chr) {
            log_line('trouble parsing chromosome', $line_num, $line);
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

sub execute { # FIXME move into Util.pm
    my $cmd = shift;
    print STDOUT "$cmd\n";
    my @cmdOut    = qx{$cmd};
    my $cmdStatus = $?;
    if ( $cmdStatus != 0 ) {
        print STDOUT "log: error: command failed with rc=$cmdStatus: $cmd\n";
        exit(-1);
    }
}

sub log_line {
    my ( $msg, $line_num, $line ) = @_;
    print STDOUT "log: error at line $line_num: $msg\n", "log: ", substr($line, 0, 100), "\n";    
}
