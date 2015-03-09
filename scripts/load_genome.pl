#!/usr/bin/perl -w

use strict;
use CoGeX;
use CoGe::Accessory::Web;
use CoGe::Core::Storage qw( index_genome_file get_tiered_path );
use CoGe::Core::Genome qw( fix_chromosome_id );
use CoGe::Accessory::Utils qw( commify units print_fasta );
use CoGe::Accessory::IRODS qw( irods_imeta $IRODS_METADATA_PREFIX );
use Data::Dumper;
use Getopt::Long;
use File::Path;
use File::Touch;
use File::Basename qw( basename dirname );
use File::Spec::Functions qw( catdir catfile );
use URI::Escape::JavaScript qw(unescape);
use POSIX qw(ceil);
use Benchmark;

use vars qw($staging_dir $install_dir $fasta_files $irods_files
  $name $description $link $version $type_id $restricted $message
  $organism_id $source_id $source_name $source_desc $user_id $user_name
  $keep_headers $split $compress $result_dir $creator_id $ignore_chr_limit
  $host $port $db $user $pass $config
  $P $MAX_CHROMOSOMES $MAX_PRINT $MAX_SEQUENCE_SIZE $MAX_CHR_NAME_LENGTH );

$MAX_CHROMOSOMES     = 200*1000;    # max number of chromosomes or contigs
$MAX_PRINT           = 50;
$MAX_SEQUENCE_SIZE   = 5 * 1024 * 1024 * 1024;    # 5 gig
$MAX_CHR_NAME_LENGTH = 255;

GetOptions(
    "staging_dir=s" => \$staging_dir,
    "install_dir=s" => \$install_dir,    # optional, for debug
    "result_dir=s"  => \$result_dir,     # results path
    "fasta_files=s" => \$fasta_files,    # comma-separated list (JS escaped) of files to load
    "irods_files=s" => \$irods_files,    # optional comma-separated list (JS escaped) of files to set metadata
    "name=s"        => \$name,           # genome name (JS escaped)
    "desc=s"        => \$description,    # genome description (JS escaped)
    "message=s"		=> \$message,		 # message (JS escaped)
    "link=s"        => \$link,           # link (JS escaped)
    "version=s"     => \$version,        # genome version (JS escaped)
    "type_id=i"     => \$type_id,        # genomic_sequence_type_id
    "restricted=i"  => \$restricted,     # genome restricted flag
    "organism_id=i" => \$organism_id,    # organism ID
    "source_id=i"	=> \$source_id,		 # data source id
    "source_name=s" => \$source_name,    # data source name (JS escaped)
    "source_desc=s" => \$source_desc,    # data source description (JS escaped)
    "user_id=i"		=> \$user_id,		 # user ID
    "creator_id=i"	=> \$creator_id,     # user ID to set genome creator
    "user_name=s"   => \$user_name,      # user name
    "keep_headers=i" => \$keep_headers,  # flag to keep original headers (no parsing)
    "ignore_chr_limit=i" => \$ignore_chr_limit, # flag to ignore chromosome/contig limit # mdb added 3/9/15 COGE-595
    "split=i"    => \$split,       # split fasta into chr directory
    "compress=i" => \$compress,    # compress fasta into RAZF before indexing
    "config=s"   => \$config       # configuration file
);

# Open log file
$| = 1;
die unless ($staging_dir);
mkpath($staging_dir); # make sure this exists
#my $logfile = "$staging_dir/log.txt";
#open( my $log, ">>$logfile" ) or die "Error opening log file: $logfile: $!";
#$log->autoflush(1);
print STDOUT "Starting $0 (pid $$)\n", qx/ps -o args $$/;

# Prevent loading again (issue #417)
my $logdonefile = "$staging_dir/log.done";
if (-e $logdonefile) {
    print STDOUT "log: error: done file already exists: $logdonefile\n";
    exit(-1);
}

# Process and verify parameters
$fasta_files = unescape($fasta_files) if ($fasta_files);
$irods_files = unescape($irods_files) if ($irods_files);
$name        = unescape($name) if ($name);
$description = unescape($description) if ($description);
$link        = unescape($link) if ($link);
$version     = unescape($version) if ($version);
$message     = unescape($message) if ($message);
$source_name = unescape($source_name) if ($source_name);
$source_desc = unescape($source_desc) if ($source_desc);
$restricted  = '0' if ( not defined $restricted );
$split       = 0 if ( not defined $split ); 	# split fasta into chr/ files (legacy method)
$compress    = 0 if ( not defined $compress ); 	# RAZF compress the fasta file
$type_id     = 1 if ( not defined $type_id );   # default genomic seq type to "unmasked"

if (not $source_id and not $source_name) {
	print STDOUT "log: error: source not specified, use source_id or source_name\n";
	exit(-1);
}
if (not $user_id and not $user_name) {
	print STDOUT "log: error: user not specified, use user_id or user_name\n";
	exit(-1);
}
if ($user_name and $user_name eq 'public') {
	print STDOUT "log: error: not logged in\n";
    exit(-1);
}
unless ($organism_id) {
    print STDOUT "log: error: organism_id not specified\n";
    exit(-1);
}

# Load config file
unless ($config) {
    print STDOUT "log: error: can't find config file\n";
    print STDERR "can't find config file\n";
    exit(-1);
}
$P    = CoGe::Accessory::Web::get_defaults($config);
$db   = $P->{DBNAME};
$host = $P->{DBHOST};
$port = $P->{DBPORT};
$user = $P->{DBUSER};
$pass = $P->{DBPASS};
my $GUNZIP = $P->{GUNZIP};
my $TAR = $P->{TAR};

# Process each file into staging area
my %sequences;
my $seqLength;
my $numSequences;
my @files = split( ',', $fasta_files );
foreach my $file (@files) {
    my $filename = basename($file);# =~ /^.+\/([^\/]+)$/;
    my $path = dirname($file);

    # Decompress file if necessary
    if ( $file =~ /\.tgz|\.tar\.gz$/ ) {
        print STDOUT "log: Unarchiving/decompressing '$filename'\n";
        my $orig = $file;
        my $filelist = execute( $TAR . ' -xvf ' . $file ); # mdb added 7/31/14 issue 438
        chomp(@$filelist);
        push @files, map { catfile($path, $_) } @$filelist;
        next;
    }
    if ( $file =~ /\.gz$/ ) {
        print STDOUT "log: Decompressing '$filename'\n";
        #execute($GUNZIP . ' ' . $file); # mdb removed 7/31/14 issue 438
        $file =~ s/\.gz$//;
        execute( $GUNZIP . ' -c ' . $file . '.gz' . ' > ' . $file ); # mdb added 7/31/14 issue 438
    }

    # Ensure text file
    if ( -B $file ) {
        print STDOUT "log: error: '$filename' is a binary file\n";
        exit(-1);
    }

    # Load file
    $seqLength += process_fasta_file( \%sequences, $file, $staging_dir );
    $numSequences = keys %sequences;

    if ( $seqLength > $MAX_SEQUENCE_SIZE ) {
        print STDOUT "log: error: total sequence size exceeds limit of "
          . units($MAX_SEQUENCE_SIZE) . "\n";
        exit(-1);
    }
    if ( !$ignore_chr_limit && $numSequences > $MAX_CHROMOSOMES ) {
        print STDOUT
          "log: error: too many sequences, limit is $MAX_CHROMOSOMES\n";
        exit(-1);
    }
}

if ( $numSequences == 0 or $seqLength == 0 ) {
    print STDOUT "log: error: couldn't parse sequences\n";
    exit(-1);
}

print STDOUT "log: Processed " . commify($numSequences) . " sequences total\n";

# Index the overall fasta file
print STDOUT "Indexing genome file\n";
my $rc = CoGe::Core::Storage::index_genome_file(
    file_path => "$staging_dir/genome.faa",
    compress  => $compress
);
if ( $rc != 0 ) {
    print STDOUT "log: error: couldn't index fasta file\n";
    exit(-1) if (!$split); # need to abort if not doing legacy method
}

################################################################################
# If we've made it this far without error then we can feel confident about our
# ability to parse all of the input files.  Now we can go ahead and
# create the db entities and install the files.
################################################################################

# Connect to database
my $connstr = "dbi:mysql:dbname=$db;host=$host;port=$port;";
my $coge = CoGeX->connect( $connstr, $user, $pass );
unless ($coge) {
    print STDOUT "log: error: couldn't connect to database\n";
    exit(-1);
}

# Retrieve organism
my $organism = $coge->resultset('Organism')->find($organism_id);
unless ($organism) {
    print STDOUT "log: error finding organism id$organism_id\n";
    exit(-1);
}

# Create datasource
print STDOUT "log: Updating database ...\n";
my $datasource;
if ($source_id) {
	$datasource = $coge->resultset('DataSource')->find($source_id);
}
else {
	$datasource = $coge->resultset('DataSource')->find_or_create({ name => $source_name, description => $source_desc });
}
die "Error creating/finding data source" unless $datasource;
print STDOUT "datasource id: " . $datasource->id . "\n";

# Retrieve user
my ($user, $creator);
if ($user_id) {
	$user = $coge->resultset('User')->find($user_id);
}
else {
	$user = $coge->resultset('User')->find( { user_name => $user_name } );
}

unless ($user) {
    print STDOUT "log: error finding user '$user_name'\n";
    exit(-1);
}

# Retreive creator
if ($creator_id) {
    $creator = $coge->resultset('User')->find($creator_id);
}

$creator = $user unless $creator;


# Create genome
my $genome = $coge->resultset('Genome')->create(
    {
        name                     => $name,
        description              => $description,
        message					 => $message,
        link                     => $link,
        version                  => $version,
        organism_id              => $organism->id,
        genomic_sequence_type_id => $type_id,
        creator_id               => $creator->id,
        restricted               => $restricted,
    }
);
unless ($genome) {
    print STDOUT "log: error creating genome\n";
    exit(-1);
}
print STDOUT "log: Created genome id" . $genome->id . "\n";

# Determine installation path
unless ($install_dir) {
    unless ($P) {
        print STDOUT "log: error: can't determine install directory, set 'install_dir' or 'config' params\n";
        exit(-1);
    }
    $install_dir = $P->{SEQDIR};
}
$install_dir = "$install_dir/"
  . CoGe::Core::Storage::get_tiered_path( $genome->id ) . "/";
print STDOUT "install path: $install_dir\n";

# mdb removed 7/29/13, issue 77
#$genome->file_path( $install_dir . $genome->id . ".faa" );
#$genome->update;

# This is a check for dev server which may be out-of-sync with prod
if ( -e $install_dir ) {
    print STDOUT "log: error: install path already exists\n";
    exit(-1);
}

# Make user owner of new genome
my $node_types = CoGeX::node_types();

# Add owner connector
my $conn       = $coge->resultset('UserConnector')->create(
    {
        parent_id   => $user->id,
        parent_type => $node_types->{user},
        child_id    => $genome->id,
        child_type  => $node_types->{genome},
        role_id     => 2                        # FIXME hardcoded
    }
);
unless ($conn) {
    print STDOUT "log: error creating user connector\n";
    exit(-1);
}

# Create datasets
my %datasets;
foreach my $file (@files) {
    my $filename = basename($file);#$file =~ /^.+\/([^\/]+)$/;
    my $dataset = $coge->resultset('Dataset')->create(
        {
            data_source_id => $datasource->id,
            name           => $filename,
            description    => $description,
            version        => $version,
            restricted     => $restricted,
        }
    );
    unless ($dataset) {
        print STDOUT "log: error creating dataset\n";
        exit(-1);
    }

    #TODO set link field if loaded from FTP
    print STDOUT "dataset id: " . $dataset->id . "\n";
    $datasets{$file} = $dataset->id;

    my $dsconn =
      $coge->resultset('DatasetConnector')->create( {
          dataset_id => $dataset->id,
          genome_id => $genome->id
      } );
    unless ($dsconn) {
        print STDOUT "log: error creating dataset connector\n";
        exit(-1);
    }
}

# Create genomic_sequence/feature/location for ea. chromosome
foreach my $chr ( sort keys %sequences ) {
    my $seqlen = $sequences{$chr}{size};
    my $dsid   = $datasets{ $sequences{$chr}{file} };

    $genome->add_to_genomic_sequences(
        { sequence_length => $seqlen, chromosome => $chr } );

	# Must add a feature of type chromosome to the dataset so the dataset
	# "knows" its chromosomes
    my $feat_type =
      $coge->resultset('FeatureType')
      ->find_or_create( { name => 'chromosome' } );
    my $feat = $coge->resultset('Feature')->find_or_create(
        {
            dataset_id      => $dsid,
            feature_type_id => $feat_type->id,
            start           => 1,
            stop            => $seqlen,
            chromosome      => $chr,
            strand          => 1
        }
    );
    my $feat_name =
      $coge->resultset('FeatureName')
      ->find_or_create(
        { name => "chromosome $chr", feature_id => $feat->id } );

    my $loc = $coge->resultset('Location')->find_or_create(
        {
            feature_id => $feat->id,
            start      => 1,
            stop       => $seqlen,
            strand     => 1,
            chromosome => $chr
        }
    );
}

# Copy files from staging directory to installation directory
my $t1 = new Benchmark;
print STDOUT "log: Copying files ...\n";
unless ( mkpath($install_dir) ) {
    print STDOUT "log: error in mkpath\n";
    exit(-1);
}
if ($split) {
    unless ( mkpath( $install_dir . "/chr" ) ) {
        print STDOUT "log: error in mkpath\n";
        exit(-1);
    }
    execute("cp -r $staging_dir/chr $install_dir");
}

# mdb changed 7/31/13, issue 77 - keep filename as "genome.faa" instead of "<gid>.faa"
#my $genome_filename = $genome->id . ".faa";
#execute( "cp $staging_dir/genome.faa $install_dir/$genome_filename" );
execute("cp $staging_dir/genome.faa $install_dir/"); #FIXME use perl copy and detect failure
execute("cp $staging_dir/genome.faa.fai $install_dir/");
if ($compress) {
    execute("cp $staging_dir/genome.faa.razf $install_dir/");
    execute("cp $staging_dir/genome.faa.razf.fai $install_dir/");
}
my $time = timestr( timediff( new Benchmark, $t1 ) );
print STDOUT "Took $time to copy\n";

print STDOUT "log: "
  . commify($numSequences)
  . " sequences loaded totaling "
  . commify($seqLength) . " nt\n";

# Update IRODS metadata
if ($irods_files) {
	my @irods = split( ',', $irods_files );
	my %metadata = (
		$IRODS_METADATA_PREFIX.'link' => $P->{SERVER} . 'GenomeInfo.pl?gid=' . $genome->id
		#TODO need to add fields that match GenomeInfo.pl
	);
	foreach my $file (@irods) {
		CoGe::Accessory::IRODS::irods_imeta($file, \%metadata);
	}
}

# Copy log file into installation directory
#print STDOUT "log: All done!\n";
#close($log);
#`cp $staging_dir/log.txt $install_dir/`; #FIXME use perl copy and detect failure

# Save result document
if ($result_dir) {
    mkpath($result_dir);
    CoGe::Accessory::TDS::write(
        catfile($result_dir, '1'),
        {
            genome_id => int($genome->id)
        }
    );
}

# Create "log.done" file to indicate completion to JEX
touch($logdonefile);

exit;

#-------------------------------------------------------------------------------
sub process_fasta_file {
    my $pSeq       = shift;
    my $filepath   = shift;
    my $target_dir = shift;
    print STDOUT "process_fasta_file: $filepath\n";
    
    my $fileSize    = -s $filepath;
    my $lineNum     = 0;
    my $totalLength = 0;
    
    # Open fasta file
    my $in; # file handle
    unless (open( $in, $filepath )) {
        print STDOUT "log: error: Error opening file for reading: $!\n";
        exit(-1);
    }
    
    # Process fasta file by sections
    $/ = '>'; # set file parsing delimiter
    while (my $section = <$in>) {
        $lineNum++;
        
        # Process the section in chunks.  There is a known problem where
        # Perl substitions fail on strings larger than 1GB.
        my $sectionName;
        my $sectionLen = length($section);
        my $filteredSeq;
        my $CHUNK_LEN = 1*1000;
        my $ofs = 0;
        while ($ofs < $sectionLen) {
            my $chunk = substr($section, $ofs, $CHUNK_LEN);
            $ofs += $CHUNK_LEN;
            
            $chunk =~ s/>//g;
            $chunk =~ s/^\n+//m;
            $chunk =~ s/\n+$//m;
            next unless $chunk;
            
            if ($ofs == $CHUNK_LEN) { # first chunk
                ( $sectionName, $chunk ) = split(/\n/, $chunk, 2);
            }
            elsif ($sectionLen < $ofs) { # last chunk
                $chunk =~ s/\s+$//; # trim trailing whitespace
            }
            $chunk =~ s/\n//g;
            $chunk =~ s/\r//g;
            $filteredSeq .= $chunk;
        }
        next unless $filteredSeq;
    
        # Convert refseq (chromosome) name
        my $chr;
        if ($keep_headers) { # Don't modify refseq name
            $chr = $sectionName;
        }
        else {
            ($chr) = split(/\s+/, $sectionName);
            $chr = fix_chromosome_id($chr);
        }

        # Check validity of chr name and sequence
        if ( not defined $chr ) {
            print STDOUT "log: error parsing section header, line $lineNum, name='$name'\n";
            exit(-1);
        }
        if ( length $filteredSeq == 0 ) {
            print STDOUT "log: warning: skipping zero-length section '$chr'\n";
            next;
        }
        if ( length($chr) > $MAX_CHR_NAME_LENGTH ) {
            print STDOUT "log: error: section header name '$chr' is too long (>$MAX_CHR_NAME_LENGTH characters)\n";
            exit(-1);
        }
        if ( defined $pSeq->{$chr} ) {
            print STDOUT "log: error: Duplicate section name '$chr'\n";
            exit(-1);
        }
        if ( $filteredSeq =~ /\W/ ) {
            print STDOUT "log: error: sequence on line $lineNum contains non-alphanumeric characters, perhaps this is not a FASTA file?\n";
            exit(-1);
        }

        # Append sequence to master file
	my $out;
        unless (open( $out, ">>$target_dir/genome.faa" )) {
            print STDOUT "log: error: Couldn't open genome.faa\n";
            exit(-1);
        }
        my $head = $chr =~ /^\d+$/ ? "gi" : "lcl";
        $head .= "|" . $chr;
        print_fasta($out, $head, \$filteredSeq);
        close($out);

        # Create individual file for chromosome
        if ($split) {
            mkpath("$target_dir/chr");
            open( $out, ">$target_dir/chr/$chr" );
            print $out $filteredSeq;
            close($out);
        }

        $pSeq->{$chr} = { size => length $filteredSeq, file => $filepath };

        # Print log message
        my $count = keys %$pSeq;
        $totalLength += length $filteredSeq;
        if ( !$ignore_chr_limit && ($count > $MAX_CHROMOSOMES or $totalLength > $MAX_SEQUENCE_SIZE) ) {
            return $totalLength;
        }
        if ( $count <= $MAX_PRINT ) {
            my $filename = basename($filepath);
            print STDOUT "log: Processed chr '$chr' in $filename (".commify(length($filteredSeq))." bp)\n";
        }
        elsif ( $count == $MAX_PRINT + 1 ) {
            print STDOUT "log: (only showing first $MAX_PRINT chromosomes)\n";
        }
        elsif ( ( $count % 10000 ) == 0 ) {
            print STDOUT "log: Processed "
              . commify($count) . " ("
              . units($totalLength) . ", "
              . int( 100 * $totalLength / $fileSize )
              . "%) sequences so far ...\n";
        }
    }
    close($in);

    return $totalLength;
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
    
    return \@cmdOut;
}
