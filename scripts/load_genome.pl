#!/usr/bin/perl -w

use strict;
use CoGeX;
use CoGe::Accessory::Web;
use CoGe::Accessory::Storage qw( index_genome_file get_tiered_path );
use CoGe::Accessory::Utils qw( commify units print_fasta );
use Roman;
use Data::Dumper;
use Getopt::Long;
use File::Path;
use URI::Escape::JavaScript qw(unescape);
use POSIX qw(ceil);
use Benchmark;

use vars qw($staging_dir $install_dir $fasta_files
  $name $description $link $version $type_id $restricted $message
  $organism_id $source_id $source_name $source_desc $user_id $user_name 
  $keep_headers $split $compress
  $host $port $db $user $pass $config
  $P $MAX_CHROMOSOMES $MAX_PRINT $MAX_SEQUENCE_SIZE $MAX_CHR_NAME_LENGTH );

$MAX_CHROMOSOMES = 200000;    # max number of chromosomes or contigs
$MAX_PRINT       = 5;
$MAX_SEQUENCE_SIZE   = 5 * 1024 * 1024 * 1024;    # 5 gig
$MAX_CHR_NAME_LENGTH = 255;

GetOptions(
    "staging_dir=s" => \$staging_dir,
    "install_dir=s" => \$install_dir,    # optional
    "fasta_files=s" => \$fasta_files,    # comma-separated list (JS escaped)
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
    "user_name=s"   => \$user_name,      # user name
    "keep_headers=i" =>
      \$keep_headers,    # flag to keep original headers (no parsing)
    "split=i"    => \$split,       # split fasta into chr directory
    "compress=i" => \$compress,    # compress fasta into RAZF before indexing

    # Database params
    "host|h=s"      => \$host,
    "port|p=s"      => \$port,
    "database|db=s" => \$db,
    "user|u=s"      => \$user,
    "password|pw=s" => \$pass,

    # Or use config file
    "config=s" => \$config
);

# Open log file
$| = 1;
die unless ($staging_dir);
mkpath($staging_dir); # make sure this exists
my $logfile = "$staging_dir/log.txt";
open( my $log, ">>$logfile" ) or die "Error opening log file";
$log->autoflush(1);
print $log "Starting $0 (pid $$)\n";

# Process and verify parameters
$fasta_files = unescape($fasta_files) if ($fasta_files);
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

if (not $source_id and not $source_name) {
	print $log "log: error: source not specified, use source_id or source_name\n";
	exit(-1);
}
if (not $user_id and not $user_name) {
	print $log "log: error: user not specified, use user_id or user_name\n";
	exit(-1);
}
if ($user_name and $user_name eq 'public') {
	print $log "log: error: not logged in\n";
    exit(-1);
}

# Load config file ... only needed to get the DB params if user specified
# a config file instead of specifying them on command-line.
if ($config) {
    $P    = CoGe::Accessory::Web::get_defaults($config);
    $db   = $P->{DBNAME};
    $host = $P->{DBHOST};
    $port = $P->{DBPORT};
    $user = $P->{DBUSER};
    $pass = $P->{DBPASS};
}

# Process each file into staging area
my @files = split( ',', $fasta_files );
my %sequences;
my $seqLength;
my $numSequences;
foreach my $file (@files) {
    if ( -B $file ) {
        my ($filename) = $file =~ /^.+\/([^\/]+)$/;
        print $log "log: error: '$filename' is a binary file\n";
        exit(-1);
    }

    $seqLength += process_fasta_file( \%sequences, $file, $staging_dir );
    $numSequences = keys %sequences;

    if ( $seqLength > $MAX_SEQUENCE_SIZE ) {
        print $log "log: error: total sequence size exceeds limit of "
          . units($MAX_SEQUENCE_SIZE) . "\n";
        exit(-1);
    }
    if ( $numSequences > $MAX_CHROMOSOMES ) {
        print $log
          "log: error: too many sequences, limit is $MAX_CHROMOSOMES\n";
        exit(-1);
    }
}

if ( $numSequences == 0 or $seqLength == 0 ) {
    print $log "log: error: couldn't parse sequences\n";
    exit(-1);
}

print $log "log: Processed " . commify($numSequences) . " sequences total\n";

# Index the overall fasta file
print $log "Indexing genome file\n";
my $rc = CoGe::Accessory::Storage::index_genome_file(
    file_path => "$staging_dir/genome.faa",
    compress  => $compress
);
if ( $rc != 0 ) {
    print $log "log: error: couldn't index fasta file\n";
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
    print $log "log: error: couldn't connect to database\n";
    exit(-1);
}

# Retrieve organism
my $organism = $coge->resultset('Organism')->find($organism_id);
unless ($organism) {
    print $log "log: error finding organism id$organism_id\n";
    exit(-1);
}

# Create datasource
print $log "log: Updating database ...\n";
my $datasource;
if ($source_id) {
	$datasource = $coge->resultset('DataSource')->find($source_id);
}
else {
	$datasource = $coge->resultset('DataSource')->find_or_create({ name => $source_name, description => $source_desc });
}
die "Error creating/finding data source" unless $datasource;
print $log "datasource id: " . $datasource->id . "\n";

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
        restricted               => $restricted
    }
);
unless ($genome) {
    print $log "log: error creating genome\n";
    exit(-1);
}
print $log "genome id: " . $genome->id . "\n";

# Determine installation path
unless ($install_dir) {
    unless ($P) {
        print $log
"log: error: can't determine install directory, set 'install_dir' or 'config' params\n";
        exit(-1);
    }
    $install_dir = $P->{SEQDIR};
}
$install_dir = "$install_dir/"
  . CoGe::Accessory::Storage::get_tiered_path( $genome->id ) . "/";
print $log "install path: $install_dir\n";

# mdb removed 7/29/13, issue 77
#$genome->file_path( $install_dir . $genome->id . ".faa" );
#$genome->update;

# This is a check for dev server which may be out-of-sync with prod
if ( -e $install_dir ) {
    print $log "log: error: install path already exists\n";
    exit(-1);
}

# Make user owner of new genome
my $user;
if ($user_id) {
	$user = $coge->resultset('User')->find($user_id);
}
else {
	$user = $coge->resultset('User')->find( { user_name => $user_name } );
}
unless ($user) {
    print $log "log: error finding user '$user_name'\n";
    exit(-1);
}
my $node_types = CoGeX::node_types();
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
    print $log "log: error creating user connector\n";
    exit(-1);
}

# Create datasets
my %datasets;
foreach my $file (@files) {
    my ($filename) = $file =~ /^.+\/([^\/]+)$/;
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
        print $log "log: error creating dataset\n";
        exit(-1);
    }

    #TODO set link field if loaded from FTP
    print $log "dataset id: " . $dataset->id . "\n";
    $datasets{$file} = $dataset->id;

    my $dsconn =
      $coge->resultset('DatasetConnector')
      ->create( { dataset_id => $dataset->id, genome_id => $genome->id } );
    unless ($dsconn) {
        print $log "log: error creating dataset connector\n";
        exit(-1);
    }
}

# Create genomic_sequence/feature/location for ea. chromosome
foreach my $chr ( sort keys %sequences ) {
    my $seqlen = $sequences{$chr}{size};
    my $dsid   = $datasets{ $sequences{$chr}{file} };

    $genome->add_to_genomic_sequences(
        { sequence_length => $seqlen, chromosome => $chr } );

# Must add a feature of type chromosome to the dataset so the dataset "knows" its chromosomes
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

print $log "log: Added genome id"
  . $genome->id
  . "\n";    # !!!! don't change, gets parsed by calling code

# Copy files from staging directory to installation directory
my $t1 = new Benchmark;
print $log "log: Copying files ...\n";
unless ( mkpath($install_dir) ) {
    print $log "log: error in mkpath\n";
    exit(-1);
}
if ($split) {
    unless ( mkpath( $install_dir . "/chr" ) ) {
        print $log "log: error in mkpath\n";
        exit(-1);
    }
    execute("cp -r $staging_dir/chr $install_dir");
}

# mdb changed 7/31/13, issue 77 - keep filename as "genome.faa" instead of "<gid>.faa"
#my $genome_filename = $genome->id . ".faa";
#execute( "cp $staging_dir/genome.faa $install_dir/$genome_filename" );
execute("cp $staging_dir/genome.faa $install_dir/");
execute("cp $staging_dir/genome.faa.fai $install_dir/");
if ($compress) {
    execute("cp $staging_dir/genome.faa.razf $install_dir/");
    execute("cp $staging_dir/genome.faa.razf.fai $install_dir/");
}

# Yay, log success!
CoGe::Accessory::Web::log_history(
    db          => $coge,
    user_id     => $user->id,
    page        => "LoadGenome",
    description => 'load genome id' . $genome->id,
    link        => 'GenomeInfo.pl?gid=' . $genome->id
);
print $log "log: "
  . commify($numSequences)
  . " sequences loaded totaling "
  . commify($seqLength) . " nt\n";

my $t2 = new Benchmark;
my $time = timestr( timediff( $t2, $t1 ) );
print $log "Took $time to copy\n";
print $log "log: Finished loading genome!\n";

close($log);

# Copy log file from staging directory to installation directory
`cp $staging_dir/log.txt $install_dir/`;

exit;

#-------------------------------------------------------------------------------
sub process_fasta_file {
    my $pSeq       = shift;
    my $filepath   = shift;
    my $target_dir = shift;

    print $log "process_fasta_file: $filepath\n";
    $/ = "\n>";
    open( my $in, $filepath ) || die "can't open $filepath for reading: $!";

    my $fileSize    = -s $filepath;
    my $lineNum     = 0;
    my $totalLength = 0;
    while (<$in>) {
        s/\n*\>//g;
        next unless $_;
        my ( $name, $seq ) = split /\n/, $_, 2;
        $seq =~ s/\n//g;
        $lineNum++;

        my $chr;
        if ($keep_headers) {
            $chr = $name;
        }
        else {
            ($chr) = split( /\s+/, $name );
            $chr =~ s/^lcl\|//;
            $chr =~ s/^gi\|//;
            $chr =~ s/chromosome//i;
            $chr =~ s/^chr//i;
            $chr =~ s/^0+// unless $chr eq '0';
            $chr =~ s/^_+//;
            $chr =~ s/\s+/ /;
            $chr =~ s/^\s//;
            $chr =~ s/\s$//;
        }

        # Check validity of chr name and sequence
        if ( not defined $chr ) {
            print $log
"log: error parsing section header, line $lineNum, name='$name'\n";
            exit(-1);
        }
        if ( length $seq == 0 ) {
            print $log "log: warning: skipping zero-length section '$chr'\n";
            next;
        }
        if ( length($chr) > $MAX_CHR_NAME_LENGTH ) {
            print $log
"log: error: section header name '$chr' is too long (>$MAX_CHR_NAME_LENGTH characters)\n";
            exit(-1);
        }
        if ( defined $pSeq->{$chr} ) {
            print $log "log: error: Duplicate section name '$chr'";
            exit(-1);
        }
        if ( $seq =~ /\W/ ) {
            print $log
"log: error: sequence on line $lineNum contains non-alphanumeric characters, perhaps this is not a FASTA file?\n";
            exit(-1);
        }

        my ($filename) = $filepath =~ /^.+\/([^\/]+)$/;

        # Append sequence to master file
        open( my $out, ">>$target_dir/genome.faa" );
        my $head = $chr =~ /^\d+$/ ? "gi" : "lcl";
        $head .= "|" . $chr;
        print_fasta($out, $head, \$seq); #print $out "$head\n$seq\n";
        close($out);

        # Create individual file for chromosome
        if ($split) {
            mkpath("$target_dir/chr");
            open( $out, ">$target_dir/chr/$chr" );
            print $out $seq;
            close($out);
        }

        $pSeq->{$chr} = { size => length $seq, file => $filepath };

        # Print log message
        my $count = keys %$pSeq;
        $totalLength += length $seq;
        if ( $count > $MAX_CHROMOSOMES or $totalLength > $MAX_SEQUENCE_SIZE ) {
            return $totalLength;
        }
        if ( $count <= $MAX_PRINT ) {
            print $log "log: Processed '$chr' in $filename\n";
        }
        elsif ( $count == $MAX_PRINT + 1 ) {
            print $log "log: (only showing first $MAX_PRINT)\n";
        }
        elsif ( ( $count % 10000 ) == 0 ) {
            print $log "log: Processed "
              . commify($count) . " ("
              . units($totalLength) . ", "
              . int( 100 * $totalLength / $fileSize )
              . "%) sequences so far ...\n";
        }
    }
    close($in);

    return $totalLength;
}

sub execute {
    my $cmd = shift;
    print $log "$cmd\n";
    my @cmdOut    = qx{$cmd};
    my $cmdStatus = $?;
    if ( $cmdStatus != 0 ) {
        print $log "log: error: command failed with rc=$cmdStatus: $cmd\n";
        exit(-1);
    }
}
