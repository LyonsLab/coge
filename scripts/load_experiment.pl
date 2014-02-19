#!/usr/bin/perl -w

use strict;
use Data::Dumper;
use CoGeX;
use Getopt::Long;
use File::Path;
use File::Touch;
use URI::Escape::JavaScript qw(unescape);
use CoGe::Accessory::Web qw(get_defaults);
use CoGe::Accessory::Utils qw( commify );

use vars qw($staging_dir $install_dir $data_file $file_type
  $name $description $version $restricted $ignore_missing_chr
  $gid $source_name $user_name $config $allow_negative
  $host $port $db $user $pass $P);

#FIXME: use these from Storage.pm instead of redeclaring them
my $DATA_TYPE_QUANT = 1;	# Quantitative data
my $DATA_TYPE_POLY	= 2;	# Polymorphism data
my $DATA_TYPE_ALIGN = 3;	# Alignments

#my $MIN_QUANT_COLUMNS = 5;
#my $MAX_QUANT_COLUMNS = 6;

my $MIN_VCF_COLUMNS = 8;
#my $MAX_VCF_COLUMNS = 10;

GetOptions(
    "staging_dir=s" => \$staging_dir,
    "install_dir=s" => \$install_dir,    # optional
    "data_file=s"   => \$data_file,      # data file (JS escape)
    "file_type=s"   => \$file_type,		 # file type
    "name=s"        => \$name,           # experiment name (JS escaped)
    "desc=s"        => \$description,    # experiment description (JS escaped)
    "version=s"     => \$version,        # experiment version (JS escaped)
    "restricted=i"  => \$restricted,     # experiment restricted flag
    "gid=s"         => \$gid,            # genome id
    "source_name=s" => \$source_name,    # data source name (JS escaped)
    "user_name=s"   => \$user_name,      # user name
    
    # Flags
    "ignore-missing-chr=i" => \$ignore_missing_chr,

    # Database params
    "host|h=s"      => \$host,
    "port|p=s"      => \$port,
    "database|db=s" => \$db,
    "user|u=s"      => \$user,
    "password|pw=s" => \$pass,

    # Or use config file
    "config=s" => \$config,

    # Optional features for debug and bulk loader
    "allow_negative=i" => \$allow_negative
);

# Open log file
$| = 1;
my $logfile = "$staging_dir/log.txt";
mkpath($staging_dir, 0, 0777) unless -r $staging_dir;
open( my $log, ">>$logfile" ) or die "Error opening log file $logfile";
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
if ($config) {
    $P = CoGe::Accessory::Web::get_defaults($config);
    $db   = $P->{DBNAME};
    $host = $P->{DBHOST};
    $port = $P->{DBPORT};
    $user = $P->{DBUSER};
    $pass = $P->{DBPASS};
}
else {
	$P = CoGe::Accessory::Web::get_defaults();
}
my $FASTBIT_LOAD  = $P->{FASTBIT_LOAD};
my $FASTBIT_QUERY = $P->{FASTBIT_QUERY};
my $SAMTOOLS = $P->{SAMTOOLS};
if (   not $FASTBIT_LOAD
    or not $FASTBIT_QUERY
    or not $SAMTOOLS
    or not -e $FASTBIT_LOAD
    or not -e $FASTBIT_QUERY
    or not -e $SAMTOOLS )
{
    print $log "FASTBIT_LOAD: $FASTBIT_LOAD\n";
    print $log "FASTBIT_QUERY: $FASTBIT_QUERY\n";
    print $log "SAMTOOLS: $SAMTOOLS\n";
    print $log "log: error: can't find required command(s)\n";
    exit(-1);
}

# Determine file type
my ($file_type, $data_type) = detect_data_type($file_type, $data_file);
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

# Copy input data file to staging area
my $cmd = "cp -f $data_file $staging_dir";
`$cmd`;
my ($filename) = $data_file =~ /^.+\/([^\/]+)$/;
my $staged_data_file = $staging_dir . '/' . $filename;

# Validate the data file
print $log "log: Validating data file\n";
my $count = 0;
my $pChromosomes;
my $data_spec;
if ( $data_type == $DATA_TYPE_QUANT ) {
    ( $staged_data_file, $data_spec, $count, $pChromosomes ) =
      validate_quant_data_file(
        file       => $staged_data_file,
        file_type  => $file_type,
        genome_chr => \%genome_chr
      );
}
elsif ( $data_type == $DATA_TYPE_POLY ) {
    ( $staged_data_file, $data_spec, $count, $pChromosomes ) =
      validate_vcf_data_file( file => $staged_data_file, genome_chr => \%genome_chr );
}
elsif ( $data_type == $DATA_TYPE_ALIGN ) {
	( $staged_data_file, $data_spec, $count, $pChromosomes ) =
      validate_bam_data_file( file => $staged_data_file, genome_chr => \%genome_chr );
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

if ( $data_type == $DATA_TYPE_QUANT or $data_type == $DATA_TYPE_POLY ) {
	# Generate fastbit database/index
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
  $coge->resultset('DataSource')
  ->find_or_create( { name => $source_name, description => "" } )
  ;    #, description => "Loaded into CoGe via LoadExperiment" } );
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

# Determine installation path
unless ($install_dir) {
    unless ($P) {
        print $log
"log: error: can't determine install directory, set 'install_dir' or 'config' params\n";
        exit(-1);
    }
    $install_dir = $P->{EXPDIR};
}
my $storage_path = "$install_dir/" . CoGe::Accessory::Storage::get_tiered_path( $experiment->id ) . '/';
print $log 'Storage path: ', $storage_path, "\n";
# mdb removed 8/7/13, issue 77
#$experiment->storage_path($storage_path);
#$experiment->update;

# This is a check for dev server which may be out-of-sync with prod
if ( -e $storage_path ) {
    print $log "log: error: install path already exists\n";
    exit(-1);
}

# Don't change, gets parsed by calling code
print $log "experiment id: " . $experiment->id . "\n";
print "experiment id: " . $experiment->id . "\n";

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
$cmd = "cp -r $staging_dir/* $storage_path";
print $log "$cmd\n";
`$cmd`;

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

#print STDERR "experiment id: ".$experiment->id,"\n";
print $log "log: All done!";
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
    
    if ( grep { $_ eq $filetype } ('csv', 'tsv', 'bed') ) { #TODO add 'bigbed', 'wig', 'bigwig', 'gff', 'gtf'
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

    print $log "validate_quant_data_file: $filepath\n";
    open( my $in, $filepath ) || die "can't open $filepath for reading: $!";
    my $outfile = $filepath . ".processed";
    open( my $out, ">$outfile" );
    while ( my $line = <$in> ) {
        $line_num++;
        next if ( $line =~ /^\s*#/ );
        chomp $line;
        
        # Interpret tokens according to file type
        my @tok;
        my ( $chr, $start, $stop, $strand, $val1, $val2 );
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
        	( $chr, $start, $stop, undef, $val1, $strand ) = @tok;
        }
        else {
        	die; # sanity check
        }

        # mdb added 2/19/14 for bulk loading
        if ($allow_negative and $val1 < 0) {
	    $val1 = -1 * $val1; 
	}

        # Validate values and set defaults
        if (   not defined $chr
            or not defined $start
            or not defined $stop
            or not defined $strand )
        {
            print $log
              "log: error at line $line_num: missing value in a column\n";
            return;
        }
        if ( not defined $val1 or $val1 < 0 or $val1 > 1 ) {
            print $log
              "log: error at line $line_num: value 1 not between 0 and 1\n";
            return;
        }
        if ( not defined $val2 ) {
            $val2 = 0;
        }

		$chr = fix_chromosome_id($chr, $genome_chr);
        if (!$chr) {
            print $log
              "log: error at line $line_num: trouble parsing chromosome\n";
            return;
        }
        $strand = $strand =~ /-/ ? -1 : 1;
        print $out join( ",", $chr, $start, $stop, $strand, $val1, $val2 ),
          "\n";
        $chromosomes{$chr}++;
        $count++;
    }
    close($in);
    close($out);
    my $format = "chr:key, start:unsigned long, stop:unsigned long, strand:byte, value1:double, value2:double";
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
            print $log
"log: error at line $line_num: missing required value in a column\n";
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
    my $format = "chr:key, start:unsigned long, stop:unsigned long, type:key, id:text, ref:key, alt:key, qual:double, info:text";
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
