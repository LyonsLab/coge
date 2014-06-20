#!/usr/bin/perl -w

use strict;
no warnings 'redefine';

use CoGeX;
use CoGe::Accessory::Web;
use CoGe::Accessory::Utils;
use Data::Dumper;
use Getopt::Long;
use File::Spec::Functions;
use File::Path;

use vars qw($conf_file $coge $gid $restricted $P $mask $uid $staging_dir $result_dir $seq_only);

GetOptions(
	"conf_file|cf=s" => \$conf_file,
	"staging_dir=s"  => \$staging_dir,
	"result_dir=s"   => \$result_dir,
	"gid=i"          => \$gid,
	"restricted|r=i" => \$restricted,
	"mask|m=i"       => \$mask,
	"uid=i"          => \$uid, # user id to assign the new genome to
	"sequence_only"  => \$seq_only, #flag for only copying the sequences -- no structural annotations
);

$staging_dir //= ".";
# Open log file
$| = 1;
unless ($staging_dir) {
	die;
	#$staging_dir = $P->{SECTEMPDIR} . "CopyMaskGenome/";
	#$staging_dir .= $uid ? $uid."/" : "public/";
	#my $uuid = CoGe::Accessory::Utils->get_unique_id;
	#$staging_dir .= $uuid."/staging/";
}
mkpath( $staging_dir, 1, 0777 );    # make sure this exists
my $logfile = "$staging_dir/log.txt";
my $log = *STDOUT;
#open( my $log, ">>$logfile" ) or die "Error opening log file: $logfile: $!";
#$log->autoflush(1);
print $log "Starting $0 (pid $$)\n";

# Process and verify parameters
if (not $uid) {
	print $log "log: error: uid not specified\n";
	exit(-1);
}

# Load config file
die "No config file specified\n" unless ($conf_file);
$P = CoGe::Accessory::Web::get_defaults($conf_file);
my $DBNAME = $P->{DBNAME};
my $DBHOST = $P->{DBHOST};
my $DBPORT = $P->{DBPORT};
my $DBUSER = $P->{DBUSER};
my $DBPASS = $P->{DBPASS};

# Connect to the database
my $connstr =
  "dbi:mysql:dbname=" . $DBNAME . ";host=" . $DBHOST . ";port=" . $DBPORT;
$coge = CoGeX->connect( $connstr, $DBUSER, $DBPASS );

# Get paths to external scripts
my $fasta_genome_loader   = $P->{SCRIPTDIR} . '/load_genome.pl';
my $replicate_annotations =
  $P->{SCRIPTDIR} . '/copy_genome/replicate_annotations.pl';
my $windowmasker = $P->{SCRIPTDIR} . '/copy_genome/windowmasker';
my $hard_mask    = $P->{SCRIPTDIR} . '/copy_genome/hard_mask.pl';

unless ( -r $fasta_genome_loader ) {
	print $log "log: error: Unable to find fasta_genome_loader\n";
	exit(-1);
}

# Load genome from db
my ($genome) = $coge->resultset('Genome')->find($gid);
$restricted = $genome->restricted
  unless defined $restricted; #use whatever is specified by the original genome unless user specifies otherwise
unless ($genome) {
	print $log "log: error: Unable to find genome $gid\n";
	exit(-1);
}

# Mask genome sequence
my $new_faa;
my $stid;
if ($mask) {
	print $log "log: Masking genome sequence (may take a few minutes):\n";
	$new_faa = mask_genome( genome => $genome );
	$stid = get_sequence_type_id_for_windowmasker();
}
else {
	$new_faa = $genome->file_path;
	$stid    = $genome->genomic_sequence_type_id;
}

# Load genome
print $log "log: Loading new genome:\n";
print $log "Calling load_genome.pl\n";
my $new_gid = load_genome( faa => $new_faa, genome => $genome, stid => $stid );
unless ($new_gid) {
	print $log "log: error: Unable to load new gnome\n";
	exit(-1);
}

# Copy annotations
unless ($seq_only)
  {
    print $log "log: Copying annotations (may take a few minutes)\n";
    add_annotations( gid1 => $gid, gid2 => $new_gid );
  }


if ($result_dir) {
    mkpath($result_dir);
    CoGe::Accessory::TDS::write(
        catfile($result_dir, '1'),
        {
            genome_id => int($new_gid)
        }
    );
}

# Done
print $log "log: Finished copying genome!\n";
exit(0);

#-------------------------------------------------------------------------------

sub add_annotations {
	my %opts   = @_;
	my $dsgid1 = $opts{gid1};
	my $dsgid2 = $opts{gid2};
	my $cmd    = $replicate_annotations;
	$cmd .= " -db " . $DBNAME;
	$cmd .= " -u " . $DBUSER;
	$cmd .= " -pw " . $DBPASS;
	$cmd .= " -dsgid1 " . $dsgid1;
	$cmd .= " -dsgid2 " . $dsgid2;
	execute($cmd);
}

sub load_genome {
	my %opts   = @_;
	my $faa    = $opts{faa};
	my $genome = $opts{genome};
	my $stid   = $opts{stid};
	my $cmd    = $fasta_genome_loader;
	$cmd .= " -staging_dir " . $staging_dir;
	$cmd .= " -organism_id " . $genome->organism->id;
	$cmd .= " -name '" . $genome->name . "'" if $genome->name;
	$cmd .= " -message 'Copy of genome gid:" . $genome->id . "'";
	$cmd .= " -link 'OrganismView.pl?gid=" . $genome->id . "'";

	my ($ds) = $genome->datasets;

	$cmd .= " -source_id " . $ds->data_source->id;
	$cmd .= " -user_id " . $uid if $uid;
	$cmd .= " -version '" . $genome->version . "'";
	#$cmd .= " -ds_name '" . $ds->name . "'";
	#$cmd .= " -ds_desc '" . $ds->description . "'" if $ds->description;
	$cmd .= " -type_id " . $stid;
	$cmd .= " -config " . $conf_file;
	$cmd .= " -restricted " . $restricted if $restricted;
	$cmd .= " -fasta_files " . $faa;

	#my $seq_dir = $genome->file_path;
	#$seq_dir =~ s/[^\/]*$//;
	#$seq_dir =~ s/\d+\///g;
	#$cmd .= " -sd " . $seq_dir;

	execute($cmd);

	my $genomeid;
	open( IN, catfile($staging_dir, "/log.txt"))
	  || die "can't find log file: " . $staging_dir . "/log.txt";
	while (<IN>) {
        print $log $_;
		if (/genome id: (\d+)/) {
			$genomeid = $1;
			print $log "Captured gid: $genomeid\n";
			last;
		}
	}
	close IN;
	
	unless ($genomeid) {
		print $log "Unable able to find gid in log file: " . $staging_dir . "/log.txt";
		exit(-1);
	}

	return $genomeid;
}

sub get_and_clean_sequence {
	my $faa     = shift;
	my $rnd     = int( rand(1000000000000) );
	my $new_faa = $staging_dir . "/" . $rnd . ".clean.faa";
	open( OUT, ">$new_faa" );
	open( IN,  $faa );
	while (<IN>) {
		if (/^>/) {
			s/lcl\|//g;
			s/gi\|//g;
			print OUT $_;
		}
		else {
			print OUT $_;
		}
	}
	close OUT;
	close IN;
	return $new_faa;
}

sub mask_genome {
	my %opts   = @_;
	my $genome = $opts{genome};
	my $rnd    = int( rand(1000000000000) );
	my $counts = $staging_dir . "/" . $rnd . ".counts";
	my $masked = $staging_dir . "/" . $rnd . ".masked";
	my $hard   = $masked . ".hard";
	my $faa    = $genome->file_path;
	
	print $log "log: Running WindowMasker (stage 1)\n";
	my $cmd    = $windowmasker . " -in " . $faa . " -mk_counts -out " . $counts;
	execute($cmd);
	
	print $log "log: Running WindowMasker (stage 2)\n";
	$cmd =
	    $windowmasker . " -in " . $faa
	  . " -ustat "
	  . $counts
	  . " -outfmt fasta -dust T -out "
	  . $masked;
	execute($cmd);
	
	print $log "log: Running hard masking\n";
	$cmd = $hard_mask . " < $masked > $hard";
	execute($cmd);
	return $hard;
}

sub get_sequence_type_id_for_windowmasker {
	my $st =
	  $coge->resultset('GenomicSequenceType')
	  ->find_or_create( { name => "NCBI WindowMasker (Hard)" } );
	return $st->id;
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

#sub usage {
#    print STDERR qq{
#Purpose:  take a genome ID from coge, generate a masked version of the genome, load it into CoGe, and map over the annotations
#Usage:  $0  -conf_file <coge configuration file> -gid <coge database id for genome to be masked and copies> -uid <coge user id> -go 1
#
#Options:
#   -go 1                  |      Make the database calls.  Default 0
#   -conf_file | cf        |      CoGe conf file
#   -gid                   |      CoGe genome id
#   -mask                  |      Mask the genome with NCBI's windowmasker
#   -restricted            |      mark genome as restricted (Defaults to whatever is set for gid's genome)
#   -uid                   |      CoGe user id to whom the new genome will be assigned
#};
#    exit;
#}
