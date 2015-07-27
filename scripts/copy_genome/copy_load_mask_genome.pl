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
use File::Slurp;
use URI::Escape::JavaScript qw(escape);

use vars qw($conf_file $coge $gid $restricted $P $mask $uid $staging_dir $seq_only $wid);

GetOptions(
    "conf_file|cf=s"   => \$conf_file,
    "staging_dir=s"    => \$staging_dir,
    "gid=i"            => \$gid,
    "restricted|r=i"   => \$restricted,
    "mask|m=i"         => \$mask,
    "uid=i"            => \$uid, # user id to assign the new genome to
    "wid=i"            => \$wid, # workflow id
    "sequence_only=i"  => \$seq_only, #flag for only copying the sequences -- no structural annotations
);

$staging_dir //= "."; #/
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
unless ($uid) {
    print $log "log: error: uid not specified\n";
    exit(-1);
}

unless ($wid) {
    print $log "log: error: wid not specified\n";
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
my $replicate_annotations = $P->{SCRIPTDIR} . '/copy_genome/replicate_annotations.pl';
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
#print $log "log: Loading new genome:\n";
print $log "Calling load_genome.pl\n";
my $new_gid = load_genome( faa => $new_faa, genome => $genome, stid => $stid );
unless ($new_gid) {
    print $log "log: error: Unable to load new gnome\n";
    exit(-1);
}

say STDOUT "SEQ_ONLY $seq_only";

# Copy annotations
unless ($seq_only and $seq_only == 1) {
    say STDERR "GID: $gid new GID: $new_gid";
    print $log "log: Copying annotations (may take a few minutes)\n";
    add_annotations( gid1 => $gid, gid2 => $new_gid );
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
    $cmd .= " -ignore_chr_limit 1"; # mdb added 3/9/15 COGE-595
    $cmd .= " -staging_dir " . $staging_dir;
    $cmd .= " -organism_id " . $genome->organism->id;
    $cmd .= " -name '" . escape($genome->name) . "'" if $genome->name;
    $cmd .= " -message 'Copy of genome gid:" . $genome->id . "'";
    $cmd .= " -link 'OrganismView.pl?gid=" . $genome->id . "'";

    my ($ds) = $genome->datasets;
    my $creator = $genome->creator;

    $cmd .= ($ds ? " -source_id " . $ds->data_source->id : " -source_name Unknown");
    $cmd .= " -user_id " . $uid if $uid;
    $cmd .= " -wid " . $wid;
    $cmd .= " -creator_id " . $creator->id if $creator;
    $cmd .= " -version '" . escape($genome->version) . "'";
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

    # Get id for newly loaded genome from log file output of load_genome.pl
    my $file = catfile($staging_dir, "log.txt");
    my $log = read_file($file);
    my ($genome_id) = $log =~ /genome id:\s+(\d+)/;
    unless ($genome_id) {
        print $log "Unable able to find genome id in log file: $file\n";
        return;
    }

    return $genome_id;
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
