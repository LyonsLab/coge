#!/usr/bin/perl -w

use strict;
use CoGeX;
use CoGe::Accessory::Web;
use CoGe::Accessory::Utils;
use Data::Dumper;
use Getopt::Long;
use File::Path;

use vars qw($DEBUG $GO $conf_file $coge $gid $restricted $P $mask $uid);

GetOptions(
    "debug=s"        => \$DEBUG,
    "go=s"           => \$GO,
    "conf_file|cf=s" => \$conf_file,
    "gid=i"          => \$gid,
    "restricted|r=i" => \$restricted,
    "mask|m"         => \$mask,
    "uid=i" => \$uid,    #coge user id to whom the genome will be assigned
);

$P = CoGe::Accessory::Web::get_defaults($conf_file);

unless ( $P && $P->{DBNAME} ) { usage(); }

my $SECTEMPDIR = $P->{TEMPDIR};
$SECTEMPDIR .= "CopyMaskGenome/";
$SECTEMPDIR .= $uid ? $uid."/" : "public/";
my $uuid = CoGe::Accessory::Utils->get_unique_id;
$SECTEMPDIR .= $uuid."/staging/";
mkpath( $SECTEMPDIR, 1, 0777 );

my $DBNAME = $P->{DBNAME};
my $DBHOST = $P->{DBHOST};
my $DBPORT = $P->{DBPORT};
my $DBUSER = $P->{DBUSER};
my $DBPASS = $P->{DBPASS};

my $connstr =
  "dbi:mysql:dbname=" . $DBNAME . ";host=" . $DBHOST . ";port=" . $DBPORT;
$coge = CoGeX->connect($connstr, $DBUSER, $DBPASS);

#$coge->storage->debugobj(new DBIxProfiler());
#$coge->storage->debug(1);

#EL: update to use Matt's load_genome.pl
my $fasta_genome_loader   = $P->{LOAD_GENOME};
my $replicate_annotations = $P->{REPLICATE_ANNOTATIONS};
my $windowmasker          = $P->{WINDOWMASKER};
my $hard_mask             = $P->{HARD_MASK};

unless ($gid) {
    usage();
}

unless ( -r $fasta_genome_loader ) {
    print STDERR qq{
Error:  unable to read fasta_genome_loader: $fasta_genome_loader
};
    exit;
}

my ($genome) = $coge->resultset('Genome')->find($gid);
$restricted = $genome->restricted
  unless defined $restricted
; #use whatever is specified by the original genome unless user specifies otherwise
unless ($genome) {
    print STDERR qq{
Error: Unable to create a datasetgroup object for $gid
};
    exit;
}

my $new_faa;
my $stid;
if ($mask) {
    $new_faa = mask_genome( genome => $genome );
    $stid = get_sequence_type_id_for_windowmasker();
}
else {
    $new_faa = $genome->file_path;
    $stid    = $genome->genomic_sequence_type_id;
}

print STDERR "Loading genome\n";
my $new_gid = load_genome( faa => $new_faa, genome => $genome, stid => $stid );
print STDERR "NEW GENOMEID: $new_gid\n" if $DEBUG;
unless ($new_gid) {
    print STDERR
"Error: unable to capture a new genome id for the masked genome.  Exiting.\n";
    exit;
}

#assign genome to user
if ($uid) {
    my $node_types = CoGeX::node_types();
    my $uc         = $coge->resultset('UserConnector')->find_or_create(
        {
            parent_id   => $uid,
            parent_type => $node_types->{user},
            child_id    => $new_gid,
            child_type  => $node_types->{genome},
            role_id     => 2
        }
    );
}

#copy_permissions (gid1=>$gid, gid2=>$new_gid);
add_annotations( gid1 => $gid, gid2 => $new_gid );

sub copy_permissions {
    my %opts = @_;
    my $gid1 = $opts{gid1};
    my $gid2 = $opts{gid2};

    #permissions will be copied from g1 to g2;
    my $g1 = $coge->resultset('Genome')->find($gid1);
    my $g2 = $coge->resultset('Genome')->find($gid2);
}

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
    $cmd .= " -go " . $GO;
    $cmd .= " -debug " . $DEBUG if $DEBUG;
    print STDERR "Running: $cmd\n";
    open( CMD, "$cmd |" );

    while (<CMD>) {
        print STDERR $_ if $DEBUG;
    }
    close CMD;
}

sub load_genome {
    my %opts = @_;
    my $faa  = $opts{faa};
    my $genome  = $opts{genome};
    my $stid = $opts{stid};
    my $cmd  = $fasta_genome_loader;
    
    $cmd .= " -staging_dir " . $SECTEMPDIR;
    $cmd .= " -organism_id " . $genome->organism->id;
    $cmd .= " -name '" . $genome->name . "'" if $genome->name;
    #NOT AN OPTION IN MATT's SCRIPT YET
    $cmd .= " -message 'Copy of genome gid:" . $genome->id . "'";
    $cmd .= " -link 'OrganismView.pl?gid=" . $genome->id . "'";

    my ($ds) = $genome->datasets;
    #NOT AN OPTION IN MATT's SCRIPT
    $cmd .= " -source_id " . $ds->data_source->id;
    #NOT AN OPTION IN MATT's SCRIPT
    $cmd .= " -user_id " . $uid if $uid;

    $cmd .= " -version '" . $genome->version . "'";
#    $cmd .= " -ds_name '" . $ds->name . "'";
#    $cmd .= " -ds_desc '" . $ds->description . "'" if $ds->description;
    $cmd .= " -type_id " . $stid;
    $cmd .= " -config " . $conf_file;
    $cmd .= " -restricted " . $restricted if $restricted;
    $cmd .= " -fasta_files " . $faa;

#    my $seq_dir = $genome->file_path;
#    $seq_dir =~ s/[^\/]*$//;
#    $seq_dir =~ s/\d+\///g;
#    $cmd .= " -sd " . $seq_dir;

    print STDERR "Running: ", $cmd, "\n";
    my $genomeid;

    if ($GO) {
      `$cmd`;
      open (IN , $SECTEMPDIR."log.txt") || die "can't find log file: ".$SECTEMPDIR."log.txt";
      while (<IN>)
	{
	  if (/genome id:\s+(\d+)/) {
	    $genomeid = $1;
	    print STDERR "Captured genomeid: $genomeid\n" if $DEBUG;
	    last;
	  }
	}
      close IN;
      die "Unable able to find genomeid in log file: ".$SECTEMPDIR."log.txt" unless $genomeid;
    }
    return $genomeid;
}

sub get_and_clean_sequence {
    my $faa     = shift;
    my $rnd     = int( rand(1000000000000) );
    my $new_faa = $SECTEMPDIR . "/" . $rnd . ".clean.faa";
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
    my $genome    = $opts{genome};
    my $rnd    = int( rand(1000000000000) );
    #EL: change this to use the SECURE TEMP DIR
    my $counts = $SECTEMPDIR . "/" . $rnd . ".counts";
    my $masked = $SECTEMPDIR . "/" . $rnd . ".masked";
    my $hard   = $masked . ".hard";
    my $faa    = $genome->file_path;
    my $cmd    = $windowmasker . " -in " . $faa . " -mk_counts -out " . $counts;
    print STDERR "running $cmd\n";
    `$cmd` if $GO;
    $cmd =
        $windowmasker . " -in " 
      . $faa
      . " -ustat "
      . $counts
      . " -outfmt fasta -dust T -out "
      . $masked;
    print STDERR "running $cmd\n";
    `$cmd` if $GO;
    $cmd = $hard_mask . " < $masked > $hard";
    print STDERR "running $cmd\n";
    `$cmd` if $GO;
    return $hard;
}

sub get_sequence_type_id_for_windowmasker {
    my $st =
      $coge->resultset('GenomicSequenceType')
      ->find_or_create( { name => "NCBI WindowMasker (Hard)" } );
    return $st->id;
}

sub usage {
    print STDERR qq{
Welcome to $0

Purpose:  take a genome ID from coge, generate a masked version of the genome, load it into CoGe, and map over the annotations

Usage:  $0  -conf_file <coge configuration file> -gid <coge database id for genome to be masked and copies> -uid <coge user id> -go 1

Options:
   -go 1                  |      Make the database calls.  Default 0

   -conf_file | cf        |      CoGe conf file

   -gid                   |      CoGe genome id

   -mask                  |      Mask the genome with NCBI's windowmasker

   -restricted            |      mark genome as restricted (Defaults to whatever is set for gid's genome)

   -uid                   |      CoGe user id to whom the new genome will be assigned


Copyright Eric Lyons 2013

};
    exit;
}
