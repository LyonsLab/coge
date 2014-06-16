#!/usr/bin/perl
#-------------------------------------------------------------------------------
# Usage:
#    ./migrate_genomes.pl
#
# Copy /storage/data/genomic_sequence/ but exclude chr/ subdirectories.
# Fix and index fasta files.
#
#-------------------------------------------------------------------------------

use DBI;
use strict;
use CoGeX;
use CoGe::Core::Storage;
use CoGe::Accessory::Utils;
use Getopt::Long;
use File::Path;

my $old_path = '/storage/coge/data/genomic_sequence';
my $new_path = '/storage/coge/data/genomic_sequence2';

my ($db, $user, $pass, $start_id, $stop_id);
GetOptions (
	"database|db=s"		=> \$db,
	"user|u=s"			=> \$user,
	"password|pw|p=s"	=> \$pass,
	"start_id=i"		=> \$start_id,	# optional
	"stop_id=i"			=> \$stop_id 	# optional
);
die "Missing DB params\n" unless ($db and $user and $pass);

$| = 1;

#-------------------------------------------------------------------------------
# Connect to database
#-------------------------------------------------------------------------------

my $connstr = "dbi:mysql:dbname=$db;host=localhost;port=3307";
my $coge = CoGeX->connect($connstr, $user, $pass);

#-------------------------------------------------------------------------------
# Copy/fix/index fasta file for each genome
#-------------------------------------------------------------------------------

my $count = 0;
foreach my $genome (sort {$a->id <=> $b->id} $coge->resultset('Genome')->all) {
	my $gid = $genome->id;
	next if (defined $start_id and $gid < $start_id);
	next if (defined $stop_id and $gid > $stop_id);
	
	print STDERR "Genome " . $genome->id . "\n";
	
	my $fasta_file = $genome->file_path();
	print STDERR "   old fasta: $fasta_file\n";
	if (not -e $fasta_file) {
		print STDERR "   fasta not found, skipping\n";
		next;
	}
	
	my $tier_path = CoGe::Core::Storage::get_tiered_path( $genome->id );
	my $new_fasta_path = $new_path . '/' . $tier_path;
	my $new_fasta_file = $new_fasta_path . '/genome.faa';
	print STDERR "   new fasta: $new_fasta_path/genome.faa\n";
	
	die "   New path already exists, please delete first\n" if (-e $new_fasta_path);
	mkpath($new_fasta_path);
	process_fasta_file($fasta_file, $new_fasta_path);
	
	print STDERR "   indexing ...\n";
	my $rc = CoGe::Core::Storage::index_genome_file( file_path => $new_fasta_file );
	die if ($rc != 0);
	
	die if (not (-r $new_fasta_file and -r "$new_fasta_file.fai"));
	
	$count++;
}

#-------------------------------------------------------------------------------
# All done!
#-------------------------------------------------------------------------------

print STDERR "$count genomes migrated\n";
print STDERR "All done!\n";
exit;

#-------------------------------------------------------------------------------

sub process_fasta_file {
    my $filepath   = shift;
    my $target_dir = shift;

    $/ = "\n>";
    open( my $in, $filepath ) || die "can't open $filepath for reading: $!";

    my $fileSize    = -s $filepath;
    my $lineNum     = 0;
    my $totalLength = 0;
    while (<$in>) {
        s/\n*\>//g;
        next unless $_;
        my ( $name, $seq ) = split /\n/, $_, 2;
        #$seq =~ s/\n//g;
        $seq =~ s/\s//g; # because of Windows end-of-line
        $lineNum++;

        my ($chr) = split( /\s+/, $name );
        die "Blank fasta header" if ( not defined $chr );

        my ($filename) = $filepath =~ /^.+\/([^\/]+)$/;

        # Append sequence to master file
        open( my $out, ">>$target_dir/genome.faa" );
        print_fasta($out, $chr, \$seq); #print $out "$head\n$seq\n";
        close($out);
    }
    close($in);

    return $totalLength;
}
