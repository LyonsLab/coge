#!/usr/bin/perl
#-------------------------------------------------------------------------------
# Usage:
#    ./index_genomes.pl
#
# Index all fasta files in the genomic_sequences directory.
#
#-------------------------------------------------------------------------------

use DBI;
use strict;
use CoGeX;
use Getopt::Long;
use File::Path;

my ($db, $user, $pass);
GetOptions (
	"database|db=s"		=> \$db,
	"user|u=s"			=> \$user,
	"password|pw|p=s"	=> \$pass
);
die "Missing DB params\n" unless ($db and $user and $pass);

$| = 1;

#-------------------------------------------------------------------------------
# Connect to database
#-------------------------------------------------------------------------------

my $connstr = "dbi:mysql:dbname=$db;host=localhost;port=3307";
my $coge = CoGeX->connect($connstr, $user, $pass);
#$coge->storage->debugobj(new DBIxProfiler());
#$coge->storage->debug(1);

#-------------------------------------------------------------------------------
# Index fasta file for each genome
#-------------------------------------------------------------------------------

my ($count, $skipped);
foreach my $genome ($coge->resultset('Genome')->all) {
	if ($genome->is_indexed) {
		print STDERR "Genome " . $genome->id . " is already indexed, skipping\n";
		next;
	}

	print STDERR "Indexing genome " . $genome->id . "\n";
	my $rc = $genome->create_index();
	if ($rc != 0) {
		print STDERR "rc=$rc, skipping\n";
		$skipped++;
		next;
	}
	$count++;
}

#-------------------------------------------------------------------------------
# All done!
#-------------------------------------------------------------------------------

print STDERR "$count genomes indexed, $skipped skipped\n";
print STDERR "All done!\n";
exit;
