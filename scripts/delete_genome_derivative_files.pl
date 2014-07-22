#!/usr/bin/perl -w

use strict;

my $genome_id = shift;

unless ($genome_id)
  {
    help();
    exit;
  }

my $bed_dir = "/storage/coge/data/bed";
my $blast_dir = "/storage/coge/data/blast/db";
my $last_dir = "/storage/coge/data/last/db";
my $cache_dir = "/storage/coge/data/cache";
my $fasta_dir = "/storage/coge/data/fasta";
my $diags_dir = "/storage/coge/diags";

my $cmd;
$cmd = "rm -f $bed_dir/$genome_id.*;";
print "$cmd","\n";

$cmd = "rm -rf $blast_dir/$genome_id;";
print "$cmd","\n";

$cmd = "rm -f $blast_dir/$genome_id-*;";
print "$cmd","\n";

$cmd = "rm -rf $last_dir/$genome_id;";
print "$cmd","\n";

$cmd = "rm -f $last_dir/$genome_id-*;";
print "$cmd","\n";

$cmd = "rm -rf $cache_dir/$genome_id;";
print "$cmd","\n";

$cmd = "rm -f $fasta_dir/$genome_id-*;";
print "$cmd","\n";

$cmd = "rm -rf $diags_dir/$genome_id;";
print "$cmd","\n";

opendir (DIR, $diags_dir) || die "Can't open $diags_dir: $!";
while (my $dir = readdir (DIR))
  {
    next if $dir =~ /^\.\.?$/;
    next unless -d "$diags_dir/$dir";
    opendir (DIR2, "$diags_dir/$dir");
    while (my $id = readdir(DIR2))
      {
	if ($id eq $genome_id)
	  {
	    $cmd = "rm -rf $diags_dir/$dir/$id;";
	    print $cmd,"\n";
	  }
      }
  }



sub help
  {

print qq{
This program generates the command to delete all derivative files for a genome in CoGe.  IT DOES NOT ACTUALLY RUN THE COMMANDS.  You will need to copy and paste these to perform the deletes.  This should not delete any primary data (e.g., genomic_sequence).

Files to be deleted:

/storage/coge/data/bed/
  file format: genome_id.bed

/storage/coge/data/blast/db/
  subdirectory with genome_id needs to be deleted
  Note: Currently an inconsistency where blastdb files are being deposited in the main directory and not a subdirectory

/storage/coge/data/last/db/
  subdirectory with genome_id needs to be deleted
  Note: Currently an inconsistency where blastdb files are being deposited in the main directory and not a subdirectory

/storage/coge/data/cache/
  subdirectory with genome_id needs to be deleted

/storage/coge/data/fasta/
  files are created with genome_id as the predicate
  need to delete: genome_id*

/storage/coge/diags/
  data are storage in subdirectories: /genome_id_1/genome_id_2/
  All directories need to be crawled to find cases where genome_id_2 exists
  Entire contents of directory need to be deleted
};
  }
