#!/usr/bin/perl -w

use strict;
use Getopt::Long;
use CoGe::Accessory::Web;

use vars qw($config $gid1 $gid2 $P $diags $fasta $blastdb $lastdb);

GetOptions(
    "gid1=i"  =>   \$gid1,   #genome id 1
    "gid2=i"  =>   \$gid2,   #genome id 2
    "config=s"=>   \$config, #coge config file
);

unless ($gid1 && $gid2 && $config) {
print STDERR qq{
Usage: $0 -gid1 <genome id 1> -gid2 <genome id 2> -config <coge configuration file>

This program deletes all the files used to generate a synmap analysis including:

  CDS fasta files for each genome
  blast databases
  SynMap Analysis results
};
exit;
}

#synmap orders the genome ids with the smaller one first
if ($gid1 > $gid2) {
    my $tmp = $gid1;
    $gid1 = $gid2;
    $gid2 = $tmp;
}

# Load config file
unless ($config) {
    print STDOUT "log: error: can't find config file\n";
    print STDERR "can't find config file\n";
    exit(-1);
}
$P    = CoGe::Accessory::Web::get_defaults($config);
$diags   = $P->{DIAGSDIR};
$fasta = $P->{FASTADIR};
$blastdb = $P->{BLASTDB};
$lastdb = $P->{LASTDB};

print STDERR join "\n", "GID1: $gid1", "GID2: $gid2", "DIAGSDIR: $diags", "FASTADIR: $fasta", "BLASTDB_DIR: $blastdb", "LASTDB_DIR: $lastdb";
print STDERR "\n\n";

my $dir;
my $file;
my $cmd;

print STDERR "Running rm commands:\n\n";

#Deleting SynMap analysis directory (everything goes)
$dir = $diags.$gid1."/".$gid2;
$cmd = "rm -rf $dir";
print STDERR $cmd,"\n";
system $cmd if -r $dir;

#Deleting fatsa files
#files may be <gid>-CDS.fasta and <gid>-protein.fasta;
$file = $fasta.$gid1."-CDS.fasta";
$cmd = "rm -f $file";
print STDERR $cmd,"\n";
system $cmd if -r $file;
$file = $fasta.$gid1."-protein.fasta";
$cmd = "rm -f $file";
print STDERR $cmd,"\n";
system $cmd if -r $file;
$file = $fasta.$gid2."-CDS.fasta";
$cmd = "rm -f $file";
print STDERR $cmd,"\n";
system $cmd if -r $file;
$file = $fasta.$gid2."-protein.fasta";
$cmd = "rm -f $file";
print STDERR $cmd,"\n";
system $cmd if -r $file;

#Delete blastdb files
#for some reason there are directories and files here, need to clean all of them up
$dir = $blastdb.$gid1;
$cmd = "rm -rf $dir";
print STDERR $cmd,"\n";
system $cmd if -r $dir;
$dir = $blastdb.$gid2;
$cmd = "rm -rf $dir";
print STDERR $cmd,"\n";
system $cmd if -r $dir;
$file = $blastdb.$gid1."-*";
$cmd = "rm -f $file";
print STDERR $cmd,"\n";
$file = $blastdb.$gid2."-*";
$cmd = "rm -f $file";
print STDERR $cmd,"\n";
 
#Deleting lastdb files
#These are all in directories;
$dir = $lastdb.$gid1;
$cmd = "rm -rf $dir";
print STDERR $cmd,"\n";
system $cmd if -r $dir;
$dir = $lastdb.$gid2;
$cmd = "rm -rf $dir";
print STDERR $cmd,"\n";
system $cmd if -r $dir;

