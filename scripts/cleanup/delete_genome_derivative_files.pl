#!/usr/bin/perl -w

use strict;
use warnings;

use DBI;
use File::Spec;

use CoGe::Accessory::Web;

open STDERR, '>', File::Spec->devnull() or die "could not open STDERR: $!\n";

my ($db, $user, $conf) = CoGe::Accessory::Web->init();

my $bed_dir = $conf->{BEDDIR};
my $blast_dir = $conf->{BLASTDB};
my $last_dir = $conf->{LASTDB};
my $cache_dir = $conf->{CACHEDIR} . '/genomes/';
my $fasta_dir = $conf->{FASTADIR};
my $diags_dir = $conf->{DIAGSDIR};

my $genome_id = shift;
if ($genome_id eq 'rm') {
	$genome_id = shift;
	if ($genome_id) {
		rm($genome_id);
	} else {
		my $ids = $db->storage->dbh->selectcol_arrayref('SELECT genome_id FROM genome WHERE deleted');
		for (@$ids) {
			rm($_);
		}
	}
} else {
	if ($genome_id) {
		ls($genome_id);
	} else {
		my $ids = $db->storage->dbh->selectcol_arrayref('SELECT genome_id FROM genome WHERE deleted');
		for (@$ids) {
			ls($_);
		}
	}
}

sub run {
	my $cmd = shift;
	my @results = `$cmd`;
	if (@results) {
		print "\n----- ",$cmd,"\n";
		for (@results) {
			chomp;
			print $_,' ';
		}
	}
}

sub ls {
	my $genome_id = shift;
	run "ls $bed_dir$genome_id.*" if -e "$bed_dir$genome_id.*";
	run "ls $blast_dir$genome_id" if -e "$blast_dir$genome_id";
	run "ls $last_dir$genome_id" if -e "$last_dir$genome_id";
	run "ls $cache_dir$genome_id" if -e "$cache_dir$genome_id";
	run "ls $fasta_dir$genome_id-*";
	run "ls $diags_dir$genome_id" if -e "$diags_dir$genome_id";
	run "ls -d $diags_dir*/$genome_id";
}

sub rm {
	my $genome_id = shift;
	system("rm -f $bed_dir$genome_id.*");
	system("rm -rf $blast_dir$genome_id");
	system("rm -f $blast_dir$genome_id-*");
	system("rm -rf $last_dir$genome_id");
	system("rm -rf $cache_dir$genome_id");
	system("rm -f $fasta_dir$genome_id-*");
	system("rm -rf $diags_dir$genome_id");
	my @dirs = `ls -d $diags_dir*/$genome_id`;
	for (@dirs) {
		system("rm -rf $_");
	}
}
