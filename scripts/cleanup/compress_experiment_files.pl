#!/usr/bin/perl -w

use strict;
use warnings;

use DBI;
use File::Basename;
use File::Spec;

use CoGe::Accessory::Web;
use CoGe::Core::Storage qw( get_experiment_path );

#open STDERR, '>', File::Spec->devnull() or die "could not open STDERR: $!\n";

my ($db, $user, $conf) = CoGe::Accessory::Web->init();

my $experiment_id = shift;
if ($experiment_id && $experiment_id eq 'c') {
	$experiment_id = shift;
	if ($experiment_id) {
		compress($experiment_id);
	} else {
		my $ids = $db->storage->dbh->selectcol_arrayref('SELECT experiment_id FROM experiment WHERE deleted');
		for (@$ids) {
			compress($_);
		}
	}
} else {
	if ($experiment_id) {
		ls($experiment_id);
	} else {
		my $ids = $db->storage->dbh->selectcol_arrayref('SELECT experiment_id FROM experiment WHERE deleted');
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
	my $experiment_id = shift;
	my $path = get_experiment_path($experiment_id);
	if ($path) {
		run "ls $path" if -e $path;
	}
}

sub compress {
	my $experiment_id = shift;
	my $path = get_experiment_path($experiment_id);
	return unless $path && -e $path;
	system('tar czf ' . dirname($path) . '/' . $experiment_id . '.tar.gz --anchored ' . $path);
	system('rm -rf ' . $path);
}
