#!/usr/bin/perl
#-------------------------------------------------------------------------------
# Verify completeness and correctness of backup.  Uses md5sum to check file
# contents.
#
# Usage:
#    ./audit_backup.pl <source_dir> <backup_dir>
#	 e.g., ./audit_backup.pl /storage/coge/data/genomic_sequence /iplant/home/coge/backup/genomic_sequence
#    source_dir is on local filesystem
#    backup_dir is in IRODS
#
# Created 9/27/13 by mdb
#-------------------------------------------------------------------------------

use strict;
use warnings;
use File::Find;
$| = 1;

my $source_dir = shift;
my $backup_dir = shift;
die "Missing arguments\n" if (!$source_dir or !$backup_dir);

# Traverse source directory tree to build list of files
my @files;
find(
	sub {
		push @files, $File::Find::name if (not -d $File::Find::name);
	},
	( $source_dir )
);

# Compare md5sum for each file to backup copy
foreach my $file (sort @files) {
	print STDERR "Comparing $file\n";
	my $backup_file = $file;
	$backup_file =~ s/$source_dir//;
	$backup_file = $backup_dir . $backup_file;
	print STDERR "       to $backup_file\n";

	my $md5_source = local_md5sum($file);
	print STDERR "   source: $md5_source\n" if ($md5_source);

	my $md5_backup = irods_md5sum($backup_file);
	print STDERR "   backup: $md5_backup\n" if ($md5_backup);

	die "Mismatch detected!\n" if (!$md5_backup or !$md5_source or $md5_backup ne $md5_source);
}

# All done!
print STDERR @files . " files verified\n";
print STDERR "All done!\n";
exit;

#-------------------------------------------------------------------------------
sub local_md5sum {
	my $filepath = shift;

	my $cmd = "cat $filepath | md5sum";
    my $out = qx{ $cmd };
    if ( $? != 0 ) {
        print STDERR "   ichksum command failed with rc=$?: $cmd\n";
        return;
    }

	my ($md5sum) = $out =~ /^(\w+)/;
	return $md5sum;
}

sub irods_md5sum { #FIXME dup of CoGe::Accessory::IRODS::irods_checksum()
	my $filepath = shift;

	my $cmd = "ichksum $filepath";
    my @out = qx{ $cmd };
    if ( $? != 0 ) {
        print STDERR "   ichksum command failed with rc=$?: $cmd\n";
        return;
    }

	my ($md5sum) = $out[0] =~ /^\s+\S+\s+(\w+)/;
	return $md5sum;
}
