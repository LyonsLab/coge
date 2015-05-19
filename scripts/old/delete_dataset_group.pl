#!/usr/bin/perl -w

use strict;
use CoGeX;
use Getopt::Long;

my ( $dsgid, $delete_seqs, $db, $user, $pass );

GetOptions(
	"dsgid=i"       => \$dsgid,
	"delete_seqs"   => \$delete_seqs,
	"database|db=s" => \$db,
	"user|u=s"      => \$user,
	"password|pw=s" => \$pass,
);

unless ( $dsgid && $db && $user && $pass ) {
	print qq{
welcome to $0

Usage:  $0 -dsgid <database id for genome> -delete_seqs -db <database name> -u <database user name> -pw <database password>

This will delete the database entry for the genome, associated entries in the database_connector table, and the datasets if the datasets do not belong to another genome!

To delete the genomic sequences associted with the genome, add flag -delete_seqs
};
	exit;
}

my $connstr = "dbi:mysql:dbname=$db;host=localhost;port=PORT";
my $coge    = CoGeX->connect( $connstr, $user, $pass );

my $dsg = $coge->resultset('Genome')->find($dsgid);
unless ($dsg) {
	print "unable to find entry for $dsgid\n";
	exit;
}

delete_dsg($dsg);
print "Finished deleting genome\n";

sub delete_dsg {
	my $dsg = shift;
	foreach my $ds ( $dsg->datasets ) {
		my %dsgs;
		foreach my $item ( $ds->genomes ) {
			next if $item->id eq $dsg->id;
			$dsgs{ $item->id } = $item;
		}
		next if keys %dsgs;
		print "Deleting dataset " . $ds->name, " (" . $ds->id . ")" . "\n";
		$ds->delete();
	}
	if ($delete_seqs) {
		my $path = $dsg->file_path;
		$path =~ s/[^\/]*$//;
		my $cmd = "rm -rf $path";
		print "Deleting genomic sequence: $cmd\n";
		`$cmd`;
	}
	else {
		print "DID NOT DELETE GENOMIC SEQUENCES: ", $dsg->file_path, "\n";
	}
	$dsg->delete();
}
