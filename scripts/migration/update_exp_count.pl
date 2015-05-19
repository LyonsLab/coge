#!/usr/bin/perl
#-------------------------------------------------------------------------------
# Usage:
#    ./update_exp_count.pl -database XXXXXXX -user XXXXXXX -password XXXXXXX
#
# Populate new row_count column for existing experiments.
#
#-------------------------------------------------------------------------------

use DBI;
use strict;
use CoGeX;
use Data::Dumper;
use Getopt::Long;
use File::Path;

use vars qw($DEBUG $coge);

my ($db, $user, $pass, $go, $cleanup);
GetOptions (
	"debug=s"			=> \$DEBUG,
	"database|db=s"		=> \$db,
	"user|u=s"			=> \$user,
	"password|pw|p=s"	=> \$pass,
	"go=i"				=> \$go
);
die "Missing DB params\n" unless ($db and $user and $pass);

$| = 1;
$DEBUG = 1 unless defined $DEBUG; # set to '1' to get updates on what's going on

print STDERR "Running $0\n";
print STDERR "PERL5LIB=" . $ENV{PERL5LIB} . "\n";

#-------------------------------------------------------------------------------
# Connect to database
#-------------------------------------------------------------------------------

my $connstr = "dbi:mysql:dbname=$db;host=localhost;port=3307";
$coge = CoGeX->connect($connstr, $user, $pass);
#$coge->storage->debugobj(new DBIxProfiler());
#$coge->storage->debug(1);

my $dbh = $coge->storage->dbh;

#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------

foreach my $e ($coge->resultset('Experiment')->all) {#search({row_count => undef})) {
	my $cmd = "ibis -d " . $e->storage_path . " -q \"select count(*) where 0.0=0.0\" 2>&1";
	my @cmdOut = qx{$cmd};

	my $row_count = 0;
	foreach (@cmdOut) {
		chomp;
		if (/^(\d+)/) {
			$row_count = $1;
			last;
		}
	}

	$e->row_count($row_count);
	$e->update;

	print STDERR 'Experiment id' . $e->id . " row_count set to $row_count\n";
}

#-------------------------------------------------------------------------------
# All done!
#-------------------------------------------------------------------------------

$coge->resultset('Log')->create( { user_id => 0, page => $0, description => 'All experiment row_counts updated' } );
print STDERR "All done!\n";
exit;
