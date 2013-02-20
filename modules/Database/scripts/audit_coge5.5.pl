#!/usr/bin/perl
#-------------------------------------------------------------------------------
# Usage:
#    ./audit_coge5.5.pl -database XXXXXXX -user XXXXXXX -password XXXXXXX
#
#	perl -I /home/mbomhoff/perl5 scripts/audit_coge5.5.pl -db coge_matt -u XXXXXXX -pw XXXXXXX
#
# Scan the database for discrepancies:
#    - 
#
#-------------------------------------------------------------------------------

use DBI;
use strict;
use CoGeX;
use Roman;
use Data::Dumper;
use Getopt::Long;
use File::Path;

use vars qw($DEBUG $coge);

my ($db, $user, $pass, $go, $cleanup);
GetOptions (
	"debug=s"		=> \$DEBUG,
	"database|db=s"	=> \$db,
	"user|u=s"		=> \$user,
	"password|pw=s"	=> \$pass,
	"go=i"			=> \$go
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
# Audit 
#-------------------------------------------------------------------------------

foreach my $group ($coge->resultset('UserGroup')->all) {
	# Check that each group has a creator
	unless ($group->creator) {
		print STDERR "Group ".$group->id." ('" . $group->name . "') is missing creator user\n";	
	}	
	
	# Check that each group has an owner user connection
	unless ($group->owner) {
		print STDERR "Group ".$group->id." ('" . $group->name . "') is missing owner user\n";	
	}
	
	# Check that non-owner users' roles match group role
#	foreach my $user ($group->users) {
#		print STDERR "User ".$user->id."\n";
#	}
}

my (%missing_parent, %missing_child);
foreach my $conn ($coge->resultset('ListConnector')->all) {
	# Check that each list connector has nodes
	unless ($conn->parent_list) {
		push @{$missing_parent{$conn->parent_id}}, $conn->id;
	}
	unless ($conn->child) {
		push @{$missing_child{$conn->child_id}}, $conn->id;
	}
}

foreach my $lid (keys %missing_parent) {
	print STDERR "Parent list " . $lid . " is missing for list connectors: " . join(', ', @{$missing_parent{$lid}}) . "\n";	
}
foreach my $child_id (keys %missing_child) {
	print STDERR "Child  " . $child_id . " is missing for list connectors: " . join(', ', @{$missing_child{$child_id}}) . "\n";	
}

#-------------------------------------------------------------------------------
# All done!
#-------------------------------------------------------------------------------

$coge->resultset('Log')->create( { user_id => 0, page => $0, description => 'database migrated to coge5.5' } );
print STDERR "All done!\n";
exit;

