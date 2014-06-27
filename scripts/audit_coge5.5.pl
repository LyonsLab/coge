#!/usr/bin/perl
#-------------------------------------------------------------------------------
# Usage:
#    ./audit_coge5.5.pl -database XXXXXXX -user XXXXXXX -password XXXXXXX
#
#	perl -I /home/mbomhoff/perl5 scripts/audit_coge5.5.pl -db coge_matt -u XXXXXXX -p XXXXXXX
#
# Scan the database and report inconsistencies.
#
#-------------------------------------------------------------------------------

use strict;
use DBI;
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
# Audit
#-------------------------------------------------------------------------------

# Verify users
print STDERR "Verifying users ---------------------------------------------\n";
foreach my $user ($coge->resultset('User')->all) {
	my %children;
	foreach my $conn ($user->child_connectors) {
		if (++$children{$conn->child_type}{$conn->child_id} == 2) {
			print STDERR "User " . $user->id . " has redundant user_connector for child " . $conn->child_id . "\n";
		}
	}
}

# Verify groups
print STDERR "Verifying groups ---------------------------------------------\n";
foreach my $group ($coge->resultset('UserGroup')->all) {
	my $group_desc = $group->id." ('" . $group->name . "')";

	# Find orphan groups
	unless ($group->user_connectors) {
		print STDERR "Group $group_desc is orphaned\n";
		next;
	}

	# Check that each group has a creator
	unless ($group->creator) {
		print STDERR "Group $group_desc is missing creator user\n";
	}

	# Check that each group has an owner user connection
	unless ($group->owner) {
		print STDERR "Group $group_desc is missing owner user\n";
	}

	# Check that non-owner users' roles match group role
	foreach my $conn ($group->user_connectors) {
		next if ($conn->role_id == 2); # skip owner
		my $user = $conn->parent;
		if ($conn->role_id != $group->role_id) {
			print STDERR "Group $group_desc has role mismatch with user " . $user->id . "\n";
		}
	}
}

print STDERR "Verifying lists ----------------------------------------------\n";
# Find orphan lists
foreach my $list ($coge->resultset('List')->all) {
	unless ($list->user_connectors) {
		print STDERR "List ".$list->id." is orphaned\n";
	}
}

# Find orphan list connectors -- FIXME: THIS IS SLOW
my (%missing_parent, %missing_child);
foreach my $conn ($coge->resultset('ListConnector')->all) {
	# Check that each list connector has nodes
	unless ($conn->parent_list) {
		push @{$missing_parent{$conn->parent_id}}, $conn->id;
		next;
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

print STDERR "Verifying images ---------------------------------------------\n";
# Find orphan images
foreach my $image ($coge->resultset('Image')->all) {
	unless ($image->user or $image->list_annotation) {
		print STDERR "Image ".$image->id." is orphaned\n";
	}
}

print STDERR "Verifying genomes --------------------------------------------\n";
foreach my $genome ($coge->resultset('Genome')->all) {
	#next if ($genome->deleted);

	my $deleted = ($genome->deleted ? '(deleted)' : '');
	unless ($genome->list_connectors or $genome->user_connectors or $genome->group_connectors) {
		print STDERR "Genome ".$genome->id." is orphaned $deleted\n";
		next;
	}
	my $file_path = $genome->file_path;
	unless ($file_path and -r $file_path) {
		print STDERR "Genome ".$genome->id." is missing sequence file $deleted\n   ".$genome->info."\n";
		next;
	}
	my @datasets;
	foreach ($genome->dataset_connectors) {
		unless ($_->dataset) {
			print STDERR "Dataset connector ".$_->id." for genome " .$genome->id. " is missing dataset\n";
			next;
		}
		push @datasets, $_->dataset;
	}
	print STDERR "Genome ".$genome->id." has no datasets\n" if (not @datasets);
}

print STDERR "Verifying experiments ----------------------------------------\n";
foreach my $experiment ($coge->resultset('Experiment')->all) {
	unless ($experiment->list_connectors or $experiment->user_connectors or $experiment->group_connectors) {
		print STDERR "Experiment ".$experiment->id." is orphaned\n";
		next;
	}
}

#-------------------------------------------------------------------------------
# All done!
#-------------------------------------------------------------------------------

$coge->resultset('Log')->create( { user_id => 0, page => $0, description => 'database audited' } );
print STDERR "All done!\n";
exit;
