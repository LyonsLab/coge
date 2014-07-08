#!/usr/bin/perl
#-------------------------------------------------------------------------------
# Usage:
#    ./migrate_coge5.pl -database XXXXXXX -user XXXXXXX -password XXXXXXX
#
# NOTE: CoGeX code must be updated to v5 BEFORE running this script!
#
# This script starts by migrating existing table entries, particularly for
# user_group and user_group_data_connector (whose contents are moving into the
# new list tables).
#
#ÊAfter the hard stuff is finished, minor modifications are made,
# mostly just renaming tables and extended indexes from 10 to 11 digits.
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

my ($db, $user, $pass);
GetOptions (
	"debug=s" => \$DEBUG,
	"database|db=s"=>\$db,
	"user|u=s"=>\$user,
	"password|pw=s"=>\$pass
);

$| = 1;
$DEBUG = 1 unless defined $DEBUG; # set to '1' to get updates on what's going on

print STDERR "Running $0\n";

#-------------------------------------------------------------------------------
# Connect to database
#-------------------------------------------------------------------------------

my $connstr = "dbi:mysql:dbname=$db;host=localhost;port=3307";
$coge = CoGeX->connect($connstr, $user, $pass);
#$coge->storage->debugobj(new DBIxProfiler());
#$coge->storage->debug(1);

my $dbh = $coge->storage->dbh;

#-------------------------------------------------------------------------------
# Migrate data:  remap genome-to-user into new list tables
#-------------------------------------------------------------------------------

# Update user_group table
drop_column('user_group', 'creator_user_id');
add_column('user_group', 'creator_user_id INT(11) NOT NULL AFTER user_group_id'); # add creator_user_id
sql('update user_group set creator_user_id=0');
drop_column('user_group', 'locked');
add_column('user_group', 'locked INT(1) NOT NULL DEFAULT 0'); # add locked
sql('alter table user_group change name name VARCHAR(255) NOT NULL'); # change from VARCHAR(50) to VARCHAR(255)
sql('alter table user_group change description description VARCHAR(1024) NOT NULL'); # change from VARCHAR(255) to VARCHAR(1024)

# Create list_connector table
sql('drop table if exists list_connector');
sql(<<SQL);
create table list_connector
(
  list_connector_id int(11) NOT NULL AUTO_INCREMENT,
  parent_id int(11) DEFAULT NULL,
  child_id int(11) NOT NULL,
  child_type tinyint(1) NOT NULL,
  PRIMARY KEY (list_connector_id),
  INDEX (parent_id, child_id)
);
SQL

# Create list_type table
sql('drop table if exists list_type');
sql(<<SQL);
create table list_type
(
  list_type_id int(11) NOT NULL AUTO_INCREMENT,
  name varchar(255) NOT NULL,
  description varchar(1024) DEFAULT NULL,
  PRIMARY KEY (list_type_id)
);
SQL
sql("insert list_type set name = 'Genome'");
sql("insert list_type set name = 'Experiment'");
sql("insert list_type set name = 'Owner'");
sql("insert list_type set name = 'Feature'");
sql("insert list_type set name = 'Mixed'");
sql("insert list_type set name = 'Other'");

# Create list table
sql('drop table if exists list');
sql(<<SQL);
CREATE TABLE list (
	list_id int(11) NOT NULL AUTO_INCREMENT,
	name varchar(255) NOT NULL,
	description varchar(1024) DEFAULT NULL,
	list_type_id int(11) NOT NULL,
	restricted tinyint(1) NOT NULL DEFAULT '0',
	locked int(1) NOT NULL DEFAULT '0',
	user_group_id int(11) NOT NULL,
	PRIMARY KEY (list_id),
	INDEX (name),
	INDEX (user_group_id)
);
SQL

# For each existing group, create a list with its assigned genomes
foreach my $group ($coge->resultset('UserGroup')->all) {
	# We changed the meaning of "owner" group
	if ($group->role->id == 2) { # is "old" owner group?
		$group->role_id(3); # set to editor
		$group->update;
	}

	my $list = $coge->resultset('List')->create(
	  { name => $group->name,
		description => $group->description,
		list_type_id => 1, # genome type
		user_group_id => $group->id,
		restricted => 1,
		locked => 0
	  } );
	die unless $list;

	my $sth = $dbh->prepare('SELECT user_group_dataset_connector_id, dataset_group_id FROM user_group_data_connector where user_group_id = ' . $group->id);
	$sth->execute();
	while (my ($ugdcid, $dsgid) = $sth->fetchrow_array()) {
		my $conn = $coge->resultset('ListConnector')->create(
		  { parent_id => $list->id,
			child_id => $dsgid,
			child_type => 2, # genome type
		  } );
		die unless $conn;
	}
	$sth->finish();
}

# Create an owner group & list for all existing users
foreach my $user ($coge->resultset('User')->all) {
	my $group;

	unless ($user->owner_group) {
		$group = $coge->resultset('UserGroup')->create(
		  { creator_user_id => $user->id,
		  	name => $user->name,
			description => 'Owner group',
			role_id => 2, # okay to hardcode owner role?
			locked => 1
		  } );
		die unless $group;

		my $conn = $coge->resultset('UserGroupConnector')->create(
		  { user_id => $user->id,
		  	user_group_id => $group->id,
		  } );
		die unless $conn;
	}

	my $list;
	if ($list = $user->owner_list) {
		$list->user_group_id($group);
		$list->update();
	}
	else {
		$list = $coge->resultset('List')->create(
		  { name => $user->name,
			description => 'Owner list',
			list_type_id => 3, # owner type
			user_group_id => $group->id,
			restricted => 1,
			locked => 1
		  } );
		die unless $list;
	}
}

#-------------------------------------------------------------------------------
# Remove obsolete tables
#-------------------------------------------------------------------------------

sql('drop table if exists list_group');
sql('drop table if exists list_group_image_connector');
sql('drop table if exists user_group_feature_list_permission_connector');
sql('drop table if exists quantitation');
sql('drop table if exists quantitation_experiment');
sql('drop table if exists sequence');
sql('drop table if exists sequence_type');
#sql('drop table if exists user_group_data_connector');

#-------------------------------------------------------------------------------
# Create new tables
#-------------------------------------------------------------------------------

# Create list_annotation table
sql('drop table if exists list_annotation');
sql(<<SQL);
create table list_annotation
(
	list_annotation_id INT(11) NOT NULL AUTO_INCREMENT,
	list_id INT(11) NOT NULL,
	annotation_type_id INT(11) NOT NULL,
	annotation TEXT NOT NULL,
	link TEXT DEFAULT NULL,
	image_id INT(11) DEFAULT NULL,
	PRIMARY KEY (list_annotation_id),
	INDEX (list_id),
	INDEX (annotation_type_id),
	FULLTEXT (annotation)
);
SQL

# image table
sql('drop table if exists image');
sql(<<SQL);
create table image
(
	image_id int(11) NOT NULL AUTO_INCREMENT,
    filename varchar(255) NOT NULL,
    image longblob NOT NULL,
    PRIMARY KEY (image_id)
);
SQL

# log table
sql('drop table if exists log');
sql(<<SQL);
create table log
(
	log_id int(11) NOT NULL AUTO_INCREMENT,
	time timestamp NOT NULL DEFAULT CURRENT_TIMESTAMP,
    user_id int(11) NOT NULL DEFAULT 0,
    page varchar(255) NOT NULL,
    description varchar(255) DEFAULT NULL,
    link varchar(255) DEFAULT NULL,
    status int(1) DEFAULT 0,
    comment varchar(255) DEFAULT NULL,
    PRIMARY (log_id),
    INDEX (user_id)
);
SQL

#-------------------------------------------------------------------------------
# Rename and retype some tables/columns
#-------------------------------------------------------------------------------

# dataset_group table
sql('rename table dataset_group to genome'); # rename table
sql('alter table genome change dataset_group_id genome_id INT(11) NOT NULL AUTO_INCREMENT'); # rename dataset_group_id to genome_id

# feature_annotation table
sql('rename table annotation to feature_annotation');
sql('alter table feature_annotation change annotation_id feature_annotation_id INT(11) NOT NULL AUTO_INCREMENT'); # rename annotation_id to feature_annotation_id
#add_column('feature_annotation', 'image_id INT(11) DEFAULT NULL');

# data_source table
sql('alter table data_source change data_source_id data_source_id INT(11) NOT NULL AUTO_INCREMENT'); # change from INT(10) to INT(11)

# dataset table
sql('alter table dataset change name name VARCHAR(255) NOT NULL'); # change from VARCHAR(100) to VARCHAR(255)
sql('alter table dataset change description description VARCHAR(1024)'); # change from VARCHAR(255) to VARCHAR(1024)

# dataset_connector table
sql('alter table dataset_connector change dataset_group_id genome_id INT(11) NOT NULL'); # rename dataset_group_id to genome_id

# experiment table
sql('alter table experiment change dataset_group_id genome_id INT(11) NOT NULL'); # rename dataset_group_id to genome_id
drop_column('experiment', 'link'); #sql('alter table experiment drop column link'); # remove link
add_column('date timestamp NOT NULL DEFAULT CURRENT_TIMESTAMP,');

# experiment_annotation table
add_column('experiment_annotation', 'link TEXT'); #sql('alter table experiment_annotation add link TEXT');
add_column('experiment_annotation', 'image_id INT(11) DEFAULT NULL');

# feature table
sql('alter table feature change feature_id feature_id INT(11) NOT NULL AUTO_INCREMENT'); # change from INT(10) to INT(11)
sql('alter table feature change feature_type_id feature_type_id INT(11) NOT NULL'); # change from INT(10) to INT(11)
sql('alter table feature change dataset_id dataset_id INT(11) NOT NULL'); # change from INT(10) to INT(11)

# feature_type table
sql('alter table feature_type change feature_type_id feature_type_id INT(11) NOT NULL AUTO_INCREMENT'); # change from INT(10) to INT(11)

# genomic_sequence table
sql('alter table genomic_sequence change dataset_group_id genome_id INT(11) NOT NULL'); # rename dataset_group_id to genome_id

# genomic_sequence_type table
sql('alter table genomic_sequence_type change genomic_sequence_type_id genomic_sequence_type_id INT(11) NOT NULL AUTO_INCREMENT'); # change from INT(10) to INT(11)

# organism table
sql('alter table organism change organism_id organism_id INT(11) NOT NULL AUTO_INCREMENT'); # change from INT(10) to INT(11)
sql('alter table organism change name name VARCHAR(255) NOT NULL'); # change from VARCHAR(200) to VARCHAR(255)

# user table
sql('alter table user change user_id user_id INT(11) NOT NULL AUTO_INCREMENT'); # change from INT(10) to INT(11)
add_column('user', 'image_id INT(11) DEFAULT NULL');
add_column('date timestamp NOT NULL DEFAULT CURRENT_TIMESTAMP,');

# user_group_connector table
sql('alter table user_group_connector change user_group_connector_id user_group_connector_id INT(11) NOT NULL AUTO_INCREMENT'); # change from INT(10) to INT(11)
sql('alter table user_group_connector change user_id user_id INT(11) NOT NULL'); # change from INT(10) to INT(11)
sql('alter table user_group_connector change user_group_id user_group_id INT(11) NOT NULL'); # change from INT(10) to INT(11)

# user_session table
sql('alter table user_session change user_session_id user_session_id INT(11) NOT NULL AUTO_INCREMENT'); # change from INT(10) to INT(11)
sql('alter table user_session change user_id user_id INT(11) NOT NULL'); # change from INT(10) to INT(11)

# web_preferences table
sql('alter table web_preferences change id id INT(11) NOT NULL AUTO_INCREMENT'); # change from INT(10) to INT(11)
sql('alter table web_preferences change user_id user_id INT(11) NOT NULL'); # change from INT(10) to INT(11)

# work table
sql('alter table work change work_id work_id INT(11) NOT NULL AUTO_INCREMENT'); # change from INT(10) to INT(11)
sql('alter table work change user_id user_id INT(11) NOT NULL'); # change from INT(10) to INT(11)
sql('alter table work change page page VARCHAR(255) NOT NULL'); # change from VARCHAR(100) to VARCHAR(255)
sql('alter table work change name name VARCHAR(255)'); # change from VARCHAR(256) to VARCHAR(255)
sql('alter table work change image_id image_id INT(11)'); # change from INT(10) to INT(11)
sql('alter table work change link link TEXT'); # change from VARCHAR(1024) to TEXT

# work_order table
sql('alter table work_order change work_order_id work_order_id INT(11) NOT NULL AUTO_INCREMENT'); # change from INT(10) to INT(11)
sql('alter table work_order change workflow_id workflow_id INT(11) NOT NULL'); # change from INT(10) to INT(11)
sql('alter table work_order change work_id work_id INT(11) NOT NULL'); # change from INT(10) to INT(11)

# workflow table
sql('alter table workflow change workflow_id workflow_id INT(11) NOT NULL AUTO_INCREMENT'); # change from INT(10) to INT(11)
sql('alter table workflow change name name VARCHAR(255) NOT NULL'); # change from VARCHAR(256) to VARCHAR(255)
sql('alter table workflow change user_id user_id INT(11) NOT NULL'); # change from INT(10) to INT(11)
sql('alter table workflow change link link TEXT'); # change from VARCHAR(1024) to TEXT

#-------------------------------------------------------------------------------
# Assign orphaned restricted genomes to admin owner list
#-------------------------------------------------------------------------------
my $admin_owner_list = $coge->resultset('List')->find({name => 'admin'});
die unless ($admin_owner_list);

# Fix owner list fields
$admin_owner_list->locked(1);
$admin_owner_list->list_type_id(3); # FIXME hardcoded owner list type
$admin_owner_list->update();

# Assign orphan restricted genomes
foreach my $genome ($coge->resultset('Genome')->all) {
	next if (not $genome->restricted);

	my $genome_owner_list;
	foreach my $lc ($genome->list_connectors) {
		if ($lc->parent_list->list_type_id == 3) {
			$genome_owner_list = $lc->parent_list;
			#print STDERR "skipping genome id" . $genome->id . "\n";
			last;
		}
	}

	unless ($genome_owner_list) {
		print STDERR "assigning genome id" . $genome->id . "\n";
		my $conn = $coge->resultset('ListConnector')->create(
		  { parent_id => $admin_owner_list->id,
			child_id => $genome->id,
			child_type => 2, # genome type
		  } );
		die unless $conn;
	}
}

#-------------------------------------------------------------------------------
# All done!
#-------------------------------------------------------------------------------

$coge->resultset('Log')->create( { user_id => 0, page => $0, description => 'database migrated to coge5' } );
print STDERR "All done!\n";
exit;

#-------------------------------------------------------------------------------
sub sql {
	my $cmd = shift;
	print $cmd, "\n" if ($DEBUG);
	$dbh->do($cmd);
}

sub drop_column {
	my $table = shift;
	my $column = shift;

	my $sth = $dbh->prepare("SELECT * FROM INFORMATION_SCHEMA.COLUMNS WHERE TABLE_SCHEMA='$db' AND TABLE_NAME='$table' AND COLUMN_NAME='$column'");
	my $count = $sth->execute();
	if ($count > 0) {
		sql("alter table $table drop $column");
	}
	$sth->finish();
}

sub add_column {
	my $table = shift;
	my $column = shift;

	my $sth = $dbh->prepare("SELECT * FROM INFORMATION_SCHEMA.COLUMNS WHERE TABLE_SCHEMA='$db' AND TABLE_NAME='$table' AND COLUMN_NAME='$column'");
	my $count = $sth->execute();
	if ($count == 0) {
		sql("alter table $table add $column");
	}
	$sth->finish();
}
