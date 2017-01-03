#!/usr/bin/env perl
use v5.14;
use strict;
use warnings;
use Data::Dumper qw(Dumper);

use File::Path qw(mkpath);
use Getopt::Long qw(GetOptions);

use CoGe::Accessory::Web qw(get_defaults);
use CoGe::Accessory::Utils qw(to_pathname);
use CoGe::Core::Notebook qw(add_items_to_notebook get_notebook);
use CoGe::Core::Storage qw(get_workflow_results add_workflow_result);
use CoGeX;

our ($userid, $wid, $log_file, $config_file, $notebook_id);

GetOptions(
    # Required workflow params
    "userid|uid=s"      => \$userid, # User creating the notebook
    "wid=s"             => \$wid,    # Workflow ID

    # Input/output files
    "log=s"             => \$log_file,
    "config|cfg=s"      => \$config_file,

    # Notebook metadata
    "notebook_id=s"     => \$notebook_id,
);

$| = 1;
print STDOUT "Starting $0 (pid $$)\n", qx/ps -o args $$/;

# Check required parameters
die "ERROR: user id not specified, use userid argument" unless $userid;
die "ERROR: notebook id not specified" unless $notebook_id;

# Connect to DB
my $conf = get_defaults($config_file);
my $db = CoGeX->dbconnect($conf);
die "ERROR: couldn't connect to the database" unless $db;

# Get user
my $user = $db->resultset("User")->find($userid);
die "ERROR: user could not be found for id: $userid" unless $user;

# Get items to add to the notebook
my $results = get_workflow_results($user->name, $wid);
exit unless $results;

# load notebook
my $notebook = get_notebook(
	db		=> $db,
	id		=> $notebook_id,
	user	=> $user
);

my @items = map { [ $_->{id}, $_->{type} ] } @$results;
print STDERR Dumper \@items, "\n";
add_items_to_notebook(
	db			=> $db,
	user		=> $user,
	notebook	=> $notebook,
	item_list	=> \@items
);

# Add notebook to workflow result file
add_workflow_result($user->name, $wid,
    {
        type => 'notebook',
        id => int($notebook->id),
        name => $notebook->name,
        description => $notebook->description,
        #type_id => $notebook->type->id, # mdb removed 1/3/17 COGE-800
        restricted => $notebook->restricted
    }
);
    
# Save notebook ID in log file -- signals task completion to JEX
my $log_path = to_pathname($log_file);
mkpath($log_path);
open(my $fh, ">>", $log_file);
say $fh "notebook id: " . $notebook->id;
close($fh);

exit;
