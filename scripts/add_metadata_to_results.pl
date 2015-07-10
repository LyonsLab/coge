#!/usr/bin/env perl
use v5.14;
use strict;
use warnings;

use File::Path qw(mkpath);
use Getopt::Long qw(GetOptions);

use CoGe::Accessory::Utils qw(to_pathname);
use CoGe::Accessory::Web qw(get_defaults);
use CoGe::Core::Metadata qw(create_annotations);
use CoGe::Core::Storage qw(get_workflow_results);
use CoGeX;

our ($userid, $wid, $log_file, $config_file, $annotations);

GetOptions(
    # Required workflow params
    "userid|uid=s"      => \$userid, # User creating the notebook
    "wid=s"             => \$wid,    # Workflow ID

    # Input/output files
    "log=s"             => \$log_file,
    "config|cfg=s"      => \$config_file,

    # Metadata (aka "annotations")
    "annotations=s"     => \$annotations # optional: semicolon-separated list of locked annotations (link:group:type:text;...)
);

$| = 1;
print STDOUT "Starting $0 (pid $$)\n", qx/ps -o args $$/;

# Check required parameters
die "ERROR: userid not specified, use userid argument" unless $userid;
die "ERROR: wid not specified" unless $wid;
die "ERROR: annotations not specified" unless $annotations;

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

# Add annotations to workflow results (genomes, experiments, notebooks)
foreach my $result (@$results) {
    CoGe::Core::Metadata::create_annotations(
        db => $db, 
        target_id => $result->{id}, 
        target_type => $result->{type}, 
        annotations => $annotations, 
        locked => 1
    );
}

# Create log file -- signals task completion to JEX
my $log_path = to_pathname($log_file);
mkpath($log_path);
open(my $fh, ">>", $log_file);
say $fh "All done!";
close($fh);

exit;
