#!/usr/bin/env perl
use v5.14;
use strict;
use warnings;

use File::Path qw(mkpath);
use Getopt::Long qw(GetOptions);
use Data::Dumper;

use CoGe::Accessory::Utils qw(to_pathname);
use CoGe::Accessory::Web qw(get_defaults);
use CoGe::Core::Metadata qw(create_annotations);
use CoGe::Core::Storage qw(get_workflow_results);
use CoGeX;

our ($userid, $wid, $log_file, $config_file, $annotations, $item_id, $item_type, $locked);

GetOptions(
    # Required workflow params
    "userid|uid=s"      => \$userid, # User creating the notebook
    "wid=s"             => \$wid,    # Workflow ID

    # Input/output files
    "log=s"             => \$log_file,
    "config|cfg=s"      => \$config_file,

    # Metadata (aka "annotations")
    "annotations=s"     => \$annotations, # optional: semicolon-separated list of locked annotations (link:group:type:text;...)
    
    # Experiment/Genome/Notebook ID -- optional alternative when no .results file exists
    "item_id=i"         => \$item_id,
    "item_type=s"       => \$item_type,
    
    # Locked or unlocked:  user cannot remove locked metadata items
    "locked=i"          => \$locked
);

$| = 1;
print STDOUT "Starting $0 (pid $$)\n", qx/ps -o args $$/;

# Check required parameters
die "ERROR: userid not specified, use userid argument" unless $userid;
die "ERROR: wid not specified" unless $wid;
die "ERROR: annotations not specified" unless $annotations;
$locked = 1 unless (defined $locked);

# Connect to DB
my $conf = get_defaults($config_file);
my $db = CoGeX->dbconnect($conf);
die "ERROR: couldn't connect to the database" unless $db;

# Get user
my $user = $db->resultset("User")->find($userid);
die "ERROR: user could not be found for id: $userid" unless $user;

# Get items to add to the notebook
my @results;
my $workflow_results = get_workflow_results($user->name, $wid);
push @results, @$workflow_results if ($workflow_results);

# Add additional item specified on command-line
push @results, { id => $item_id, type => $item_type } if ($item_id && $item_type);

# Add annotations to workflow results (genomes, experiments, notebooks)
#print STDERR Dumper $results, "\n";
foreach my $result (@results) {
    next unless ($result->{id} && $result->{type} 
        && ($result->{type} eq 'genome' || $result->{type} eq 'experiment' || $result->{type} eq 'notebook' ));
    CoGe::Core::Metadata::create_annotations(
        db => $db,
        user => $user,
        target_id => $result->{id}, 
        target_type => $result->{type}, 
        annotations => $annotations, 
        locked => $locked
    );
}

# Create log file -- signals task completion to JEX
my $log_path = to_pathname($log_file);
mkpath($log_path);
open(my $fh, ">>", $log_file);
say $fh "All done!";
close($fh);

exit;
