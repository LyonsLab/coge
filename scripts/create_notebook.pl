#!/usr/bin/env perl
use v5.14;
use strict;
use warnings;
use Data::Dumper qw(Dumper);

use File::Path qw(mkpath);
use Getopt::Long qw(GetOptions);

use CoGe::Accessory::Web qw(get_defaults);
use CoGe::Accessory::Utils qw(to_pathname);
use CoGe::Core::Notebook qw(create_notebook);
use CoGe::Core::Metadata qw(create_annotations);
use CoGe::Core::Storage qw(get_workflow_results add_workflow_result);
use CoGeX;

our ($log_file, $config_file, $name, $description, $type, $userid, 
    $restricted, $wid, $annotations);

GetOptions(
    # Required workflow params
    "userid|uid=s"      => \$userid, # User creating the notebook
    "wid=s"             => \$wid,    # Workflow ID

    # Input/output files
    "log=s"             => \$log_file,
    "config|cfg=s"      => \$config_file,

    # Notebook metadata
    "name|n=s"          => \$name,
    "type|t=s"          => \$type,
    "description|d=s"   => \$description,
    "restricted|r=s"    => \$restricted,
    "annotations=s"     => \$annotations, # optional: semicolon-separated list of locked annotations (link:group:type:text;...)
);

$| = 1;
print STDOUT "Starting $0 (pid $$)\n", qx/ps -o args $$/;

# Check required parameters
die "ERROR: user id not specified, use userid argument" unless $userid;
die "ERROR: notebook name not specified" unless $name;

# Set default parameters
$restricted  = '1' unless (defined $restricted && (lc($restricted) eq 'false' || $restricted eq '0'));

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

# Create notebook
my $notebook;
if (@$results > 1) {
    my @items = map { [ $_->{id}, $_->{type} ] } @$results;
    print STDOUT Dumper \@items, "\n";
    $notebook = create_notebook(
        db         => $db,
        user       => $user,
        name       => $name,
        desc       => $description,
        type_id    => $type,
        item_list  => \@items,
        restricted => $restricted
    );
    unless ($notebook) {
        print STDERR "ERROR: couldn't create notebook\n";
        exit(-1);
    }
}

# Add annotations to notebook and results
if ($annotations) {
    if ($notebook) {
        CoGe::Core::Metadata::create_annotations(
            db => $db, 
            target => $notebook, 
            annotations => $annotations, 
            locked => 1
        );
    }
    foreach my $result (@$results) {
        CoGe::Core::Metadata::create_annotations(
            db => $db, 
            target_id => $result->{id}, 
            target_type => $result->{type}, 
            annotations => $annotations, 
            locked => 1
        );
    }
}

# Add notebook to workflow result file
add_workflow_result($user->name, $wid,
    {
        type => 'notebook',
        id => int($notebook->id),
        name => $name,
        description => $description,
        type_id => $type,
        restricted => $restricted
    }
);
    
# Save notebook ID in log file -- signals task completion to JEX
my $log_path = to_pathname($log_file);
mkpath($log_path);
open(my $fh, ">>", $log_file);
say $fh "notebook id: " . $notebook->id;
close($fh);

exit;
