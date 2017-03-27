#!/usr/bin/perl
#-------------------------------------------------------------------------------
# Workflow staging path cleanup script
#
# Intended to be run as cron process to monitor staging area disk usage and
# remove old files/directories as needed.  Files/directories in use by currently
# running jobs will not be removed.
#
# Usage:
#    perl cleanup_staging.pl <conf_file> <max_size> <older_than>
#
#    Delete staging subdirectories older than 7 days:
#    perl ./scripts/cleanup/cleanup_staging.pl -conf ./coge.conf -older_than 7
#-------------------------------------------------------------------------------

use strict;
use warnings;

use Data::Dumper;
use Getopt::Long;
use File::Find;
use File::Path qw(remove_tree);
use File::Spec::Functions qw(catdir catfile);
use CoGe::Accessory::Web;
use CoGe::JEX::Jex;

$| = 1;

our ($conf_file, $max_size, $older_than, $job_id, $debug);
GetOptions(
    "conf=s"       => \$conf_file,  # CoGe configuration file #TODO autofind this based on script path
    "max_size=i"   => \$max_size,   # max size of staging directory (in bytes)
    "older_than=i" => \$older_than, # maximum age of subdirectories (in days)
    "job_id=i"     => \$job_id,     # only delete directory for given job ID
    "debug"        => \$debug       # set to 1 to only show output without actually deleting anything
);
die "Missing 'conf' argument\n" unless $conf_file;

# Load config file
unless ($conf_file) {
    print STDERR "Can't find config file: $conf_file\n";
    exit(-1);
}
my $conf = CoGe::Accessory::Web::get_defaults($conf_file);
my $staging_root_path = catdir($conf->{SECTEMPDIR}, 'staging');
my $jex_host = $conf->{JOBSERVER};
my $jex_port = $conf->{JOBPORT};

# Traverse source directory tree to build list of subdirectories
my @dirs;
my $min_depth = () = $staging_root_path =~ /(\/)/g;
my $max_depth = $min_depth + 3;
find(
	sub {
        if (-d $File::Find::name) {
            my $depth = () = $File::Find::name =~ /(\/)/g;
            push @dirs, $File::Find::name if ($depth > $min_depth+1 && $depth < $max_depth);
        }
	},
	( $staging_root_path )
);
unless (@dirs) {
    print STDERR "No staging directories to process.";
    exit;
}

# Connect to JEX
my $jex = CoGe::JEX::Jex->new( host => $jex_host, port => $jex_port );
unless ($jex) {
    print STDERR "Couldn't connect to JEX\n";
    exit(-1);
}

# Delete directories that exceed time and space constraints
my $staging_size = get_dir_size($staging_root_path);
print "Staging: $staging_root_path, size=$staging_size bytes\n";

my $num_deleted = 0;
foreach my $d (sort { (stat $a)[9] <=> (stat $b)[9] } @dirs) { # sort directories by time
    my $dir_size = get_dir_size($d);
    my $dir_age = (stat $d)[9];
    my $id = $d =~ /\/(\d+)$/;
    next if ($job_id && $id ne $job_id);
    print "$d $dir_age $dir_size\n";

    # Check if directory should be deleted
    my $delete = 0;
    if ($max_size && $staging_size > $max_size) {
        print "   DELETE due to size\n";
        $delete = 1;
    }
    if ($older_than && (time - $dir_age) > $older_than*24*60*60) {
        print "   DELETE due to age\n";
        $delete = 1;
    }
    next unless ($delete || $job_id);

    # Check if job is still running
    my $job = $jex->get_job($id);
    unless ($job) {
        print STDERR "Couldn't retrieve job $id\n";
        exit(-1);
    }
    print "   STATUS ", $job->{status}, "\n";
    if (lc($job->{status}) eq 'running') {
        print "Skipping running job $id\n";
        next;
    }

    # Delete directory
    #
    print "   DELETING\n";
    next if $debug;
    remove_tree($d);
    $staging_size -= $dir_size;
    $num_deleted ++;
}

# All done!
print "All done: deleted $num_deleted out of ", scalar(@dirs), " directories" . ($debug ? ' (debug mode)' : '') . "\n";
exit;

#-------------------------------------------------------------------------------

# Return directory size in bytes
sub get_dir_size {
    my $dir = shift;
    my $size = `du -bs $dir | awk '{print \$1}'`;
    chomp($size);
    return $size;
}
