#!/usr/bin/env perl
use v5.14;
use strict;
use warnings;
use Data::Dumper qw(Dumper);

use File::Basename qw(fileparse basename);
use File::Path qw(mkpath);
use File::Spec::Functions qw(catdir catfile);
use Getopt::Long qw(GetOptions);
use JSON qw(decode_json);
use URI::Escape::JavaScript qw(unescape);

use CoGe::Accessory::Web qw(get_defaults);
use CoGe::Accessory::Utils qw(to_pathname);
use CoGe::Core::Notebook qw(create_notebook);
use CoGe::Core::Metadata qw(create_annotations);
use CoGeX;

our ($LOG, $DEBUG, $PAGE, $P, $db, $host, $port, $user, $pass, $config_file,
     $name, $description, $version, $type, $userid, $restricted, $result_file,
     $annotations, @ITEMS);

GetOptions(
    "userid|uid=s"      => \$userid, # User creating the notebook

    # General configuration options
    "log=s"             => \$LOG,
    "config|cfg=s"      => \$config_file,
    "result_file=s"     => \$result_file, # results file

    # Notebook options
    "page=s"            => \$PAGE, # The reference page
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
die "ERROR: no items found" unless @ARGV;
die "ERROR: notebook name not specified" unless $name;

# Set default parameters
$restricted  = '1' unless (defined $restricted && (lc($restricted) eq 'false' || $restricted eq '0'));

sub parse_files {
    my @logs = @_;
    my @ids;
    my $node_types = CoGeX::node_types();

    for my $log (@logs) {
        open( my $fh, $log );

        while (<$fh>) {
            if ( $_ =~ /experiment id: (\d+)/i ) {
                push @ids, [$1, $node_types->{experiment}];
            }
            elsif ( $_ =~ /log: error/i ) {
                return;
            }
        }

        close($fh);
    }

    return @ids;
}

sub main {
    # Items being added to the notebook
    @ITEMS = @ARGV;
    unless (@ITEMS) {
        print STDOUT "No log paths specified\n";
        exit(-1);
    }
    my @ids = parse_files(@ITEMS);
    
    # Connect to DB
    my $conf = get_defaults($config_file);
    my $db = CoGeX->dbconnect($conf);
    die "ERROR: couldn't connect to the database" unless $db;

    # Get user
    my $user = $db->resultset("User")->find($userid);
    die "ERROR: user could not be found for id: $userid" unless $user;

    # Create notebook
    my $notebook = create_notebook(
        db         => $db,
        user       => $user,
        name       => $name,
        desc       => $description,
        type_id    => $type,
        item_list  => \@ids,
        restricted => $restricted
    );
    unless ($notebook) {
        print STDERR "ERROR: couldn't create notebook\n";
        exit(-1);
    }

    # Create annotations
    if ($annotations) {
        CoGe::Core::Metadata::create_annotations(db => $db, target => $notebook, annotations => $annotations, locked => 1);
    }

    # Save notebook ID in log file
    open(my $fh, ">>", $LOG);
    say $fh "notebook id: " . $notebook->id;
    close($fh);

    # Save result document
    if ($result_file) {
        my $result_dir = to_pathname($result_file);
        mkpath($result_dir);
        
        my $rc = CoGe::Accessory::TDS::write(
            $result_file,
            {
                type => 'notebook',
                id => int($notebook->id),
                name => $name,
                description => $description,
                type_id => $type,
                restricted => $restricted
            }
        );
        unless ($rc) {
            print STDOUT "Error: couldn't write result file\n";
            exit(-1);
        }
    }
}

main;
