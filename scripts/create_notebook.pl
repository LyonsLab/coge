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

use CoGe::Accessory::Utils qw(to_pathname);
use CoGe::Core::Notebook qw(create_notebook);
use CoGe::Core::Metadata qw(create_annotations);
use CoGeX;

our ($LOG, $DEBUG, $PAGE, $P, $db, $host, $port, $user, $pass, $config,
     $name, $description, $version, $type, $userid, $restricted, $result_file,
     $annotations, @ITEMS);

GetOptions(
    "userid|uid=s"      => \$userid, # User creating the notebook

    # General configuration options
    "log=s"             => \$LOG,
    "config|cfg=s"      => \$config,
    "result_file=s"      => \$result_file, # results file

    # Notebook options
    "page=s"            => \$PAGE, # The reference page
    "name|n=s"          => \$name,
    "type|t=s"          => \$type,
    "description|d=s"   => \$description,
    "restricted|r=s"    => \$restricted,
    "annotations=s"     => \$annotations, # optional: semicolon-separated list of locked annotations (link:group:type:text;...)
);

$| = 1;

# Check required parameters
die "ERROR: user id not specified use userid" unless $userid;
die "ERROR: no items found" unless @ARGV;
die "ERROR: notebook name not specified" unless $name;

# Set default parameters
$restricted  = '0' unless (defined $restricted && lc($restricted) eq 'true');

sub setup {
    $P    = CoGe::Accessory::Web::get_defaults($config);
    $db   = $P->{DBNAME};
    $host = $P->{DBHOST};
    $port = $P->{DBPORT};
    $user = $P->{DBUSER};
    $pass = $P->{DBPASS};

    # Items being added to the notebook
    @ITEMS = @ARGV;

    my $connstr = "dbi:mysql:dbname=$db;host=$host;port=$port;";
    return CoGeX->connect( $connstr, $user, $pass );
}

sub parse_files {
    my @logs = @_;
    my @ids;
    my $node_types = CoGeX::node_types();

    for my $log (@logs) {
        open( my $fh, $log );

        while (<$fh>) {
            if ( $_ =~ /experiment id: (\d+)/i ) {
                #XXX: Hard code to experiment
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
    my $coge = setup();
    die "ERROR: couldn't connect to the database" unless $coge;

    my $user = $coge->resultset("User")->find($userid);
    die "ERROR: user could not be found for id: $userid" unless $user;

    my @ids = parse_files(@ITEMS);

    # Create notebook
    my $notebook = create_notebook(
        db         => $coge,
        user       => $user,
        name       => $name,
        desc       => $description,
        type_id    => $type,
        item_list  => \@ids,
        restricted => $restricted
    );
    exit(-1) unless $notebook;

    # Create annotations
    if ($annotations) {
        CoGe::Core::Metadata::create_annotations(db => $coge, target => $notebook, annotations => $annotations, locked => 1);
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
