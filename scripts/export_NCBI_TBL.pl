#!/usr/bin/perl -w

use strict;
use CoGeX;
use Getopt::Long;
use File::Copy;
use File::Spec::Functions;
use Sort::Versions;
use URI::Escape::JavaScript qw(unescape);

our ($DEBUG, $db, $user, $pass, $gid, $config, $host, $port, $P,
     $filename, $download_dir);

GetOptions(
    "debug=s"           => \$DEBUG,
    "gid=i"             => \$gid,
    "download_dir=s"    => \$download_dir,
    "filename|f=s"      => \$filename,

    # Database params
    "host|h=s"          => \$host,
    "port|p=s"          => \$port,
    "database|db=s"     => \$db,
    "user|u=s"          => \$user,
    "password|pw=s"     => \$pass,

    # Or use config file
    "config=s"          => \$config,
);

$| = 1;

# Set default download directory
$download_dir //= ".";

my $logfile = catfile($download_dir, "$filename.log");
open (my $logh, ">", $logfile) or die "Error opening log file";

if ($config) {
    $P    = CoGe::Accessory::Web::get_defaults($config);
    $db   = $P->{DBNAME};
    $host = $P->{DBHOST};
    $port = $P->{DBPORT};
    $user = $P->{DBUSER};
    $pass = $P->{DBPASS};
}

# Verify parameters
$gid = unescape($gid) if $gid;
$filename = unescape($filename) if $filename;
my $file_temp = $filename . ".tmp";

if (not $gid) {
    say $logh "log: error: genome not specified use gid";
    exit(-1);
}

if (not $filename) {
    say $logh "log: error: output file not specified use output";
    exit(-1);
}

my $connstr = "dbi:mysql:dbname=$db;host=$host;port=$port;";
my $coge = CoGeX->connect( $connstr, $user, $pass );
#$coge->storage->debugobj(new DBIxProfiler());
#$coge->storage->debug(1);

unless ($coge) {
    say $logh "log: error: couldn't connect to database";
    exit(-1);
}

my $file = catfile($download_dir, $filename);

return if -r $file;

sub main {
    my $dsg = $coge->resultset('Genome')->find($gid);

    open(my $fh, ">", $file_temp);

    my %chr2ds;
    foreach my $ds ( $dsg->datasets ) {
        map { $chr2ds{$_} = $ds } $ds->chromosomes;
    }
    foreach my $chr ( sort { versioncmp( $a, $b ) } $dsg->chromosomes ) {
        print $fh ">Features $chr\n";
        foreach my $feat (
            sort {
                $a->start <=> $b->start
                || $a->feature_type_id <=> $b->feature_type_id
            } $chr2ds{$chr}
            ->features( { chromosome => $chr }, { "order_by" => "start ASC" } )
        )
        {
            next if $feat->type->name eq "chromosome";
            next if $feat->type->name =~ /utr/i;

            my $locs = get_locs($feat);

            my $item = shift @$locs;
            print $fh join( "\t", @$item, $feat->type->name ), "\n";
            foreach my $loc (@$locs) {
                print $fh join( "\t", @$loc ), "\n";
            }
            my $name_tag = get_name_tag($feat);
            my ( $pri_name, @names ) = $feat->names;
            print $fh "\t\t\t", $name_tag, "\t";
            if ( $name_tag =~ /_id/ ) {
                print $fh "gnl|dbname|";
            }
            print $fh $pri_name, "\n";

            foreach my $name (@names) {
                print $fh "\t\t\t", join( "\t", "alt_name", $name ), "\n";
            }
            foreach my $anno ( $feat->annotations ) {
                print $fh "\t\t\t", join( "\t", $anno->type->name, $anno->annotation ),
                "\n";
            }
        }
    }
    close($fh);
    exit 1 unless move($file_temp, $file);
}

sub get_name_tag {
    my $feat = shift;
    my $tag;
    if    ( $feat->type->name =~ /mRNA/ ) { $tag = "transcript_id" }
    elsif ( $feat->type->name =~ /CDS/ )  { $tag = "protein_id" }
    else                                  { $tag = "locus_tag" }
    return $tag;
}

sub get_locs {
    my $feat = shift;
    my @locs =
      map { [ $_->start, $_->stop ] }
      sort { $a->start <=> $b->start } $feat->locations;
    if ( $feat->strand =~ /-/ ) {
        @locs = reverse @locs;
        foreach my $item (@locs) {
            ( $item->[0], $item->[1] ) = ( $item->[1], $item->[0] );
        }
    }
    return \@locs;
}

main;
