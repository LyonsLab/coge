#!/usr/bin/perl -w

use strict;
use CoGeX;
use Getopt::Long;
use File::Copy;
use File::Path;
use File::Spec::Functions;
use URI::Escape::JavaScript qw(unescape);

our ($DEBUG, $db, $user, $pass, $id, $config, $host, $port, $P,
     $filename, $annos, $cds, $name_unique, $staging_dir,
     $id_type, $upa, $coge);

GetOptions(
    "debug=s"                         => \$DEBUG,
    "id=i"                            => \$id,
    "staging_dir=s"                  => \$staging_dir,
    "filename|f=s"                    => \$filename,
    "annos=i"                         => \$annos,
    "cds=i"                           => \$cds,
    "name_unique|nu=i"                => \$name_unique,
    "id_type|type=s"                  => \$id_type,
    "unique_parent_annotations|upa=i" => \$upa,

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
#open (my $logh, ">", $logfile) or die "Error opening log file";
$staging_dir //= ".";
$filename = unescape($filename) if $filename;
$id = unescape($id) if $id;

if (not $filename) {
    say STDERR "log: error: output file not specified use output";
    exit(-1);
}

my $file = catfile($staging_dir, $filename);
my $file_temp = $file . ".tmp";

# Check if file already exists
return if -r $file;

# Verify parameters
if (not $id) {
    say STDERR "log: error: genome not specified use id";
    exit(-1);
}

if ($config) {
    $P    = CoGe::Accessory::Web::get_defaults($config);
    $db   = $P->{DBNAME};
    $host = $P->{DBHOST};
    $port = $P->{DBPORT};
    $user = $P->{DBUSER};
    $pass = $P->{DBPASS};

    my $connstr = "dbi:mysql:dbname=$db;host=$host;port=$port;";
    $coge = CoGeX->connect( $connstr, $user, $pass );
    #$coge->storage->debugobj(new DBIxProfiler());
    #$coge->storage->debug(1);
}

unless ($coge) {
    say STDERR "log: error: couldn't connect to database";
    exit(-1);
}

my ($item, $org);
$item = $coge->resultset('Genome')->find($id);
$org = $item->organism->name . "id";

open(my $fh, ">", $file_temp) or die "Error creating gff file";

print $fh $item->gff(
    print                     => 0,
    annos                     => $annos,
    cds                       => $cds,
    name_unique               => $name_unique,
    id_type                   => $id_type,
    unique_parent_annotations => $upa,
    base_url => $P->{SERVER},
    debug => $DEBUG,
);

close($fh);
exit 1 unless move($file_temp, $file);
