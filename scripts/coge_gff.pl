#!/usr/bin/perl -w

use strict;
use CoGeX;
use Getopt::Long;
use File::Path;
use File::Spec;
use URI::Escape::JavaScript qw(unescape);

our ($DEBUG, $db, $user, $pass, $gid, $config, $host, $port, $P,
     $filename, $download_dir, $dsid, $annos, $cds, $name_unique,
     $id_type, $upa);

GetOptions(
    "debug=s"                         => \$DEBUG,
    "gid=i"                           => \$gid,
    "dsid=i"                          => \$dsid,
    "download_dir=s"                  => \$download_dir,
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

mkpath($download_dir, 0, 0777) unless -r $download_dir;

my $logfile = File::Spec->catdir($download_dir, "$filename.log");
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

my $file = File::Spec->catdir($download_dir, $filename);


unless ( -r $file and -r "$file.finished") {
    open(my $fh, ">", $file) or die "Error creating gff file";

    my $item;
    my $id;
    my $org;

    if ($gid) {
        $item = $coge->resultset('Genome')->find($gid);
        $id  = $gid;
        $org = $item->organism->name . "_dsgid";
    }
    elsif ($dsid) {
        $item = $coge->resultset('Dataset')->find($dsid);
        $id   = $dsid;
        $org  = $item->organism->name . "_dsid";
    }

    print $fh $item->gff(
        print                     => 0,
        annos                     => $annos,
        cds                       => $cds,
        name_unique               => $name_unique,
        id_type                   => $id_type,
        unique_parent_annotations => $upa
    );

    close($fh);
    system("touch $file.finished");
}
