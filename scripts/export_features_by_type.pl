#!/usr/bin/perl -w

use strict;
use CoGeX;
use CoGe::Accessory::Web;
use Data::Dumper;
use Getopt::Long;
use File::Copy;
use File::Spec::Functions;
use Text::Wrap;

$|++;

our ($P, $db, $host, $port, $user, $pass, $dsid, $gid, $ftid, $prot,
     $filename, $workdir, $config);

GetOptions(
    "dir=s"    => \$workdir,
    "gid=n"    => \$gid,
    "dsid"     => \$dsid,
    "ftid=n"   => \$ftid,
    "prot=n"   => \$prot,
    "file|f=s" => \$filename,

    # Database params
    "host|h=s"          => \$host,
    "port|p=s"          => \$port,
    "database|db=s"     => \$db,
    "user|u=s"          => \$user,
    "password|pw=s"     => \$pass,

    # Or use config file
    "config=s"          => \$config,
);

$workdir //= ".";

my $logfile = catfile($workdir, "export.log");

open (my $logh, ">", $logfile) or die "Error opening log file";

if ($config) {
    $P    = CoGe::Accessory::Web::get_defaults($config);
    $db   = $P->{DBNAME};
    $host = $P->{DBHOST};
    $port = $P->{DBPORT};
    $user = $P->{DBUSER};
    $pass = $P->{DBPASS};
}

my $connstr = "dbi:mysql:dbname=$db;host=$host;port=$port;";
my $coge = CoGeX->connect( $connstr, $user, $pass );

unless ($coge) {
    say $logh "Failed to connect to database";
    exit;
}

unless ( $dsid or $gid ) {
    say $logh "No genome or dataset id specified.";
    exit;
}
unless ($ftid) {
    say $logh "No feature type id specified.";
    exit;
}

my ( $dsg, $ds );
($dsg) = $coge->resultset('Genome')->search( { "me.genome_id" => $gid },
    { join => 'genomic_sequences', prefetch => 'genomic_sequences' } )
  if $gid;

$ds = $coge->resultset('Dataset')->find($dsid) if $dsid;

($dsg) = $ds->genomes if $ds;

my $ft = $coge->resultset('FeatureType')->find($ftid);
my $file = catfile($workdir, $filename);
my $file_temp = $file . ".tmp";

# Exit if the file already exists
exit if -r $file;

my $count = 1;
my @feats;
if ($ds) {
    @feats = $ds->features( { feature_type_id => $ftid } );
}
else {
    @feats = $dsg->features( { feature_type_id => $ftid } );
}

my $org = $dsg->organism->name;

open(my $fh, ">", $file_temp);

foreach my $feat (@feats) {
    my ($chr) = $feat->chromosome;    #=~/(\d+)/;
    my $name;

    eval {
        foreach my $n ( $feat->names ) {
            $name = $n;
            last unless $name =~ /\s/;
        }
        $name =~ s/\s+/_/g;
        my $title = join( "||",
            $org,              $chr,      $feat->start,
            $feat->stop,       $name,     $feat->strand,
            $feat->type->name, $feat->id, $count );
        my $seq;
        if ($prot) {
            my (@seqs) = $feat->protein_sequence( dsgid => $dsg->id );
            next unless scalar @seqs;
            next if scalar @seqs > 1;    #didn't find the correct reading frame;
            next unless $seqs[0] =~ /[^x]/i;
            $seq = $seqs[0];
        }
        else {
            $seq = $feat->genomic_sequence;
        }
        $title = ">" . $title . "\n";

        say $fh $title, $seq, "\n";

        $count++;
    };

    if ($@) {
        say $logh "DOWNLOAD ABORTED: $filename";
        last;
    }
}

close($fh);
close($logh);

# ensure the file has been copied successfully
exit 1 unless move($file_temp, $file);
