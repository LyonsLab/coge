#!/usr/bin/env perl
use v5.14;
use strict;
use warnings;
use CoGe::Accessory::Web qw(get_defaults);
use Data::Dumper;
use File::Basename;
use File::Spec::Functions;

my %opts = @ARGV;
my $gid = $opts{'gid'};
my $chr = $opts{'chr'};
my $config = get_defaults(catfile(File::Basename::dirname(File::Basename::dirname($0)), 'coge.conf'));
my $path = catfile($config->{SECTEMPDIR}, "downloads/genome", $gid);
mkpath( $path, 0, 0777 ) unless -d $path;
my $file = catfile($path, $gid . "_" . $chr . ".faa");
unless (-e $file) {
    my $db = CoGeX->dbconnect($config);
    die "ERROR: couldn't connect to the database" unless $db;
    my $genome = $db->resultset('Genome')->find($gid);
    die "ERROR: unable to create genome object using id $gid" unless ($genome);

    # Get sequence from file
    my $seq = $genome->get_genomic_sequence(
        chr   => $chr,
        start => 1,
        stop  => $genome->get_chromosome_length($chr)
    );
    open(my $fh, '>', $file) or die "Could not open file '$file' $!";
    print $fh '>chromosome ' . $chr . "\n";
    for (my $i=0; $i<length($seq); $i+=70) {
        print $fh substr($seq, $i, 70);
        print $fh "\n";
    }
    close $fh;
}
