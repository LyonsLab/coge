#!/usr/bin/env perl
use v5.14;
use strict;
use warnings;
use CoGe::Accessory::Web qw(get_defaults);

our ($gid, $chr);
GetOptions(
    # Required workflow params
    "gid=s"        => \$gid,    # genome ID
    "chr=s"        => \$chr     # chromosome
);

my $config = get_defaults();
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
