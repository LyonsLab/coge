#!/usr/bin/perl

use strict;
use Data::Dumper;
use CoGeX;
use File::Path;
use File::Spec::Functions;
use POSIX;

my $BASEDIR = "/opt/apache/CoGe/data/genomic_sequence";
my $connstr = 'dbi:mysql:dbname=DB;host=HOST;port=PORT';
my $coge = CoGeX->connect($connstr, 'USER', 'PASSWORD' );

foreach my $item ( $coge->resultset('Dataset')->all() ) {
    next unless $item->sequence_type;
    my $dsid      = "dataset_" . $item->id;
    my $type_id   = "seqtype_" . $item->sequence_type->id;
    my $chr_total = $item->chromosomes;
    $chr_total = scalar @$chr_total;
    my $orgid = "org_" . $item->organism->id;
    my $dir   = catfile( $BASEDIR, ceil( $item->organism->id / 1000 ),
        $orgid, $type_id, $dsid );
    my $chr_level = 1;
    mkpath( $dir . "/" . $chr_level );
    my $chr_count = 0;

    foreach my $chr ( sort $item->chromosomes ) {
        my $seq = $item->genomic_sequence(
            chr               => $chr,
            start             => 1,
            stop              => $item->last_chromosome_position($chr),
            skip_length_check => 1
        );
        open( OUT, ">$dir/$chr_level/$chr" );
        print OUT $seq;
        close OUT;
        $chr_count++;
        unless ( $chr_count % 1000 ) {
            $chr_level++;
            mkpath( $dir . "/" . $chr_level );
        }
    }

#    print $dir."/".$chr_level,"\t",$chr_count, "\t", $chr_total,"\n" if $chr_total > 1000;
}
