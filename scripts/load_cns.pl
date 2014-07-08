#! /usr/bin/perl -w

use strict;
use Data::Dumper;
use CoGeX;
use Getopt::Long;

my $connstr = 'dbi:mysql:dbname=DB;host=HOST;port=PORT';
my $coge = CoGeX->connect($connstr, 'USER', 'PASSWORD' );

my $GO = 0;
my $ds1;
my $ds2;
GetOptions(
    "ds1=i" => \$ds1,
    "ds2=i" => \$ds2,
    "go=s"  => \$GO,
);

# perl scripts/load_cns.pl -ds1 34465 -ds2 34978 < /tmp/rice_sorghum_for_paper.csv

# to clear out db:
#DELETE from annotation where annotation.feature_id IN (SELECT f.feature_id  FROM feature f WHERE f.feature_type_id =202 AND f.dataset_id IN (34465, 34978));
#DELETE from location where location.feature_id IN (SELECT f.feature_id  FROM feature f WHERE f.feature_type_id =202 AND f.dataset_id IN (34465, 34978));

#DELETE from feature_name where feature_name.feature_id IN (SELECT f.feature_id  FROM feature f WHERE f.feature_type_id =202 AND f.dataset_id IN (34465, 34978));
#DELETE FROM feature WHERE feature.feature_id IN (
# SELECT fid from (
#    SELECT f.feature_id as fid
#    FROM feature f
#    WHERE f.feature_type_id =202
#    AND f.dataset_id
#    IN ( 34465, 34978 )
# ) as foo
#);

exit() unless $ds1 && $ds2;

my $dsq = $coge->resultset('Dataset')->find($ds1);
my $dss = $coge->resultset('Dataset')->find($ds2);
my ($anno_type) =
  $coge->resultset('AnnotationType')->search( { name => "note" } );
my $cns_id  = 202;     #feature type id
my $anno_id = 16156;

print STDERR $dsq->dataset_id . "\n";
print STDERR $dss->dataset_id . "\n";

warn "-go flag is not true, nothing will be added to the database.\n"
  unless $GO;

my %data;
my %annos;
my %seen;

# line like:
#qname   qchr    qstart  qstop   qstrand sname   schr    sstart  sstop   sstrand
my $outfile = "cns_fids_$ds1" . "_" . "$ds2.txt";
print "printing ids to: " . $outfile . "\n";
open( SAVEFIDS, ">", $outfile );

while ( my $line = <> ) {
    next if $line =~ /^qname/ || $line =~ /^#/;
    chomp $line;
    $line =~ s/,\s/,/g;

    #my ($qname, $sname, $qchr, $schr, @line)= split(/[,|\|]/, $line);
    my (
        $qname, $qchr, $qstart, $qstop, $qstrand,
        $sname, $schr, $sstart, $sstop, $sstrand
    ) = split( /\t/, $line );
    $qchr =~ s/^0+//;
    $schr =~ s/^0+//;

    if (
        $seen{
                $qchr . ", "
              . $qstart . ","
              . $qstop . ","
              . $schr . ","
              . $sstart . ","
              . $sstop
        }
      )
    {
        next;
    }
    $seen{ $qstart . "," . $qstop . "," . $sstart . "," . $sstop } = 1;

    my $qfeat = $dsq->add_to_features(
        {
            feature_type_id => $cns_id,
            start           => $qstart,
            stop            => $qstop,
            chromosome      => $qchr,
            strand          => $qstrand,
        }
    ) if $GO;
    my $sfeat = $dss->add_to_features(
        {
            feature_type_id => $cns_id,
            start           => $sstart,
            stop            => $sstop,
            chromosome      => $schr,
            strand          => $sstrand,
        }
    ) if $GO;

    my $qloc = $qfeat->add_to_locations(
        {
            start      => $qstart,
            stop       => $qstop,
            strand     => $qstrand,
            chromosome => $qchr
        }
    ) if $GO;
    my $sloc = $sfeat->add_to_locations(
        {
            start      => $sstart,
            stop       => $sstop,
            strand     => $sstrand,
            chromosome => $schr
        }
    ) if $GO;

    my $qfeat_name =
      $qfeat->add_to_feature_names( { name => "CNS $qname with $sname", } )
      if $GO;
    my $sfeat_name =
      $sfeat->add_to_feature_names( { name => "CNS $sname with $qname", } )
      if $GO;

    my $qannos = $qfeat->add_to_annotations(
        {
            annotation         => "CNS for $qname, $sname",
            annotation_type_id => $anno_id,
        }
    ) if $GO;
    my $sannos = $sfeat->add_to_annotations(
        {
            annotation         => "CNS for $qname, $sname",
            annotation_type_id => $anno_id,
        }
    ) if $GO;

    print SAVEFIDS $qfeat->feature_id . "\n" if $qfeat;
    print SAVEFIDS $sfeat->feature_id . "\n" if $sfeat;

}
