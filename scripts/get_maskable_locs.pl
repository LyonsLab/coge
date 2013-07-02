# -*- perl -*-
use strict;

use CoGeX;
my $connstr = 'dbi:mysql:genomes:biocon:3306';
my $s = CoGeX->connect($connstr, 'bpederse', 'brent_cnr');

my $datasets = [@ARGV] ;
my %seen;
foreach my $feat ( $s->resultset('Feature')->search( {
                    'me.dataset_id'   => { 'IN' => $datasets }
                    ,'me.feature_type_id' => {'IN' => [ 3, 191 ] } # CDS, exon
                }, {
                    prefetch => [ 'locations']
                    ,order_by => ['me.chromosome', 'locations.start']
                } )) {
                my $chr = $feat->chromosome;
                next if $chr =~ /super/;
                foreach my $loc ( $feat->locations() ) {
                    my $locstr = $chr . "," . $loc->start . "," . $loc->stop;
                    if($seen{$locstr}){ next; }
                    $seen{$locstr} = 1;
                    print $locstr . "\n";
                }
}
