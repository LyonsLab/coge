use strict;
use Data::Dumper;
use CoGeX;

my $connstr = 'dbi:mysql:genomes:biocon:3306';
my $s = CoGeX->connect($connstr, 'bpederse', 'brent_cnr');

my $datasets = \@ARGV;


get_locs($datasets);

sub get_locs {
    my $datasets = shift;
    my $rs = $s->resultset('Feature')->search( {
              'me.dataset_id' => { 'IN' => $datasets }
            , 'feature_type.name'  => { 'NOT LIKE' => '%contig%' }
            } , { 
               'prefetch'           => [ 'feature_type', 'locations'] 
              ,'order_by'           => [ 'me.chromosome', 'me.start']
            });

    
    while(my $g = $rs->next()){
        my $f = {'locs' => [map { {'start' => $_->start, 'stop' => $_->stop, 'strand' => $_->strand, 'chr' => $_->chromosome }
                                            } $g->locations({} , {'order' => 'start'})]
                            , 'type' => $g->feature_type->name, 'names' => [sort map { $_->name } $g->feature_names() ]};
        my $locstr = "";
        map { $locstr .= $_->{start} . "|" . $_->{stop} . "|" } @{$f->{'locs'}};
        chop $locstr;
        
        print join("|", @{$f->{names}}) . "," . $f->{'type'} . "," . $f->{'locs'}[0]{'chr'} . "," . $f->{'locs'}[0]{'strand'} . "," . $locstr . "\n";
    }
}


