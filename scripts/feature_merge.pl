use strict;
use Data::Dumper;
use CoGeX;

# THOUGH this is called feature_merge, it just gets all features, and
# all chromosomal sequences for a given organism.

my $connstr = 'dbi:mysql:dbname=DB;host=HOST;port=PORT';
my $s = CoGeX->connect($connstr, 'USER', 'PASSWORD' );

my $organism = shift;
my $datasets = \@ARGV;

get_locs($datasets);
get_sequence($datasets);

sub get_locs {
    my $SEP      = "\t";
    my $datasets = shift;
    my $rs       = $s->resultset('Feature')->search(
        {
            'me.dataset_id'     => { 'IN'       => $datasets },
            'feature_type.name' => { 'NOT LIKE' => '%contig%' }
        },
        {
            'prefetch' => [ 'feature_type',  'locations' ],
            'order_by' => [ 'me.chromosome', 'me.start' ]
        }
    );

    while ( my $g = $rs->next() ) {
        if ( $g->feature_type->name eq 'CNS' ) { next; }
        my $f = {
            'locs' => [
                map {
                    {
                        'start'  => $_->start,
                        'stop'   => $_->stop,
                        'strand' => $_->strand,
                        'chr'    => $_->chromosome
                    }
                  } $g->locations( {}, { 'order' => 'start' } )
            ],
            'type'  => uc( $g->feature_type->name ),
            'names' => [ map { uc( $_->name ) } $g->feature_names() ]
        };
        my $locstr = "";
        map { $locstr .= $_->{start} . "|" . $_->{stop} . "|" }
          @{ $f->{'locs'} };
        chop $locstr;

        print join( "|", @{ $f->{names} } )
          . $SEP
          . $f->{'type'}
          . $SEP
          . $f->{'locs'}[0]{'chr'}
          . $SEP
          . $f->{'locs'}[0]{'strand'}
          . $SEP
          . $locstr . "\n";
    }
}

sub get_sequence {
    open( FA, ">", $organism . ".fasta" );
    foreach my $ds (@$datasets) {
        my $ds = $s->resultset('Dataset')->resolve($ds);
        my %seen;
        foreach my $chr ( $ds->get_chromosomes ) {

            # TODO: this will break with contigs.
            #next if $chr =~ /^contig/;
            next if $chr =~ /random/;

       #$chr =~ s/scaffold/super/g; print STDERR "CHANGING scaffold => super\n";

            #next if $chr =~ /scaffold/;

            print STDERR $chr . "\n" unless $seen{$chr};
            $seen{$chr} = 1;

            #  rice/chr01.fasta
            print FA "> $chr\n";
            print FA $ds->get_genomic_sequence( chromosome => $chr ) . "\n";
        }
    }
    close FA;
}
