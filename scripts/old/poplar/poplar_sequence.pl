#! /usr/bin/perl -w

use strict;
use CoGeX;
use Data::Dumper;

my @chrs =
  qw/I II III IV V VI VII VIII IX X XI XII XIII XIV XV XVI XVII XVIII XIX/;
my %cmap;
my $i = 1;
map { $cmap{ $chrs[ $_ - 1 ] } = ( $_ < 10 ) ? '0' . $_ : $_ } 1 .. 19;

my %FH;

foreach my $chr ( values %cmap ) {
    open( $FH{$chr}, ">", "/tmp/poplar_/poplarchr" . $chr . ".fasta" );
}

# 9min 9 sec
my $connstr = 'dbi:mysql:DB:HOST:PORT';
my $s = CoGeX->connect( $connstr, 'USER', 'PASSWORD' );

$s->storage->debug(0);

#my $org = $s->resultset('Dataset')->search(
#    { 'organism.name' => {like => "%Poplar%" },
#      'version'       => '1.1'
#    },
#    { prefetch => 'organism' }
#);
#
#my $did = $org->next()->dataset_id;

my $rs = $s->resultset('Feature')->search(
    {
        'me.dataset_id'             => 505,
        'feature_names.description' => 'poplar_gene_name',
    },
    {
        prefetch => [ 'feature_names', 'dataset' ],
        order_by => 'me.feature_id'
    }
);

my $currentpid  = '';
my $sequence    = "";
my $currentgene = 0;
my $pid         = 1;

my %seen;
while ( my $feat = $rs->next() ) {
    my $chrnum = $cmap{ $feat->chromosome };
    my ($name) = map { $_->name }
      grep { $_->description eq 'poplar_gene_name' } $feat->feature_names();
    next if $seen{$name};
    $seen{$name} = 1;

    #print $name . "\t" . $chrnum . "\n";
    my $featset = $s->resultset('Feature')->esearch(
        {
            'me.dataset_id'      => 505,
            'feature_names.name' => $name
        },
        { join => 'feature_names' }
    );
    my @thisfeats;
    while ( my $f = $featset->next() ) {
        push( @thisfeats, $f );
    }
    @thisfeats = sort { $a->start <=> $b->start } @thisfeats;
    my $sequence = '';
    map { $sequence .= $_->genomic_sequence() } @thisfeats;
    my $fh = $FH{$chrnum};

    print $fh ">", $name, "\n";

    #print $feat->feature_type->name . "\n";
    print $fh $sequence, "\n\n";

}
