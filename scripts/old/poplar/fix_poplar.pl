#! /usr/bin/perl -w

use strict;
use CoGeX;
use Data::Dumper;

my @chrs =
  qw/I II III IV V VI VII VIII IX X XI XII XIII XIV XV XVI XVII XVIII XIX/;
my %cmap;
my $i = 1;
map { $cmap{ $chrs[ $_ - 1 ] } = ( $_ < 10 ) ? '0' . $_ : $_ } 1 .. 19;

my $connstr = 'dbi:mysql:DB:HOST:PORT';
my $s = CoGeX->connect( $connstr, 'USER', 'PASSWORD' );

$s->storage->debug(0);

=begin
my $rs = $s->resultset('Feature')->esearch(
        {
            'me.dataset_id' => 505,
            'feature_type.name' => 'CDS'
        },
        {
            join => 'feature_names'
        });

while( my $feat = $rs->next()){
    my $name = ($feat->names())[0];
    print $name . "\n";
    #my ($newname1) = $name =~ /(proteinId\s\d+)/;
    my ($newname) = $name =~ /name\s([^;]+)/;
    my $feat_id = $feat->feature_id;
    print "$newname\t$feat_id\n\n";

    if( !$newname ){ next; }

    my $name_obj = $s->resultset('FeatureName')->new( {
            name => $newname,
            description=>'poplar_group',
            feature_id=>$feat->feature_id
        });
    $name_obj->insert();
}
=cut

my $rs = $s->resultset('Feature')->search(
    {
        'me.dataset_id'             => 505,
        'feature_names.description' => 'poplar_group'
    },
    { prefetch => 'feature_names' }
);

my %feats;
while ( my $feat = $rs->next() ) {
    my $chr = $feat->chromosome;
    next if $chr =~ /scaffold/;
    $chr = $cmap{$chr};

    my @names = $feat->feature_names;
    my ($d) = map { $_->name }
      grep { $_->description eq 'poplar_group' } @names;

    #print $d . "\n";
    if ( !$feats{$chr}{$d} ) { $feats{$chr}{$d} = $feat->start; }
    else {
        $feats{$chr}{$d} =
          $feats{$chr}{$d} < $feat->start ? $feats{$chr}{$d} : $feat->start;
    }
}

foreach my $chr ( keys %feats ) {
    my %chrhash = %{ $feats{$chr} };

    sub valsort {
        return $chrhash{$a} <=> $chrhash{$b};
    }
    my $gene_num = 1;
    foreach my $key ( sort valsort ( keys %chrhash ) ) {
        my $gn = "0" x ( 5 - length($gene_num) ) . $gene_num;

        my $name = 'P' . $chr . 'G' . $gn;
        print $name . "\t" . $key . "\t" . $chrhash{$key} . "\n";
        ++$gene_num;
        my $fids =
          $s->resultset('Feature')->search( { 'feature_names.name' => $key },
            { 'join' => 'feature_names' } );
        while ( my $fid = $fids->next() ) {
            my $feat_id  = $fid->feature_id;
            my $name_obj = $s->resultset('FeatureName')->new(
                {
                    name        => $name,
                    description => 'poplar_gene_name',
                    feature_id  => $feat_id
                }
            );
            print "#insert commented out as this has already been run\n";

            # $name_obj->insert();
        }
    }

    #exit();
}
