#!/usr/bin/perl -w

use CoGeX;
use strict;
use Data::Dumper;
$| = 1;

my $connstr = 'dbi:mysql:dbname=DB;host=HOST;port=PORT';
my $s= CoGeX->connect($connstr, 'USER', 'PASSWORD' );

$s->storage->debug(1);

my $rs = $s->resultset('Feature')->search(
    {
        'feature_type.name' => 'gene',
        'dataset.version'   => 6,
        'feature_names.name' =>
          [ '-and', { '-ilike' => 'At%' }, { '-not_like' => ['%.%'] } ]
    },
    {
        join     => [ 'feature_names', 'feature_type', 'dataset' ],
        prefetch => [ 'feature_type',  'dataset' ],
    },
);

#print "got resultset\n";

my %seen;
$s->storage->debug(0);
my %annotations;
while ( my $feat = $rs->next() ) {
    my $fn   = $feat->feature_names;
    my $type = $feat->feature_type->name;

    my $annos = $feat->annotations;
    my @names = grep { $_->name =~ /^At\dg\d{5}$/i } $fn->next();
    print join( "|", @names ), "\n";
    next unless @names;

    #print join("|", sort map { $_->name } @names) . "\n";
    my @gos = ();
    while ( my $a = $annos->next() ) {
        next unless $a->annotation_type;
        next unless $a->annotation_type->annotation_type_group;
        next
          unless $a->annotation_type->annotation_type_group->name =~
              /^[Gg][Oo]/;

        #print $a->annotation_type->annotation_type_group->name  . "\n";
        push(
            @{ $annotations{ $names[0]->name } },
            'G' . $a->annotation_type->name
        );
    }

    #map { print $_->annotation .  "\n" }  @gos;
}
foreach my $accn ( keys %annotations ) {
    print $accn . "\t" . join( ";", @{ $annotations{$accn} } ) . "\n";
}

#TODO: fix me
