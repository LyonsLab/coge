#!/usr/bin/perl -w

use CoGeX;
use strict;
use Data::Dumper;
use Benchmark qw/timethese/;
$| = 1;

my $connstr = 'dbi:mysql:dbname=DB;host=HOST;port=PORT';
my $s = CoGeX->connect($connstr, 'USER', 'PASSWORD' );

$s->storage->debug(0);

timethese( 10, { 'search' => \&search, 'esearch' => \&esearch } );

sub esearch {

    my @results;
    my $rs = $s->resultset('Feature')->esearch(
        {
            'feature_type.name' => 'gene',
            'feature_names.name' =>
              [ '-and', { '-like' => 'At2g2%' }, { '-not_like' => ['%.%'] } ]
        }
    );

    #print "got resultset\n";
    while ( my $feat = $rs->next() ) {
        my $fn = $feat->feature_names;
        foreach my $name ( $fn->next() ) {
            push( @results, $name->name );
        }
        push( @results, $feat->dataset->organism->name );
        map { push( @results, $_->annotation ) } $feat->annotations;
    }
}

sub search {

    my @results;
    my $rs = $s->resultset('Feature')->search(
        {
            'feature_type.name' => 'gene',
            'feature_names.name' =>
              [ '-and', { '-like' => 'At2g2%' }, { '-not_like' => ['%.%'] } ]
        },
        {
            'join' =>
              [ 'feature_names', 'feature_type', { 'dataset' => 'organism' } ]
        }

    );

    #print "got resultset\n";
    while ( my $feat = $rs->next() ) {
        my $fn = $feat->feature_names;
        foreach my $name ( $fn->next() ) {
            push( @results, $name->name );
        }
        push( @results, $feat->dataset->organism->name );
        map { push( @results, $_->annotation ) } $feat->annotations;
    }
}
