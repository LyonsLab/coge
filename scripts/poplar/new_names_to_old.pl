#! /usr/bin/perl -w

use strict;
use CoGeX;
use Data::Dumper;

my $connstr = 'dbi:mysql:DB:HOST:PORT';
my $s = CoGeX->connect( $connstr, 'USER', 'PASSWORD' );

$s->storage->debug(0);

my $rs =
  $s->resultset('FeatureName')
  ->search( { 'description' => 'poplar_gene_name', }, {} );

my %seen;
while ( my $name = $rs->next() ) {
    next if $seen{ $name->name };
    $seen{ $name->name } = 1;
    my $ns = $s->resultset('FeatureName')->search(
        {
            'feature_id'  => $name->feature_id,
            'description' => 'poplar_group'
        },
        {}
    );
    my $oldname = $ns->next();
    print $name->name . "," . $oldname->name . "\n";

}
