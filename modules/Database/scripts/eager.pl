#!/usr/bin/perl -w

use CoGeX;
use strict;
use Data::Dumper;
$| = 1;


my $connstr = 'dbi:mysql:genomes:biocon:3306';
my $s = CoGeX->connect($connstr, 'cnssys', 'CnS' );

$s->storage->debug(1);


my $rs = $s->resultset('Feature')->esearch( { 
                    'feature_type.name' =>   'gene' ,
                    'feature_names.name' => ['-and' , 
                            { '-like' => 'At2g26%'} ,
                            { '-not_like' => ['%.%'] }
                        ]
                    } 
                    );

#print "got resultset\n";
while (my $feat = $rs->next()){
    my $fn = $feat->feature_names;
    foreach my $name ($fn->next()){
        print $name->name . "\t";
    }
    print "\nprefetched: " . $feat->dataset->organism->name . "\n";
    print "\nnot prefetched: ";
    map { print $_->annotation . "\t" } $feat->annotations;
}

