#! /usr/bin/perl -w

use strict;
use CoGeX;

#time
#mysql before update: real    0m54.935s
#mysql after  update: real    0m44.592s
#postgresql         : real    0m30.214s

my $connstr = 'dbi:mysql:dbname=DB;host=HOST;port=PORT';
my $s = CoGeX->connect($connstr, 'USER', 'PASSWORD' );

my $rs = $s->resultset('Feature')->search(
    {
        'feature_type.name'  => 'CDS',
        'feature_names.name' => { like => 'At%g%', '-not_like' => "%.%" },
        'me.dataset_id'      => { 'IN' => [ 6, 7, 8, 9, 10 ] }
    },
    {
        join     => [ 'feature_names', 'feature_type' ],
        prefetch => [ 'feature_names', 'feature_type' ],
        order_by => ['feature_names.name']
    }

);

while ( my $feat = $rs->next() ) {
    my $fn   = $feat->feature_names;
    my $type = $feat->feature_type->name;
    print ">";
    map { print $_->name . ":" . $type . "\t" } $fn->next();
    print "\n";

    # this prefetch avoids n calls where n is number of sequences #
    foreach my $seq ( $feat->sequences ) {
        print $seq->sequence_data;
    }
    print "\n\n";
}
