# -*- perl -*-

use CoGeX;
my $connstr = 'dbi:mysql:dbname=DB;host=HOST;port=PORT';
my $s = CoGeX->connect($connstr, 'USER', 'PASSWORD' );
use strict;
use DBIxProfiler;

my $profiler;

#my $profiler = new DBIxProfiler();
#$s->storage->debugobj($profiler);
#$s->storage->debug(1);

# Thu Oct 11 15:52:34 2007
# number of database queries: 3001
# total time for database queries: 30.8494563102722

my $rs = $s->resultset('Feature')->search(
    {
        'feature_names.name' => { 'like' => 'At3g3%' },
        'feature_type.name'  => 'CDS'
    },
    {
        join     => [ 'feature_names',      'feature_type' ],
        order_by => [ 'feature_names.name', 'me.start', 'feature_type.name' ]
    }
);

open( OUT, ">", "scripts/seq_new.txt" );
while ( my $feat = $rs->next() ) {
    my $s = $feat->genomic_sequence();
    print OUT $s . "\n";
}
if ($profiler) {
    if ( length(`diff scripts/seq.txt scripts/seq_new.txt`) ) {
        print "*" x 80 . "\n";
        print "OH NOOOOOOO... not same as test output\n";
        print "*" x 80 . "\n";
    }

    my ( $qcount, $qtime ) = $profiler->get_query_count_and_time();
    print scalar localtime() . "\n";
    print "number of database queries: " . $qcount . "\n";
    print "total time for database queries: " . $qtime . "\n";
}
