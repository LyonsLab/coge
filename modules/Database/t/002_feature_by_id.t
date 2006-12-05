# -*- perl -*-

# t/002_feature_by_id.t - check feature object creation works by using
# a feature_id query

use Test::More tests => 4;

BEGIN { use_ok( 'CoGeX' ); }

my $connstr = 'dbi:mysql:genomes:biocon:3306';
my $s = CoGeX->connect($connstr, 'cnssys', 'CnS' );

isa_ok ($s, 'CoGeX');

my $rs = $s->resultset('Feature')->search( feature_id => 42 );

my $f = $rs->next();

is( $f->feature_type_id(), 3 );
is( $f->dataset_id(),      2 );
