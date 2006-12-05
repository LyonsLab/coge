# -*- perl -*-

# t/001_load.t - check module loading and create testing directory,
# including connecting to the DB on biocon

use Test::More tests => 2;

BEGIN {
  use_ok( 'CoGeX' );
}

my $connstr = 'dbi:mysql:genomes:biocon:3306';
my $s = CoGeX->connect($connstr, 'cnssys', 'CnS' );

isa_ok ($s, 'CoGeX');
