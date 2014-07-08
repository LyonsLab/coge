# -*- perl -*-

# t/005_feature_by_organsim - get the feature objects for an organism,
# given by name, "Nostoc%"

use Test::More tests => 4;

BEGIN {
  use_ok( 'CoGeX' );
  use_ok( 'DBIxProfiler' );
}

my $connstr = 'dbi:mysql:genomes:biocon:3306';
my $s = CoGeX->connect($connstr, 'cnssys', 'CnS' );

$s->storage->debugobj(new DBIxProfiler());
$s->storage->debug(1);

isa_ok ($s, 'CoGeX');

my $rs = $s->resultset('Feature')->search(
                    {
                      'organism.name' => { 'like' => 'Nostoc%' }
                    },
                    {
                      join => { dataset => 'organism' },
                      prefetch => qw/dataset/,
                    }
                );

my @features = $rs->all();
is( scalar(@features), 10863 );
