# -*- perl -*-

# t/009_feature_from_dataset_repeated - get the feature objects

use Test::More tests => 5;

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
                      prefetch => { dataset => 'organism' },
                    }
                );

my @features = $rs->all();
is( scalar(@features), 10863 );

$rs = $s->resultset('Feature')->search(
                    {
                      'organism.name' => { 'like' => 'Drosophila%' },
                      'me.feature_id' => 1367613
                    },
                    {
                      join => { dataset => 'organism' },
                      prefetch => qw/locations/
                    }
                );

my @features2 = $rs->all();
is( scalar(@features2), 1 );
