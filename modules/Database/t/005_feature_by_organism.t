# -*- perl -*-

# t/005_feature_by_organsim - get the feature objects for an organism,
# given by name, "Nostoc%"

use Test::More tests => 2;

BEGIN { use_ok( 'CoGeX' ); }

my $connstr = 'dbi:mysql:genomes:biocon:3306';
my $s = CoGeX->connect($connstr, 'cnssys', 'CnS' );

$s->storage->debug(1);

isa_ok ($s, 'CoGeX');

my $rs = $s->resultset('Feature')->search(
                                      { 
                                        'organism.name' => { 'like' => 'Nostoc%' }
                                      },
                                      {
                                        prefetch => { dataset => 'organism' },
                                      }
                                    );

my @features = $rs->all();
is( scalar(@features), 10863 );
