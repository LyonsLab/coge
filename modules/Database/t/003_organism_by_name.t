# -*- perl -*-

# t/003_organsim_by_name.t - get an organism object by name

use Test::More tests => 4;

BEGIN { use_ok( 'CoGeX' ); }

my $connstr = 'dbi:mysql:genomes:biocon:3306';
my $s = CoGeX->connect($connstr, 'cnssys', 'CnS' );

isa_ok ($s, 'CoGeX');

my $rs = $s->resultset('Organism')->search(
                                      {
                                        name => { 'like' => 'Nostoc%' }
                                      }
                                    );

my @orgs = $rs->all();

is( scalar(@orgs),           1  );
is( $orgs[0]->organism_id(), 24 );
