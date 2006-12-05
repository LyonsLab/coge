# -*- perl -*-

# t/029_feature_by_organsim - get the feature objects for an organism,
# given by name, "Nostoc%"

use Test::More tests => 4;

BEGIN { use_ok( 'CoGeX' ); }

my $connstr = 'dbi:mysql:genomes:biocon:3306';
my $s = CoGeX->connect($connstr, 'cnssys', 'CnS' );

isa_ok ($s, 'CoGeX');

my $rs = $s->resultset('Feature')->search(
                                      { 
                                        organism.name => { 'like' => 'Nostoc%' }
                                      },
                                      {
                                        join => [ qw/dataset organism/ ],
                                      }
                                    );

my $f = $rs->next();
print STDERR "\n", $f->feature_id(), "\n";
