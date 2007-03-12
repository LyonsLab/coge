# -*- perl -*-

# t/004_features_from_organsim - get the feature objects for an
# organism using the built in features() method

use Data::Dumper;
use Test::More tests => 5;

BEGIN { use_ok( 'CoGeX' ); }

my $connstr = 'dbi:mysql:genomes:biocon:3306';
my $s = CoGeX->connect($connstr, 'cnssys', 'CnS' );

#$s->storage->debug(1);

isa_ok ($s, 'CoGeX');

my $rs = $s->resultset('Dataset')->search( { name => 'CHR1.R5v01212004.xml' } );

my $d = $rs->next();
#diag( "\n\tDescription: " . $d->description() );

my @features = $d->features( { feature_id => 2377 } );

#diag( "Feature array size is " . scalar(@features) );
is( scalar(@features), 1);

my $f = $features[0];
my @fnames = $f->feature_names();
is( scalar(@fnames), 5);

#diag( $fnames[0]->name() );
is( $fnames[0]->name(), "At1g55470" );
