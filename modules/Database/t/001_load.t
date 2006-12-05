# -*- perl -*-

# t/001_load.t - check module loading and create testing directory,
# including connecting to the DB on biocon

use Test::More tests => 27;

BEGIN {
  use_ok( 'CoGeX' );
  use_ok( 'CoGeX::Annotation' );
  use_ok( 'CoGeX::FeatureType' );
  use_ok( 'CoGeX::AnnotationType' );
  use_ok( 'CoGeX::GenomicSequence' );
  use_ok( 'CoGeX::AnnotationTypeGroup' );
  use_ok( 'CoGeX::Image' );
  use_ok( 'CoGeX::DataInformation' );
  use_ok( 'CoGeX::Location' );
  use_ok( 'CoGeX::DataSource' );
  use_ok( 'CoGeX::Organism' );
  use_ok( 'CoGeX::Dataset' );
  use_ok( 'CoGeX::Permission' );
  use_ok( 'CoGeX::Sequence' );
  use_ok( 'CoGeX::SequenceType' );
  use_ok( 'CoGeX::Feature' );
  use_ok( 'CoGeX::User' );
  use_ok( 'CoGeX::FeatureList' );
  use_ok( 'CoGeX::UserGroup' );
  use_ok( 'CoGeX::FeatureListConnector' );
  use_ok( 'CoGeX::UserGroupConnector' );
  use_ok( 'CoGeX::FeatureListGroup' );
  use_ok( 'CoGeX::UserGroupFeatureListPermissionConnector' );
  use_ok( 'CoGeX::FeatureListGroupImageConnector' );
  use_ok( 'CoGeX::UserSession' );
  use_ok( 'CoGeX::FeatureName' );
}

my $connstr = 'dbi:mysql:genomes:biocon:3306';
my $s = CoGeX->connect($connstr, 'cnssys', 'CnS' );

isa_ok ($s, 'CoGeX');
