package CoGe::Core::Experiment;
use strict;

use Sort::Versions;

our ( @EXPORT, @EXPORT_OK, @ISA, $VERSION );

BEGIN {
    require Exporter;
    $VERSION = 0.1;
    @ISA = qw ( Exporter );
    @EXPORT = qw ( );
    @EXPORT_OK = qw(experimentcmp);
}

sub experimentcmp($$) {
    my ($a, $b) = @_;

    versioncmp( $b->version, $a->version )
      || $a->name cmp $b->name
      || $b->id cmp $a->id;
}

1;
