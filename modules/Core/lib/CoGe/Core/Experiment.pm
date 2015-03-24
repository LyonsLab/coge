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

sub experimentcmp($$) { # for sorting DBI-X objects or DBI hashes
    my ($a, $b) = @_;

    if ( ref($a) eq 'HASH' && ref($b) eq 'HASH' ) { # DBI
        versioncmp( $b->{version}, $a->{version} )
          || $a->{name} cmp $b->{name}
          || $b->{id} cmp $a->{id};
    }
    else { # DBI-X
        versioncmp( $b->version, $a->version )
          || $a->name cmp $b->name
          || $b->id cmp $a->id;
    }
}

1;
