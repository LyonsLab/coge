package CoGe::Core::Genome;

use strict;
use warnings;

use File::Spec;

use CoGe::Accessory::Storage qw(get_genome_path);
use CoGe::Accessory::TDS qw(write read);

BEGIN {
    our ( @EXPORT, @ISA, $VERSION );
    require Exporter;

    $VERSION = 0.1;
    @ISA = qw ( Exporter );
    @EXPORT = qw ( );
}


1;

