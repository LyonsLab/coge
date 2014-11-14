package CoGe::Pipelines::Misc::IPut;

use Moose;

use File::Basename qw(basename);
use File::Spec::Functions;
use URI::Escape::JavaScript qw(escape);

use CoGe::Accessory::IRODS qw(irods_iput);
use CoGe::Accessory::Utils;

BEGIN {
    our ( @EXPORT, @ISA, $VERSION );
    require Exporter;

    $VERSION = 0.1;
    @ISA = qw ( Exporter );
    @EXPORT = qw ( export_to_irods );
}

sub export_to_irods {
    my ($src, $dest, $overwrite, $done_file) = @_;

    $overwrite = 0 unless defined $overwrite;

    my $cmd = irods_iput($src, $dest, { no_execute => 1, overwrite => $overwrite });

    my $filename = basename($done_file);

   return (
        cmd => qq[$cmd && touch $filename],
        description => "Exporting file to IRODS",
        args => [],
        inputs => [$src],
        outputs => [$done_file]
    );
}

1;
