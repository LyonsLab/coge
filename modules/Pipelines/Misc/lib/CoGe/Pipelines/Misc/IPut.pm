package CoGe::Pipelines::Misc::IPut;

use Moose;

use File::Basename qw(basename);
use File::Spec::Functions;
use URI::Escape::JavaScript qw(escape);

use CoGe::Accessory::IRODS qw(irods_get_base_path irods_iput);
use CoGe::Accessory::Utils;

BEGIN {
    our ( @EXPORT, @ISA, $VERSION );
    require Exporter;

    $VERSION = 0.1;
    @ISA = qw ( Exporter );
    @EXPORT = qw ( export_to_irods );
}

sub export_to_irods {
    my ($src, $options, $user) = @_;

    my $base = $options->{dest_path};
    $base = irods_get_base_path($user->name) unless $base;
    my $output = catfile($base, basename($src));

    return $output, (
        cmd => irods_iput($src, $output, { no_execute => 1 }),
        description => "Exporting file to IRODS",
        args => [],
        inputs => [$src],
        outputs => []
    );
}

1;
