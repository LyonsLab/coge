package CoGe::Pipelines::Misc::IPut;

use Moose;

use File::Spec::Functions;
use URI::Escape::JavaScript qw(escape);

use CoGe::Accessory::IRODS;
use CoGe::Accessory::Utils;
use CoGe::Core::Genome qw(get_download_path);

BEGIN {
    our ( @EXPORT, @ISA, $VERSION );

    $VERSION = 0.01;
    require Exporter;

    $VERSION = 0.1;
    @ISA = qw ( Exporter );
    @EXPORT = qw ( export_to_irods );
}

sub export_to_irods {
    my ($src, $options, $user) = @_;

    my $env_file = CoGe::Accessory::IRODS::_irods_get_env_file();

    my $dest = $options->{dest_path};
    $dest = get_default_path($user->name) unless $dest;

    return (
        cmd => "export irodsEnvFile='$env_file'; iput",
        description => "Exporting file to IRODS",
        args => [
            ["-fT", '', 0],
            [$src, '', 0],
            [$dest, '', 0],
        ],
        inputs => [$src],
        outputs => []
    );
}

sub get_default_path {
    return "/iplant/home/" . shift . "/coge_data/";
}

1;
