package CoGe::Pipelines::Experiment;

use File::Basename qw(basename);
use File::Spec::Functions;
use URI::Escape::JavaScript qw(escape);

use CoGe::Accessory::Utils;

BEGIN {
    our ( @EXPORT, @ISA, $VERSION );
    require Exporter;

    $VERSION = 0.1;
    @ISA = qw ( Exporter );
    @EXPORT = qw ( export_experiment );
}

sub export_experiment {
    my ($params, $output, $conf) = @_;

    return (
        cmd => catdir($conf->{SCRIPTDIR}, "export_experiment.pl"),
        description => "Generating experiment files",
        args => [
            ["-eid", $params->{eid}, 0],
            ["-output", $output, 1],
            ["-conf", $conf->{_CONFIG_PATH}, 0],
            ["-dir", ".", ""]
        ],
        inputs => [],
        outputs => [$output]
    );
}

1;
