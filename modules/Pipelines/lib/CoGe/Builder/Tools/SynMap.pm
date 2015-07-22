package CoGe::Builder::Tools::SynMap;

use v5.14;
use strict;
use warnings;

use File::Spec::Functions;

our (@EXPORT_OK, @ISA, $VERSION);

BEGIN {
    require Exporter;

    $VERSION = 0.1;
    @ISA = qw(Exporter);
    @EXPORT_OK = qw(generate_pseudo_assembly);
}

sub generate_pseudo_assembly {
    my ($JEX, $config, $input, $output, $flip) = @_;
    $flip = 0 unless $flip;

    my $cmd = "synmap/order_contigs_to_chromosome.pl";

    my $workflow = $JEX->create_workflow(
        name => "Generate Pseudo Assembly"
    );

    $workflow->add_job({
        cmd  => catfile($config->{SCRIPTDIR}, $cmd),
        args => [
            ["-cfg", $config->{_CONFIG_PATH}, 1],
            ["-input", $input, 1],
            ["-output", $output, 1],
	    ["-flip", $flip,1],
        ],
        inputs    => [$input, $config->{_CONFIG_PATH}],
        outputs   => [$output],
        description => "Generating pseudo assembly"
    });

    my $response = $JEX->submit_workflow($workflow);

    return {
        id => $response->{id},
        success => $JEX->is_successful($response),
        output => $output
    };
}

1;
