package CoGe::Pipelines::SNP::Samtools;

use v5.14;
use warnings;
use strict;

use Carp;
use File::Basename qw(basename);

use CoGe::Accessory::Web;

require Exporter;
our @ISA = qw(Exporter);
our @EXPORT_OK = qw(build);
our $CONFIG = CoGe::Accessory::Web::get_defaults();

sub build {
}

sub _find_snps {
    my $options = shift;

    my @subcommands =  (
        _subcommand($options->{samtools}),
        _subcommand($options->{bcf}),
    );

    # Get the name of the output file
    my $output = basename($options->{output});

    # Pipe commands together
    my $command = join "|", @subcommands;

    return {
        command => qq[$command - > $options],
        inputs => $options->{inputs},
        outputs => $options->{outputs},
    };
}

sub _filter_snps {
    my $options = shift;

    my @subcommands =  (
        _subcommand($options->{bcf}),
        _subcommand($options->{vcf}),
    );

    # Get the name of the output file
    my $output = basename($options->{output});

    # Pipe commands together
    my $command = join "|", @subcommands;

    return {
        command => qq[$command > $output],
        inputs  => $options->{inputs},
        outputs => $options->{outputs},
    };
}

sub _subcommand {
    my $opts = shift;
    my $cmd = $opts->{command};
    my $subtask = $opts->{subtask};
    my $args = $opts->{args};

    # Filter parameters and create parameter string
    my @params = map { qq[-$_ @{$args->{$_}}] } grep { @{$args->{$_}} } keys $args;

    return qq[$cmd $subtask @params @{$opts->{inputs}}];
}
