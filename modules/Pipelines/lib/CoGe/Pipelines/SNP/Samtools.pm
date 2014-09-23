package CoGe::Pipelines::SNP::Samtools;

use v5.14;
use warnings;
use strict;

use Carp;

require Exporter;
our @ISA = qw(Exporter);
our @EXPORT_OK = qw(build);

sub build {
    carp "Not Implemented";
}

sub _find_snps {
    my $options = shift;

    my @subcommands =  (
        _samtools($options->{samtools}),
        _bcftools($options->{bcf}),
    );

    # Pipe commands together
    my $command = join "|", @subcommands;

    return {
        command => $command,
        inputs => $options->{inputs},
        outputs => $options->{outputs},
    };
}

sub _filter_snps {
    my $options = shift;

    my @subcommands =  (
        _bcftools($options->{bcf}),
        _vcfutils($options->{vcf}),
    );

    # Pipe commands together
    my $command = join "|", @subcommands;

    return {
        command => $command,
        inputs => $options->{inputs},
        outputs => $options->{outputs},
    };
}

sub _samtools {
    return "";
}

sub _bcftools {
    return "";
}

sub _vcfutils {
    return "";
}
