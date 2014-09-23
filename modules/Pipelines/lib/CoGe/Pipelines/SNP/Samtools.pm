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

sub find_snps {
    my $options = shift;

    my @subcommands =  (
        pileup($options),
        bcftools($options),
    );

    # Pipe commands together
    my $command = join "|", @subcommands;

    return {
        command => $command,
    };
}

sub filter_snps {
    my $options = shift;

    my @subcommands =  (
        bcftools($options),
        vcfutils($options),
    );

    # Pipe commands together
    my $command = join "|", @subcommands;

    return {
        command => $command,
    };
}

sub samtools {
    return "";
}

sub bcftools {
    return "";
}

sub vcfutils {
    return "";
}
