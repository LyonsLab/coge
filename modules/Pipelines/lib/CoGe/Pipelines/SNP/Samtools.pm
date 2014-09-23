package CoGe::Pipelines::SNP::Samtools;

use v5.14;
use warnings;
use strict;

use Carp;

use CoGe::Accessory::Web;

require Exporter;
our @ISA = qw(Exporter);
our @EXPORT_OK = qw(build);
our $CONFIG = CoGe::Accessory::Web::get_defaults();

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
    my $cmd = $CONFIG->{SAMTOOLS};
    return "";
}

sub _bcftools {
    my $cmd = $CONFIG->{BCFTOOLS};
    return "";
}

sub _vcfutils {
    my $cmd = $CONFIG->{VCFUTILS};
    return "";
}
