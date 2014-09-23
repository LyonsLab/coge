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
    my $opts = shift;
    my $cmd = $CONFIG->{SAMTOOLS};
    my $subtask = $opts->{subtask};

    return qq[$cmd $subtask];
}

sub _bcftools {
    my $opts = shift;
    my $cmd = $CONFIG->{BCFTOOLS};
    my $subtask = $opts->{subtask};

    return qq[$cmd $subtask];
}

sub _vcfutils {
    my $opts = shift;
    my $cmd = $CONFIG->{VCFUTILS};
    my $subtask = $opts->{subtask};

    return qq[$cmd $subtask];
}
