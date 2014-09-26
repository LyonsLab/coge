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
    my $opts = shift;

    my @subopts = [
        {
            command => $CONFIG->{SAMTOOLS},
            subtask => "mpileup",
            args    => {
                u => [],
                f => [],
            },
            inputs => [
                $opts->{reference},
                $opts->{align1},
                $opts->{align2},
            ],
        },
        {
            command => $CONFIG->{BCFTOOLS},
            subtask => "view",
            args    => {
                b => [],
                v => [],
                c => [],
                g => [],
            },
            inputs => [
            ],
        },
    ];

    # Generate subcommand strings
    my @subcommands =  map { _subcommand($_) } @subopts;

    # Pipe commands together
    my $command = join "|", @subcommands;

    # Get the name of the output file
    my $output = qq[$opts->{basename}.raw.bcf];

    return {
        command => qq[$command - > $output],
        inputs => [
            catfile($opts->{result_dir}, $opts->{reference}),
            catfile($opts->{result_dir}, $opts->{align1}),
            catfile($opts->{result_dir}, $opts->{align2}),
        ],
        outputs => [
            catfile($opts->{result_dir}, $output),
        ],
    };
}

sub _filter_snps {
    my $opts = shift;
    my $depth = $opts->{depth} || 100;

    my $subopts = [
        bcf => {
            command => $CONFIG->{BCFTOOLS},
            subtask => "view",
            args    => {
            },
            inputs => [
                $opts->{snps}
            ],
        },

        vcf => {
            command => $CONFIG->{VCFTOOLS},
            subtask => "varFilter",
            args    => {
                D => [$depth]
            },
            inputs => [
            ],
        },
    ];

    my @subcommands =  (
        _subcommand($subopts->{bcf}),
        _subcommand($subopts->{vcf}),
    );

    # Pipe commands together
    my $command = join "|", @subcommands;

    # Get the name of the output file
    my $output = qq[$opts->{basename}.flt.vcf];

    return {
        cmd => qq[$command > $output],
        inputs  => [
            catfile($opts->{result_dir}, $opts->{snps}),
        ],
        outputs => [
            catfile($opts->{result_dir}, $output),
        ],
    };
}

sub _subcommand {
    my $opts = shift;
    my $cmd = $opts->{command};
    my $subtask = $opts->{subtask};
    my $args = $opts->{args};

    # create parameter string
    my @params = map { @{$args->{$_}} ? qq[-$_ @{$args->{$_}}] : qq[-$_] } keys $args;

    return qq[$cmd $subtask @params @{$opts->{inputs}}];
}
