package CoGe::Builder::SNP::Analyzer;

use Moose;
extends 'CoGe::Builder::Buildable';

use Switch;
use Data::Dumper qw(Dumper);
use File::Spec::Functions qw(catdir catfile);
use Clone qw(clone);
use File::Basename qw(basename dirname);
use File::Path qw(make_path);
use String::ShellQuote qw(shell_quote);

use CoGe::Accessory::Web qw(get_defaults get_command_path);
use CoGe::Accessory::Utils;
use CoGe::Accessory::TDS;
use CoGe::Core::Storage;
use CoGe::Core::Metadata;
use CoGe::Builder::SNP::GATK;
use CoGe::Builder::SNP::Platypus;
use CoGe::Builder::SNP::Samtools;
use CoGe::Builder::SNP::CoGeSNPs;
use CoGe::Exception::Generic;

# Settings
has NUM_CPUS => (is => 'ro', isa => 'Int', default => 16); # number of CPUs to use for SNP analysis tasks

# Outputs
has vcf => (is => 'rw', isa => 'Str'); # VCF file

sub build {
    my $self = shift;
    my %opts = @_;

    my ($bam_file, $gid);
    if ($opts{data_files}) {
        ($bam_file) = @{$opts{data_files}};
        $gid = $self->request->genome->id;
    }
    else {
        my $experiment = $self->request->experiment;
        $bam_file = get_experiment_files($experiment->id, $experiment->data_type)->[0];
        $gid = $experiment->genome->id;
    }

    unless ($self->params->{metadata}) {
        my $experiment = $self->request->experiment;
        $self->params->{metadata} = { # could almost use experiment->to_hash here except for source_name
            name       => $experiment->name,
            version    => $experiment->version,
            source     => $experiment->source->name,
            restricted => $experiment->restricted
        };
    }

    my $snp_params  = $self->params->{snp_params};

    # Reheader and index FASTA file
    my ($reheader_fasta) = $self->add(
        $self->reheader_fasta($gid)
    );
    $self->add(
        $self->index_fasta($reheader_fasta)
    );

    # Add SNP analysis workflow
    my $snp;
    switch( lc($snp_params->{method}) ) {
        case 'coge'     { $snp = CoGe::Builder::SNP::CoGeSNPs->new($self) }
        case 'samtools' { $snp = CoGe::Builder::SNP::Samtools->new($self) }
        case 'platypus' { $snp = CoGe::Builder::SNP::Platypus->new($self) }
        case /gatk|gatk-haplotype-vcf|gatk-haplotype-gvcf/ {
            $snp = CoGe::Builder::SNP::GATK->new($self)
        }
        default {
            CoGe::Exception::Generic->throw(message => 'Invalid SNP method');
        }
    }
    $snp->build(fasta_file => $reheader_fasta, data_files => [$bam_file]);
    $self->add($snp);
    $self->vcf($snp->vcf);
}

__PACKAGE__->meta->make_immutable;

1;
