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
use CoGe::Core::Storage;
use CoGe::Core::Metadata;
use CoGe::Builder::SNP::GATK;
use CoGe::Builder::SNP::Platypus;
use CoGe::Builder::SNP::Samtools;
use CoGe::Builder::SNP::CoGeSNPs;
use CoGe::Exception::Generic;

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
        case 'gatk'     { $snp = CoGe::Builder::SNP::GATK->new($self) }
        default {
            CoGe::Exception::Generic->throw(message => 'Invalid SNP method');
        }
    }
    $snp->build(fasta_file => $reheader_fasta, data_files => [$bam_file], is_sorted => $opts{is_sorted});
    $self->add($snp);
    $self->vcf($snp->vcf);
}

sub load_vcf {
    my $self = shift;
    my %opts = @_;
    my $annotations = $opts{annotations};
    my $gid = $opts{gid};
    my $vcf = $opts{vcf};

    my $metadata = $self->params->{metadata};
    my $method   = $self->params->{snp_params}->{method};

    my $desc = 'Single nucleotide polymorphisms' . ($method ? " (determined by $method method)" : '');

    my $output_path = catdir($self->staging_dir, "load_vcf");

    my $annotations_str = '';
    $annotations_str = join(';', @$annotations) if (defined $annotations && @$annotations);

    my @tags = ( 'VCF' ); # add VCF tag
    push @tags, @{$metadata->{tags}} if $metadata->{tags};
    my $tags_str = tags_to_string(\@tags);

    return {
        cmd => catfile($self->conf->{SCRIPTDIR}, "load_experiment.pl"),
        args => [
            ['-user_name',   shell_quote($self->user->name), 0],
            ['-name',        ($metadata->{name} ? shell_quote($metadata->{name}." (SNPs)") : '""'), 0],
            ['-desc',        ($desc ? shell_quote($desc) : '""'), 0],
            ['-version',     ($metadata->{version} ? shell_quote($metadata->{version}) : '""'), 0],
            ['-link',        ($metadata->{link} ? shell_quote($metadata->{link}) : '""'), 0],
            ['-restricted',  shell_quote($metadata->{restricted}), 0],
            ['-exit_without_error_for_empty_input', 1, 0],
            ['-gid',         $gid, 0],
            ['-wid',         $self->workflow->id, 0],
            ['-source_name', ($metadata->{source_name} ? shell_quote($metadata->{source_name}) : '""'), 0],
            ['-tags',        shell_quote($tags_str), 0],
            ['-annotations', shell_quote($annotations_str), 0],
            ['-staging_dir', "./load_vcf", 0],
            ['-file_type',   "vcf", 0],
            ['-data_file',   $vcf, 0],
            ['-config',      $self->conf->{_CONFIG_PATH}, 0]
        ],
        inputs => [
            $vcf,
        ],
        outputs => [
            [$output_path, '1'],
            catfile($output_path, "log.done"),
        ],
        description => "Loading SNPs as new experiment"
    };
}

__PACKAGE__->meta->make_immutable;

1;
