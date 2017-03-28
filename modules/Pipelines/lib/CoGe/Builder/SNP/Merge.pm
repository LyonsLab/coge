package CoGe::Builder::SNP::Merge;

use Moose;
extends 'CoGe::Builder::Buildable';

use Switch;
use Data::Dumper qw(Dumper);
use File::Spec::Functions qw(catdir catfile);
use Clone qw(clone);
use File::Basename qw(basename dirname);
use File::Path qw(make_path);
use String::ShellQuote qw(shell_quote);

use CoGe::Accessory::Web qw(get_defaults get_command_path url_for);
use CoGe::Accessory::Utils;
use CoGe::Accessory::TDS;
use CoGe::Core::Storage;
use CoGe::Core::Metadata;
use CoGe::Exception::Generic;

# Outputs
has gvcf => (is => 'rw', isa => 'Str'); # GVCF file

my $JAVA_MAX_MEM = '200g'; # suggested memory limit for GATK GenotypeGVCFs by P. Ozersky from TERRA-REF
my ($PICARD, $GATK);
sub BUILD { # called immediately after constructor
    my $self = shift;

    $PICARD = $self->conf->{PICARD};
    unless ($PICARD) {
        CoGe::Exception::Generic->throw(message => 'Missing PICARD in config file');
    }

    $GATK = $self->conf->{GATK};
    unless ($GATK) {
        CoGe::Exception::Generic->throw(message => 'Missing GATK in config file');
    }
}

sub build {
    my $self = shift;
    #my %opts = @_;

    # Validate inputs not already checked in Request
    my $metadata = $self->params->{metadata};
    unless ($metadata && $metadata->{name}) {
        CoGe::Exception::MissingField->throw(message => "Missing metadata name field");
    }
    $metadata->{source_name} = $self->user->display_name unless $metadata->{source_name};
    $metadata->{version} = 1 unless $metadata->{version};

    # Verify experiments and build list of VCF files for all experiments
    my $experiments = $self->request->experiments;
    my $gid = $experiments->[0]->genome->id;
    my @vcf_files;
    foreach my $experiment (@$experiments) {
        if ($experiment->genome->id != $gid) {
            CoGe::Exception::Generic->throw(message => 'Inconsistent genome for experiment ' . $experiment->id);
        }
        if ($experiment->data_type != $DATA_TYPE_POLY) {
            CoGe::Exception::Generic->throw(message => 'Type is not polymorphism for experiment ' . $experiment->id);
        }
        my $vcf = get_experiment_files($experiment->id, $experiment->data_type)->[0];
        push @vcf_files, $vcf;
    }

    # Reheader and index FASTA file
    my ($reheader_fasta) = $self->add(
        $self->reheader_fasta($gid)
    );
    $self->add(
        $self->index_fasta($reheader_fasta)
    );
    $self->add(
        $self->fasta_dict($reheader_fasta, $gid)
    );

    # There are two ways to merge
    my $method = lc($self->params->{method});
    if ($method eq 'combine') { # Combine individual gVCFs into multisample GVCF
        $self->add_to_previous(
            $self->gatk_CombineGVCFs(
                input_fasta => $reheader_fasta,
                input_vcfs  => \@vcf_files
            )
        );
    }
    elsif ($method eq 'genotype') { # Merge individual gVCFs into a single-sample GVCF
        $self->add_to_previous(
            $self->gatk_GenotypeGVCFs(
                input_fasta => $reheader_fasta,
                input_vcfs  => \@vcf_files
            )
        );
    }

    # Add metadata describing this pipeline to the resulting experiment
    my $md = $self->generate_additional_metadata();
    push @{$self->params->{additional_metadata}}, @$md;

    # Load GVCF experiment
    $self->add_to_previous(
        $self->load_vcf(
            gid => $gid,
            vcf => $self->previous_output
        )
    );
}

sub gatk_CombineGVCFs {
    my $self = shift;
    my %opts = @_;
    my $input_fasta = $opts{input_fasta};
    my $input_vcfs  = $opts{input_vcfs}; # ref to array of VCF files

    my ($first_vcf) = @$input_vcfs;
    my $output_vcf = to_filename_base($first_vcf) . '.vcf';

    my $args = [
        ['-R', qq[$input_fasta.fa], 0],
        ['-o', $output_vcf,         1]
    ];
    push @$args, map { ['-V', $_, 1 ] } @$input_vcfs;

    return {
        cmd => "ln -sf $input_fasta $input_fasta.fa && " . # GATK expects the filename to end in .fa or .fasta
               "ln -sf $input_fasta.fai $input_fasta.fa.fai && " .
               "ln -sf $input_fasta.dict $input_fasta.fa.dict && " .
                qq[java -Xmx$JAVA_MAX_MEM -jar $GATK -T CombineGVCFs],
        args => $args,
        inputs => [
            @$input_vcfs,
            $input_fasta,
            qq[$input_fasta.fai],
            qq[$input_fasta.dict]
        ],
        outputs => [
            catfile($self->staging_dir, $output_vcf)
        ],
        description => "Combining VCFs using GATK CombineGVCFs"
    };
}

sub gatk_GenotypeGVCFs {
    my $self = shift;
    my %opts = @_;
    my $input_fasta = $opts{input_fasta};
    my $input_vcfs  = $opts{input_vcfs}; # ref to array of VCF files

    my ($first_vcf) = @$input_vcfs;
    my $output_vcf = to_filename_base($first_vcf) . '.vcf';

    my $args = [
        ['-R', qq[$input_fasta.fa], 0],
        ['-o', $output_vcf,         1]
    ];
    push @$args, map { ['-V', $_, 1 ] } @$input_vcfs;

    return {
        cmd => "ln -sf $input_fasta $input_fasta.fa && " . # GATK expects the filename to end in .fa or .fasta
               "ln -sf $input_fasta.fai $input_fasta.fa.fai && " .
               "ln -sf $input_fasta.dict $input_fasta.fa.dict && " .
                qq[java -Xmx$JAVA_MAX_MEM -jar $GATK -T GenotypeGVCFs],
        args => $args,
        inputs => [
            @$input_vcfs,
            $input_fasta,
            qq[$input_fasta.fai],
            qq[$input_fasta.dict]
        ],
        outputs => [
            catfile($self->staging_dir, $output_vcf)
        ],
        description => "Genotyping VCFs using GATK GenotypeGVCFs"
    };
}

sub generate_additional_metadata {
    my $self = shift;
    my $experiments = $self->request->experiments;

    my @md = ({
        type => 'SNP merge method',
        text => lc($self->params->{method})
    });

    foreach my $experiment (@$experiments) {
        push @md, {
            group => 'Source experiments',
            type => 'Experiment ID',
            text => $experiment->id,
            link => url_for('ExperimentView.pl', eid => $experiment->id) #FIXME hardcoded page url
        };
    }

    return \@md;
}

__PACKAGE__->meta->make_immutable;

1;
