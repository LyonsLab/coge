package CoGe::Builder::PopGen::SummaryStats;

use Moose;
extends 'CoGe::Builder::Buildable';

use Data::Dumper qw(Dumper);
use File::Spec::Functions qw(catdir catfile);

use CoGe::Accessory::Utils;
use CoGe::Accessory::Web qw(get_command_path);
use CoGe::Core::Storage;
use CoGe::Core::Metadata;
use CoGe::Exception::Generic;

sub get_name {
    #my $self = shift;
    return 'Compute Summary Stats'; #TODO add source experiment info
}

sub build {
    my %opts = @_;
    my $experiment = $opts{experiment};

    my $genome = $experiment->genome;

    #
    # Build workflow
    #

    # Reheader the fasta file
    $self->add_task(
        $self->reheader_fasta($genome->id)
    );
    my $reheader_fasta = $self->previous_output;
    
    # Check if genome has annotations
    my $isAnnotated = $genome->has_gene_features;
    unless ($isAnnotated) {
        CoGe::Exception::Generic->throw(message => 'Genome must be annotated to compute summary stats');
    }
    
    # Generate cached gff if genome is annotated
    $self->add_task(
        $self->create_gff( #FIXME simplify this
            gid => $genome->id,
            output_file => get_gff_cache_path(
                gid => $genome->id,
                genome_name => sanitize_name($genome->organism->name),
                output_type => 'gff',
                params => {}
            )
        )
    );
    my $gff_file = $self->previous_output;

    # Get experiment VCF file
    my $vcf_file = get_experiment_files($experiment->id, $experiment->data_type)->[0];
    
    # Compress file using bgzip
    $self->add_task(
        $self->bgzip($vcf_file)
    );
    $vcf_file = $self->previous_output;
    
    # Create a Tabix index
    $self->add_task(
        $self->tabix_index($vcf_file, 'vcf')
    );

    # Determine output path for result files
    my $result_path = get_popgen_result_path($experiment->id);
    unless ($result_path) {
        CoGe::Exception::Generic->throw(message => 'Cannot determine sumstat result path');
    }
    
    # Compute summary stats
    $self->add_task(
        $self->sumstats(
            vcf => $vcf_file,
            gff => $gff_file,
            fasta => $fasta_file,
            output_path => $result_path
        )
    );
    
    # Add workflow result
    $self->add_task_chain(
        $self->add_result(
            result => {
                type => "popgen",
                experiment_id => $experiment->id,
                name => "Diversity analysis results"
            }
        )
    );
}

sub sumstats {
    my $self = shift;
    my %opts = @_;
    my $vcf = $opts{vcf};
    my $gff = $opts{gff};
    my $fasta = $opts{fasta};
    my $output_path = $opts{output_path};

    return {
        cmd => catfile($self->conf->{SCRIPTDIR}, "popgen/sumstats.pl"),
        args => [
            ['-vcf',    $vcf,         0],
            ['-gff',    $gff,         0],
            ['-fasta',  $fasta,       0],
            ['-output', $output_path, 0],
            ['-debug',  '',           0]
        ],
        inputs => [
            $vcf,
            $gff,
            $fasta
        ],
        outputs => [
            catfile($output_path, "sumstats.done"),
        ],
        description => "Calculating summary statistics"
    };
}

1;