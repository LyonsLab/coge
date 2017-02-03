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
    my $self = shift;

    my $experiment = $self->request->experiment;
    my $genome = $experiment->genome;

    #
    # Build workflow
    #

    # Reheader the fasta file
    my ($reheader_fasta) = $self->add(
        $self->reheader_fasta($genome->id)
    );

    # Check if genome has annotations
    my $isAnnotated = $genome->has_gene_features;
    unless ($isAnnotated) {
        CoGe::Exception::Generic->throw(message => 'Genome must be annotated to compute summary stats');
    }
    
    # Generate cached gff if genome is annotated
    $self->add(
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
    $self->add(
        $self->bgzip($vcf_file)
    );
    $vcf_file = $self->previous_output;
    
    # Create a Tabix index
    $self->add(
        $self->tabix_index($vcf_file, 'vcf')
    );

    # Determine output path for result files
    my $result_path = get_popgen_result_path($experiment->id);
    unless ($result_path) {
        CoGe::Exception::Generic->throw(message => 'Cannot determine sumstat result path');
    }
    
    # Compute summary stats
    $self->add(
        $self->sumstats(
            vcf => $vcf_file,
            gff => $gff_file,
            fasta => $reheader_fasta,
            output_path => $result_path
        )
    );
    
    # Add workflow result
    $self->add_to_previous(
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

sub bgzip {
    my $self = shift;
    my $input_file = shift;
    my $output_file = $input_file . '.bgz';

    my $cmd = get_command_path('BGZIP');

    return {
        cmd => "$cmd -c $input_file > $output_file && touch $output_file.done",
        args => [],
        inputs => [
            $input_file
        ],
        outputs => [
            $output_file,
            "$output_file.done"
        ],
        description => "Compressing " . basename($input_file) . " with bgzip"
    };
}

sub tabix_index {
    my $self = shift;
    my $input_file = shift;
    my $index_type = shift;
    my $output_file = $input_file . '.tbi';

    my $cmd = $self->conf->{TABIX} || 'tabix';

    return {
        cmd => "$cmd -p $index_type $input_file && touch $output_file.done",
        args => [],
        inputs => [
            $input_file,
            $input_file . '.done'
        ],
        outputs => [
            $output_file,
            "$output_file.done"
        ],
        description => "Indexing " . basename($input_file)
    };
}

1;