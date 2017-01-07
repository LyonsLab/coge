package CoGe::Builder::PopGen::SummaryStats;

use Moose;
extends 'CoGe::Builder::Buildable';

use Data::Dumper qw(Dumper);
use File::Spec::Functions qw(catdir catfile);
use CoGe::Core::Storage qw(get_genome_file get_experiment_files get_popgen_result_path get_genome_cache_path);
use CoGe::Builder::CommonTasks;
use CoGe::Accessory::Utils qw(is_gzipped to_filename);

sub get_name {
    #my $self = shift;
    return 'Compute Summary Stats';
}

sub build {
    my %opts = @_;
    my $experiment = $opts{experiment};
    my $input_file = $opts{input_file}; # path to vcf file
    my $user = $opts{user};
    my $wid = $opts{wid};
    my $metadata = $opts{metadata};
    
    # Setup paths
    my $genome = $experiment->genome;
    my $gid = $genome->id;
    my $CONF = CoGe::Accessory::Web::get_defaults();
    my $FASTA_CACHE_DIR = get_genome_cache_path($gid);
    die "ERROR: CACHEDIR not specified in config" unless $FASTA_CACHE_DIR;
    my $fasta_file = get_genome_file($gid);
    my $reheader_fasta = to_filename($fasta_file) . ".reheader.faa";
    
    # Build the workflow's tasks
    my ($task, @tasks, @done_files);
    $task = create_fasta_reheader_job(
        fasta => $fasta_file,
        reheader_fasta => $reheader_fasta,
        cache_dir => $FASTA_CACHE_DIR
    );
    push @tasks, $task;
    $fasta_file = $task->{outputs}->[0];
    
    # Check if genome has annotations
    my $isAnnotated = $genome->has_gene_features;
    
    # Generate cached gff if genome is annotated
    my $gff_file;
    if ($isAnnotated) {
        $task = create_gff_generation_job(gid => $gid, organism_name => $genome->organism->name);
        push @tasks, $task;
        $gff_file = $task->{outputs}->[0];
    }
    else {
        # TODO fail here, GFF is required by sumstats.pl
    }
    
    # Get experiment VCF file
    my $vcf_file = get_experiment_files($experiment->id, $experiment->data_type)->[0];
    
    # Compress file using bgzip
    $task = create_bgzip_job( $vcf_file );
    push @tasks, $task;
    $vcf_file = $task->{outputs}->[0];
    
    # Create a Tabix index
    $task = create_tabix_index_job( $vcf_file, 'vcf' );
    push @tasks, $task;
    
    # Determine output path for result files
    my $result_path = get_popgen_result_path($experiment->id);
    die "Cannot determine sumstat result path" unless $result_path;
    
    # Compute summary stats
    $task = create_sumstats_job(
        vcf => $vcf_file,
        gff => $gff_file,
        fasta => $fasta_file,
        output_path => $result_path
    );
    push @tasks, $task;
    
    # Add workflow result
    $task = add_workflow_result(
        username => $user->name,
        wid => $wid,
        result => {
            type => "popgen",
            experiment_id => $experiment->id,
            name => "Diversity analysis results"
        },
        dependency => catfile($result_path, 'sumstats.done')
    );
    push @tasks, $task;    
    
    return {
        tasks => \@tasks,
        done_files => \@done_files
    };
}

1;