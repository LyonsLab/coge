package CoGe::Builder::PopGen::Diversity;

use Moose;
with qw(CoGe::Builder::Buildable);

use Data::Dumper qw(Dumper);
use CoGe::Builder::CommonTasks;
use CoGe::Accessory::Utils qw(is_gzipped);

sub get_name {
    #my $self = shift;
    return 'Compute Summary Stats';
}

sub build {
    my $self = shift;
    
    # Required arguments
    my $genome = $opts->{genome};
    my $input_file = $opts->{input_file}; # path to vcf file
    my $user = $opts->{user};
    my $wid = $opts->{wid};
    my $metadata = $opts->{metadata};
    
    # Setup paths
    my $gid = $genome->id;
    my $FASTA_CACHE_DIR = catdir($CONF->{CACHEDIR}, $gid, "fasta");
    die "ERROR: CACHEDIR not specified in config" unless $FASTA_CACHE_DIR;
    my ($staging_dir, $result_dir) = get_workflow_paths($user->name, $wid);
    my $fasta_file = get_genome_file($gid);
    my $reheader_fasta =  to_filename($fasta_file) . ".reheader.faa";
    
    # Build the workflow's tasks
    my ($task, @tasks, @done_files);
    $task = create_fasta_reheader_job(
        fasta => $fasta_file,
        reheader_fasta => $reheader_fasta,
        cache_dir => $FASTA_CACHE_DIR
    );
    push @tasks, $task;
    my $fasta_file = $task->{outputs}->[0];
    
    # Check if genome has annotations
    my $isAnnotated = $genome->has_gene_features;
    
    # Generate cached gff if genome is annotated
    my $gff_file;
    if ($isAnnotated) {
        my $gff = create_gff_generation_job(gid => $gid, organism_name => $genome->organism->name);
        $gff_file = $gff->{outputs}->[0];
        push @tasks, $gff;
    }
    else {
        # TODO fail here, GFF is required by sumstats.pl
    }
    
    if (is_gzipped($input_file)) {
        $task = create_gunzip_job( input_file => $input_file );
        push @tasks, $task;
        $input_file = $task->{outputs}->[0];
    }
    
    $task = create_bgzip_job( input_file => $input_file );
    $input_file = $task->{outputs}->[0];
    
    $task = create_tabix_index_job(
        input_file => $input_file,
        index_type => 'vcf'
    );
    
    my $result_path = '';
    
    my $sumstats_task = create_sumstats_job(
        vcf => $input_file,
        gff => $gff_file
        fasta => $fasta_file,
        output_path => $result_path
    );
    push @tasks, $sumstats_task;

    return {
        tasks => \@tasks,
        done_files => \@done_files
    };
}

1;