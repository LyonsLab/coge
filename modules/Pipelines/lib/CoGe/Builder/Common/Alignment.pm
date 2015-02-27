package CoGe::Builder::Common::Alignment;

use strict;
use warnings;

use Data::Dumper qw(Dumper);
use File::Spec::Functions qw(catdir catfile);
use CoGe::Core::Storage qw(get_genome_file get_workflow_paths get_upload_path get_genome_cache_path);
use CoGe::Accessory::Utils qw(is_fastq_file to_filename);
use CoGe::Builder::CommonTasks;

our $CONF = CoGe::Accessory::Web::get_defaults();

BEGIN {
    use vars qw ($VERSION @ISA @EXPORT @EXPORT_OK);
    require Exporter;

    $VERSION = 0.1;
    @ISA     = qw (Exporter);
    @EXPORT_OK = qw(build);
}

sub build {
    my %opts = @_;
    my $user        = $opts{user};
    my $wid         = $opts{wid};
    my $input_files = $opts{input_files};
    my $genome      = $opts{genome};
    my $metadata    = $opts{metadata};
    my $options     = $opts{options};
    my $alignment_params = $opts{alignment_params};
    my $cutadapt_params  = $opts{cutadapt_params};
    
    my @tasks;
    
    # Setup paths to data files
    #FIXME this is for LoadExperiment, also need to handle IRODS/FTP/HTTP data from API
    my ($staging_dir, $result_dir) = get_workflow_paths($user->name, $wid);
    my $upload_dir = get_upload_path($user->name, $options->{load_id});
    my @files = map { catfile($upload_dir, $_->{path}) } @$input_files;

    # Check multiple files (if more than one file then all should be FASTQ)
    my $numFastq = 0;
    foreach (@files) {
        $numFastq++ if (is_fastq_file($_));
    }
    return return encode_json({ error => 'Unsupported combination of file types' }) if ($numFastq > 0 and $numFastq != @files);
    return return encode_json({ error => 'Too many files' }) if ($numFastq == 0 and @files > 1);
    
    # Process the fastq input files
    my (@validated, @trimmed);
    foreach my $file (@files) {
        # Decompress (if necessary)
        my $done_file;
        if ( $file =~ /\.gz$/ ) {
            push @tasks, create_gunzip_job($file);
            $file =~ s/\.gz$//;
            $done_file = "$file.decompressed";
        }
        
        # Validate
        my $validate_task = create_validate_fastq_job($file, $done_file);
        push @validated, @{$validate_task->{outputs}}[0];
        push @tasks, $validate_task;
    
        # Trim
        my $trim_reads = 1; #FIXME hook this as option in LoadExperiment interface
        if ($trim_reads) {
            my $trim_task = create_cutadapt_job(
                fastq => $file, 
                validated => "$file.validated", 
                staging_dir => $staging_dir,
                params => $cutadapt_params
            );
            push @trimmed, @{$trim_task->{outputs}}[0];
            push @tasks, $trim_task;
        }
        else {
            push @trimmed, $file;
        }
    }

    # Get genome cache path
    my $gid = $genome->id;
    my $fasta_cache_dir = get_genome_cache_path($gid);

    # Reheader the fasta file
    my $fasta = get_genome_file($gid);
    my $reheader_fasta = to_filename($fasta) . ".reheader.faa";
    push @tasks, create_fasta_reheader_job(
        fasta => $fasta,
        reheader_fasta => $reheader_fasta,
        cache_dir => $fasta_cache_dir
    );

    # Index the fasta file
    push @tasks, create_fasta_index_job(
        fasta => catfile($fasta_cache_dir, $reheader_fasta),
        cache_dir => $fasta_cache_dir
    );
    
    # Generate gff if genome annotated
    my $gff_file;
    if ( $genome->has_gene_features ) {
        my $gff_task = create_gff_generation_job(
            gid => $gid,
            organism_name => $genome->organism->name,
            #validated => "$fastq.validated"
        );
        $gff_file = @{$gff_task->{outputs}}[0];
        push @tasks, $gff_task;
    }

    # Add aligner workflow
    my ($bam, @alignment_tasks);
    if ($alignment_params->{tool} eq 'tophat') {
        ($bam, @alignment_tasks) = create_tophat_workflow(
            gid => $gid,
            fasta => catfile($fasta_cache_dir, $reheader_fasta),
            fastq => \@trimmed,
            validated => \@validated,
            read_type => $alignment_params->{read_type},
            gff => $gff_file,
            staging_dir => $staging_dir,
            params => $alignment_params,
        );
    }
    elsif ($alignment_params->{tool} eq 'gsnap') {
        ($bam, @alignment_tasks) = create_gsnap_workflow(
            gid => $gid,
            fasta => $reheader_fasta,
            fastq => \@trimmed,
            validated => \@validated,
            read_type => $alignment_params->{read_type},
            staging_dir => $staging_dir,
            params => $alignment_params,
        );
    }
    else {
        print STDERR "Error: unrecognized alignment tool '", $alignment_params->{tool}, "'\n";
        return;
    }
    push @tasks, @alignment_tasks;

    # Load alignment
    push @tasks, create_load_bam_job(
        user => $user,
        metadata => $metadata,
        staging_dir => $staging_dir,
        result_dir => $result_dir,
        annotations => '', #FIXME 12/12/14
        wid => $wid,
        gid => $gid,
        bam_file => $bam
    );
    
    # Save outputs for retrieval by downstream tasks
    my %outputs;
    $outputs{reheader_fasta} = $reheader_fasta;
    $outputs{bam_file} = $bam;
    $outputs{gff_file} = $gff_file;
    
    return (\@tasks, \%outputs);
}

1;