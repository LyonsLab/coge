package CoGe::Builder::Load::BatchExperiment;

use Moose;
extends 'CoGe::Builder::Buildable';

use Data::Dumper qw(Dumper);
use Switch;
use File::Spec::Functions qw(catfile);

use CoGe::Accessory::Utils qw(get_unique_id);
use CoGe::Core::Storage qw(get_upload_path);
use CoGe::Builder::CommonTasks;

sub get_name {
    my $self = shift;
    my $metadata = $self->params->{metadata};
    my $info = '"' . $metadata->{name};
    $info .= ": " . $metadata->{description} if $metadata->{description};
    $info .= " (v" . $metadata->{version} . ")";
    $info .= '"';
    return "Load batch " . $info;
}

sub build {
    my $self = shift;
    
    # Validate inputs
    my $genome_id = $self->params->{genome_id};
    return unless $genome_id;
    my $data = $self->params->{source_data};
    return unless (defined $data && @$data);
    my $metadata = $self->params->{metadata};
    return unless $metadata;
    my $load_id = $self->params->{load_id} || get_unique_id();
    
    # mdb added 2/25/15 - convert from Mojolicious boolean: bless( do{\\(my $o = 1)}, 'Mojo::JSON::_Bool' )
    $metadata->{restricted} = $metadata->{restricted} ? 1 : 0;

    # Get organism
    my $genome = $self->db->resultset('Genome')->find($genome_id);
    return unless $genome;
    
    # Determine file type if not set
    my $file_type = $data->[0]->{file_type}; # type of first data file
    ($file_type) = detect_data_type($file_type, $data->[0]->{path}) unless $file_type;
    
    #
    # Build workflow
    #
    my (@tasks, @input_files, @done_files);
    
    # Create tasks to retrieve files
    my $upload_dir = get_upload_path($self->user->name, $load_id);
    my $data_workflow = create_data_retrieval_workflow(upload_dir => $upload_dir, data => $data);
    push @tasks, @{$data_workflow->{tasks}} if ($data_workflow->{tasks});
    push @input_files, @{$data_workflow->{outputs}} if ($data_workflow->{outputs});
    
    # Add load batch task
    my $task = create_load_batch_job(
        user => $self->user,
        staging_dir => $self->staging_dir,
        wid => $self->workflow->id,
        gid => $genome->id,
        nid => $self->params->{notebook_id},
        input_files => \@input_files,
        metadata => $metadata,
    );
    push @tasks, $task;
    push @done_files, $task->{outputs}->[1];
    
    #print STDERR Dumper \@tasks, "\n";
    $self->workflow->add_jobs(\@tasks);
    
    return 1;
}

1;
