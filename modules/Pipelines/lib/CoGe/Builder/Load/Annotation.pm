package CoGe::Builder::Load::Annotation;

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
    my $gid = $self->params->{genome_id};
    my $metadata = $self->params->{metadata};

    my $genome = $self->db->resultset('Genome')->find($gid);
    my $info = '"';
    $info .= $genome->organism->name;
    $info .= " (" . $metadata->{name} . ")"  if $metadata->{name};
    $info .= ": " . $metadata->{description} if $metadata->{description};
    $info .= " (v" . $metadata->{version} . ")";
    $info .= '"';
    return "Load Annotation " . $info;
}

sub build {
    my $self = shift;
    
    # Validate inputs
    my $gid = $self->params->{genome_id};
    unless ($gid) {
        Mojo::Exception->throw("Missing genome_id");
    }
    my $data = $self->params->{source_data};
    unless (defined $data && @$data) {
        Mojo::Exception->throw("Missing source_data");
    }
    my $metadata = $self->params->{metadata};
    unless ($metadata) {
        Mojo::Exception->throw("Missing metadata");
    }
    my $load_id = $self->params->{load_id} || get_unique_id();

    # Get genome
    my $genome = $self->db->resultset('Genome')->find($gid);
    unless ($genome) {
        Mojo::Exception->throw("Genome $gid not found");
    }
    
    #
    # Build workflow
    #
    my (@tasks, @input_files);
    
    # Create tasks to retrieve files
    my $upload_dir = get_upload_path($self->user->name, $load_id);
    my $data_workflow = create_data_retrieval_workflow(upload_dir => $upload_dir, data => $data);
    push @tasks, @{$data_workflow->{tasks}} if ($data_workflow->{tasks});
    push @input_files, @{$data_workflow->{outputs}} if ($data_workflow->{outputs});
    
    # Submit workflow to generate annotation
    my $load_task = create_load_annotation_job(
        user => $self->user,
        staging_dir => $self->staging_dir,
        wid => $self->workflow->id,
        gid => $genome->id,
        input_file => $input_files[0],
        metadata => $metadata
    );
    push @tasks, $load_task;
    
    #print STDERR Dumper \@tasks, "\n";
    $self->workflow->add_jobs(\@tasks);
    
    return 1;
}

__PACKAGE__->meta->make_immutable;

1;
