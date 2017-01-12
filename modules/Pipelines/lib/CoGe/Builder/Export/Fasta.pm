package CoGe::Builder::Export::Fasta;

use Moose;
extends 'CoGe::Builder::Buildable';

use File::Spec::Functions qw(catfile);

use CoGe::Accessory::IRODS qw(irods_get_base_path);
use CoGe::Accessory::Utils qw(sanitize_name);
use CoGe::Accessory::TDS;
use CoGe::Accessory::Web qw(download_url_for);
use CoGe::Core::Storage qw(get_genome_file);
use CoGe::Core::Genome qw(get_irods_metadata);
use CoGe::Exception::Generic;
use CoGe::Exception::MissingField;

sub get_name {
    return "Export FASTA"; #TODO add genome id
}

sub build {
    my $self = shift;

    # Get genome data file path
    my $gid = $self->params->{gid} || $self->params->{genome_id};
    unless ($gid) {
        CoGe::Exception::MissingField->throw(message => "Missing genome_id");
    }
    my $genome = $self->db->resultset("Genome")->find($gid);
    unless ($genome) {
        CoGe::Exception::Generic->throw(message => "Genome $gid not found");
    }

    # Determine name of exported file
    my $genome_file = get_genome_file($gid);
    my $genome_name = sanitize_name($genome->organism->name);#escape($genome->organism->name);
       $genome_name = 'genome_'.$gid unless $genome_name;
    my $output_file = $genome_name.'.faa';

    # Setup tasks to export/download the file
    my $dest_type = $self->params->{dest_type} // 'http';
    if ($dest_type eq "irods") { # irods export
        # Set IRODS destination path
        my $irods_base = $self->params->{dest_path};
        $irods_base = irods_get_base_path($self->user->name) unless $irods_base;
        my $irods_dest = catfile($irods_base, $output_file);

        # Export file task
        $self->add_task(
            $self->export_to_irods(
                src_file => $genome_file, 
                dest_file => catfile($irods_base, $output_file), 
                overwrite => $self->params->{overwrite} 
            )
        );
        
        # Set file metadata task
        my $md = get_irods_metadata($genome);
        my $md_file = catfile($self->staging_dir, 'irods_metadata.json');
        CoGe::Accessory::TDS::write($md_file, $md);
        $self->add_task_chain(
            $self->create_irods_imeta_add(
                dest_file => $irods_dest,
                metadata_file => $md_file
            )
        );
        
        # Add to results
        $self->add_task_chain(
            $self->add_result(
                result   => {
                    type => 'irods',
                    path => $irods_dest
                }
            )
        );
    }
    else { # http download
        $self->add_task_chain(
            $self->add_result(
                result   => {
                    type => 'url',
                    path => download_url_for( gid => $genome->id )
                }
            )
        );
    }
    
    return 1;
}

__PACKAGE__->meta->make_immutable;

1;
