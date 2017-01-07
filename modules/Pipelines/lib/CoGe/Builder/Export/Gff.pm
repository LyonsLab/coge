package CoGe::Builder::Export::Gff;

use Moose;
extends 'CoGe::Builder::Buildable';

use CoGe::Accessory::IRODS qw(irods_get_base_path);
use CoGe::Accessory::Utils qw(sanitize_name);
use CoGe::Accessory::Web qw(download_url_for);
use CoGe::Core::Genome qw(get_irods_metadata);
use CoGe::Core::Storage qw(get_genome_cache_path);

use File::Basename qw(basename);
use File::Spec::Functions;
use Data::Dumper;

sub get_name {
    return "Generate/export gff"; #TODO add genome id and output types
}

sub build {
    my $self = shift;
    
    # Verify required parameters and set defaults
    my $dest_type = $self->params->{dest_type};
    $dest_type = "http" unless $dest_type;

    my $gid = $self->params->{gid} || $self->params->{genome_id};
    $self->params->{gid} = $gid; # required by create_gff() call below
    return unless $gid;

    # Get genome
    my $genome = $self->db->resultset("Genome")->find($gid);
    my $genome_name = $self->params->{basename} = sanitize_name($genome->organism->name);

    # Parse output types
    my $outputs = $self->params->{output_types};
    my %outputs;
    %outputs = map { lc($_) => 1 } @$outputs if $outputs;
    $outputs{'gff'} = 1 if (keys %outputs == 0); # default

    # Generate GFF/BED/TBL generation and export tasks
    foreach my $output_type (keys %outputs) {
        # Generate output filename based on params
        my $param_string = join( "-", map { $_ . ($self->params->{$_} // 0) } qw(annos cds id_type nu upa add_chr) );
        my $output_filename = $genome_name . "_" . $param_string;
        $output_filename .= "_" . $self->params->{chr} if $self->params->{chr};
        $output_filename .= ".gid" . $gid;
        $output_filename =~ s/\s+/_/g;
        $output_filename =~ s/\)|\(/_/g;
        $output_filename .= $output_type;
        $output_filename = catfile( get_genome_cache_path($gid), $output_filename );

        # Add task to generate output file(s)
        if ($output_type eq 'gff') {
            $self->add_task(
                $self->create_gff(
                    %{$self->params},
                    output_file => $output_filename
                )
            );
        }
        if ($output_type eq 'bed') {
            $self->add_task(
                $self->create_bed(
                    %{$self->params},
                    output_file => $output_filename
                )
            );
        }
        if ($output_type eq 'tbl') {
            $self->add_task(
                $self->create_tbl(
                    %{$self->params},
                    output_file => $output_filename
                )
            );
        }

        # Add tasks to export/download output file(s)
        if ($dest_type eq "irods") { # irods export
            # Set IRODS destination path
            my $irods_base = $self->params->{dest_path};
            $irods_base = irods_get_base_path($self->user->name) unless $irods_base;
            my $irods_dest = catfile($irods_base, basename($output_filename));

            # Export file task
            $self->add_task_chain(
                $self->export_to_irods(
                    src_file  => $output_filename,
                    dest_file => $irods_dest,
                    overwrite => $self->params->{overwrite}
                )
            );

            # Set file metadata task
            my $md = get_irods_metadata($genome);
            my $md_file = catfile($self->staging_dir, 'irods_metadata.json');
            CoGe::Accessory::TDS::write($md_file, $md);
            $self->add_task_chain(
                $self->create_irods_imeta(
                    dest_file     => $irods_dest,
                    metadata_file => $md_file
                )
            );

            # Add to results
            $self->add_task_chain(
                $self->add_result(
                    result => {
                        type => 'irods',
                        path => $irods_dest
                    }
                )
            );
        }
        else { # http download
            $self->add_task_chain(
                $self->add_result(
                    result => {
                        type => 'url',
                        path => download_url_for(
                            wid  => $self->workflow->id,
                            file => $output_filename
                        )
                    }
                )
            );
        }
    }
    
    return 1;
}

__PACKAGE__->meta->make_immutable;

1;
