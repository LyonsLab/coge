package CoGe::Services::Data::Experiment;

use Mojo::Base 'Mojolicious::Controller';
use Mojo::JSON qw( decode_json );
use Data::Dumper;
use File::Path qw( mkpath );
use File::Basename qw( basename );
#use IO::Compress::Gzip 'gzip';
use CoGeX;
use CoGe::Services::Auth;
use CoGe::Accessory::Utils;
use CoGe::Accessory::Jex;
use CoGe::Accessory::Workflow;
use CoGe::Accessory::IRODS;

sub search {
    my $self = shift;
    my $search_term = $self->stash('term');

    # Validate input
    if (!$search_term or length($search_term) < 3) {
        $self->render(json => { error => { Error => 'Search term is shorter than 3 characters' } });
        return;
    }

    # Authenticate user and connect to the database
    my ($db, $user) = CoGe::Services::Auth::init($self);

    # Search experiments
    my $search_term2 = '%' . $search_term . '%';
    my @experiments = $db->resultset("Experiment")->search(
        \[
            'experiment_id = ? OR name LIKE ? OR description LIKE ?',
            [ 'experiment_id', $search_term  ],
            [ 'name',        $search_term2 ],
            [ 'description', $search_term2 ]
        ]
    );

    # Filter response
    my @filtered = grep {
        $user->has_access_to_experiment($_);
    } @experiments;

    # Format response
    my @result = map {
      {
        id => int($_->id),
        name => $_->name,
        description => $_->description,
      }
    } @filtered;

    $self->render(json => { experiments => \@result });
}

sub fetch {
    my $self = shift;
    my $id = int($self->stash('id'));

    # Authenticate user and connect to the database
    my ($db, $user) = CoGe::Services::Auth::init($self);

    # Get experiment
    my $experiment = $db->resultset("Experiment")->find($id);
    unless (defined $experiment) {
        $self->render(json => {
            error => { Error => "Item not found" }
        });
        return;
    }

    # Check permissions
    unless ($user->has_access_to_experiment($experiment)) {
        $self->render(json => {
            error => { Auth => "Access denied" }
        }, status => 401);
        return;
    }

    # Format metadata
    my @metadata = map {
        {
            text => $_->annotation,
            link => $_->link,
            type => $_->type->name,
            type_group => $_->type->group
        }
    } $experiment->annotations;

    # Format types
    my @types = map {
        {
            name => $_->name,
            description => $_->description
        }
    } $experiment->types;

    $self->render(json => {
        id => int($experiment->id),
        name => $experiment->name,
        description => $experiment->description,
        version => $experiment->version,
        genome_id  => int($experiment->genome->id),
        source => {
            name => $experiment->source->name,
            description => $experiment->source->description,
            link => $experiment->source->link
        },
        types => \@types,
        metadata => \@metadata,
        restricted => $experiment->restricted ? Mojo::JSON->true : Mojo::JSON->false,
    });
}

sub add {
    my $self = shift;
    my $data = $self->req->json;
    #print STDERR Dumper $data, "\n";
    
    # Authenticate user and connect to the database
    my ($db, $user, $conf) = CoGe::Services::Auth::init($self);
    
    # Get genome
    my $gid = $data->{genome_id};
    my $genome = $db->resultset("Genome")->find($gid);
    unless ($genome) {
        $self->render(json => {
            error => { Error => "Item not found" }
        });
        return;
    }
    
    # Check permissions to add to genome
    unless ($user && ($user->is_admin || $user->is_owner_editor( dsg => $gid ))) {
        $self->render(json => {
            error => { Auth => "Access denied" }
        });
        return;
    }
    
    # TODO validate metadata parameters
    if (!@{ $data->{items} }) {
        $self->render(json => {
            error => { Error => "No data items specified" }
        });
        return;
    }
    
    # Connect to workflow engine and get an id
    my $jex = CoGe::Accessory::Jex->new( host => $conf->{JOBSERVER}, port => $conf->{JOBPORT} );
    my $job = $db->resultset('Job')->create(
        {
            link       => '',
            page       => 'api',
            #process_id => getpid(),
            user_id    => $user->id,
            status     => 0,
        }
    );
    unless ($jex and $job) {
        $self->render(json => {
            error => { Error => "Could not connect to JEX" }
        });
        return;
    }
    
    my $load_id = CoGe::Accessory::Utils::get_unique_id();
    my $tempdir = File::Spec->catdir($conf->{SECTEMPDIR}, 'api', $user->name, $load_id);
    my $staging_dir = File::Spec->catdir($tempdir, 'staging');
    mkpath($staging_dir);
    
    # Create the workflow
    my $workflow = $jex->create_workflow(
        id => $job->id,
        name => 'blah blah blah',
        logfile => $staging_dir.'/log_main.txt'
    );
    
    my %load_params;
    my @staged_files;
    foreach my $item (@{ $data->{items} }) {
        next unless ($item->{type} eq 'irods');
        %load_params = _create_iget_job($conf, $item->{path}, $staging_dir);
        unless ( %load_params ) {
            $self->render(json => {
                error => { Error => "Could not create iget task" }
            });
            return;
        }
        $workflow->add_job(%load_params);
        push @staged_files, $load_params{outputs}[0];
    }
    
    %load_params = _create_load_job($conf, $data, $user, $staging_dir, \@staged_files);
    unless ( $workflow and %load_params ) {
        $self->render(json => {
            error => { Error => "Could not create load task" }
        });
        return;
    }
    $workflow->add_job(%load_params);
    
    
    # Submit the workflow
    my $result = $jex->submit_workflow($workflow);
    my $status = decode_json($result);
    if ($status->{error}) {
        $self->render(json => {
            error => { Error => "Could not submit workflow: " . $status->{error} }
        });
        return;
    }
    
    $self->render(json => 
        {
            success => Mojo::JSON->true,
            job_id => $job->id,
            link => ''
        }
    );
}

#-------------------------------------------------------------------------------

sub _create_iget_job {
    my ($conf, $irods_path, $staging_dir) = @_;
    
    my $dest = $irods_path;
    my $cmd = CoGe::Accessory::IRODS::irods_iget( $irods_path, undef, { no_execute => 1 } );
    
    return (
        cmd => $cmd,
        script => undef,
        args => [
            ['', $dest, 1]
        ],
        inputs => [],
        outputs => [
            File::Spec->catdir($staging_dir, 'irods', $dest)
        ],
        description => "Fetching $irods_path..."    
    );
}

sub _create_load_job {
    my ($conf, $data, $user, $staging_dir, $files) = @_;
    my $cmd = File::Spec->catdir(($conf->{SCRIPTDIR}, "load_experiment.pl"));
    return unless $cmd; # SCRIPTDIR undefined
    
    my $file_str = join(',', map { basename($_) } @$files);
    
    return (
        cmd => $cmd,
        script => undef,
        args => [
            ['-user_name', $user->name, 0],
            ['-name', '"' . $data->{name} . '"', 0],
            ['-desc', '"' . $data->{description} . '"', 0],
            ['-version', '"' . $data->{version} . '"', 0],
            ['-restricted', ( $data->{name} ? 1 : 0 ), 0],
            ['-gid', $data->{genome_id}, 0],
            ['-source_name', '"' . $data->{source_name} . '"', 0],
            #['-types', qq{"Expression"}, 0], # FIXME
            #['-annotations', $ANNOTATIONS, 0],
            ['-staging_dir', 'staging', 0],
            ['-file_type', "csv", 0], # FIXME
            ['-data_file', "$file_str", 0],
            ['-config', $conf->{_CONFIG_PATH}, 1]
        ],
        inputs => [
            ($conf->{_CONFIG_PATH}, @$files)
        ],
        outputs => [
            [$staging_dir, 1],
            File::Spec->catdir($staging_dir, 'log.done')
        ],
        description => "Loading experiment data..."
    );
}

1;
