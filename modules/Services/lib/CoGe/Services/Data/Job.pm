package CoGe::Services::Data::Job;

use Mojo::Base 'Mojolicious::Controller';
#use IO::Compress::Gzip 'gzip';
use Data::Dumper;
use File::Basename qw( basename );
use File::Spec::Functions qw( catdir catfile );
use CoGeX;
use CoGe::Services::Auth;
use CoGe::Accessory::Jex;
use CoGe::Accessory::TDS;
use CoGe::Core::Storage qw( get_workflow_paths );
use CoGe::Requests::RequestFactory;
use CoGe::Pipelines::PipelineFactory;

sub add {
    my $self = shift;
    my $payload = $self->req->json;

    # Authenticate user and connect to the database
    my ($db, $user, $conf) = CoGe::Services::Auth::init($self);

    # User authentication is required
    unless (defined $user) {
        return $self->render(json => {
            error => { Auth => "Access denied" }
        });
    }

    my $jex = CoGe::Accessory::Jex->new( host => $conf->{JOBSERVER}, port => $conf->{JOBPORT} );
    my $request_factor = CoGe::Requests::RequestFactory->new(db => $db, user => $user, jex => $jex);
    my $request_handler = $request_factor->get($payload);

    # Validate the request has all required fields
    unless ($request_handler and $request_handler->is_valid) {
        return $self->render(json => {
            error => { Invalid => "The request was not valid." }
        });
    }

    # Check users permissions to execute the request
    unless ($request_handler->has_access) {
        return $self->render(json => {
            error => { Auth => "Access denied" }
        });
    }

    my $pipeline_factory = CoGe::Pipelines::PipelineFactory->new(conf => $conf);
    my $workflow = $pipeline_factory->get($payload);

    return $self->render(json => $request_handler->execute($workflow));
}

sub fetch {
    my $self = shift;
    my $id = $self->stash('id');

    # Authenticate user and connect to the database
    my ($db, $user, $conf) = CoGe::Services::Auth::init($self);

    # User authentication is required
    unless (defined $user) {
        $self->render(json => {
            error => { Auth => "Access denied" }
        });
        return;
    }

    # Get job status from JEX
    my $jex = CoGe::Accessory::Jex->new( host => $conf->{JOBSERVER}, port => $conf->{JOBPORT} );
    my $job_status = $jex->get_job($id);
    unless ($job_status) {
        $self->render(json => {
            error => { Error => "Item not found" }
        });
        return;
    }

#    unless ($job_status->{status} =~ /completed|running/i) {
#        $self->render(json => {
#            id => int($id),
#            status => $job_status->{status}
#        });
#        return;
#    }

    # TODO COGE-472: add read from "debug.log" in results path if JEX no longer has the log

    # Add tasks (if any)
    my @tasks;
    foreach my $task (@{$job_status->{jobs}}) {
        my $t = {
            started => $task->{started},
            ended => $task->{ended},
            elapsed => $task->{elapsed},
            description => $task->{description},
            status => $task->{status},
            log => undef
        };

        if (defined $task->{output}) {
            foreach (split(/\\n/, $task->{output})) {
                print STDERR $_, "\n";
                next unless ($_ =~ /^log\: /);
                $_ =~ s/^log\: //;
                $t->{log} .= $_ . "\n";
            }
        }

        push @tasks, $t;
    }

    # Add results (if any)
    #FIXME add routine to Storage.pm to get results path
    my ( undef, $result_dir ) = get_workflow_paths( $user->name, $id ); # FIXME mdb 8/22/14 directory "experiments" used to be here so whatever puts results there is broken now
    my @results;
    if (-r $result_dir) {
        # Get list of result files in results path
        opendir(my $fh, $result_dir);
        foreach my $file ( readdir($fh) ) {
            my $fullpath = catfile($result_dir, $file);
            next unless -f $fullpath;

            my $name = basename($file);
            push @results, {
                type => 'http',
                name => $name,
                path => $conf->{SERVER}.'api/v1/jobs/'.$id.'/results/'.$name # FIXME move api path into conf file ...?
            };
        }
        closedir($fh);
    }

    $self->render(json => {
        id => int($id),
        status => $job_status->{status},
        tasks => \@tasks,
        results => \@results
    });
}

sub results {
    my $self = shift;
    my $id = $self->stash('id');
    my $name = $self->stash('name');

    # Authenticate user and connect to the database
    my ($db, $user, $conf) = CoGe::Services::Auth::init($self);

    # User authentication is required
    unless (defined $user) {
        $self->render(json => {
            error => { Auth => "Access denied" }
        });
        return;
    }

    #FIXME add routine to Storage.pm to get results path
    my $result_dir = catdir($conf->{SECTEMPDIR}, 'results', 'experiment', $user->name, $id);
    my $result_file = catfile($result_dir, $name);

#    print STDERR $result_file, "\n";

    unless (-r $result_file) {
        $self->render(json => {
            error => { Error => "Item not found" }
        });
        return;
    }

    my $pResult = CoGe::Accessory::TDS::read($result_file);

    $self->render(json => $pResult);
}

1;
