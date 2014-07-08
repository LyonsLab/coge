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

    unless ($job_status->{status} =~ /completed/i) {
        $self->render(json => {
            id => int($id),
            status => $job_status->{status}
        });
        return;
    }

    #FIXME add routine to Storage.pm to get results path
    my $result_dir = catdir($conf->{SECTEMPDIR}, 'results', 'experiment', $user->name, $id);
    my @results;
    opendir(my $fh, $result_dir);
    foreach my $file ( readdir($fh) ) {
        my $fullpath = catfile($result_dir, $file);
        next unless -f $fullpath;

        my $name = basename($file);
        push @results, {
            type => 'http',
            name => $name,
            path => $conf->{SERVER}.'api/v1/jobs/'.$id.'/results/'.$name # FIXME move api path into conf file ...?
        }
    }
    closedir($fh);

    $self->render(json => {
        id => int($id),
        status => $job_status->{status},
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

    print STDERR $result_file, "\n";

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
