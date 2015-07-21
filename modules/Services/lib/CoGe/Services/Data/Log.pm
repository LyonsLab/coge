package CoGe::Services::Data::Log;

use Mojo::Base 'Mojolicious::Controller';
#use IO::Compress::Gzip 'gzip';
use Data::Dumper;
use CoGeX;
use CoGe::Services::Auth;

sub fetch {
    my $self = shift;
    my $id   = $self->stash('id');
    my $type = $self->stash('type');

    # Authenticate user and connect to the database
    my ($db, $user, $conf) = CoGe::Services::Auth::init($self);

    # User authentication is required
    unless (defined $user) {
        $self->render(json => {
            error => { Auth => "Access denied" }
        });
        return;
    }

    

    $self->render(json => {
        id => int($id),
        status => $job_status->{status},
        tasks => \@tasks,
        results => $results
    });
}

1;
