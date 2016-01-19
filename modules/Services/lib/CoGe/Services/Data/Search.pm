package CoGe::Services::Data::Search;

use Mojo::Base 'Mojolicious::Controller';
use Mojo::JSON;
use CoGeX;
use CoGe::Services::Auth;
use CoGe::Services::Data::Job;
use CoGe::Core::Search;

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
    my $show_users = 0;
    if ($user && $user->is_admin) {
		$show_users = 1;
    }

    my @results = CoGe::Core::Search::search(
        db => $db, 
        user => $user, 
        search_term => $search_term, 
        show_users => $show_users
    );

    $self->render(json => { results => \@results });
}

1;
