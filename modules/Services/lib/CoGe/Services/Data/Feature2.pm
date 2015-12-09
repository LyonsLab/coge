package CoGe::Services::Data::Feature2;

use Mojo::Base 'Mojolicious::Controller';
use Mojo::JSON;
use CoGe::Services::Auth qw(init);
use CoGe::Core::Feature qw(search_features get_feature);
use Data::Dumper;

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

    # Search features
    my @results = search_features($db, $search_term);
    
    $self->render(json => { features => \@results });
}

sub fetch {
    my $self = shift;
    my $id = int($self->stash('id'));

    # Authenticate user and connect to the database
    my ($db, $user) = CoGe::Services::Auth::init($self);

    # Get feature from DB
    my $feature = $db->resultset('Feature')->find($id);
    unless (defined $feature) {
        $self->render(json => {
            error => { Error => 'Item not found' }
        });
        return;
    }

    # Verify that user has access to all genomes associated with this feature
    foreach my $genome ($feature->genomes) {
        unless ( !$genome->restricted || (defined $user && $user->has_access_to_genome($genome)) ) {
            $self->render(json => {
                error => { Auth => 'Access denied' }
            }, status => 401);
            return;
        }
    }
    
    # Generate response
    $self->render(json => get_feature(feature => $feature));
}

sub sequence {
    my $self = shift;
    my $id = int($self->stash('id'));

    # Authenticate user and connect to the database
    my ($db, $user) = CoGe::Services::Auth::init($self);

    # Get feature from DB
    my $feature = $db->resultset('Feature')->find($id);
    unless (defined $feature) {
        $self->render(json => {
            error => { Error => 'Item not found' }
        });
        return;
    }

    # Verify that user has access to all genomes associated with this feature
    foreach my $genome ($feature->genomes) {
        unless ( !$genome->restricted || (defined $user && $user->has_access_to_genome($genome)) ) {
            $self->render(json => {
                error => { Auth => 'Access denied' }
            }, status => 401);
            return;
        }
    }
    
    # Generate response
    $self->render(text => $feature->genomic_sequence);
}

1;
