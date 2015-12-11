package CoGe::Services::Data::Feature2;

use Mojo::Base 'Mojolicious::Controller';
use Mojo::JSON;
use CoGe::Services::Auth qw(init);
use CoGe::Core::Feature qw(search_features get_feature get_genome_for_feature);
use Data::Dumper;

sub search {
    my $self = shift;
    my $search_term = $self->stash('term');
    my $exact = $self->param('exact'); # perform exact match

    # Validate input
    if (!$search_term or length($search_term) < 3) {
        $self->render(json => { error => { Error => 'Search term is shorter than 3 characters' } });
        return;
    }

    # Authenticate user and connect to the database
    my ($db, $user) = CoGe::Services::Auth::init($self);

    # Search features
    my @results = search_features( db => $db, user => $user, search_term => $search_term, exact => $exact );
    
    $self->render(json => { features => \@results });
}

sub fetch {
    my $self = shift;
    my $id = int($self->stash('id'));

    # Authenticate user and connect to the database
    my ($db, $user) = CoGe::Services::Auth::init($self);

    # Generate response
    $self->render(json => get_feature(db => $db, user => $user, fid => $id));
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

    # Verify that user has access to a genome associated with this feature
    unless (get_genome_for_feature(feature => $feature)) {
        $self->render(json => {
            error => { Error => 'Access denied' }
        });
        return;
    }
    
    # Generate response
    $self->render(text => $feature->genomic_sequence);
}

1;
