package CoGe::Services::Data::Feature2;

use Mojo::Base 'Mojolicious::Controller';
use Mojo::JSON;
use CoGe::Services::Auth qw(init);
use CoGe::Core::Feature qw(search_features);
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
    my $feature = $db->resultset("Feature")->find($id);
    unless (defined $feature) {
        $self->render(json => {
            error => { Error => "Item not found"}
        });
        return;
    }

    # Verify that user has access to all genomes associated with this feature
    my $first_genome;
    foreach my $genome ($feature->genomes) {
        unless ( !$genome->restricted || (defined $user && $user->has_access_to_genome($genome)) ) {
            $self->render(json => {
                error => { Auth => "Access denied"}
            }, status => 401);
            return;
        }
        $first_genome = $genome unless $first_genome;
    }
    
    # Generate response
    $self->render(json => {
        id => int($feature->id),
        type => $feature->type->name,
        start => $feature->start,
        stop => $feature->stop,
        strand => $feature->strand,
        chromosome => $feature->chromosome,
        genome => {
            id => $first_genome->id,
            version => $first_genome->version,
            summary => $first_genome->info
        },
        sequence => $feature->genomic_sequence,
        site_url => 
    });
}

1;
