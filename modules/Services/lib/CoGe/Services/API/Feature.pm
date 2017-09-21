package CoGe::Services::API::Feature;

use Mojo::Base 'Mojolicious::Controller';
use Mojo::JSON qw(decode_json);
use Data::Dumper;

use CoGe::Services::Auth qw(init);
use CoGe::Core::Feature qw(search_features get_feature get_genome_for_feature);
use CoGe::Services::API::Job;
use CoGe::Services::Error;

sub search {
    my $self = shift;
    my $search_term = $self->stash('term');
    my $exact = $self->param('exact'); # perform exact match

    # Validate input
    if (!$search_term or length($search_term) < 3) {
        $self->render(API_STATUS_SEARCHTERM);
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
    
    # Validate input
    unless ($id) {
        $self->render(API_STATUS_MISSING_ID);
        return;
    }

    # Authenticate user and connect to the database
    my ($db, $user) = CoGe::Services::Auth::init($self);

    # Generate response
    $self->render(json => get_feature(db => $db, user => $user, fid => $id));
}

sub fasta {
    my $self = shift;
    my $id = int($self->stash('id'));
    
    # Validate input
    unless ($id) {
        $self->render(API_STATUS_MISSING_ID);
        return;
    }

    # Authenticate user and connect to the database
    my ($db, $user) = CoGe::Services::Auth::init($self);

    # Get feature from DB
    my $feature = $db->resultset('Feature')->find($id);
    unless (defined $feature) {
        $self->render(API_STATUS_NOTFOUND);
        return;
    }

    # Verify that user has access to a genome associated with this feature
    unless (get_genome_for_feature(feature => $feature, user => $user)) {
        $self->render(API_STATUS_UNAUTHORIZED);
        return;
    }
    
    # Generate response
    $self->render(text => $feature->fasta(
        col        => 100
    ));
}

sub protein {
    my $self = shift;
    my $id = int($self->stash('id'));
    
    # Validate input
    unless ($id) {
        $self->render(API_STATUS_MISSING_ID);
        return;
    }

    # Authenticate user and connect to the database
    my ($db, $user) = CoGe::Services::Auth::init($self);

    # Get feature from DB
    my $feature = $db->resultset('Feature')->find($id);
    unless (defined $feature) {
        $self->render(API_STATUS_NOTFOUND);
        return;
    }

    # Verify that user has access to a genome associated with this feature
    unless (get_genome_for_feature(feature => $feature, user => $user)) {
        $self->render(API_STATUS_UNAUTHORIZED);
        return;
    }
    
    # Generate response
    $self->render(text => $feature->fasta(
        col        => 100,
        protein    => 1
    ));
}

sub sequence {
    my $self = shift;
    my $id = int($self->stash('id'));
    
    # Validate input
    unless ($id) {
        $self->render(API_STATUS_MISSING_ID);
        return;
    }

    # Authenticate user and connect to the database
    my ($db, $user) = CoGe::Services::Auth::init($self);

    # Get feature from DB
    my $feature = $db->resultset('Feature')->find($id);
    unless (defined $feature) {
        $self->render(API_STATUS_NOTFOUND);
        return;
    }

    # Verify that user has access to a genome associated with this feature
    unless (get_genome_for_feature(feature => $feature, user => $user)) {
        $self->render(API_STATUS_UNAUTHORIZED);
        return;
    }
    
    # Generate response
    $self->render(text => $feature->genomic_sequence);
}

sub add {
    my $self = shift;
    my $data = $self->req->body; #$self->req->json; # mdb replaced 11/30/16 -- req->json hides JSON errors, doing conversion manually prints them to STDERR
    unless ($data) {
        $self->render(API_STATUS_MISSING_BODY);
        return;
    }
    $data = decode_json($data);
    #print STDERR "CoGe::Services::Data::Genome::add\n", Dumper $data, "\n";

    # Valid data items # TODO move into request validation
    unless ($data->{source_data} && @{$data->{source_data}}) {
        $self->render(API_STATUS_MISSING_DATA);
        return;
    }

    # Alias to Job Submit -- is there a better way to do this using Mojolicious routing?
    my $request = {
        type => 'load_annotation',
        parameters => $data
    };

    return CoGe::Services::API::Job::add($self, $request);
}

1;
