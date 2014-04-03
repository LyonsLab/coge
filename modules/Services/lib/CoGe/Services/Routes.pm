package CoGe::Services::Routes;
use Mojo::Base "Mojolicious";

sub startup {
    my $self = shift;

    my $r = $self->routes->namespaces(["CoGe::Services::Data"]);

    # Organism routes
    $r->get("/organisms/search/#term")
        ->name("organisms-search")
        ->to("organism#search", term => undef);

    $r->get("/organisms/fetch/:id" => [id => qr/\d+/])
        ->name("organisms-fetch")
        ->to("organism#fetch", id => undef);

    $r->post("/organisms")
        ->name("organisms-add")
        ->to("organism#add");

    # Genome routes
    $r->get("/genomes/search/#term")
        ->name("genomes-search")
        ->to("genome2#search", term => undef);

    $r->get("/genomes/fetch/:id" => [id => qr/\d+/])
        ->name("genomes-fetch")
        ->to("genome2#fetch", id => undef);

    # Experiment routes
    $r->get("/experiments/search/#term")
        ->name("experiments-search")
        ->to("experiment#search", term => undef);

    $r->get("/experiments/fetch/:id" => [id => qr/\d+/])
        ->name("experiments-fetch")
        ->to("experiment#fetch", id => undef);

    # User routes
    $r->get("/users/search/#term")
        ->name("users-search")
        ->to("user#search", term => undef);

    $r->get("/users/fetch/:id" => [id => qr/\d+/])
        ->name("users-fetch")
        ->to("user#fetch", id => undef);

    $r->get("/users/fetch/:id/items" => [id => qr/\d+/])
        ->name("users-fetch-items")
        ->to("user#items", id => undef);

    # User group routes
    $r->get("/groups/search/#term")
        ->name("groups-search")
        ->to("group#search", term => undef);

    $r->get("/groups/fetch/:id" => [id => qr/\d+/])
        ->name("groups-fetch")
        ->to("group#fetch", id => undef);

    $r->get("/groups/fetch/:id/items" => [id => qr/\d+/])
        ->name("groups-fetch-items")
        ->to("group#items", id => undef);

    # Job routes
}

1;
