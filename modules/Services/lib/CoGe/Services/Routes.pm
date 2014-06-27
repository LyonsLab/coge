package CoGe::Services::Routes;
use Mojo::Base "Mojolicious";

sub startup {
    my $self = shift;

    my $r = $self->routes->namespaces(["CoGe::Services::Data"]);

    # Organism routes
    $r->get("/organisms/search/#term")
        ->name("organisms-search")
        ->to("organism#search", term => undef);

    $r->get("/organisms/:id" => [id => qr/\d+/])
        ->name("organisms-fetch")
        ->to("organism#fetch", id => undef);

    $r->post("/organisms")
        ->name("organisms-add")
        ->to("organism#add");

    # Genome routes
    $r->get("/genomes/search/#term")
        ->name("genomes-search")
        ->to("genome2#search", term => undef);

    $r->get("/genomes/:id" => [id => qr/\d+/])
        ->name("genomes-fetch")
        ->to("genome2#fetch", id => undef);

    # Experiment routes
    $r->get("/experiments/search/#term")
        ->name("experiments-search")
        ->to("experiment#search", term => undef);

    $r->get("/experiments/:id" => [id => qr/\d+/])
        ->name("experiments-fetch")
        ->to("experiment#fetch", id => undef);

    $r->put("/experiments")
        ->name("experiments-add")
        ->to("experiment#add");

    # Notebook routes
    $r->get("/notebooks/search/#term")
        ->name("notebooks-search")
        ->to("notebook#search", term => undef);

    $r->get("/notebooks/:id" => [id => qr/\d+/])
        ->name("notebooks-fetch")
        ->to("notebook#fetch", id => undef);

    # User routes
    # mdb 4/10/14: Removing search & fetch because they are available via the
    # iPlant Trellis API.
#    $r->get("/users/search/#term")
#        ->name("users-search")
#        ->to("user#search", term => undef);

#    $r->get("/users/:id" => [id => qr/\w+/])
#        ->name("users-fetch")
#        ->to("user#fetch", id => undef);

    $r->get("/users/:id/items" => [id => qr/\w+/])
        ->name("users-items")
        ->to("user#items", id => undef);

    # User group routes
    $r->get("/groups/search/#term")
        ->name("groups-search")
        ->to("group#search", term => undef);

    $r->get("/groups/:id" => [id => qr/\d+/])
        ->name("groups-fetch")
        ->to("group#fetch", id => undef);

    $r->get("/groups/:id/items" => [id => qr/\d+/])
        ->name("groups-items")
        ->to("group#items", id => undef);

    # Job routes
    $r->get("/jobs/:id" => [id => qr/\d+/])
        ->name("jobs-fetch")
        ->to("job#fetch", id => undef);

    $r->get("/jobs/:id/results/:name" => { id => qr/\d+/, name => qr/\w+/ })
        ->name("jobs-results")
        ->to("job#results", id => undef, name => undef);
}

1;
