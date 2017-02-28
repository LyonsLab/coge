package CoGe::Services::API::Genome;

use Mojo::Base 'Mojolicious::Controller';
use Mojo::JSON qw(decode_json);
use Data::Dumper;

use CoGe::Services::Auth qw(init);
use CoGe::Services::API::Job;
use CoGe::Core::Genome qw(genomecmp search_genomes);
use CoGe::Core::Storage qw(get_genome_seq);
use CoGe::Core::Favorites;
use CoGe::Accessory::Utils qw(sanitize_name);
use CoGeDBI qw(get_feature_counts);

sub search {
    my $self = shift;
    my $search_term = $self->stash('term');
    my $fast = $self->param('fast');
    $fast = (defined $fast && ($fast eq '1' || $fast eq 'true')); # default to false
    my $sort = $self->param('sort');
    $sort = (defined $sort && ($sort eq '1' || $sort eq 'true')); # default to false

    # Validate input
    if (!$search_term or length($search_term) < 3) {
        $self->render(status => 400, json => { error => { Error => 'Search term is shorter than 3 characters' } });
        return;
    }

    # Authenticate user and connect to the database
    my ($db, $user) = CoGe::Services::Auth::init($self);
    
    # Search notebooks and filter based on user permissions
    my $filtered = search_genomes(db => $db, search_term => $search_term, user => $user, sort => $sort);

    # Get user's favorites
    my $favorites = CoGe::Core::Favorites->new(user => $user);

    # Format response
    my @result;
    if ($fast) {
        @result = map {
          {
            id => int($_->id),
            info => $_->info,
            certified  => $_->certified ? Mojo::JSON->true : Mojo::JSON->false,
            restricted => $_->restricted ? Mojo::JSON->true : Mojo::JSON->false,
            favorited  => $favorites->is_favorite($_) ? Mojo::JSON->true : Mojo::JSON->false
          }
        } @$filtered;
    }
    else {
        @result = map {
          {
            id => int($_->id),
            name => $_->name,
            description => $_->description,
            link => $_->link,
            version => $_->version,
            info => $_->info,
            organism_id  => int($_->organism->id),
            sequence_type => {
                id => ($_->type ? $_->type->id : 0),
                name => ($_->type ? $_->type->name : ''),
                description => ($_->type ? $_->type->description : ''),
            },
            restricted => $_->restricted ? Mojo::JSON->true : Mojo::JSON->false,
            certified => $_->certified ? Mojo::JSON->true : Mojo::JSON->false,
            chromosome_count => int($_->chromosome_count),
            organism => {
                id => int($_->organism->id),
                name => $_->organism->name,
                description => $_->organism->description
            }
          }
        } @$filtered;
    }
    
    $self->render(json => { genomes => \@result });
}

sub fetch {
    my $self = shift;
    my $id = int($self->stash('id'));
    
    # Validate input
    unless ($id) {
        $self->render(status => 400, json => {
            error => { Error => "Invalid input"}
        });
        return;
    }

    # Authenticate user and connect to the database
    my ($db, $user) = CoGe::Services::Auth::init($self);

    # Get genome
    my $genome = $db->resultset("Genome")->find($id);
    unless (defined $genome) {
        $self->render(status => 404, json => {
            error => { Error => "Resource not found"}
        });
        return;
    }

    # Verify permission
    unless ( !$genome->restricted || (defined $user && $user->has_access_to_genome($genome)) ) {
        $self->render(json => {
            error => { Auth => "Access denied"}
        }, status => 401);
        return;
    }

    # Format metadata
    my @metadata = map {
        {
            text => $_->annotation,
            link => $_->link,
            type => $_->type->name,
            type_group => $_->type->group
        }
    } $genome->annotations;
    
    # Build chromosome list
    my $chromosomes = $genome->chromosomes_all;
    my $feature_counts = get_feature_counts($db->storage->dbh, $genome->id);
    foreach (@$chromosomes) {
        my $name = $_->{name};
        $_->{gene_count} = int($feature_counts->{$name}{1}{count} // 0);
        $_->{CDS_count} = int($feature_counts->{$name}{3}{count} // 0);
    }
    
    # Generate response
    $self->render(json => {
        id => int($genome->id),
        name => $genome->name,
        description => $genome->description,
        link => $genome->link,
        version => $genome->version,
        restricted => $genome->restricted ? Mojo::JSON->true : Mojo::JSON->false,
        certified => $genome->certified ? Mojo::JSON->true : Mojo::JSON->false,
        organism => {
            id => int($genome->organism->id),
            name => $genome->organism->name,
            description => $genome->organism->description
        },
        sequence_type => {
            name => $genome->type->name,
            description => $genome->type->description,
        },
        chromosome_count => int($genome->chromosome_count),
        chromosomes => $chromosomes,
        experiments => [ map { int($_->id) } $genome->experiments ],
        additional_metadata => \@metadata
    });
}

sub sequence {
    my $self   = shift;
    my $gid    = $self->stash('id');
    return unless $gid;
    my $chr    = $self->stash('chr');    # optional
    my $start  = $self->param('start');  # optional
    my $stop   = $self->param('stop') || $self->param('end'); # optional
    my $strand = $self->param('strand'); # optional
    print STDERR "Data::Genome::fetch_sequence gid=$gid ",
        (defined $chr ? "chr=$chr " : ''),
        (defined $start ? "start=$start " : ''),
        (defined $stop ? "stop=$stop " : ''), "\n";

    # Connect to the database
    my ($db, $user, $conf) = CoGe::Services::Auth::init($self);

    # Retrieve genome
    my $genome = $db->resultset('Genome')->find($gid);
    unless ($genome) {
        print STDERR "Data::Sequence::get genome $gid not found in db\n";
        return;
    }

    # Check permissions
    if ( $genome->restricted
        and ( not defined $user or not $user->has_access_to_genome($genome) ) )
    {
        print STDERR "Data::Sequence::get access denied to genome $gid\n";
        return;
    }

    # Force browser to download whole genome as attachment
    my $format;
    my $genome_name = sanitize_name($genome->organism->name);
    $genome_name = 'genome_'.$gid unless $genome_name;
    if ( (!defined($chr) || $chr eq '') ) {
        $self->res->headers->content_disposition("attachment; filename=$genome_name.faa;");
    }
    elsif (defined($chr) && !defined($start) && !defined($stop)) {
        $genome_name .= '_' . $chr;
        $stop = $genome->get_chromosome_length($chr);
        $format = 'fasta';
        $self->res->headers->content_disposition("attachment; filename=$genome_name.faa;");
    }

    # Get sequence from file
    $self->render(text => get_genome_seq(
        gid   => $gid,
        chr   => $chr,
        start => $start,
        stop  => $stop,
        strand => $strand,
        format => $format
    ));
}

sub export {
    my $self = shift;
    my $gid  = $self->stash('id');
    return unless $gid;
    my $data = $self->req->json;
    $data->{gid} = $gid;

    # Alias to Job Submit -- is there a better way to do this using Mojolicious routing?
    my $request = {
        type => 'export_genome',
        parameters => $data
    };
    
    return CoGe::Services::API::Job::add($self, $request);
}

sub add {
    my $self = shift;
    my $data = $self->req->body; #$self->req->json; # mdb replaced 11/30/16 -- req->json hides JSON errors, doing conversion manually prints them to STDERR
    unless ($data) {
        $self->render(status => 400, json => {
            error => { Error => "No request body specified" }
        });
        return;
    }
    $data = decode_json($data);
    #print STDERR "CoGe::Services::Data::Genome::add\n", Dumper $data, "\n";

    # Valid data items # TODO move into request validation
    unless ($data->{source_data} && @{$data->{source_data}}) {
        $self->render(status => 400, json => {
            error => { Error => "No data items specified" }
        });
        return;
    }
    
    # Alias to Job Submit -- is there a better way to do this using Mojolicious routing?
    my $request = {
        type => 'load_genome',
        parameters => $data
    };
    
    return CoGe::Services::API::Job::add($self, $request);
}

1;
