package CoGe::Services::API::Notebook;
use Mojo::Base 'Mojolicious::Controller';
use Mojo::JSON qw(encode_json);
use CoGeX;
use CoGe::Services::Auth;
use CoGe::Core::Notebook;
use CoGe::Accessory::Web qw(log_history);
use Data::Dumper;

sub search {
    my $self = shift;
    my $search_term = $self->stash('term');

    # Validate input
    if (!$search_term or length($search_term) < 3) {
        $self->render(status => 400, json => { error => { Error => 'Search term is shorter than 3 characters' } });
        return;
    }

    # Authenticate user and connect to the database
    my ($db, $user) = CoGe::Services::Auth::init($self);

    # Search notebooks and filter based on user permissions
    my $aNotebooks = search_notebooks(db => $db, search_term => $search_term, user => $user);

    # Format response
    my @result = map {
      {
        id => int($_->id),
        name => $_->name,
        description => $_->description
      }
    } @$aNotebooks;

    $self->render(json => { notebooks => \@result });
}

sub fetch {
    my $self = shift;
    my $id = int($self->stash('id'));
#    print STDERR "CoGe::Services::API::Notebook fetch id=$id\n";
    
    # Validate input
    unless ($id) {
        $self->render(status => 400, json => {
            error => { Error => "Invalid input" }
        });
        return;
    }

    # Authenticate user and connect to the database
    my ($db, $user) = CoGe::Services::Auth::init($self);

    # Get notebook from DB
    my $notebook = $db->resultset("List")->find($id);
    unless (defined $notebook) {
        $self->render(status => 404, json => {
            error => { Error => "Resource not found" }
        });
        return;
    }

    # Verify that user has read access to the notebook
    unless ( !$notebook->restricted || (defined $user && $user->has_access_to_list($notebook)) ) {
        $self->render(json => {
            error => { Auth => "Access denied" }
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
    } $notebook->annotations;

    # Format items
    my @items = map {
        { type => 'genome', id => $_->id }
    } $notebook->genomes;

    push @items, map {
        { type => 'experiment', id => $_->id }
    } $notebook->experiments;

    # Format response
    $self->render(json => {
        id => int($notebook->id),
        name => $notebook->name,
        description => $notebook->description,
#        type => $notebook->type->name, # mdb removed 2/22/17 COGE-800
        restricted => $notebook->restricted ? Mojo::JSON->true : Mojo::JSON->false,
        additional_metadata => \@metadata,
        items => \@items
    });
}

sub add {
    my $self = shift;
    my $data = $self->req->json;
    print STDERR "CoGe::Services::Data::Notebook::add\n", Dumper $data, "\n";

    # Authenticate user and connect to the database
    my ($db, $user) = CoGe::Services::Auth::init($self);

    # User authentication is required to add notebook
    unless (defined $user) {
        $self->render(status => 401, json => {
            error => { Auth => "Access denied" }
        });
        return;
    }
    
    # Validate parameters
    my $metadata = $data->{metadata};
    unless ($metadata->{name}) {
        $self->render(status => 400, json => {
            error => { Invalid => "Notebook name not specified" }
        });
        return;
    }
    
    # Create the notebook
    my $notebook = create_notebook(
        db => $db,
        user => $user,
        name => $metadata->{name},
        desc => $metadata->{description} || '',
        type => $metadata->{type} || '',
        page => 'API'
    );
    unless ($notebook) {
        $self->render(json => {
            error => { Error => "Could not create notebook" }
        });
        return;
    }    

    $self->render(status => 201, json => { 
        success => Mojo::JSON->true,
        id => $notebook->id
    });
}

sub add_items {
    my $self = shift;
    my $id = int($self->stash('id'));
    
    # Validate input
    unless ($id) {
        $self->render(status => 400, json => {
            error => { Error => "Notebook ID missing"}
        });
        return;
    }
    
    # Authenticate user and connect to the database
    my ($db, $user) = CoGe::Services::Auth::init($self);
    unless ($user) {
        $self->render(status => 404, json => {
            error => { Error => "User not logged in"}
        });
        return;
    }

    # Get notebook from DB
    my $notebook = $db->resultset("List")->find($id);
    unless (defined $notebook) {
        $self->render(status => 404, json => {
            error => { Error => "Notebook $id not found: $notebook"}
        });
        return;
    }

    my $data = $self->req->json;
    my @items;
    foreach my $item (@{$data->{items}}) {
    	push @items, [ $item->{id}, $item->{type} ];
    }

	my $error = add_items_to_notebook(
		db => $db,
		user => $user,
		notebook => $notebook,
		item_list => \@items
	);
	if ($error) {
        $self->render(status => 401, json => {
            error => { Error => $error}
        });
        return;
	}

    log_history(
        db          => $db,
        user_id     => $user->id,
        page        => 'API',
    	description => 'added items ' . encode_json($data) . ' to notebook ' . $notebook->info_html,
        link        => 'NotebookView.pl?nid=' . $notebook->list_id,
        parent_id   => $notebook->list_id,
        parent_type => 1 #FIXME magic number
    );

	$self->render(json => {
		success => Mojo::JSON->true
	});
}

sub remove {
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

    # Get notebook from DB
    my $notebook = $db->resultset("List")->find($id);
    unless (defined $notebook) {
        $self->render(status => 404, json => {
            error => { Error => "Resource not found"}
        });
        return;
    }

    # Attempt to delete/undelete the notebook
    my $success;
    if ($notebook->deleted) {
        $success = undelete_notebook(
            db => $db, 
            user => $user, 
            notebook_id => $id, page => 'API'
        );
    }
    else {
        $success = delete_notebook(
            db => $db,
            user => $user,
            notebook_id => $id, page => 'API'
        );
    }
    unless ($success) {
        $self->render(status => 401, json => {
            error => { Error => "Access denied"}
        });
        return;
    }
    
    $self->render(json => { 
        success => Mojo::JSON->true,
        id => $notebook->id
    });    
}

sub remove_items {
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
    unless ($user) {
        $self->render(status => 404, json => {
            error => { Error => "User not logged in"}
        });
        return;
    }

    # Get notebook
    my $notebook = $db->resultset("List")->find($id);
    unless (defined $notebook) {
        $self->render(status => 404, json => {
            error => { Error => "Resource not found" }
        });
        return;
    }

    # Check permissions
    unless ($user->is_admin || $user->is_owner_editor(list => $id)) {
        $self->render(json => {
            error => { Auth => "Access denied" }
        }, status => 401);
        return;
    }

    my $data = $self->req->json;
    my @items;
    foreach my $item (@{$data->{items}}) {
    	push @items, [ $item->{id}, $item->{type} ];
    }
    my $error = remove_items_from_notebook(
    	db => $db,
    	user => $user,
    	notebook => $notebook,
    	item_list => \@items
    );
	if ($error) {
        $self->render(status => 401, json => {
            error => { Error => $error}
        });
        return;
	}
    
	$self->render(json => {
		success => Mojo::JSON->true
	});	
}

sub update {
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

    # Get notebook
    my $notebook = $db->resultset("List")->find($id);
    unless (defined $notebook) {
        $self->render(status => 404, json => {
            error => { Error => "Resource not found" }
        });
        return;
    }

    # Check permissions
    unless ($user->is_admin || $user->is_owner_editor(list => $id)) {
        $self->render(json => {
            error => { Auth => "Access denied" }
        }, status => 401);
        return;
    }

    my $data = $self->req->json;
    if (exists($data->{metadata}->{id})) {
	    delete $data->{metadata}->{id};
    }
	$notebook->update($data->{metadata});
	$self->render(json => {
		success => Mojo::JSON->true
	});
}

1;
