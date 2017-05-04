package CoGe::Services::API::Notebook;
use Mojo::Base 'Mojolicious::Controller';
use Mojo::JSON qw(encode_json);
use Data::Dumper;
use CoGeX;
use CoGe::Accessory::Web qw(log_history);
use CoGe::Core::Metadata qw( create_annotation delete_annotation get_annotation get_annotations );
use CoGe::Core::Notebook;
use CoGe::Services::Auth;
use CoGe::Services::Error;

sub search {
    my $self = shift;
    my $search_term = $self->stash('term');

    # Validate input
    if (!$search_term or length($search_term) < 3) {
        $self->render(API_STATUS_SEARCHTERM);
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

    my ($db, $user) = CoGe::Services::Auth::init($self);
    my $notebook = $self->_get_notebook($id, 0, $db, $user);
    return unless $notebook;

    # Format metadata
    my @metadata = map {
        {
            text => $_->annotation,
            link => $_->link,
            type => $_->type->name,
            type_group => $_->type->group
        }
    } $notebook->annotations;

    # Check permissions and format results
    my @items = map {
        { type => 'genome', id => $_->id, info => $_->info, date => ($_->date eq '0000-00-00 00:00:00' ? undef : $_->date) }
    } grep { !$_->restricted || (defined $user && $user->has_access_to_genome($_)) } $notebook->genomes;
    push @items, map {
        { type => 'experiment', id => $_->id, info => $_->info, date => ($_->date eq '0000-00-00 00:00:00' ? undef : $_->date), data_type => $_->data_type_desc }
    } grep { !$_->restricted || (defined $user && $user->has_access_to_experiment($_)) } $notebook->experiments;
    foreach my $f ($notebook->features) {
        foreach my $g ($f->genomes) {
            if  (!$g->restricted || (defined $user && $user->has_access_to_genome($g))) {
                push @items, { type => 'feature', id => $f->id, info => $f->info };
                last;
            }
        }
    }

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

sub fetch_annotations {
    my $self = shift;
    my $id = int($self->stash('id'));
    my ($db) = CoGe::Services::Auth::init($self);
    #TODO add error checking on ID param

    $self->render(json => get_annotations($id, 'List', $db, 1));
}

sub fetch_annotation {
    my $self = shift;
    my $id = int($self->stash('id'));
    my $aid = int($self->stash('aid'));
    #TODO add error checking on ID params

    my ($db, $user) = CoGe::Services::Auth::init($self);
    my $notebook = $self->_get_notebook($id, 0, $db, $user);
    unless ($notebook) {
        $self->render(API_STATUS_NOTFOUND);
        return;
    }

    my $annotation = get_annotation($aid, 'List', $db);
    unless ($annotation) {
        $self->render(API_STATUS_NOTFOUND);
        return;
    }
    $self->render(json => $annotation);
}

sub add {
    my $self = shift;
    my $data = $self->req->json;
    print STDERR "CoGe::Services::Data::Notebook::add\n", Dumper $data, "\n";

    # Authenticate user and connect to the database
    my ($db, $user) = CoGe::Services::Auth::init($self);

    # User authentication is required to add notebook
    unless (defined $user) {
        $self->render(API_STATUS_UNAUTHORIZED);
        return;
    }
    
    # Validate parameters
    my $metadata = $data->{metadata};
    unless ($metadata->{name}) {
        $self->render(API_STATUS_BAD_REQUEST("Notebook name not specified"));
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
        $self->render(API_STATUS_CUSTOM(200, "Could not create notebook"));
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
    
    my ($db, $user) = CoGe::Services::Auth::init($self);
    my $notebook = $self->_get_notebook($id, 1, $db, $user);
    unless ($notebook) {
        $self->render(API_STATUS_NOTFOUND);
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
        $self->render(API_STATUS_CUSTOM(401, $error));
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

sub add_annotation {
    my $self = shift;
    my ($db, $user, $conf) = CoGe::Services::Auth::init($self);
    #TODO add error checking on params
    create_annotation(
        conf => $conf,
        db => $db,
        filename => $self->param('filename'),
        group_name => $self->param('group_name'),
        image => $self->param('image'),
        link => $self->param('link'),
        target_id => int($self->stash('id')),
        target_type => 'notebook',
        text => $self->param('annotation'),
        type_name => $self->param('type_name'),
        user => $user
    );
    $self->render(json => { success => Mojo::JSON->true });
}

sub remove {
    my $self = shift;
    my $id = int($self->stash('id'));

    my ($db, $user) = CoGe::Services::Auth::init($self);
    my $notebook = $self->_get_notebook($id, 1, $db, $user);
    unless ($notebook) {
        $self->render(API_STATUS_NOTFOUND);
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
        $self->render(API_STATUS_UNAUTHORIZED);
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

    my ($db, $user) = CoGe::Services::Auth::init($self);
    my $notebook = $self->_get_notebook($id, 1, $db, $user);
    return unless $notebook;

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
        $self->render(API_STATUS_CUSTOM(401, $error));
        return;
	}
    
	$self->render(json => {
		success => Mojo::JSON->true
	});	
}

sub delete_annotation {
    my $self = shift;
    my $id = int($self->stash('id'));
    my $aid = int($self->stash('aid'));
    #TODO add error checking on params

    my ($db, $user) = CoGe::Services::Auth::init($self);
    my $error = CoGe::Core::Metadata::delete_annotation($aid, $id, 'List', $db, $user);
    if ($error) {
        $self->render(API_STATUS_BAD_REQUEST($error));
        return;
    }
    $self->render(json => { success => Mojo::JSON->true });
}

sub update {
	my $self = shift;
    my $id = int($self->stash('id'));

    my ($db, $user) = CoGe::Services::Auth::init($self);
    my $notebook = $self->_get_notebook($id, 1, $db, $user);
    return unless $notebook;

    my $data = $self->req->json;
    if (exists($data->{metadata}->{id})) {
	    delete $data->{metadata}->{id};
    }
	$notebook->update($data->{metadata});
	$self->render(json => { success => Mojo::JSON->true });
}

sub update_annotation {
    my $self = shift;
    my ($db, $user) = CoGe::Services::Auth::init($self);
    #TODO add error checking on params
    CoGe::Core::Metadata::update_annotation(
        annotation_id => int($self->stash('aid')),
        db => $db,
        delete_bisque_image => $self->param('delete_bisque_image'),
        filename => $self->param('filename'),
        group_name => $self->param('group_name'),
        image => $self->param('image'),
        link => $self->param('link'),
        target_id => $self->stash('id'),
        target_type => 'notebook',
        text => $self->param('annotation'),
        type_name => $self->param('type_name'),
        user => $user
    );
    $self->render(json => { success => Mojo::JSON->true });
}

sub _get_notebook {
    my ($self, $id, $own_or_edit, $db, $user) = @_;
    
    unless ($id) {
        $self->render(API_STATUS_MISSING_ID);
        return;
    }

    my $notebook = $db->resultset("List")->find($id);
    unless (defined $notebook) {
        $self->render(API_STATUS_NOTFOUND);
        return;
    }

    if ($own_or_edit) {
        unless ($user) {
            $self->render(API_STATUS_CUSTOM(404, "User not logged in"));
            return;
        }
        return $notebook if $user->is_admin;
        unless ($user->is_owner_editor(list => $id)) {
            $self->render(API_STATUS_UNAUTHORIZED);
            return;
        }
    }

    unless ( !$notebook->restricted || (defined $user && $user->has_access_to_notebook($notebook)) ) {
        $self->render(API_STATUS_UNAUTHORIZED);
        return;
    }

    return $notebook;
}

1;
