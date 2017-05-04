package CoGe::Core::Notebook;
use v5.14;
use strict;
use warnings;
use Data::Dumper qw(Dumper);
use Switch;

use CoGeX;

our %ITEM_TYPE;

BEGIN {
    our (@ISA, $VERSION, @EXPORT);
    require Exporter;

    $VERSION = 0.0.1;
    @ISA = qw(Exporter);
    @EXPORT = qw( create_notebook search_notebooks add_items_to_notebook remove_items_from_notebook
                  get_notebook notebookcmp %ITEM_TYPE 
                  delete_notebook undelete_notebook 
              );

    my $node_types = CoGeX::node_types();

    %ITEM_TYPE = (    # content/toc types   # FIXME doesn't belong here, and also view constants and database constants
        all                 => 100,
        mine                => 101,
        shared              => 102,
        activity            => 103,
        trash               => 104,
        activity_viz        => 105,
        activity_analyses   => 106,
        user                => $node_types->{user},
        group               => $node_types->{group},
        notebook            => $node_types->{list},
        genome              => $node_types->{genome},
        experiment          => $node_types->{experiment}
    );
}

sub get_notebook {
    my %opts = @_;
    my $db = $opts{db};
    my $id = $opts{id};
    my $user = $opts{user}; # optional database user object
    return unless $db and $id;

	my $notebook = $db->resultset('List')->find($id);
	unless ($notebook) {
    	print STDERR "error reading notebook from db in CoGe::Core::Notebook::get_notebook\n";
		return undef;
	}

	if ($user) { # check permissions if user specified
		if ($notebook->restricted && !$user->has_access_to_notebook($notebook)) {
	    	print STDERR "attempt to load notebook without user permissions in CoGe::Core::Notebook::get_notebook\n";
			return undef;
		}
	}

	return $notebook;	
}

sub search_notebooks {
    my %opts = @_;
    my $db = $opts{db};
    my $search_term = $opts{search_term}; # id value, or keyword in name/description
    my $user = $opts{user}; # optional database user object
    return unless $db and $search_term;
    my $include_deleted = $opts{include_deleted}; # optional boolean flag

    # Search genomes
    my $search_term2 = '%' . $search_term . '%';
    my @notebooks = $db->resultset("List")->search(
        \[
            'list_id = ? OR name LIKE ? OR description LIKE ?',
            [ 'list_id',     $search_term  ],
            [ 'name',        $search_term2 ],
            [ 'description', $search_term2 ]
        ]
    );

    # Filter result by permissions
    my @filtered = grep {
        (!$_->deleted || $include_deleted) &&
        (!$_->restricted || (defined $user && $user->has_access_to_notebook($_)))
    } @notebooks;

    return \@filtered;
}

sub create_notebook {
    my %opts = @_;
    my $db      = $opts{db}; #FIXME use add_to_* functions to create new connectors and remove this param
    my $user    = $opts{user};
    my $name    = $opts{name};
    my $desc    = $opts{desc};
       $desc    = $opts{description} unless $desc;
    my $type    = $opts{type};    # string type name (or use type_id)
    my $type_id = $opts{type_id}; # type id (or use type)
    my $page    = $opts{page};
    return unless ($name and ($type or $type_id) and $db and $user);
    my $item_list = $opts{item_list}; # optional
    return if ( $user->is_public );
    
    # Convert type string to id #FIXME move this into function in DBIX module ...?
    unless ($type_id) {
        $type_id = 5; #FIXME hardcoded to "mixed"
        switch($type) {
            case 'genome'     { $type_id = 1; }
            case 'experiment' { $type_id = 2; }
            case 'owner'      { $type_id = 3; } # deprecated
            case 'feature'    { $type_id = 4; }
            case 'mixed'      { $type_id = 5; }
            case 'other'      { $type_id = 6; }
        }
    }

    # Create the new list
    my $notebook = $db->resultset('List')->create({
        name         => $name,
        description  => $desc,
        #list_type_id => $type_id, # mdb removed 12/14/16 COGE-800
        creator_id   => $user->id,
        restricted   => 1
    });
    return unless $notebook;

    # Set user as owner
    my $conn = $db->resultset('UserConnector')->create({
        parent_id   => $user->id,
        parent_type => 5,           #FIXME hardcoded to "user"
        child_id    => $notebook->id,
        child_type  => 1,           #FIXME hardcoded to "list"
        role_id     => 2            #FIXME hardcoded to "owner"
    });
    return unless $conn;

    # Add selected items to new notebook
    my $rc = add_items_to_notebook( user => $user, db => $db, notebook => $notebook, item_list => $item_list) if ($item_list);
    return if $rc;

    # Record in log
    if ($page) {
        CoGe::Accessory::Web::log_history(
            db          => $db,
            user_id     => $user->id,
            page        => "$page",
            description => 'created notebook ' . $notebook->info,
            parent_id   => $notebook->id,
            parent_type => 1 #FIXME magic number
        );
    }

    return $notebook;
}

# this returns a string on error or 0 otherwise
sub add_items_to_notebook {
    my %opts = @_;
    my $db        = $opts{db}; #FIXME use add_to_* functions to create new connectors and remove this param
    my $user      = $opts{user};      # user object
    my $notebook  = $opts{notebook};  # notebook object
    my $item_list = $opts{item_list}; # ref to array of [item_id, item_type(string)]
    my $items     = $opts{items};     # ref to array of genome/experiment DBIX objects
    return "Missing arguments: db, notebook, user" unless ($db and $notebook and $user);

    # Check permissions
    return 'Access denied' unless ( $user->admin || $user->is_owner_editor(list => $notebook) );

    # Convert item list to item array
    if ($item_list) {
        $items = [] unless $items;
        foreach (@$item_list) {
            my ( $item_id, $item_type ) = @$_;
            next unless ($item_id and $item_type);
            next if ($item_type eq 'notebook'); # adding notebooks to a notebook is not allowed

            my $item;
            if ($item_type eq 'experiment') {
                $item = $db->resultset('Experiment')->find($item_id);
            }
            elsif ($item_type eq 'genome') {
                $item = $db->resultset('Genome')->find($item_id);
            }
            return 'Invalid item: id=' . $item_id . ' type=' . $item_type unless $item;

            push @$items, $item;
        }
    }

    if ($items) {
        foreach my $item (@$items) {
            #TODO add permission check for each item
            my $conn = $db->resultset('ListConnector')->find_or_create({
                parent_id   => $notebook->id,
                child_id    => $item->id,
                child_type  => $item->item_type # numeric type ID
            });
            return 'Error adding ' . $item->type . ' ' . $item->id . ' to notebook' unless $conn;
        }
    }

    return 0;
}

sub remove_items_from_notebook {
    my %opts = @_;
    my $db        = $opts{db};
    my $user      = $opts{user};      # user object
    my $notebook  = $opts{notebook};  # notebook object
    my $item_list = $opts{item_list}; # ref to array refs of item_id, item_type
    my $items     = $opts{items};     # ref to array of genomes and/or experiments
    return "Missing arguments: $db and $notebook and $user" unless ($db and $notebook and $user);
    
    # Check permissions
    return 'User does not have permission to remove items from this notebook' unless ($user->admin || $user->is_owner_editor(list => $notebook));

    # Create connections for each item
    if ($item_list) {
        foreach (@$item_list) {
            my ( $item_id, $item_type ) = @$_;
            return 'Item id or type not specified' unless ( $item_id and $item_type );
            $item_type = $ITEM_TYPE{$item_type} if ($item_type eq 'genome' or $item_type eq 'experiment');
    
            #TODO check access permission on each item
    
            $db->resultset('ListConnector')->search({
                    parent_id   => $notebook->id,
                    child_id    => $item_id,
                    child_type  => $item_type
            })->delete;
        }
    }
    if ($items) {
        foreach my $item (@$items) {
            #TODO check access permission on each item
    
            $db->resultset('ListConnector')->search({
                    parent_id   => $notebook->id,
                    child_id    => $item->id,
                    child_type  => $item->item_type
            })->delete;
        }
    }

    return 0;
}

sub delete_notebook {
    my %opts = @_;
    my $db       = $opts{db}; #FIXME use add_to_* functions to create new connectors and remove this param
    my $user     = $opts{user};     # user object
    my $notebook_id = $opts{notebook_id}; # notebook id
    my $page     = $opts{page};
    return unless ($db and $notebook_id and $user);
    #print STDERR "delete_notebook\n";
    
    # Get notebook from DB
    my $notebook = $db->resultset('List')->find($notebook_id);
    return 0 unless $notebook;
    
    # Verify that user has editor access to the notebook
    if (!$user->is_admin && ($notebook->locked || !$user->is_owner( list => $notebook_id )) ) {
        return 0;
    }
    
    # Delete notebook
    $notebook->deleted(1);
    $notebook->update;
    
    # Record in log
    if ($page) {
        CoGe::Accessory::Web::log_history(
            db          => $db,
            user_id     => $user->id,
            page        => "$page",
            description => 'deleted notebook "' . $notebook->info . '"',
            parent_id   => $notebook->id,
            parent_type => 1 #FIXME magic number
        );
    }
    
    return 1;
}

sub undelete_notebook {
    my %opts = @_;
    my $db       = $opts{db}; #FIXME use add_to_* functions to create new connectors and remove this param
    my $user     = $opts{user};     # user object
    my $notebook_id = $opts{notebook_id}; # notebook id
    my $page     = $opts{page};
    return unless ($db and $notebook_id and $user);
    #print STDERR "undelete_notebook\n";
    
    # Get notebook from DB
    my $notebook = $db->resultset('List')->find($notebook_id);
    return 0 unless $notebook;
    
    # Verify that user has editor access to the notebook
    if (!$user->is_admin && ($notebook->locked || !$user->is_owner( list => $notebook_id )) ) {
        return 0;
    }
    
    # Undelete notebook
    $notebook->deleted(0);
    $notebook->update;
    
    # Record in log
    if ($page) {
        CoGe::Accessory::Web::log_history(
            db          => $db,
            user_id     => $user->id,
            page        => "$page",
            description => 'undeleted notebook "' . $notebook->info . '"',
            parent_id   => $notebook->id,
            parent_type => 1 #FIXME magic number
        );
    }
    
    return 1;
}

sub notebookcmp($$) {
    my ($a, $b) = $_;

    my $namea = "";
    my $nameb = "";

    $namea = $a->name if defined $a and $a->name;
    $nameb = $b->name if defined $b and $b->name;

    $namea cmp $nameb;
}

1;
