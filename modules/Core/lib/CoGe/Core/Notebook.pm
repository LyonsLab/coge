package CoGe::Core::Notebook;
use v5.14;
use strict;
use warnings;
use Data::Dumper qw(Dumper);

use CoGeX;

our %ITEM_TYPE;

BEGIN {
    our (@ISA, $VERSION, @EXPORT);
    require Exporter;

    $VERSION = 0.0.1;
    @ISA = qw(Exporter);
    @EXPORT = qw( create_notebook add_item_to_notebook %ITEM_TYPE );

    my $node_types = CoGeX::node_types();

    %ITEM_TYPE = (    # content/toc types
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

sub create_notebook {
    my %opts = @_;
    my $db = $opts{db};
    my $user = $opts{user};
    my $name    = $opts{name};
    my $desc    = $opts{desc};
    my $type_id = $opts{type_id};
    my $PAGE = $opts{PAGE};
    $PAGE //= "";

    return unless $name and $type_id and $db and $user;
    my $items = $opts{item_list}; # optional
    return if ( $user->user_name eq "public" );

    # Create the new list
    my $list = $db->resultset('List')->create(
        {
            name         => $name,
            description  => $desc,
            list_type_id => $type_id,

            # user_group_id => $owner->id,
            restricted => 1
        }
    );
    return unless $list;

    # Set user as owner
    my $conn = $db->resultset('UserConnector')->create(
        {
            parent_id   => $user->id,
            parent_type => 5,           #FIXME hardcoded to "user"
            child_id    => $list->id,
            child_type  => 1,           #FIXME hardcoded to "list"
            role_id     => 2            #FIXME hardcoded to "owner"
        }
    );
    return unless $conn;

    # Add selected items to new notebook
    add_items_to_notebook( user => $user, db => $db, nid => $list->id, item_list => $items)
      if ($items);

    # Record in log
    CoGe::Accessory::Web::log_history(
        db          => $db,
        user_id     => $user->id,
        page        => "$PAGE",
        description => 'create notebook id' . $list->id
    );

    return $list;
}

sub add_items_to_notebook {
    my %opts = @_;
    my $db = $opts{db};
    my $user = $opts{user};
    my $nid  = $opts{nid};

    return unless $nid and $db and $user;

    my $items = $opts{item_list};
    return unless $items;

    # print STDERR "add_items_to_notebook $nid $item_list\n";

    my $notebook = $db->resultset('List')->find($nid);
    return unless $user->has_access_to_list($notebook);

    foreach (@$items) {
        my ( $item_id, $item_type ) = @$_;
        say STDERR "ITEMS: $item_id $item_type";
        next unless ( $item_id and $item_type );
        next
          unless ( $item_type eq $ITEM_TYPE{notebook}
            or $item_type eq $ITEM_TYPE{genome}
            or $item_type eq $ITEM_TYPE{experiment} );

        #TODO check access permission on each item

        # print STDERR "add_item_to_notebook $item_id $item_type\n";

        my $conn = $db->resultset('ListConnector')->find_or_create(
            {
                parent_id  => $nid,
                child_id   => $item_id,
                child_type => $item_type
            }
        );
        return unless $conn;
    }

    return 1;
}

1;
