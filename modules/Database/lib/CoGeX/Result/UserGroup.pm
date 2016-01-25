package CoGeX::Result::UserGroup;

use strict;
use warnings;
use Data::Dumper;

use base 'DBIx::Class::Core';

=head1 NAME

CoGeX::UserGroup

=head1 SYNOPSIS

This object uses the DBIx::Class to define an interface to the C<user_group> table in the CoGe database.

=head1 DESCRIPTION

=head1 AUTHORS

 Eric Lyons
 Brent Pedersen

=head1 COPYRIGHT

The full text of the license can be found in the
LICENSE file included with this module.

=head1 SEE ALSO

=cut

my $node_types = CoGeX::node_types();

__PACKAGE__->table("user_group");
__PACKAGE__->add_columns(
    "user_group_id",
    {
        data_type     => "INT",
        default_value => undef,
        is_nullable   => 0,
        size          => 11
    },
    "creator_user_id",
    {
        data_type     => "INT",
        default_value => undef,
        is_nullable   => 0,
        size          => 11
    },
    "date",
    { data_type => "TIMESTAMP", default_value => undef, is_nullable => 0 },
    "name",
    {
        data_type     => "VARCHAR",
        default_value => "",
        is_nullable   => 0,
        size          => 250
    },
    "description",
    {
        data_type     => "TEXT",
        default_value => undef,
        is_nullable   => 1,
        size          => 255
    },
    "role_id",
    {
        data_type     => "INT",
        default_value => undef,
        is_nullable   => 0,
        size          => 11
    },
    "locked",
    { data_type => "BOOLEAN", default_value => 0, is_nullable => 0, size => 1 },
    "deleted",
    { data_type => "BOOLEAN", default_value => 0, is_nullable => 0, size => 1 }
);
__PACKAGE__->set_primary_key("user_group_id");
__PACKAGE__->belongs_to( 'role' => "CoGeX::Result::Role", 'role_id' );
__PACKAGE__->belongs_to(
    "creator" => "CoGeX::Result::User", 
    { 'foreign.user_id' => 'self.creator_user_id' }
);
__PACKAGE__->has_many(    # parent users
    'user_connectors' => 'CoGeX::Result::UserConnector',
    { "foreign.child_id" => "self.user_group_id" },
    {
        where => [
            -and => [
                parent_type => $node_types->{user},
                child_type  => $node_types->{group}
            ]
        ]
    }
);
__PACKAGE__->has_many(    # all children (genomes/experiments/lists)
    'child_connectors' => 'CoGeX::Result::UserConnector',
    { "foreign.parent_id" => "self.user_group_id" },
    { where               => { parent_type => $node_types->{group} } }
);
__PACKAGE__->has_many(    # child genomes
    'genome_connectors' => "CoGeX::Result::UserConnector",
    { "foreign.parent_id" => "self.user_group_id" },
    {
        where => [
            -and => [
                parent_type => $node_types->{group},
                child_type  => $node_types->{genome}
            ]
        ]
    }
);
__PACKAGE__->has_many(    # child experiments
    'experiment_connectors' => "CoGeX::Result::UserConnector",
    { "foreign.parent_id" => "self.user_group_id" },
    {
        where => [
            -and => [
                parent_type => $node_types->{group},
                child_type  => $node_types->{experiment}
            ]
        ]
    }
);
__PACKAGE__->has_many(    # child features
    'feature_connectors' => "CoGeX::Result::UserConnector",
    { "foreign.parent_id" => "self.user_group_id" },
    {
        where => [
            -and => [
                parent_type => $node_types->{group},
                child_type  => $node_types->{feature}
            ]
        ]
    }
);
__PACKAGE__->has_many(    # child lists
    'list_connectors' => "CoGeX::Result::UserConnector",
    { "foreign.parent_id" => "self.user_group_id" },
    {
        where => [
            -and => [
                parent_type => $node_types->{group},
                child_type  => $node_types->{list}
            ]
        ]
    }
);

sub item_type {
    return $node_types->{group};   
}

################################################ subroutine header begin ##

=head2 users

 Usage     :
 Purpose   : Returns users objects
 Returns   : wantArray of users objects
 Argument  : None
 Throws    : None
 Comments  :

=cut

################################################## subroutine header end ##

sub users {
    my $self = shift;

    my @users;
    foreach ( $self->user_connectors ) {
        push @users, $_->parent;
    }

    return wantarray ? @users : \@users;
}

################################################ subroutine header begin ##

=head2 has_member

 Usage     :
 Purpose   : Check if group contains specified user
 Returns   : 0 or 1
 Argument  : user id or user object
 Throws    : None
 Comments  :

=cut

################################################## subroutine header end ##

sub has_member {
    my $self = shift;
    my $user = shift;
    my $uid  = $user =~ /^\d+$/ ? $user : $user->id;

    return 1 if $self->user_connectors( { parent_id => $uid } );

    return 0;
}

################################################ subroutine header begin ##

=head2 owner

 Usage     :
 Purpose   : Returns user object
 Returns   : user object
 Argument  : None
 Throws    : None
 Comments  :

=cut

################################################## subroutine header end ##

sub owner {
    my $self = shift;

    foreach ( $self->user_connectors( { role_id => 2 } ) ) {    #FIXME hardcoded
        return $_->parent;
    }
}

################################################ subroutine header begin ##

=head2 is_editable

 Usage     : is this group editable by the specified user?
 Purpose   :
 Returns   : 0 or 1
 Argument  :
 Throws    : None
 Comments  :

=cut

################################################## subroutine header end ##
sub is_editable {
    my $self = shift;
    my $user = shift;

    return (
             $user->is_admin
          or $user->is_owner_editor( group => $self )
          or $user->id == $self->creator_user_id
          or $self->is_editor
    );
}

################################################ subroutine header begin ##

=head2 is_<ROLE>

 Usage     : does this group have the specified role?
 Purpose   :
 Returns   : 0 or 1
 Argument  :
 Throws    : None
 Comments  :

=cut

################################################## subroutine header end ##

sub is_owner {
    my $self = shift;
    return $self->role->name =~ /owner/i;
}

sub is_editor {
    my $self = shift;
    return $self->role->name =~ /editor/i;
}

################################################ subroutine header begin ##

=head2 lists

 Usage     :
 Purpose   : Returns the set of lists associated with the user group
 Returns   : wantArray of lists
 Argument  : None
 Throws    : None
 Comments  :

=cut

################################################## subroutine header end ##

sub lists {
    my $self = shift;

    my @lists;
    foreach ( $self->list_connectors ) {
        push @lists, $_->child;
    }

    return wantarray ? @lists : \@lists;
}

################################################ subroutine header begin ##

=head2 experiments

 Usage     :
 Purpose   : Returns the set of experiments associated with the user group
 Returns   : wantArray of experiments
 Argument  : None
 Throws    : None
 Comments  :

=cut

################################################## subroutine header end ##

sub experiments {
    my $self = shift;

    my @experiments;
    foreach ( $self->experiment_connectors ) {
        push @experiments, $_->child;
    }
    return wantarray ? @experiments : \@experiments;
}

################################################ subroutine header begin ##

=head2 features

 Usage     :
 Purpose   : Returns the set of features associate with the user group
 Returns   : wantArray of features
 Argument  : None
 Throws    : None
 Comments  :

=cut

################################################## subroutine header end ##

sub features {
    my $self = shift;

    my @features;
    foreach ( $self->feature_connectors ) {
        push @features, $_->child;
    }
    return wantarray ? @features : \@features;
}

################################################ subroutine header begin ##

=head2 genomes

 Usage     :
 Purpose   : Returns the set of genomes associated with the user group
 Returns   : wantArray of genomes
 Argument  : None
 Throws    : None
 Comments  :

=cut

################################################## subroutine header end ##

sub genomes {
    my $self = shift;

    my @genomes;
    foreach ( $self->genome_connectors ) {
        push @genomes, $_->child;
    }
    return wantarray ? @genomes : \@genomes;
}

################################################ subroutine header begin ##

=head2 datasets

 Usage     :
 Purpose   : Returns the set of datasets associated with the user group
 Returns   : wantArray of datasets
 Argument  : None
 Throws    : None
 Comments  :

=cut

################################################## subroutine header end ##

sub datasets {
    my $self = shift;
    my @ds;
    foreach my $genome ( $self->genomes ) {
        push @ds, $genome->datasets;
    }
    return wantarray ? @ds : \@ds;
}

################################################ subroutine header begin ##

=head2 restricted_datasets

 Usage     : $self->restricted_datasets
 Purpose   : Returns the set of restricted datasets associated with the user group
 Returns   : wantArray of datasets
 Argument  : None
 Throws    : None
 Comments  :

=cut

################################################## subroutine header end ##

sub restricted_datasets {
    my $self = shift;
    my @ds;
    foreach my $genome ( $self->genomes ) {
        push @ds, $genome->datasets( restricted => 1 );
    }
    return wantarray ? @ds : \@ds;
}

################################################ subroutine header begin ##

=head2 info

 Usage     : $self->info
 Purpose   : generate a string of information about the user group
 Returns   : a string
 Argument  : None
 Throws    : None
 Comments  : uses name, description and role

=cut

################################################## subroutine header end ##

sub info {
    my $self = shift;
    my $info = $self->name;
    $info .= ": " . $self->description if $self->description;
    # $info .= " (" . $self->role->name . ")";
    $info .= ' (id' . $self->id . ')';
    return $info;
}

################################################ subroutine header begin ##

=head2 info

 Usage     : $self->info_html
 Purpose   : generate a string of information about the user group
 Returns   : a string wrapped with a linnk to GroupView
 Argument  : None
 Throws    : None
 Comments  : uses name, description and role

=cut

################################################## subroutine header end ##

sub info_html {
    my $self = shift;
    my $info = $self->info;
    $info =
        qq{<span class=link onclick='window.open("GroupView.pl?ugid=}
      . $self->id
      . qq{")'>}
      . $info
      . "</span>";
    return $info;
}

################################################ subroutine header begin ##

=head2 annotation_pretty_print_html

 Usage     : my $pretty_annotation_html = $feat->annotation_pretty_print_html
 Purpose   : returns a string with information and annotations about a user group
             in a nice html format with breaks and class tags (called "annotation")
 Returns   : returns a string
 Argument  : none
 Throws    :
 Comments  : uses Coge::Genome::Accessory::Annotation to build the annotations,
           : specifying delimters, and printing to string.   Pretty cool object.

See Also   : CoGe::Accessory::Annotation

=cut

################################################## subroutine header end ##

sub annotation_pretty_print_html {
    my $self         = shift;
    my %opts         = @_;
    my $minimal      = $opts{minimal};
    my $allow_delete = $opts{allow_delete};

    my $anno_obj = new CoGe::Accessory::Annotation( Type => "anno" );
    $anno_obj->Val_delimit("\n");
    $anno_obj->Add_type(0);
    $anno_obj->String_end("\n");

    my $anno_type =
      new CoGe::Accessory::Annotation(
        Type => "<tr><td nowrap='true'><span class=\"title5\">" . "Name"
          . "</span>" );
    $anno_type->Type_delimit(": <td class=\"data5\">");
    $anno_type->add_Annot( $self->name . "</td>" );
    $anno_obj->add_Annot($anno_type);

    $anno_type =
      new CoGe::Accessory::Annotation(
            Type => "<tr><td nowrap='true'><span class=\"title5\">"
          . "Description"
          . "</span>" );
    $anno_type->Type_delimit(": <td class=\"data5\">");
    $anno_type->add_Annot( $self->description . "</td>" );
    $anno_obj->add_Annot($anno_type);

    $anno_type =
      new CoGe::Accessory::Annotation(
            Type => "<tr><td valign='top' nowrap='true'><span class=\"title5\">"
          . "Role"
          . "</span>" );
    $anno_type->Type_delimit(": <td class=\"data5\">");
    $anno_type->Val_delimit("<br>");
    $anno_type->add_Annot(
        $self->role->name
          . (
            $self->role->description
            ? ' (' . $self->role->description . ')'
            : ''
          )
    );
    $anno_type->add_Annot(
"<span style='color:red;font-style:italic;'>Note: this group is locked and cannot be edited.</span>"
    ) if ( $self->locked );
    $anno_obj->add_Annot($anno_type);

    $anno_type =
      new CoGe::Accessory::Annotation(
            Type => "<tr><td valign='top' nowrap='true'><span class=\"title5\">"
          . "Users"
          . "</span>" );
    $anno_type->Type_delimit(": <td class=\"data5\">");
    $anno_type->Val_delimit("<br>");
    foreach my $user ( sort { $a->info cmp $b->info } $self->users ) {
        my $user_info = $user->info
          ; #qq{<span class="link" onclick="window.open('User.pl?uid=} . $user->id . qq{');">} . $user->info . '</span>';
        my @special_roles = ();
        if ( $user->id == $self->creator_user_id ) {
            push @special_roles, '<b>creator</b>';
        }
        foreach ( $self->user_connectors( { parent_id => $user->id } ) ) {
            push @special_roles, '<b>owner</b>' if ( $_->role->is_owner );
        }
        $anno_type->add_Annot( $user_info
              . ( @special_roles ? ' - ' . join( ', ', @special_roles ) : '' )
        );
    }
    $anno_obj->add_Annot($anno_type);

    my %genomes     = map { $_->id, $_ } $self->genomes;
    my %experiments = map { $_->id, $_ } $self->experiments;
    my %lists;
    my @lists = $self->lists;
    foreach my $list (@lists) {
        $lists{ $list->id } = $list;
        my $items = $self->_process_list(
            list      => $list,
            anno_type => $anno_type,
            anno_obj  => $anno_obj
        );
        foreach my $exp ( @{ $items->{experiments} } ) {
            $experiments{ $exp->id } = $exp;
        }
        foreach my $genome ( @{ $items->{genomes} } ) {
            $genomes{ $genome->id } = $genome;
        }
        foreach my $item ( @{ $items->{lists} } ) {
            $lists{ $item->id } = $item;
            my $data = $self->_process_list( list => $item );
            map { $genomes{ $_->id }     = $_ } @{ $data->{genomes} };
            map { $experiments{ $_->id } = $_ } @{ $data->{experiments} };
            map { $lists{ $_->id }       = $_ } @{ $data->{lists} };
        }
    }
    if ( values %lists ) {
        $anno_type =
          new CoGe::Accessory::Annotation(
            Type => "<tr valign='top'><td nowrap='true'><span class=\"title5\">"
              . "Notebooks"
              . "</span>" );
        $anno_type->Type_delimit(": <td class=\"data5\">");
        $anno_type->Val_delimit("<br>");
        foreach my $list ( sort { $b->id <=> $a->id } values %lists ) {
            my $a =
              $list->info_html
              . ( $allow_delete
                ? "<span class='link ui-icon ui-icon-trash' onclick=\"remove_list_from_group({ugid: '"
                  . $self->id
                  . "', lid: '"
                  . $list->id
                  . "'});\"></span>"
                : '' );
            $anno_type->add_Annot($a);
        }
        $anno_obj->add_Annot($anno_type);
    }

    if ( values %genomes ) {
        $anno_type =
          new CoGe::Accessory::Annotation(
            Type => "<tr valign='top'><td nowrap='true'><span class=\"title5\">"
              . "Genomes"
              . "</span>" );
        $anno_type->Type_delimit(": <td class=\"data5\">");
        $anno_type->Val_delimit("<br>");
        foreach my $item ( sort { $b->id <=> $a->id } values %genomes ) {
            my $a = $item->info_html;
            $anno_type->add_Annot($a);
        }
        $anno_obj->add_Annot($anno_type);
    }

    if ( values %experiments ) {
        $anno_type =
          new CoGe::Accessory::Annotation(
            Type => "<tr valign='top'><td nowrap='true'><span class=\"title5\">"
              . "Experiments"
              . "</span>" );
        $anno_type->Type_delimit(": <td class=\"data5\">");
        $anno_type->Val_delimit("<br>");
        foreach my $item ( sort { $b->id <=> $a->id } values %experiments ) {
            my $a = $item->info_html;
            $anno_type->add_Annot($a);
        }
        $anno_obj->add_Annot($anno_type);
    }

    if ( $self->deleted ) {
        $anno_type =
          new CoGe::Accessory::Annotation(
            Type => "<tr><td nowrap='true'><span class=\"alert\">" . "Note"
              . "</span>" );
        $anno_type->Type_delimit(": <td class=\"alert\">");
        $anno_type->add_Annot( "This group is deleted" . "</td>" );
        $anno_obj->add_Annot($anno_type);
    }

    return
        "<table cellpadding=0 class='small'>"
      . $anno_obj->to_String
      . "</table>";
}

sub _process_list {
    my $self = shift;
    my %opts = @_;
    my $list = $opts{list};

    my %experiments = map { $_->id, $_ } $list->experiments;
    my %genomes     = map { $_->id, $_ } $list->genomes;
    my %lists       = map { $_->id, $_ } $list->lists;
    return {
        experiments => [ $list->experiments ],
        genomes     => [ $list->genomes ],
        lists       => [ $list->lists ]
    };
}

1;
