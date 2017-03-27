package CoGeX::Result::List;

use strict;
use warnings;
use base 'DBIx::Class::Core';
use CoGe::Accessory::Annotation;
use CoGe::Accessory::BisQue qw(set_bisque_visiblity);

=head1 NAME

CoGeX::List

=head1 SYNOPSIS

This object uses the DBIx::Class to define an interface to the C<list> table in the CoGe database.

=head1 DESCRIPTION

=head1 AUTHORS

 Eric Lyons

=head1 COPYRIGHT

The full text of the license can be found in the
LICENSE file included with this module.

=head1 SEE ALSO

=cut

my $node_types = CoGeX::node_types();

__PACKAGE__->table("list");
__PACKAGE__->add_columns(
    "list_id",
    {
        data_type     => "INT",
        default_value => undef,
        is_nullable   => 0,
        size          => 11
    },
    "name",
    {
        data_type     => "VARCHAR",
        default_value => "",
        is_nullable   => 0,
        size          => 255
    },
    "description",
    {
        data_type     => "VARCHAR",
        default_value => undef,
        is_nullable   => 1,
        size          => 1024
    },
# mdb removed 12/14/16 COGE-800
#    "list_type_id",
#    {
#        data_type     => "INT",
#        default_value => undef,
#        is_nullable   => 1,
#        size          => 11
#    },
    "restricted",
    { data_type => "BOOLEAN", default_value => 0, is_nullable => 0, size => 1 },
    "locked",
    { data_type => "BOOLEAN", default_value => 0, is_nullable => 0, size => 1 },
    "deleted",
    { data_type => "BOOLEAN", default_value => 0, is_nullable => 0, size => 1 },
    "creator_id",
    { data_type => "INT", default_value => 0, is_nullable => 0, size => 11 },
    "date",
    { data_type => "TIMESTAMP", default_value => undef, is_nullable => 0 },
);
__PACKAGE__->set_primary_key("list_id");
# mdb removed 12/14/16 COGE-800
#__PACKAGE__->belongs_to(
#    "list_type" => "CoGeX::Result::ListType",
#    'list_type_id'
#);
__PACKAGE__->belongs_to(
    "creator" => "CoGeX::Result::User", 
    { 'foreign.user_id' => 'self.creator_id' }
);
__PACKAGE__->has_many(
    "list_annotations" => "CoGeX::Result::ListAnnotation",
    'list_id'
);
__PACKAGE__->has_many(    # parent users
    "user_connectors" => "CoGeX::Result::UserConnector",
    { 'foreign.child_id' => 'self.list_id' },
    {
        where => [
            -and => [
                parent_type => $node_types->{user},
                child_type  => $node_types->{list}
            ]
        ]
    }
);
__PACKAGE__->has_many(    # parent groups
    "group_connectors" => "CoGeX::Result::UserConnector",
    { 'foreign.child_id' => 'self.list_id' },
    {
        where => [
            -and => [
                parent_type => $node_types->{group},
                child_type  => $node_types->{list}
            ]
        ]
    }
);
__PACKAGE__->has_many(    # all children (genomes/experiments/features/lists)
    'child_connectors' => 'CoGeX::Result::ListConnector',
    { "foreign.parent_id" => "self.list_id" }
);
__PACKAGE__->has_many(    # child genomes
    'genome_connectors' => "CoGeX::Result::ListConnector",
    { "foreign.parent_id" => "self.list_id" },
    { where               => { child_type => $node_types->{genome} } }
);
__PACKAGE__->has_many(    # child experiments
    'experiment_connectors' => "CoGeX::Result::ListConnector",
    { "foreign.parent_id" => "self.list_id" },
    { where               => { child_type => $node_types->{experiment} } }
);
__PACKAGE__->has_many(    # child features
    'feature_connectors' => "CoGeX::Result::ListConnector",
    { "foreign.parent_id" => "self.list_id" },
    { where               => { child_type => $node_types->{feature} } }
);
__PACKAGE__->has_many(    # child lists
    'list_connectors' => "CoGeX::Result::ListConnector",
    { "foreign.parent_id" => "self.list_id" },
    { where               => { child_type => $node_types->{list} } }
);
__PACKAGE__->has_many(
    "favorite_connectors" => "CoGeX::Result::FavoriteConnector",
    { "foreign.child_id" => "self.list_id" },
    { where => [ -and => [ child_type  => $node_types->{list} ] ] }
);

sub item_type {
    return $node_types->{list};   
}

sub owner {
    my $self = shift;

    foreach ( $self->user_connectors( { role_id => 2 } ) ) { #FIXME hardcoded role ID
        return $_->parent;
    }
}

sub lists # return child lists within this list -- DEPRECATED
{
    my $self       = shift;
    my %opts       = @_;
    my $restricted = $opts{restricted};    # limit result to restricted lists
    my $count      = $opts{count};         # return count;

    my @lists;
    foreach my $conn ( $self->list_connectors ) {
        next if ( $restricted and not $conn->child->restricted );
        push @lists, $conn->child;
    }

    if ($count) {
        return scalar @lists;
    }

    return wantarray ? @lists : \@lists;
}

sub features {
    my $self       = shift;
    my %opts       = @_;
    my $restricted = $opts{restricted};    # limit result to restricted features
    my $count      = $opts{count};         # return count;

    my @features;
    foreach my $conn ( $self->feature_connectors ) {
        next if ( $restricted and not $conn->child->restricted );
        push @features, $conn->child;
    }

    if ($count) {
        return scalar @features;
    }

    return wantarray ? @features : \@features;
}

sub genomes {
    my $self = shift;
    my %opts = @_;
    my $restricted = $opts{restricted};           # option to limit result to restricted genomes
    my $include_deleted = $opts{include_deleted}; # optional flag to include deleted genomes
    my $count = $opts{count};                     # optional flag to return count only

    my @genomes;
    foreach my $conn ( $self->genome_connectors ) {
        my $genome = $conn->child;
        next if ( $genome->deleted and not $include_deleted );
        next if ( $restricted and not $genome->restricted );
        push @genomes, $genome;
    }

    if ($count) {
        return scalar @genomes;
    }

    return wantarray ? @genomes : \@genomes;
}

sub experiments {
    my $self       = shift;
    my %opts       = @_;
    my $restricted = $opts{restricted};           # limit result to restricted experiments
    my $data_type  = $opts{data_type};            # limit result to experiments of given type
    my $include_deleted = $opts{include_deleted};
    my $count           = $opts{count};           # return count;

    my @experiments;
    foreach my $conn ( $self->experiment_connectors ) {
        my $experiment = $conn->child;
        next if ( $experiment->deleted and not $include_deleted );
        next if ( $restricted and not $experiment->restricted );
        next if ( $data_type and $experiment->data_type != $data_type);
        push @experiments, $experiment;
    }
    
    if ($count) {
        return scalar @experiments;
    }

    return wantarray ? @experiments : \@experiments;
}

sub children {
    my $self       = shift;
    my %opts       = @_;
    my $restricted = $opts{restricted};    # limit result to restricted children

    my @children;
    foreach my $conn ( $self->child_connectors ) {
        next if ( $restricted and not $conn->child->restricted );
        push @children, $conn->child;
    }
    return wantarray ? @children : \@children;
}

sub children_by_type {
    my $self       = shift;
    my %opts       = @_;
    my $restricted = $opts{restricted};    # limit result to restricted children

    my %children;
    foreach my $conn ( $self->child_connectors ) {
        next if ( $restricted and not $conn->child->restricted );
        my $type = $conn->child_type;
        push @{ $children{$type} }, $conn->child;
    }
    return \%children;
}

sub groups
{ # FIXME: mdb re-added 8/6/13, but need to incorporate into User.pm groups_with_access.
    my $self = shift;

    my @groups = ();
    foreach my $conn ( $self->group_connectors ) {
        push @groups, $conn->parent;
    }

    return wantarray ? @groups : \@groups;
}

sub users
{ # FIXME: mdb re-added 8/6/13, but need to incorporate into User.pm groups_with_access.
    my $self = shift;
    my %opts = @_;

    my %users;
    foreach ( $self->user_connectors ) {
        $users{ $_->parent_id } = $_->user;
    }
    foreach my $group ( $self->groups ) {
        map { $users{ $_->id } = $_ } $group->users;
    }

    return wantarray ? values %users : [ values %users ];
}

################################################ subroutine header begin ##

=head2 annotations

 Usage     :
 Purpose   : Alias to the list_annotations() method.
 Returns   : See list_annotations()
 Argument  : None
 Throws    :
 Comments  :
           :

See Also   : ListAnnotation()

=cut

################################################## subroutine header end ##

sub annotations {
    return shift->list_annotations();
}

################################################ subroutine header begin ##

=head2 type

 Usage     :
 Purpose   : Alias to the list_type() method.
 Returns   : See list_type()
 Argument  : None
 Throws    :
 Comments  :
           :

See Also   : ListType

=cut

################################################## subroutine header end ##

# mdb removed 12/14/16 COGE-800
#sub type {
#    return shift->list_type;
#}
#
#sub is_genome {
#    return shift->list_type->is_genome;
#}
#
#sub is_experiment {
#    return shift->list_type->is_experiment;
#}
#
#sub is_owner {
#    return shift->list_type->is_owner;
#}
#
#sub is_feature {
#    return shift->list_type->is_feature;
#}
#
#sub is_mixed {
#    return shift->list_type->is_mixed;
#}
#
#sub is_other {
#    return shift->list_type->is_other;
#}

################################################ subroutine header begin ##

=head2 is_editable

 Usage     : is this notebook editable by the specified user?
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
          || ( !$self->locked && $user->is_owner_editor( list => $self->id ) )
    );
}

sub is_deletable {
    my $self = shift;
    my $user = shift;

    return (
             $user->is_admin
          || ( !$self->locked && $user->is_owner( list => $self->id ) )
    );
}

################################################ subroutine header begin ##

=head2 info

 Usage     :
 Purpose   : provides quick information about the list
 Returns   : a string
 Argument  :
 Throws    :
 Comments  : name, description, restricted, type
           :

See Also   :

=cut

################################################## subroutine header end ##

sub info { 
    my $self = shift;
    my %opts = @_;

    my $info = $self->name;
    $info = '&#x1f512; ' . $info if $self->restricted && !$opts{hideRestrictedSymbol}; #TODO move this into view code
    $info .= ': ' . $self->description if $self->description;
    #$info .= ' (' . $self->type->name . ')' if $self->type;
    $info .= ' (id' . $self->id . ')';
    return $info;
}

################################################ subroutine header begin ##

=head2 info_html

 Usage     :
 Purpose   : provides quick information about the list wrapped with a link to LIstView
 Returns   : a string
 Argument  :
 Throws    :
 Comments  : name, description, restricted, type
           :

See Also   :

=cut

################################################## subroutine header end ##

sub info_html {  # FIXME deprecate this -- don't want view code in the model
    my $self = shift;
    my $info = $self->info;
    return
        qq{<span class=link onclick='window.open("NotebookView.pl?nid=}
      . $self->id
      . qq{")'>}
      . $info
      . "</span>";
}

################################################ subroutine header begin ##

=head2 data_summary

 Usage     :
 Purpose   : provides quick summary of the data contained in the list
 Returns   : a string
 Argument  :
 Throws    :
 Comments  :
           :

See Also   :

=cut

################################################## subroutine header end ##

sub data_summary {  # FIXME deprecate this -- don't want view code in the model
    my $self = shift;
    my @stuff;
    my $exps = $self->experiments( count => 1 );
    push @stuff, "Experiments: " . $exps if $exps;
    my $feats = $self->features( count => 1 );
    push @stuff, "Features: " . $feats if $feats;
    my $genomes = $self->genomes( count => 1 );
    push @stuff, "Genomes: " . $genomes if $genomes;
    my $lists = $self->lists( count => 1 );
    push @stuff, "Notebooks: " . $lists if $lists;
    return join( "; ", @stuff );
}

sub contents_summary_html { # FIXME deprecate this -- don't want view code in the model
    my $self = shift;

    my $html;
    $html .= 'Genomes: ' . @{ $self->genomes } . '<br>'
      if ( @{ $self->genomes } );
    $html .= 'Experiments: ' . @{ $self->experiments } . '<br>'
      if ( @{ $self->experiments } );
    $html .= 'Features: ' . @{ $self->features } . '<br>'
      if ( @{ $self->features } );
    $html .= 'Notebooks: ' . @{ $self->lists } . '<br>'
      if ( @{ $self->lists } );

    return $html;
}

################################################ subroutine header begin ##

=head2 annotation_pretty_print_html

 Usage     : my $pretty_annotation_html = $feat->annotation_pretty_print_html
 Purpose   : returns a string with information and annotations about a feature
             in a nice html format with breaks and class tags (called "annotation")
 Returns   : returns a string
 Argument  : none
 Throws    :
 Comments  : uses Coge::Genome::Accessory::Annotation to build the annotations,
           : specifying delimters, and printing to string.   Pretty cool object.

See Also   : CoGe::Accessory::Annotation

=cut

################################################## subroutine header end ##

sub annotation_pretty_print_html { # FIXME deprecate this -- don't want view code in the model
    my $self    = shift;
    my %opts    = @_;

    my $anno_obj = new CoGe::Accessory::Annotation( Type => "anno" );
    $anno_obj->Val_delimit("\n");
    $anno_obj->Add_type(0);
    $anno_obj->String_end("\n");
    my $anno_type =
      new CoGe::Accessory::Annotation(
        Type => "<tr><td nowrap='true'><span class='title5'>" . "ID" . "</span>" );
    $anno_type->Type_delimit(": <td class='data5'>");
    $anno_type->add_Annot( $self->id . '</td>' );
    $anno_obj->add_Annot($anno_type);
    $anno_type =
      new CoGe::Accessory::Annotation(
        Type => "<tr><td nowrap='true'><span class='title5'>" . "Name" . "</span>" );
    $anno_type->Type_delimit(": <td class='data5'>");
    $anno_type->add_Annot( $self->name . '</td>' );
    $anno_obj->add_Annot($anno_type);
    $anno_type =
      new CoGe::Accessory::Annotation(
            Type => "<tr valign='top'><td nowrap='true'><span class='title5'>"
          . "Description"
          . "</span>" );
    $anno_type->Type_delimit(": <td class='data5' style='max-width:400px;overflow:hidden;word-wrap:break-word;'>");
    $anno_type->add_Annot( $self->description ) if $self->description;
    $anno_type->add_Annot( "</td>" );
    $anno_obj->add_Annot($anno_type);

# mdb removed 12/14/16 COGE-800
#    $anno_type =
#      new CoGe::Accessory::Annotation(
#        Type => "<tr><td nowrap='true'><span class='title5'>" . "Type"
#          . "</span>" );
#    $anno_type->Type_delimit(": <td class='data5'>");
#    if ( $self->type ) {
#        my $type = $self->type->name;
#        $type .= ": " . $self->type->description if $self->type->description;
#        $anno_type->add_Annot( $type . "</td>" );
#    }
#    $anno_obj->add_Annot($anno_type);

    if ( $self->locked ) {
        $anno_type =
            new CoGe::Accessory::Annotation(
                Type => "<tr><td valign='top' nowrap='true'><span class='title5'>"
                    ."Note"
                    ."</span>" );
        $anno_type->Type_delimit( ": <td class='data5'>" );
        $anno_type->add_Annot(
            "<span style='color:red;font-style:italic;'>this list is locked and cannot be edited.</span>"
        );
        $anno_type->add_Annot( "</td>" );
        $anno_obj->add_Annot( $anno_type );
    }

    $anno_type =
      new CoGe::Accessory::Annotation(
            Type => "<tr><td nowrap='true'><span class='title5'>"
          . "Restricted"
          . "</span>" );
    $anno_type->Type_delimit(": <td class='data5'>");
    my $restricted = $self->restricted ? "Yes" : "No";
    $anno_type->add_Annot( $restricted . "</td>" );
    $anno_obj->add_Annot($anno_type);
    
    my $creation = ($self->creator_id ? $self->creator->display_name  . ' ' : '') . ($self->date ne '0000-00-00 00:00:00' ? $self->date : '');
    if ($creation) {
        $anno_type =
          new CoGe::Accessory::Annotation(
                Type => "<tr><td nowrap='true'><span class='title5'>"
              . "Creation"
              . "</span>" );
        $anno_type->Type_delimit(": <td class='data5'>");
        $anno_type->add_Annot( $creation . "</td>" );
        $anno_obj->add_Annot($anno_type);
    }
    
    my $owner = $self->owner;
    if ($owner) {
        $anno_type =
          new CoGe::Accessory::Annotation(
                Type => "<tr><td nowrap='true'><span class='title5'>"
              . "Owner"
              . "</span>" );
        $anno_type->Type_delimit(": <td class='data5'>");
        $anno_type->add_Annot( $owner->display_name . "</td>" );
        $anno_obj->add_Annot($anno_type);
    }

    $anno_type =
      new CoGe::Accessory::Annotation(
            Type => "<tr><td valign='top' nowrap='true'><span class='title5'>"
          . "Users with access"
          . "</span>" );
    $anno_type->Type_delimit(": <td class='data5'>");
    my $users = ( $self->restricted ? 
        join( ', ', sort map { $_->display_name } grep { defined $_ } $self->users ) : # mdb changed 1/3/17 -- added check for defined user for case where NULL (not sure how that happened)
        'Everyone' );
    $anno_type->add_Annot( $users . "</td>" );
    $anno_obj->add_Annot($anno_type);
    
    my $groups = ( $self->restricted ? 
        join( ', ', sort map { $_->name } $self->groups ) :
        undef );
    if ($groups) {
        $anno_type =
          new CoGe::Accessory::Annotation(
                Type => "<tr><td valign='top' nowrap='true'><span class='title5'>"
              . "Groups with access"
              . "</span>" );
        $anno_type->Type_delimit(": <td class='data5'>");
        $anno_type->add_Annot( $groups . "</td>" );
        $anno_obj->add_Annot($anno_type);
    }

    if ( $self->deleted ) {
        $anno_type =
          new CoGe::Accessory::Annotation(
            Type => "<tr><td nowrap='true'><span class=\"alert\">" . "Note"
              . "</span>" );
        $anno_type->Type_delimit(": <td class=\"alert\">");
        $anno_type->add_Annot( "This notebook is deleted" . "</td>" );
        $anno_obj->add_Annot($anno_type);
    }

    return
        "<table cellpadding=0 class='border-top border-bottom'>"
      . $anno_obj->to_String
      . "</table>";
}

sub to_hash {
    my $self = shift;
    return {
        id          => $self->id,
        type        => 'notebook',
        name        => $self->name,
        description => $self->description,
        restricted  => $self->restricted
    };
}

1;
