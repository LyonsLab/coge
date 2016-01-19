package CoGeX::Result::Experiment;

use strict;
use warnings;
use base 'DBIx::Class::Core';
use File::Spec::Functions;
use Data::Dumper;
use POSIX;
use Switch;

use CoGe::Accessory::Annotation;

=head1 NAME

CoGeX::Result::Experiment

=head1 SYNOPSIS

This object uses the DBIx::Class to define an interface to the C<experiment> table in the CoGe database.

=head1 DESCRIPTION

=head1 AUTHORS

 Eric Lyons
 Matt Bomhoff

=head1 COPYRIGHT

The full text of the license can be found in the
LICENSE file included with this module.

=head1 SEE ALSO

=cut

my $node_types = CoGeX::node_types();

__PACKAGE__->table("experiment");
__PACKAGE__->add_columns(
    "experiment_id",
    {
        data_type     => "INT",
        default_value => undef,
        is_nullable   => 0,
        size          => 11
    },
    "genome_id",
    {
        data_type     => "INT",
        default_value => undef,
        is_nullable   => 0,
        size          => 11
    },
    "data_source_id",
    {
        data_type     => "INT",
        default_value => undef,
        is_nullable   => 0,
        size          => 11
    },
    "data_type",
    { data_type => "INT", default_value => undef, is_nullable => 1, size => 1 },
    "name",
    {
        data_type     => "VARCHAR",
        default_value => "",
        is_nullable   => 0,
        size          => 255
    },
    "description",
    { data_type => "TEXT", default_value => undef, is_nullable => 1 },
    "version",
    {
        data_type     => "VARCHAR",
        default_value => undef,
        is_nullable   => 0,
        size          => 50
    },
# mdb removed 8/7/13 issue 77
#	"storage_path",
#	{ data_type => "VARCHAR", default_value => undef, is_nullable => 0, size => 255 },
    "restricted",
    { data_type => "INT", default_value => "0", is_nullable => 0, size => 1 },
    "row_count",
    { data_type => "INT", default_value => "0", is_nullable => 1, size => 10 },
    "date",
    { data_type => "TIMESTAMP", default_value => undef, is_nullable => 0 },
    "deleted",
    { data_type => "INT", default_value => "0", is_nullable => 0, size => 1 },
    "creator_id",
    { data_type => "INT", default_value => 0, is_nullable => 0, size => 11 },
);

__PACKAGE__->set_primary_key("experiment_id");
__PACKAGE__->belongs_to(
    "data_source" => "CoGeX::Result::DataSource",
    'data_source_id'
);
__PACKAGE__->belongs_to( "genome" => "CoGeX::Result::Genome", 'genome_id' );
__PACKAGE__->belongs_to(
    "creator" => "CoGeX::Result::User", 
    { 'foreign.user_id' => 'self.creator_id' }
);
__PACKAGE__->has_many(
    "experiment_type_connectors" => "CoGeX::Result::ExperimentTypeConnector",
    'experiment_id'
);
__PACKAGE__->has_many(
    "experiment_annotations" => "CoGeX::Result::ExperimentAnnotation",
    'experiment_id'
);
__PACKAGE__->has_many(    # parent lists
    'list_connectors' => 'CoGeX::Result::ListConnector',
    { 'foreign.child_id' => 'self.experiment_id' },
    { where => [ -and => [ child_type => $node_types->{experiment} ] ] }
);
__PACKAGE__->has_many(    # parent users
    'user_connectors' => 'CoGeX::Result::UserConnector',
    { "foreign.child_id" => "self.experiment_id" },
    {
        where => [
            -and => [
                parent_type => $node_types->{user},
                child_type  => $node_types->{experiment}
            ]
        ]
    }
);
__PACKAGE__->has_many(    # parent groups
    'group_connectors' => 'CoGeX::Result::UserConnector',
    { "foreign.child_id" => "self.experiment_id" },
    {
        where => [
            -and => [
                parent_type => $node_types->{group},
                child_type  => $node_types->{experiment}
            ]
        ]
    }
);

sub item_type {
    return $node_types->{experiment};   
}

sub desc {
    return shift->description(@_);
}

sub source {
    shift->data_source(@_);
}

sub annotations {
    shift->experiment_annotations(@_);
}

sub types {
    shift->experiment_types(@_);
}

sub tags {
    shift->experiment_types(@_);
}

sub data_type_desc {
    my $self = shift;
    
    # Experiment Data Types
    switch ($self->data_type) {
        case 1 { return 'quantitative'; } # Quantitative data
        case 2 { return 'polymorphism'; } # Polymorphism data
        case 3 { return 'alignment';    } # Alignments
        case 4 { return 'marker';       } # Markers
        else   { return 'unknown';      }
    }    
}

sub owner {
    my $self = shift;

    foreach ( $self->user_connectors( { role_id => 2 } ) ) {    #FIXME hardcoded
        return $_->parent;
    }
}

################################################ subroutine header begin ##

=head2 experiment_types

 Usage     : $self->experiment_types
 Purpose   : pass through experiment_type_connector to fake a many-to-many connection with experiment_type
 Returns   :
 Argument  : none, really
 Throws    :
 Comments  :

See Also   :

=cut

################################################## subroutine header end ##

sub experiment_types {
    map { $_->experiment_type } shift->experiment_type_connectors();
}

# mdb removed 8/7/13 issue 77, moved into Accessory::Storage
#sub get_path
#{
#	my $self          = shift;
#	my $experiment_id = $self->id;
#	my $level0        = floor( $experiment_id / 1000000000 ) % 1000;
#	my $level1        = floor( $experiment_id / 1000000 ) % 1000;
#	my $level2        = floor( $experiment_id / 1000 ) % 1000;
#	my $path          = catdir( $level0, $level1, $level2, $experiment_id );    #adding experiment_id for final directory.
#	return $path;
#}

# mdb added 8/7/13 issue 77
sub storage_path {
    my $self = shift;
    return CoGe::Core::Storage::get_experiment_path( $self->id );
}

################################################ subroutine header begin ##

=head2 lists

 Usage     : $self->lists
 Purpose   : returns lists containing with experiments
 Returns   :
 Argument  :
 Throws    :
 Comments  :

See Also   :

=cut

################################################## subroutine header end ##

sub lists {
    my $self = shift;

    my %lists;
    foreach my $conn ( $self->list_connectors ) {
        my $list = $conn->parent_list;
        $lists{ $list->id } = $list if ($list);
    }

    return wantarray ? values %lists : [ values %lists ];
}

sub notebooks {
    shift->lists(@_);
}

sub notebooks_desc { #FIXME move this view code out of here
    my $self = shift;
    my $embed = shift || 0;
    my @notebooks;
    foreach my $notebook ($self->notebooks) {
        my $link = "NotebookView.pl?nid=" . $notebook->id . "&embed=$embed";
        my $desc = "<a href='$link'>" . $notebook->name . '</a>';
        $desc =~ s/&reg;\s*//;
        push @notebooks, $desc;
    }
    return join(',', @notebooks) || '';
}

#sub groups {
#	my $self = shift;
#	my %opts = @_;
#
#	my @groups = ();
#	foreach my $conn ( $self->user_connectors )
#	{
#		if ($conn->parent_type == 6) { #FIXME hardcoded type
#			push @groups, $conn->group;
#		}
#	}
#
#	return wantarray ? @groups : \@groups;
#}
#
#sub users {
#	my $self = shift;
#	my %users;
#
#	foreach	( $self->user_connectors )
#	{
#		if ($_->is_parent_user) {
#			$users{$_->parent_id} = $_->user;
#		}
#		elsif ($_->is_parent_group) {
#			map { $users{$_->id} = $_ } $_->group->users;
#		}
#	}
#
#	return wantarray ? values %users : [ values %users ];
#}

# mdb removed 4/30/15 - deprecated, don't want view code in the model
#sub annotation_pretty_print_html {
#    my $self         = shift;
#    my %opts         = @_;
#    my $minimal      = $opts{minimal};
#    my $allow_delete = $opts{allow_delete};
#
#    my $anno_obj = new CoGe::Accessory::Annotation( Type => "anno" );
#    $anno_obj->Val_delimit("\n");
#    $anno_obj->Add_type(0);
#    $anno_obj->String_end("\n");
#
#    my $anno_type =
#      new CoGe::Accessory::Annotation(
#        Type => "<tr><td nowrap='true'><span class=\"title5\">" . "ID"
#          . "</span>" );
#    $anno_type->Type_delimit(": <td class=\"data5\">");
#    $anno_type->add_Annot( $self->id . "</td>" );
#    $anno_obj->add_Annot($anno_type);
#
#    $anno_type =
#      new CoGe::Accessory::Annotation(
#        Type => "<tr><td nowrap='true'><span class=\"title5\">" . "Name"
#          . "</span>" );
#    $anno_type->Type_delimit(": <td class=\"data5\">");
#    $anno_type->add_Annot( $self->name . "</td>" );
#    $anno_obj->add_Annot($anno_type);
#
#    $anno_type =
#      new CoGe::Accessory::Annotation(
#            Type => "<tr><td nowrap='true'><span class=\"title5\">"
#          . "Description"
#          . "</span>" );
#    $anno_type->Type_delimit(": <td class=\"data5\">");
#    $anno_type->add_Annot( $self->description . "</td>" );
#    $anno_obj->add_Annot($anno_type);
#
#    $anno_type =
#      new CoGe::Accessory::Annotation(
#            Type => "<tr><td nowrap='true'><span class=\"title5\">"
#          . "Genome"
#          . "</span>" );
#    $anno_type->Type_delimit(": <td class=\"data5\">");
#    $anno_type->add_Annot( $self->genome->info_html . "</td>" );
#    $anno_obj->add_Annot($anno_type);
#
#    $anno_type =
#      new CoGe::Accessory::Annotation(
#            Type => "<tr><td nowrap='true'><span class=\"title5\">"
#          . "Source"
#          . "</span>" );
#    $anno_type->Type_delimit(": <td class=\"data5\">");
#    $anno_type->add_Annot(
#        ( $self->source ? $self->source->info_html : '<no source>' )
#        . "</td>" );
#    $anno_obj->add_Annot($anno_type);
#
#    $anno_type =
#      new CoGe::Accessory::Annotation(
#            Type => "<tr><td nowrap='true'><span class=\"title5\">"
#          . "Version"
#          . "</span>" );
#    $anno_type->Type_delimit(": <td class=\"data5\">");
#    $anno_type->add_Annot( $self->version . "</td>" );
#    $anno_obj->add_Annot($anno_type);
#
#    $anno_type =
#      new CoGe::Accessory::Annotation(
#            Type => "<tr valign='top'><td nowrap='true'><span class=\"title5\">"
#          . "Types"
#          . "</span>" );
#    $anno_type->Type_delimit(": <td class=\"data5\">");
#    $anno_type->Val_delimit("<br>");
#    foreach my $type ( $self->types ) {
#        my $info = $type->name;
#        $info .= ": " . $type->description if $type->description;
#        if ($allow_delete) {
#            # NOTE: it is undesirable to have a javascript call in a DB object, but it works
#            $info .=
#                "<span onClick=\"remove_experiment_type({eid: '"
#              . $self->id
#              . "', etid: '"
#              . $type->id
#              . "'});\" class=\"link ui-icon ui-icon-trash\"></span>";
#        }
#        $anno_type->add_Annot($info);
#    }
#    $anno_obj->add_Annot($anno_type);
#
#    $anno_type =
#      new CoGe::Accessory::Annotation(
#            Type => "<tr valign='top'><td nowrap='true'><span class=\"title5\">"
#          . "Notebooks"
#          . "</span>" );
#    $anno_type->Type_delimit(": <td class=\"data5\">");
#    $anno_type->Val_delimit("<br>");
#    foreach my $list ( $self->lists ) {
#        next if ( not $list );    #FIXME mdb 8/21/12 - why necessary?
#        my $info = $list->info_html;
#        $anno_type->add_Annot($info);
#    }
#    $anno_obj->add_Annot($anno_type);
#
#    $anno_type =
#      new CoGe::Accessory::Annotation(
#            Type => "<tr><td nowrap='true'><span class=\"title5\">"
#          . "Restricted"
#          . "</span>" );
#    $anno_type->Type_delimit(": <td class=\"data5\">");
#    my $restricted = $self->restricted ? "Yes" : "No";
#    $anno_type->add_Annot( $restricted . "</td>" );
#    $anno_obj->add_Annot($anno_type);
#
#    if ( $self->deleted ) {
#        $anno_type =
#          new CoGe::Accessory::Annotation(
#            Type => "<tr><td nowrap='true'><span class=\"alert\">" . "Note"
#              . "</span>" );
#        $anno_type->Type_delimit(": <td class=\"alert\">");
#        $anno_type->add_Annot( "This experiment is deleted" . "</td>" );
#        $anno_obj->add_Annot($anno_type);
#    }
#
#    return
#        "<table cellpadding=0 class='ui-widget-content ui-corner-all'>"
#      . $anno_obj->to_String
#      . "</table>";
#}

################################################ subroutine header begin ##

=head2 info

 Usage     : $self->info
 Purpose   : returns a string of information about the experiment.

 Returns   : returns a string
 Argument  : none
 Throws    :
 Comments  : To be used to quickly generate a string about the experiment

See Also   :

=cut

################################################## subroutine header end ##

sub info {
    my $self = shift;
    my %opts = @_;

    my $source;
    if ( $opts{source} ) {    # added for performance
        $source = $opts{source};
    }
    elsif ( $self->source ) {
        $source = $self->source->name;
    }
    else {
        $source = '<no source>';
    }

    my $info;
    $info .= "&reg; "                  if $self->restricted;
    $info .= $self->name;
    $info .= ": " . $self->description if $self->description;
    $info .= " (v" . $self->version . ", id" . $self->id . "): " . $source;
    return $info;
}

############################################### subroutine header begin ##

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

sub info_html {
    my $self = shift;
    my $info = $self->info;
    return
        qq{<span class=link onclick='window.open("ExperimentView.pl?eid=}
      . $self->id
      . qq{")'>}
      . $info
      . "</span>";
}

sub info_file {
    my $self = shift;
    
    my $restricted = ($self->restricted) ? "yes" : "no";
    my $types = join ",", map($_->name, $self->types);
    my $genome_name = $self->genome->info;
    $genome_name =~ s/&reg;\s*//;

    my @lines = (
        qq{"Name","} . $self->name . '"',
        qq{"Description","} . $self->description . '"',
        qq{"Genome","$genome_name"},
        qq{"Source","} . $self->source->info . '"',
        qq{"Version","} . $self->version . '"',
        qq{"Types","$types"},
        qq{"Notebooks","} . $self->notebooks_desc . '"',
        qq{"Restricted","$restricted"}
    );

    return join("\n", @lines);
}

sub to_hash {
    my $self = shift;
    return {
        id => $self->id,
        name => $self->name,
        description => $self->description,
        version => $self->version,
        source_name => $self->source->name,
        restricted => $self->restricted,
        data_type => $self->data_type
    };
}

1;
