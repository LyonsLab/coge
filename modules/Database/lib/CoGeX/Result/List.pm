package CoGeX::Result::List;

use strict;
use warnings;
use base 'DBIx::Class::Core';

use CoGe::Accessory::Annotation;

=head1 NAME

CoGeX::FeatureList

=head1 SYNOPSIS

This object uses the DBIx::Class to define an interface to the C<feature_list> table in the CoGe database.

=head1 DESCRIPTION

Has columns:
C<feature_list_id> (Primary Key)
Type: INT, Default: undef, Nullable: no, Size: 10

C<name>
Type: VARCHAR, Default: "", Nullable: no, Size: 50

C<description>
Type: VARCHAR, Default: undef, Nullable: yes, Size: 255

C<feature_list_group_id>
Type: INT, Default: undef, Nullable: yes, Size: 10

C<notes>
Type: VARCHAR, Default: undef, Nullable: yes, Size: 1024


=head1 USAGE

  use CoGeX;

=head1 METHODS

=cut

__PACKAGE__->table("list");
__PACKAGE__->add_columns(
	"list_id",
	{ data_type => "INT", default_value => undef, is_nullable => 0, size => 11 },
	"name",
	{ data_type => "VARCHAR", default_value => "", is_nullable => 0, size => 255 },
	"description",
	{
		data_type     => "VARCHAR",
		default_value => undef,
		is_nullable   => 1,
		size          => 1024,
	},
	"list_type_id",
	{ data_type => "INT", default_value => undef, is_nullable => 1, size => 11 },
	"user_group_id",
	{ data_type => "INT", default_value => undef, is_nullable => 1, size => 11 },
	"restricted",
	{
		data_type     => "BOOLEAN",
		default_value => 0,
		is_nullable   => 0,
		size          => 1
	},
	"locked",
	{ data_type => "INT", default_value => "0", is_nullable => 0, size => 1 }
);
__PACKAGE__->set_primary_key("list_id");

#__PACKAGE__->has_many( "experiment_list_connectors" => "CoGeX::Result::ExperimentListConnector", 'list_id' );
#__PACKAGE__->has_many( "genome_list_connectors"     => "CoGeX::Result::GenomeListConnector",     'list_id' );
#__PACKAGE__->has_many( "feature_list_connectors"    => "CoGeX::Result::FeatureListConnector",    'list_id' );
#__PACKAGE__->has_many( "list_collection_connectors" => "CoGeX::Result::ListCollectionConnector", 'list_id' );
__PACKAGE__->has_many( "list_annotations"           => "CoGeX::Result::ListAnnotation",          'list_id' );
__PACKAGE__->belongs_to( "user_group" => "CoGeX::Result::UserGroup", 'user_group_id' );
__PACKAGE__->belongs_to( "list_type"  => "CoGeX::Result::ListType",  'list_type_id' );

__PACKAGE__->has_many( "list_connectors_as_child"    => "CoGeX::Result::ListConnector",  {'foreign.child_id' => 'self.list_id' } );
__PACKAGE__->has_many( "list_connectors_as_parent"    => "CoGeX::Result::ListConnector",  {'foreign.parent_id' => 'self.list_id' } );

sub group
{
	return shift->user_group(@_);
}

sub lists
{
	my $self = shift;
	my %opts = @_;
	
	my @lists;
	foreach my $conn ( $self->list_connectors_as_parent )
	{
		push @lists, $conn->child if ($conn->is_list);
	}
	return wantarray ? @lists : \@lists;		
}

sub features
{
	my $self = shift;
	my %opts = @_;
	
	my @features;
	foreach my $conn ( $self->list_connectors_as_parent )
	{
		push @features, $conn->child if ($conn->is_feature);
	}
	return wantarray ? @features : \@features;	
}

sub genomes
{
	my $self = shift;
	my %opts = @_;
	
	my @genomes;
	foreach my $conn ( $self->list_connectors_as_parent )
	{
		push @genomes, $conn->child if ($conn->is_genome);
	}
	return wantarray ? @genomes : \@genomes;	
}

sub experiments
{
	my $self = shift;
	my %opts = @_;

	my @experiments;
	foreach my $conn ( $self->list_connectors_as_parent )
	{
		push @experiments, $conn->child if ($conn->is_experiment);
	}
	return wantarray ? @experiments : \@experiments;
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

sub annotations
{
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

sub type
{
	return shift->list_type();
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

sub info
{
	my $self = shift;
	my $info = $self->name;
	$info = "(R)" . $info if $self->restricted;
	$info .= ": " . $self->description      if $self->description;
	$info .= " (" . $self->type->name . ")" if $self->type;
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

sub info_html
{
	my $self = shift;
	my $info = $self->info;
	return qq{<span class=link onclick='window.open("ListView.pl?lid=} .$self->id. qq{")'>} . $info . "</span>";
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

sub data_summary
{
	my $self = shift;
	my @stuff;
	push @stuff, "Experiments: " . @{ $self->experiments } if @{ $self->experiments };
	push @stuff, "Features: "    . @{ $self->features }    if @{ $self->features };
	push @stuff, "Genomes: "     . @{ $self->genomes }     if @{ $self->genomes };
	push @stuff, "Lists: "       . @{ $self->lists }       if @{ $self->lists };
	return join( "; ", @stuff );
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

sub annotation_pretty_print_html
{
	my $self     = shift;
	my %opts     = @_;
	my $minimal  = $opts{minimal};
	my $anno_obj = new CoGe::Accessory::Annotation( Type => "anno" );
	$anno_obj->Val_delimit("\n");
	$anno_obj->Add_type(0);
	$anno_obj->String_end("\n");
	my $anno_type = new CoGe::Accessory::Annotation( Type => "<tr><td nowrap='true'><span class=\"title5\">" . "Name" . "</span>" );
	$anno_type->Type_delimit(": <td class=\"data5\">");
	$anno_type->add_Annot( $self->name . "</td>" );
	$anno_obj->add_Annot($anno_type);
	$anno_type = new CoGe::Accessory::Annotation( Type => "<tr><td nowrap='true'><span class=\"title5\">" . "Description" . "</span>" );
	$anno_type->Type_delimit(": <td class=\"data5\">");
	$anno_type->add_Annot( $self->description . "</td>" );
	$anno_obj->add_Annot($anno_type);

#	foreach my $anno ( sort { uc( $a->type->name ) cmp uc( $b->type->name ) } $self->annotations( {}, { prefetch => { annotation_type => 'annotation_type_group' } } ) )
#	{
#		my $type      = $anno->type();
#		my $group;
#		my $anno_name;
#		
#		if (not $type) {
#			$group = 'undefined';
#			$anno_name = 'undefined';			
#		}
#		else {
#			$group     = $type->group();
#			$anno_name = $type->name;
#			$anno_name .= ", " . $type->description if $type->description;
#			if ( ref($group) =~ /group/i && !( $type->name eq $group->name ) )
#			{
#				{
#					$anno_name .= ":" unless $anno_name =~ /:$/;
#					$anno_name = "<span class=\"title5\">" . $anno_name . "</span>";
#				}
#			}
#			else
#			{
#				if ( $anno->link )
#				{
#					$anno_name = "<tr><td nowrap='true'><span class=\"coge_link\">" . $anno_name . "</span>";
#				}
#				else
#				{
#					$anno_name = "<tr><td nowrap='true'><span class=\"title5\">" . $anno_name . "</span>";
#				}
#			}
#		}
#		
#		my $anno_type = new CoGe::Accessory::Annotation( Type => $anno_name );
#		$anno_type->Val_delimit("<br>");
#		$anno_type->Type_delimit(" ");
#		my $annotation = "<span class=\"data5";
#		$annotation .= qq{ link" onclick="window.open('} . $anno->link . qq{')} if $anno->link;
#		$annotation .= "\">" . $anno->annotation . "</span>";
#		$anno_type->add_Annot($annotation) if $anno->annotation;
#
#		if ( ref($group) =~ /group/i && !( $type->name eq $group->name ) )
#		{
#			my $class = $anno->link ? "coge_link" : "title5";
#			my $anno_g = new CoGe::Accessory::Annotation( Type => "<tr><td nowrap='true'><span class=\"$class\">" . $group->name . "</span>" );
#			$anno_g->add_Annot($anno_type);
#			$anno_g->Type_delimit(":<td>");
#			$anno_g->Val_delimit(", ");
#
#			#	    $anno_g->Val_delimit(" ");
#			$anno_obj->add_Annot($anno_g);
#		}
#		else
#		{
#			$anno_type->Type_delimit(":<td>");
#			$anno_obj->add_Annot($anno_type);
#		}
#	}
	
	$anno_type = new CoGe::Accessory::Annotation( Type => "<tr><td nowrap='true'><span class=\"title5\">" . "Type" . "</span>" );
	$anno_type->Type_delimit(": <td class=\"data5\">");
	if ($self->type)
	{
		my $type = $self->type->name;
		$type .= ": " . $self->type->description if $self->type->description;
		$anno_type->add_Annot( $type . "</td>" );
	}
	$anno_obj->add_Annot($anno_type);
  
	$anno_type = new CoGe::Accessory::Annotation( Type => "<tr><td nowrap='true'><span class=\"title5\">" . "User Group" . "</span>" );
	$anno_type->Type_delimit(": <td class=\"data5\">");
	my $group = $self->group->info_html;
	$anno_type->add_Annot( $group . "</td>" );
	$anno_obj->add_Annot($anno_type);

	$anno_type = new CoGe::Accessory::Annotation( Type => "<tr><td nowrap='true'><span class=\"title5\">" . "Restricted" . "</span>" );
	$anno_type->Type_delimit(": <td class=\"data5\">");
	my $restricted = $self->restricted ? "Yes" : "No";
	$anno_type->add_Annot( $restricted . "</td>" );
	$anno_obj->add_Annot($anno_type);
	
	if ( $self->locked )
	{
		$anno_type = new CoGe::Accessory::Annotation( Type => "<tr><td nowrap='true'><span class=\"title5\">" . "Locked" . "</span>" );
		$anno_type->Type_delimit(": <td class=\"data5\">");
		$anno_type->Val_delimit("<br>");
		$anno_type->add_Annot("Yes</td>");
		$anno_obj->add_Annot($anno_type);
	}
	
#	if ( my @cols = $self->lists )
#	{
#		foreach my $col (@cols)
#		{
#			$anno_type = new CoGe::Accessory::Annotation( Type => "<tr><td nowrap='true'><span class=\"title5\">" . "List of lists" . "</span>" );
#			$anno_type->Type_delimit(": <td class=\"data5\">");
#			$anno_type->Val_delimit("<br>");
#			my $col_name = $col->name;
#			$col_name .= ": " . $col->description if $col->description;
#			#$col_name = qq{<span class='link' onclick="window.open('} . $col->link . qq{')} . $col_name . "</span>" if $col->link;
#			$anno_type->add_Annot( $col_name . "</td>" );
#			$anno_obj->add_Annot($anno_type);
#		}
#	}

	return "<table cellpadding=0 class='ui-widget-content ui-corner-all small'>" . $anno_obj->to_String . "</table>";
}

1;

=head1 BUGS


=head1 SUPPORT


=head1 AUTHORS

 Eric Lyons

=head1 COPYRIGHT 2012

This program is free software; you can redistribute
it and/or modify it under the same terms as Perl itself.

The full text of the license can be found in the
LICENSE file included with this module.


=head1 SEE ALSO

=cut

