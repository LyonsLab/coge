package CoGeX::Result::Experiment;

use strict;
use warnings;
use base 'DBIx::Class::Core';
use File::Spec::Functions;
use Data::Dumper;
use POSIX;
use CoGe::Accessory::Annotation;

=head1 NAME

CoGeX::Result::Experiment

=head1 SYNOPSIS

This object uses the DBIx::Class to define an interface to the C<dataset_group> table in the CoGe database.

=head1 DESCRIPTION

Has columns:
C<experiment_id> (Primary Key)
Type: INT, Default: undef, Nullable: no, Size: 11

C<genome_id> 
Type: INT, Default: undef, Nullable: no, Size: 11

C<data_source_id>
Type: INT, Default: 0, Nullable: 0, Size: 11

C<name>
Type:VARCHAR, Default: "", Nullable: no, Size: 255

C<description>
Type: Text, Default: undef, Nullable: yes

C<version>
Type:VARCHAR, Default: undef, Nullable: no, Size: 50

C<storage_path>
Type: VARCHAR, Default: undef, Nullable: 0, Size: 255

C<link>
Type: TEXT, Defaullt:  undef, Nullable: 1

C<restricted>
Type: INT, Default: 0, Nullable: 0, Size: 1

Belongs to CCoGeX::Result::DataSource> via C<data_source_id>
Belongs to CCoGeX::Result::Genome> via C<genome_id>
Has many CCoGeX::Result::ExperimentTypeConnector> via C<experiment_id>
Has many CCoGeX::Result::ExperimentAnnotation> via C<experiment_id>


=head1 USAGE

 use CoGeX;

=head1 METHODS

=cut

__PACKAGE__->table("experiment");
__PACKAGE__->add_columns(
	"experiment_id",
	{ data_type => "INT", default_value => undef, is_nullable => 0, size => 11 },
	"genome_id",
	{ data_type => "INT", default_value => undef, is_nullable => 0, size => 11 },
	"data_source_id",
	{ data_type => "INT", default_value => undef, is_nullable => 0, size => 11 },
	"name",
	{ data_type => "VARCHAR", default_value => "", is_nullable => 0, size => 255 },
	"description",
	{
		data_type     => "TEXT",
		default_value => undef,
		is_nullable   => 1,
	},
	"version",
	{
		data_type     => "VARCHAR",
		default_value => undef,
		is_nullable   => 0,
		size          => 50,
	},
	"storage_path",
	{
		data_type     => "VARCHAR",
		default_value => undef,
		is_nullable   => 0,
		size          => 255,
	},
	"restricted",
	{ data_type => "int", default_value => "0", is_nullable => 0, size => 1 },
	"access_count",
	{ data_type => "int", default_value => "0", is_nullable => 1, size => 10 }
);

__PACKAGE__->set_primary_key("experiment_id");
__PACKAGE__->has_many( "experiment_type_connectors" => "CoGeX::Result::ExperimentTypeConnector", 'experiment_id' );
__PACKAGE__->has_many( "experiment_annotations"     => "CoGeX::Result::ExperimentAnnotation",    'experiment_id' );
__PACKAGE__->has_many( "list_connectors" => "CoGeX::Result::ListConnector", {'foreign.child_id' => 'self.experiment_id'} );
__PACKAGE__->belongs_to( "data_source" => "CoGeX::Result::DataSource", 'data_source_id' );
__PACKAGE__->belongs_to( "genome"      => "CoGeX::Result::Genome",     'genome_id' );


################################################ subroutine header begin ##

=head2 desc

 Usage     : 
 Purpose   : alias for $self->description
 Returns   : 
 Argument  : 
 Throws    : 
 Comments  : 

See Also   : 

=cut

################################################## subroutine header end ##

sub desc
{
	return shift->description(@_);
}

################################################ subroutine header begin ##

=head2 source

 Usage     : 
 Purpose   : alias for $self->data_source
 Returns   : 
 Argument  : 
 Throws    : 
 Comments  : 

See Also   : 

=cut

################################################## subroutine header end ##

sub source
{
	shift->data_source(@_);
}

################################################ subroutine header begin ##

=head2 annotations

 Usage     : 
 Purpose   : alias for $self->experiment_annotations
 Returns   : 
 Argument  : 
 Throws    : 
 Comments  : 

See Also   : 

=cut

################################################## subroutine header end ##

sub annotations
{
	shift->experiment_annotations(@_);
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

sub experiment_types
{
	map { $_->experiment_type } shift->experiment_type_connectors();
}

################################################ subroutine header begin ##

=head2 types

 Usage     : 
 Purpose   : alias for $self->experiment_types
 Returns   : 
 Argument  : 
 Throws    : 
 Comments  : 

See Also   : 

=cut

################################################## subroutine header end ##

sub types
{
	shift->experiment_types(@_);
}

################################################ subroutine header begin ##

=head2 user_groups

 Usage     : 
 Purpose   : alias for $self->experiment_types
 Returns   : 
 Argument  : 
 Throws    : 
 Comments  : 

See Also   : 

=cut

################################################## subroutine header end ##

sub user_groups
{
	map { $_->user_group } shift->user_group_experiment_connectors();
}

################################################ subroutine header begin ##

=head2 get_path

 Usage     : 
 Purpose   : This method determines the correct directory structure for storing
 			the sequence files for an experiment.
 Returns   : 
 Argument  : 
 Throws    : none
 Comments  : The idea is to build a dir structure that holds large amounts of
 			files, and is easy to lookup based on experiment ID number.
			The strucuture is three levels of directorys, and each dir holds
			1000 files and/or directorys.
			Thus:
			./0/0/0/ will hold files 0-999
			./0/0/1/ will hold files 1000-1999
			./0/0/2/ will hold files 2000-2999
			./0/1/0 will hold files 1000000-1000999
			./level0/level1/level2/

See Also   : 

=cut

################################################## subroutine header end ##

sub get_path
{
	my $self          = shift;
	my $experiment_id = $self->id;
	my $level0        = floor( $experiment_id / 1000000000 ) % 1000;
	my $level1        = floor( $experiment_id / 1000000 ) % 1000;
	my $level2        = floor( $experiment_id / 1000 ) % 1000;
	my $path          = catdir( $level0, $level1, $level2, $experiment_id );    #adding experiment_id for final directory.
	return $path;
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

sub lists
{
	my $self = shift;
	my %opts = @_;
	my @lists;
	foreach my $lc ($self->list_connectors)
	{
		push @lists, $lc->parent_list if ($lc->child_type == 3); # FIXME hardcoded type value
	}
	return wantarray ? @lists : \@lists;
}


sub annotation_pretty_print_html
{
	my $self     = shift;
	my %opts     = @_;
	my $minimal  = $opts{minimal};
	my $allow_delete = $opts{allow_delete};
	
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

	$anno_type = new CoGe::Accessory::Annotation( Type => "<tr><td nowrap='true'><span class=\"title5\">" . "Genome" . "</span>" );
	$anno_type->Type_delimit(": <td class=\"data5\">");
	$anno_type->add_Annot( $self->genome->info_html . "</td>" );
	$anno_obj->add_Annot($anno_type);

	$anno_type = new CoGe::Accessory::Annotation( Type => "<tr><td nowrap='true'><span class=\"title5\">" . "Source" . "</span>" );
	$anno_type->Type_delimit(": <td class=\"data5\">");
	$anno_type->add_Annot( $self->source->info_html . "</td>" );
	$anno_obj->add_Annot($anno_type);

	$anno_type = new CoGe::Accessory::Annotation( Type => "<tr><td nowrap='true'><span class=\"title5\">" . "Version" . "</span>" );
	$anno_type->Type_delimit(": <td class=\"data5\">");
	$anno_type->add_Annot( $self->version . "</td>" );
	$anno_obj->add_Annot($anno_type);
	
	$anno_type = new CoGe::Accessory::Annotation( Type => "<tr valign='top'><td nowrap='true'><span class=\"title5\">" . "Types" . "</span>" );
	$anno_type->Type_delimit(": <td class=\"data5\">");
	$anno_type->Val_delimit("<br>");
	foreach my $type ($self->types)
		{
			my $info = $type->name;
			$info .= ": ".$type->description if $type->description;
			if ($allow_delete) {
				# NOTE: it is somewhat undesirable to have a javascript call in a DB object, but it works
				$info .= "<span onClick=\"remove_experiment_type({eid: '" . $self->id . "', etid: '" . $type->id . "'});\" class=\"link ui-icon ui-icon-trash\"></span>";
			}
			$anno_type->add_Annot( $info);
		}
	$anno_obj->add_Annot($anno_type);

	$anno_type = new CoGe::Accessory::Annotation( Type => "<tr valign='top'><td nowrap='true'><span class=\"title5\">" . "Lists" . "</span>" );
	$anno_type->Type_delimit(": <td class=\"data5\">");
	$anno_type->Val_delimit("<br>");
	foreach my $list ($self->lists)
		{
			next if (not $list); #FIXME mdb 8/21/12 - why necessary?
			my $info = $list->info_html;
			$anno_type->add_Annot( $info);
		}
	$anno_obj->add_Annot($anno_type);

#	foreach my $anno ( sort { uc( $a->type->name ) cmp uc( $b->type->name ) } $self->annotations( {}, { prefetch => { annotation_type => 'annotation_type_group' } } ) )
#	{
#		my $type      = $anno->type();
#		my $group     = $type->group();
#		my $anno_name = $type->name;
#		$anno_name .= ", " . $type->description if $type->description;
#		if ( ref($group) =~ /group/i && !( $type->name eq $group->name ) )
#		{
#			{
#				$anno_name .= ":" unless $anno_name =~ /:$/;
#				$anno_name = "<span class=\"title5\">" . $anno_name . "</span>";
#			}
#		}
#		else
#		{
#			if ( $anno->link )
#			{
#				$anno_name = "<tr><td nowrap='true'><span class=\"coge_link\">" . $anno_name . "</span>";
#			}
#			else
#			{
#				$anno_name = "<tr><td nowrap='true'><span class=\"title5\">" . $anno_name . "</span>";
#			}
#		}
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
#			#$anno_g->Val_delimit(" ");
#			$anno_obj->add_Annot($anno_g);
#		}
#		else
#		{
#			$anno_type->Type_delimit(":<td>");
#			$anno_obj->add_Annot($anno_type);
#		}
#	}
	
	$anno_type = new CoGe::Accessory::Annotation( Type => "<tr><td nowrap='true'><span class=\"title5\">" . "Restricted" . "</span>" );
	$anno_type->Type_delimit(": <td class=\"data5\">");
	my $restricted = $self->restricted ? "Yes" : "No";
	$anno_type->add_Annot( $restricted . "</td>" );
	$anno_obj->add_Annot($anno_type);	
	
	return "<table cellpadding=0 class='ui-widget-content ui-corner-all small'>" . $anno_obj->to_String . "</table>";
}

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

sub info
{
	my $self = shift;
	my $info;
	$info .= "&reg; " if $self->restricted;
	$info .= $self->name;
	$info .= ": " . $self->description if $self->description;
	$info .= " (v" . $self->version . ", eid" . $self->id . "): " . $self->source->name;
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

sub info_html
{
	my $self = shift;
	my $info = $self->info;
	return qq{<span class=link onclick='window.open("ExperimentView.pl?eid=} .$self->id. qq{")'>} . $info . "</span>";
}

=head1 BUGS


=head1 SUPPORT


=head1 AUTHORS

 Eric Lyons
 Matt Bomhoff

=head1 COPYRIGHT

This program is free software; you can redistribute
it and/or modify it under the same terms as Perl itself.

The full text of the license can be found in the
LICENSE file included with this module.

=head1 SEE ALSO

=cut

1;
