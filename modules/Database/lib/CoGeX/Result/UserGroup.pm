package CoGeX::Result::UserGroup;

use strict;
use warnings;

use base 'DBIx::Class::Core';

=head1 NAME

CoGeX::UserGroup

=head1 SYNOPSIS

This object uses the DBIx::Class to define an interface to the C<user_group> table in the CoGe database.

=head1 DESCRIPTION

Has columns:
C<user_group_id> (Primary Key)
Type: INT, Default: undef, Nullable: no, Size: 10

C<name>
Type: VARCHAR, Default: "", Nullable: no, Size: 50

C<description>
Type: VARCHAR, Default: undef, Nullable: yes, Size: 255


=head1 USAGE

  use CoGeX;

=head1 METHODS

=cut

__PACKAGE__->table("user_group");
__PACKAGE__->add_columns(
	"user_group_id",
	{ data_type => "INT", default_value => undef, is_nullable => 0, size => 11 },
	"creator_user_id",
	{ data_type => "INT", default_value => undef, is_nullable => 0, size => 11 },
	"name",
	{ data_type => "VARCHAR", default_value => "", is_nullable => 0, size => 250 },
	"description",
	{ data_type => "TEXT", default_value => undef, is_nullable => 1, size => 255 },
	"role_id",
	{ data_type => "INT", default_value => undef, is_nullable => 0, size => 11 },
	"locked",
	{ data_type => "INT", default_value => "0", is_nullable => 0, size => 1 }
);
__PACKAGE__->set_primary_key("user_group_id");
__PACKAGE__->has_many( 'user_group_connectors' => "CoGeX::Result::UserGroupConnector", 'user_group_id' );
__PACKAGE__->has_many( 'lists' => "CoGeX::Result::List", 'user_group_id' );
__PACKAGE__->belongs_to( 'role' => "CoGeX::Result::Role", 'role_id' );
__PACKAGE__->belongs_to( 'creator' => "CoGeX::Result::User", 'creator_user_id' );

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

sub users
{
	my $self  = shift;
	my @users = ();

	foreach my $ugc ( $self->user_group_connectors() )
	{
		push @users, $ugc->user;
	}

	return wantarray ? @users : \@users;
}

################################################ subroutine header begin ##

=head2 creator

 Usage     : 
 Purpose   : Returns user object for creator
 Returns   : user object
 Argument  : None
 Throws    : None
 Comments  : 

=cut

################################################## subroutine header end ##

sub creator
{
	my $self  = shift;
	my @users = ();

	foreach my $ugc ( $self->user_group_connectors() )
	{
		push @users, $ugc->user;
	}

	return wantarray ? @users : \@users;
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

sub experiments
{

	my $self = shift;
	my @experiments = ();
	foreach my $list ( $self->lists )
	{
		foreach my $experiment ( $list->experiments )
		{
			push @experiments, $experiment;
		}
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

sub features
{

	my $self     = shift;
	my @features = ();
	foreach my $list ( $self->lists )
	{
		foreach my $feature ( $list->features )
		{
			push @features, $feature;
		}
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

sub genomes
{

	my $self    = shift;
	my @genomes = ();
	foreach my $list ( $self->lists )
	{
		foreach my $genome ( $list->genomes )
		{
			push @genomes, $genome;
		}
	}
	return wantarray ? @genomes : \@genomes;
}

################################################ subroutine header begin ##

=head2 private_genomes

 Usage     : $self->genomes
 Purpose   : alias for sub private_genomes
 Returns   : wantArray of genomes
 Argument  : None
 Throws    : None
 Comments  : 

=cut

################################################## subroutine header end ##

sub private_genomes { return shift->genomes(@_); }

sub private_datasets # added by Eric 9/26/2012 
{ 
	my $self = shift;
	my @ds;
	foreach my $genome ($self->genomes)
	{
		push @ds, $genome->datasets;
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

sub info
{
	my $self = shift;
	my $info = $self->name;
	$info .= ": " . $self->description if $self->description;
	$info .= " (" . $self->role->name . ")";
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

sub info_html
{
	my $self = shift;
	my $info = $self->info;
	$info = qq{<span class=link onclick='window.open("GroupView.pl?ugid=} . $self->id . qq{")'>} . $info . "</span>";
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

	$anno_type = new CoGe::Accessory::Annotation( Type => "<tr><td valign='top' nowrap='true'><span class=\"title5\">" . "Role" . "</span>" );
	$anno_type->Type_delimit(": <td class=\"data5\">");
	$anno_type->Val_delimit("<br>");
	$anno_type->add_Annot( $self->role->name . ($self->role->description ? ' (' . $self->role->description . ')' : '') );
	$anno_type->add_Annot( "<span style='color:red;font-style:italic;'>Note: this group was created automatically and cannot be edited.</span>" ) if ($self->locked);
	$anno_obj->add_Annot($anno_type);

	$anno_type = new CoGe::Accessory::Annotation( Type => "<tr><td valign='top' nowrap='true'><span class=\"title5\">" . "Users" . "</span>" );
	$anno_type->Type_delimit(": <td class=\"data5\">");
	$anno_type->Val_delimit("<br>");
	foreach my $user (sort { $a->info cmp $b->info } $self->users) {
		$anno_type->add_Annot($user->info . ($user->id == $self->creator_user_id ? ' (creator)' : ''));
	}
	$anno_obj->add_Annot($anno_type);

	$anno_type = new CoGe::Accessory::Annotation( Type => "<tr valign='top'><td nowrap='true'><span class=\"title5\">" . "Lists" . "</span>" );
	$anno_type->Type_delimit(": <td class=\"data5\">");
	$anno_type->Val_delimit("<br>");
	foreach my $list ($self->lists) {
		my $a = $list->info_html . ($allow_delete ? "<span class='link ui-icon ui-icon-trash' onclick=\"remove_list_from_group({ugid: '" . $self->id . "', lid: '" . $list->id . "'});\"></span>" : '');
		$anno_type->add_Annot($a);
	}
	$anno_obj->add_Annot($anno_type);

  return "<table cellpadding=0 class='ui-widget-content ui-corner-all small'>".$anno_obj->to_String."</table>";
}

################################################ subroutine header begin ##

=head2 owner_list

 Usage     : return group's owner list
 Purpose   : 
 Returns   : list object
 Argument  : 
 Throws    : None
 Comments  : 

=cut

################################################## subroutine header end ##

sub owner_list {
	my $self = shift;
	foreach my $list ( $self->lists ) {
		return $list if ($list->list_type_id == 3); # FIXME list type hardcoded
	}
	return;
}

=head1 BUGS


=head1 SUPPORT


=head1 AUTHORS

 Eric Lyons
 Brent Pedersen

=head1 COPYRIGHT

This program is free software; you can redistribute
it and/or modify it under the same terms as Perl itself.

The full text of the license can be found in the
LICENSE file included with this module.




=head1 SEE ALSO


=cut

1;
