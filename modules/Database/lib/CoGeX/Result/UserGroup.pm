package CoGeX::Result::UserGroup;

# Created by DBIx::Class::Schema::Loader v0.03009 @ 2006-12-01 18:13:38

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
			 "user_group_id",  { data_type => "INT", default_value => undef, is_nullable => 0, size => 11 },
			 "name",  { data_type => "VARCHAR", default_value => "", is_nullable => 0, size => 250 },
			 "description",
			 {
			  data_type => "TEXT",
			  default_value => undef,
			  is_nullable => 1,
			  size => 255,
			 },
			 "role_id",  { data_type => "INT", default_value => undef, is_nullable => 0, size => 11 },
);
__PACKAGE__->set_primary_key("user_group_id");
__PACKAGE__->has_many('user_group_connectors'=>"CoGeX::Result::UserGroupConnector",'user_group_id');
__PACKAGE__->has_many('user_group_data_connectors'=>"CoGeX::Result::UserGroupDataConnector",'user_group_id');
__PACKAGE__->belongs_to('role'=>"CoGeX::Result::Role",'role_id');


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
    my @users=();

    foreach my $ugc ($self->user_group_connectors())
      {
        push @users, $ugc->user;
    }

    return wantarray ? @users : \@users;
}



################################################ subroutine header begin ##

=head2 private_genomes

 Usage     : 
 Purpose   : Returns the set of genomes associated with a genome (dataset_group)
 Returns   : wantArray of genomes
 Argument  : None
 Throws    : None
 Comments  : 



=cut

################################################## subroutine header end ##


sub private_genomes {

    my $self = shift;
    my @private_genomes=();

    foreach my $group_dataset ($self->user_group_data_connectors())
      {
	push(@private_genomes,$group_dataset->genome()) if $group_dataset->genome;
    }

    return wantarray ? @private_genomes : \@private_genomes;
}

################################################ subroutine header begin ##

=head2 genomes

 Usage     : $self->genomes
 Purpose   : alias for sub private_genomes
 Returns   : wantArray of genomes
 Argument  : None
 Throws    : None
 Comments  : 



=cut

################################################## subroutine header end ##


sub genomes { return shift->private_genomes(@_);}


################################################ subroutine header begin ##

=head2 private_datasets

 Usage     : 
 Purpose   : Returns the set of genomes associated with a genome (dataset_group)
 Returns   : Array of Groups
 Argument  : None
 Throws    : None
 Comments  : 



=cut

################################################## subroutine header end ##


sub private_datasets {

    my $self = shift;
    my @private_genomes=();

    foreach my $group_dataset ($self->user_group_data_connectors()){
	push(@private_genomes,$group_dataset->dataset())if $group_dataset->dataset;
    }

    return wantarray ? @private_genomes : \@private_genomes;
}


################################################ subroutine header begin ##

=head2 datasets

 Usage     : $self->datasets
 Purpose   : alias for sub private_datasets
 Returns   : wantArray of datasets
 Argument  : None
 Throws    : None
 Comments  : 



=cut

################################################## subroutine header end ##


sub datasets { return shift->private_datasets(@_);}

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
