package CoGeX_dev::Result::UserGroup;

# Created by DBIx::Class::Schema::Loader v0.03009 @ 2006-12-01 18:13:38

use strict;
use warnings;

use base 'DBIx::Class::Core';


=head1 NAME

CoGeX_dev::UserGroup

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

  use CoGeX_dev;

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
			 "locked",  { data_type => "INT", default_value => "0", is_nullable => 0, size => 1 }
);
__PACKAGE__->set_primary_key("user_group_id");
__PACKAGE__->has_many('user_group_connectors'=>"CoGeX_dev::Result::UserGroupConnector",'user_group_id');
__PACKAGE__->has_many('lists'=>"CoGeX_dev::Result::List",'user_group_id');
__PACKAGE__->belongs_to('role'=>"CoGeX_dev::Result::Role",'role_id');


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
    my @experiments=();
    foreach my $list ($self->lists)
      {
	foreach my $experiment ($list->experiments)
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


sub features {

    my $self = shift;
    my @features=();
    foreach my $list ($self->lists)
      {
	foreach my $feature ($list->features)
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


sub genomes {

    my $self = shift;
    my @genomes=();
    foreach my $list ($self->lists)
      {
	foreach my $genome ($list->genomes)
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


sub private_genomes { return shift->genomes(@_);}



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
