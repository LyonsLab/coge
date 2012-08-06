  package CoGeX_dev::ResultSet::Genome;

  use strict;
  use warnings;
  use base 'DBIx::Class::ResultSet';

 ################################################ subroutine header begin ##

=head2 resolve

 Usage     : 
 Purpose   : 
 Returns   : 
 Argument  : 
 Throws    : 
 Comments  : 

See Also   : 

=cut

################################################## subroutine header end ##

sub resolve {
    my $self = shift;
    my $info = shift;
    return $info if ref($info) =~ /Genome/i;
    return $self->find($info) if $info =~ /^\d+$/;
    return $self->search({ 'name' => { '-like' => '%' . $info . '%'}},
			 ,{});
  }

  1;
