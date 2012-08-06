package CoGeX_dev::ResultSet::List;

use strict;
use warnings;
use base 'DBIx::Class::ResultSet';

################################################ subroutine header begin ##

=head2 public_lists

 Usage     : 
 Purpose   : 
 Returns   : 
 Argument  : 
 Throws    : 
 Comments  : 
           : 

See Also   : 

=cut

################################################## subroutine header end ##

sub public_lists
  {
    my $self = shift;
    my %opts = @_;
    return $self->search({restricted=>0});
  }

1;