  package CoGeX_dev::ResultSet::FeatureName;

  use strict;
  use warnings;
  use base 'DBIx::Class::ResultSet';

################################################ subroutine header begin ##

=head2 esearch

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

sub esearch {
    my $self = shift;
    my $join = $_[1]{'join'};
    my $prefetch = $_[1]{'prefetch'};

    map { push(@$prefetch, $_ ) } 
        ({ 'feature' => ['locations','feature_type'] });


    $_[1]{'join'} = $join;
    $_[1]{'prefetch'} = $prefetch;
    my $rs = $self->search(
         @_
    );
    return $rs;

}

  1;
