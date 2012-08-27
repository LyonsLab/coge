package CoGeX::ResultSet::Experiment;

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
	return $info if ref($info) =~ /experiment/i;
	return $self->find($info) if $info =~ /^\d+$/;
	return $self->search( { 'name' => { '-like' => '%' . $info . '%' } },, {} );
}

################################################ subroutine header begin ##

=head2 public

 Usage     : 
 Purpose   : 
 Returns   : 
 Argument  : 
 Throws    : 
 Comments  : 

See Also   : 

=cut

################################################## subroutine header end ##

sub public {
	my $self  = shift;
	my %opts  = @_;
	my $limit = $opts{limit};
	return $self->search( { 'restricted' => 0 } );
}

1;
