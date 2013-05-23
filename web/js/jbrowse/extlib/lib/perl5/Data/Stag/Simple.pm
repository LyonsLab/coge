package Data::Stag::Simple;

=head1 NAME

  Data::Stag::Simple

=head1 SYNOPSIS


=cut

=head1 DESCRIPTION

=head1 PUBLIC METHODS -

=cut

use strict;
use base qw(Data::Stag::Base);

use vars qw($VERSION);
$VERSION="0.11";


sub start_event { shift;print "OPEN :$_[0]\n"} 
sub end_event { shift;print "CLOSE:$_[0]\n"} 
sub evbody { shift;my $b=shift || "";print "BODY:$b\n"} 

1;
