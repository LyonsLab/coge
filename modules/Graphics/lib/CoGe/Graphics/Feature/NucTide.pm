package CoGe::Graphics::Feature::NucTide;
use strict;
use base qw(CoGe::Graphics::Feature);


BEGIN {
    use vars qw($VERSION $HEIGHT $WIDTH $ATC $GCC);
    $VERSION     = '0.1';
    $HEIGHT = 5;
    $WIDTH = 5;
    $ATC= [150,150,255];
    $GCC= [150,255,150];
    __PACKAGE__->mk_accessors(
"nt",
);
}

sub _initialize
  {
    my $self = shift;
    my %opts = @_;
    my $h = $HEIGHT; #total image height 
    my $w = $WIDTH;
    $self->image_width($w);
    $self->image_height($h);
    $self->bgcolor([255,255,255]) unless $self->bgcolor;
    $self->fill(1);
    $self->order(1);
    $self->stop($self->start + length $self->nt) unless $self->stop;
    my $at = 0;
    my $cg = 0;
    my $seq = $self->nt;
    while ($seq=~ /a|t|n/ig)
      {
	$at++;
      }
    while ($seq =~ /c|g|n/ig)
      {
	$cg++;
      }
    print STDERR "SEQ: $seq\n" unless ($at+$cg) > 0;
    my @color;
    for my $i (0..2)
      {
	push @color, $ATC->[$i]*$at/($at+$cg)+ $GCC->[$i]*$cg/($at+$cg);
      }
    $self->color(\@color);
#    if ($self->nt && $self->nt =~ /a|t/i)
#      {
#	$self->color($ATC) unless $self->color;
#      }
#    else
#      {
#	$self->color($GCC) unless $self->color;
#      }
    $self->label($self->nt) if $self->nt;
  }

sub _post_initialize
  {
    my $self = shift;
    my %opts = @_;
    my $gd = $self->gd;
    $gd->fill(0,0, $self->get_color($self->color));
  }

#################### subroutine header begin ####################

=head2 sample_function

 Usage     : How to use this function/method
 Purpose   : What it does
 Returns   : What it returns
 Argument  : What it wants to know
 Throws    : Exceptions and other anomolies
 Comment   : This is a sample subroutine header.
           : It is polite to include more pod and fewer comments.

See Also   : 

=cut

#################### subroutine header end ####################




#################### main pod documentation begin ###################
## Below is the stub of documentation for your module. 
## You better edit it!


=head1 NAME

CoGe::Graphics::Feature::Base

=head1 SYNOPSIS

  use CoGe::Graphics::Feature::Base


=head1 DESCRIPTION

=head1 USAGE



=head1 BUGS



=head1 SUPPORT



=head1 AUTHOR

	Eric Lyons
	elyons@nature.berkeley.edu

=head1 COPYRIGHT

This program is free software licensed under the...

	The Artistic License

The full text of the license can be found in the
LICENSE file included with this module.


=head1 SEE ALSO

perl(1).

=cut

#################### main pod documentation end ###################


1;
# The preceding line will help the module return a true value

