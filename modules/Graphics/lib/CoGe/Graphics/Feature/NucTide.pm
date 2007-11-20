package CoGe::Graphics::Feature::NucTide;
use strict;
use base qw(CoGe::Graphics::Feature);


BEGIN {
    use vars qw($VERSION $HEIGHT $WIDTH $ATC $GCC $NC $XC %EXTERNAL_IMAGES);
    $VERSION     = '0.1';
    $HEIGHT = 5;
    $WIDTH = 5;
    $ATC= [255,255,255]; #ATs are white
    $GCC= [150,255,150]; #GCs are green
    $NC= [255,225,150]; #Ns are orange
    $XC= [255,200,255]; #Xs are Pale Dull Magenta
    %EXTERNAL_IMAGES = (
			A=>'/opt/apache/CoGe/picts/A.png',
			T=>'/opt/apache/CoGe/picts/T.png',
			C=>'/opt/apache/CoGe/picts/C.png',
			G=>'/opt/apache/CoGe/picts/G.png',
		       );
    __PACKAGE__->mk_accessors(
"nt",
"gc",
"show_label",
"motif",
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
    $self->type('nt') unless $self->type;
    $self->stop($self->start + length $self->nt-1) unless $self->stop;
    $self->skip_overlap_search(1); #make sure to skip searching for overlap for these guys.  Search can be slow
    my $at = 0;
    my $cg = 0;
    my $n = 0;
    my $x = 0;
    my $seq = $self->nt;
    my $seq_len = length($seq) if $seq;
    $self->color(255,255,255) unless $self->color();
    if ($self->options)
      {
	($at) = $seq =~ tr/atrywmkhbvdATRYWMKHBVD/atrywmkhbvdATRYWMKHBVD/;
	($cg) = $seq =~ tr/cgrysmkhbvdCGRYWMKHBVD/cgrysmkhbvdCGRYWMKHBVD/;
	($n) = $seq =~ tr/nN/nN/;
	($x) = $seq =~ tr/xX/xX/;
	
	print STDERR "SEQ: $seq\n" unless ($at+$cg+$n+$x) > 0;
	
	my @color;
	
	for my $i (0..2)
	  {
	    my $color = 0;
	    $color += $self->options eq "gc" ? $ATC->[$i]*$at/($seq_len)+ $GCC->[$i]*$cg/($seq_len) : $ATC->[$i]*($at+$cg) /($seq_len)  ;
	    $color += $NC->[$i]*$n/($seq_len);
	    $color += $XC->[$i]*$x/($seq_len);

	    push @color, $color;
	  }
#	unless ($self->options eq "gc")
#	  {
#	    if ($color[2])
#	      {
#		$color[0] = 255;
#		$color[1] = 200;
#	      }
#	  }
	$self->color(\@color) if ($color[0] || $color[1] || $color[2]);
      }
    $self->label($self->nt) if $self->nt && $self->show_label;
  }

sub _post_initialize
  {
    my $self = shift;
    my %opts = @_;
    my $gd = $self->gd;
    $gd->fill(0,0, $self->get_color($self->color));
    if ($self->label && length ($self->label) == 1 && $EXTERNAL_IMAGES{uc($self->label)} && -r $EXTERNAL_IMAGES{uc($self->label)})
      {
	my $ei= GD::Image->new($EXTERNAL_IMAGES{uc($self->label)});
#	print STDERR Dumper ($ei);
	if ($self->strand =~ /-/ || $self->strand =~ /bot/i)
	  {
	    $ei->rotate180();
	  }
	$self->external_image($ei) if $self->use_external_image;
      }
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

