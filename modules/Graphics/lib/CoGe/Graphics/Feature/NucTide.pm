package CoGe::Graphics::Feature::NucTide;
use strict;
use base qw(CoGe::Graphics::Feature);


BEGIN {
    use vars qw($VERSION $HEIGHT $WIDTH $ATC $GCC $NC %EXTERNAL_IMAGES);
    $VERSION     = '0.1';
    $HEIGHT = 5;
    $WIDTH = 5;
    $ATC= [255,255,255]; #ATs are white
    $GCC= [100,255,100]; #GCs are green
    $NC= [255,200,125]; #Ns are orange
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
    $self->type('nt');
    $self->stop($self->start + length $self->nt-1) unless $self->stop;
    $self->skip_overlap_search(1); #make sure to skip searching for overlap for these guys.  Search can be slow
    my $at = 0;
    my $cg = 0;
    my $n = 0;
    my $seq = $self->nt;
    $self->color(255,255,255);
    if ($self->options && $self->options eq "gc")
      {
	while ($seq=~ /a|t|r|y|w|m|k|h|b|v|d/ig)
	  {
	    $at++;
	  }
	while ($seq =~/c|g|r|y|s|m|k|h|b|v|d/ig)
	  {
	    $cg++;
	  }
	while ($seq =~/n|\?/ig)
	  {
	    $n++;
	  }
	
	print STDERR "SEQ: $seq\n" unless ($at+$cg+$n) > 0;
	
	my $pat = $at/($at+$cg+$n) if $at+$cg+$n > 0;
	my $pcg = $cg/($at+$cg+$n) if $at+$cg+$n > 0;
	my $pn = $n/($at+$cg+$n) if $at+$cg+$n > 0;
	my @color;
	for my $i (0..2)
	  {
	    push @color, $ATC->[$i]*$at/($at+$cg+$n)+ $GCC->[$i]*$cg/($at+$cg+$n)+ $NC->[$i]*$n/($at+$cg+$n);
	  }
      	#    print join (":", @color),"  ", $at, ":", $cg,"\n";
	$self->color(\@color)
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

