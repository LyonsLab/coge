package CoGe::Graphics::Feature::HSP;
use strict;
use base qw(CoGe::Graphics::Feature);


BEGIN {
    use vars qw($VERSION $HEIGHT $BLOCK_HEIGHT);
    $VERSION     = '0.1';
    $HEIGHT = 25;
    $BLOCK_HEIGHT = 20;
    __PACKAGE__->mk_accessors(
"alignment",
"color_matches", #flag for whether to just fill the hsp with a color, or color each match individually
);
}

sub _initialize
  {
    my $self = shift;
    my %opts = @_;
    my $h = $HEIGHT; #total image height 
    my $s = $self->start;
    my $e = $self->stop;;
    my $w = $e-$s;
    $self->image_width($w);
    $self->image_height($h);
    $self->merge_percent(100);
    $self->bgcolor([255,255,255]) unless $self->bgcolor;
    $self->skip_duplicate_search(1);
    $self->force_label(1);
    $self->type('HSP');
    $self->mag(1);
    $self->font_size(10);
    $self->color_matches(1) unless defined $self->color_matches();
  }

sub _post_initialize
  {
    my $self = shift;
    my %opts = @_;
    my $gd = $self->gd;
    my $w = $self->image_width;
    my $h = $self->image_height;
#    unless ($self->
    unless ($self->color_matches)
      {
	$gd->fill(0,0, $self->get_color($self->color));
        unless ($w < 10)
          {
            $self->_rounded_edges(x1=>0, y1=>0);
            $self->_rounded_edges(x1=>$w, y1=>0, negx=>-1);
            $self->_rounded_edges(x1=>$w, y1=>$h, negx=>-1, negy=>-1);
            $self->_rounded_edges(x1=>0, y1=>$h, negy=>-1);
          }
	return;
      }
    my $alignment = $self->alignment;
    my $count = 0;
    #make a colored box around the hsp
    $gd->fill(0,0, $self->get_color([255,255,255]));

    foreach my $chr (split //, $alignment)
      {
	if ($chr eq "|")
	  {
	    $gd->line($count, 0,  $count, $self->ih, $self->get_color($self->color));
	  }
	elsif ($chr ne " ") #probably an amino acid match
	  {
	    $gd->filledRectangle($count*3, 0,  $count*3+2, $self->ih, $self->get_color($self->color));
	  }  
	$count++;
      }
    $gd->rectangle(0,0,$self->iw-1, $self->ih-1, $self->get_color($self->color));
  }

sub _rounded_edges
{
	my $self = shift;
	my %opts = @_;
	my $x1 = $opts{x1};
	my $y1 = $opts{y1};
	my $neg_x = $opts{negx} || 1;
	my $neg_y = $opts{negy} || 1;
	my $gd = $self->gd;
	my $poly1 = GD::Polygon->new;
        my $poly2 = GD::Polygon->new;
        my $poly3 = GD::Polygon->new;
        $poly1->addPt($x1,$y1);
        $poly1->addPt($x1+($neg_x*3),$y1);
        $poly1->addPt($x1,$y1+($neg_y*1)); 
        $poly2->addPt($x1,$y1);
        $poly2->addPt($x1+($neg_x*2),$y1);
        $poly2->addPt($x1,$y1+($neg_y*2)); 
        $poly3->addPt($x1,$y1);
        $poly3->addPt($x1+($neg_x*1),$y1);
        $poly3->addPt($x1,$y1+($neg_y*3));
        $gd->filledPolygon($poly1, $self->get_color(255,255,255));
        $gd->filledPolygon($poly2, $self->get_color(255,255,255));
        $gd->filledPolygon($poly3, $self->get_color(255,255,255));
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

