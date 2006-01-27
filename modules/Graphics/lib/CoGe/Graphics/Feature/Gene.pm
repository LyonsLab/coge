package CoGe::Graphics::Feature::Gene;
use strict;
use base qw(CoGe::Graphics::Feature);


BEGIN {
    use vars qw($VERSION $HEIGHT $BLOCK_HEIGHT);
    $VERSION     = '0.1';
    $HEIGHT = 25;
    $BLOCK_HEIGHT = 20;
    __PACKAGE__->mk_accessors(
"block_height",
"segments", 
"print_label", #flag for printing feature label in gene
);
}

sub add_segment
  {
    my $self = shift;
    my %opts = @_;
    my $start = $opts{start} || $opts{begin} || $opts{START} || $opts{BEGIN};
    my $stop = $opts{stop} || $opts{end} || $opts{STOP} || $opts{END};
    if ($start > $stop)
      {
	my $tmp = $start;
	$start = $stop;
	$stop = $tmp;
      }
    my @segs;
    push @segs,  @{$self->segments} if $self->segments;
    push @segs, [$start, $stop];
    $self->segments([sort {$a->[0]<=>$b->[0]} @segs]);
  }

sub _initialize
  {
    my $self = shift;
    my %opts = @_;
    my $h = $HEIGHT; #total image height 
    my $s;
    my $e;

    foreach my $seg (sort {$a->[0] <=> $b->[0]} @{$self->segments})
      {
	$s = $seg->[0] unless $s;
	$e = $seg->[1] unless $e;
	$s = $seg->[0] if $seg->[0] < $s;
	$e = $seg->[1] if $seg->[1] > $e;
      }
    my $w = $e-$s;
    $self->start($s);
    $self->stop($e);
    $self->image_width($w);
    $self->image_height($h);
    $self->bgcolor([255,255,255]) unless $self->bgcolor;
    $self->block_height($BLOCK_HEIGHT) unless $self->block_height;
    $self->print_label(1) unless defined $self->print_label();
  }

sub _post_initialize
  {
    my $self = shift;
    my %opts = @_;
    my $gd = $self->gd;
    $gd->fill(0,0, $self->get_color($self->bgcolor));
    $gd->transparent($self->get_color($self->bgcolor));
    my $s = $self->start;
    my $black = $self->get_color(0,0,0);
    my $color = $self->get_color($self->color);
    my $last;
    my $c = $self->ih()/2;
    my $bh = $self->block_height/2;
    my @sorted = sort {$a->[0] <=> $b->[0]} @{$self->segments};
    foreach my $seg (@sorted)
      {
	my $x1 = $seg->[0] - $s;
	my $x2 = $seg->[1] - $s;
	my $y1 = $c-$bh;
	my $y2 = $c+$bh;;
	$gd->filledRectangle($x1,$y1, $x2, $y2, $color);
	$gd->rectangle($x1,$y1, $x2, $y2, $black);
	$gd->setStyle($black, $black, $black, GD::gdTransparent, GD::gdTransparent);
	if ($last)
	  {
	    my $mid = ($x1-$last)/2+$last;
	    $gd->line($last, $y1, $mid, 0, GD::gdStyled);
	    $gd->line($mid, 0, $x1, $y1, GD::gdStyled);
	  }
	$last = $x2;
      }

    my $x = $self->draw_arrow;
    
    $self->_gd_string(y=>$c-$bh+2, x=>$x, text=>$self->label, size=>$self->block_height-4) if $self->print_label;
  }

sub draw_arrow
  {
    my $self = shift;
    my %opts = @_;
    my $seg = $self->strand =~ /-/ ? $self->segments->[0]:$self->segments->[-1];

    my $gd = $self->gd;
    my $c = $self->ih()/2;
    my $y = $self->ih-1;
    my $w = ($seg->[1] - $seg->[0]); #width of segment
    my $x = $seg->[0] - $self->start; #start of segment
    my $poly = GD::Polygon->new;
    my $arrow_width = $self->ih*4;
    $arrow_width = $w if $arrow_width > $w;
#    print STDERR "Arrowhead: X: $x, width: $arrow_width\n";
    my $arrow_end;
    if ($self->strand =~ /-/i)
      {
	$gd->filledRectangle($x,0, $x+($arrow_width), $y, $self->get_color(255,255,255));
	$poly->addPt($x+($arrow_width), 0);
	$poly->addPt($x+($arrow_width), $y);
	$poly->addPt($x, $c);
	$arrow_end = $x+($arrow_width);
      }
    else
      {
	$gd->filledRectangle($x+($w-$arrow_width),0, $x+$w, $y, $self->get_color(255,255,255));
	$poly->addPt($x+($w-$arrow_width), 0);
	$poly->addPt($x+($w-$arrow_width), $y);
	$poly->addPt($x+($w), $c);
	$arrow_end=2;
      }
    $gd->filledPolygon($poly, $self->get_color($self->color));
    $gd->openPolygon($poly, $self->get_color());
    return $arrow_end;
  }

sub _draw_join
  {
    my $self = shift;
    my $s = shift;
    my $e = shift;
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

