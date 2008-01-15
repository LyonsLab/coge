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
"add_type", #flag to set if gene type should be appended to label
"no_three_D", #switch between flat images and "3D" images, if given value, use flat images
"arrow_width", #width of arrow in pixels
);
}

sub add_segment
  {
    my $self = shift;
    my %opts = @_;
    my $start = $opts{start} || $opts{begin} || $opts{START} || $opts{BEGIN};
    my $stop = $opts{stop} || $opts{end} || $opts{STOP} || $opts{END};
#    use Data::Dumper;
#    print STDERR Dumper \%opts;
#    print STDERR "$start - $stop\n";
    if ($start > $stop)
      {
	my $tmp = $start;
	$start = $stop;
	$stop = $tmp;
      }
#    $start = 1 if $start < 1;
    return unless $start && $stop;
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
    $self->segments([]) unless $self->segments;
    foreach my $seg (sort {$a->[0] <=> $b->[0]} @{$self->segments})
      {
	$s = $seg->[0] unless $s;
	$e = $seg->[1] unless $e;
	$s = $seg->[0] if $seg->[0] < $s;
	$e = $seg->[1] if $seg->[1] > $e;
      }
    my $w = $e || $s ? $e-$s : 1;
    $self->start($s);
    $self->stop($e);
    $self->image_width($w);
    $self->image_height($h);
    $self->merge_percent(100);
    $self->bgcolor([255,255,255]) unless $self->bgcolor;
    $self->block_height($BLOCK_HEIGHT) unless $self->block_height;
    $self->print_label(0) unless defined $self->print_label();
    $self->skip_duplicate_search(1);
  }

sub _post_initialize
  {
    my $self = shift;
    my %opts = @_;
    $self->label($self->label." (".$self->type.")") if $self->add_type && $self->type;
    my $gd = $self->gd;
    $gd->fill(0,0, $self->get_color($self->bgcolor));
#    $gd->transparent($self->get_color($self->bgcolor));
    my $s = $self->start;
    my @tmp;
    foreach my $c (@{$self->color})
      {
	my $ct = $c;
	$ct-=100;
	$ct = 1 if $ct < 1;
	push @tmp, $ct;
      }
    my $border = $self->get_color(@tmp);
    my $color = $self->get_color($self->color);
    my $black = $self->get_color(1,1,1);
    my $last;
    my $c = $self->ih()/2;
    my $bh = $self->block_height/2;
    my @sorted = sort {$a->[0] <=> $b->[0]} @{$self->segments};
    foreach my $seg (@sorted)
      {
	my $x1 = $seg->[0] - $s;
	my $x2 = $seg->[1] - $s;
	my $y1 = $c-$bh;
	my $y2 = $c+$bh;
	my $gray_lvl = 15 || $opts{gray_lvl};
	if (!$self->no_three_D)
	  {
	    my @colors = $gd->rgb($color);
	    my ($r, $g, $b) = ($colors[0],$colors[1],$colors[2]);
	    if ($r -75 >= 0)
	      {$r -= 75;}
	    if ($g - 75 >= 0)
	      {$g -= 75;}
	    if ($b -75 >= 0)
	      {$b -= 75;}
	    $self->_make_3d(r=>$r, g=>$g, b=>$b, x1=>$x1, x2=>$x2, y1=>$y1, y2=>$y2, gray_lvl=>$gray_lvl);
	  }
	else 
	  {
	    $gd->filledRectangle($x1,$y1, $x2, $y2, $color);
	  }
	$gd->rectangle($x1,$y1, $x2, $y2, $border);
	#$gd->setStyle($black, $black, $black, GD::gdTransparent, GD::gdTransparent);
	#if ($last)
	#  {
	#    my $mid = ($x1-$last)/2+$last;
	#    $gd->line($last, $y1, $mid, 0, GD::gdStyled);
	#    $gd->line($mid, 0, $x1, $y1, GD::gdStyled);
	#xs  }
	#$last = $x2;
      }

    my $x = $self->draw_arrow;
    
    $self->_gd_string(y=>$c-$bh+2, x=>$x, text=>$self->label, size=>$self->block_height-4) if $self->print_label;
  }

sub draw_arrow
  {
    my $self = shift;
    my %opts = @_;
    my $seg = $self->strand =~ /-/ ? $self->segments->[0]:$self->segments->[-1];
    my $black = $self->get_color(1,1,1);
    my $gd = $self->gd;
    my $c = $self->ih()/2;
    my $y = $self->ih-1;
    my $w = ($seg->[1] - $seg->[0]); #width of segment
    use Data::Dumper;
    print STDERR Dumper $seg;
    print STDERR $self->start, "-", $self->stop,"\n";
    my $x = $seg->[0] - $self->start; #start of segment
    
    
    my $gray_lvl = 15;
    my $arrow_width = $self->arrow_width;
    $arrow_width = ($self->stop-$self->start)/10 if (!$arrow_width && $self->start && $self->stop);
    $arrow_width = $self->ih*6 unless $arrow_width;
    $arrow_width = $w if $arrow_width > $w;
#    print STDERR "Arrowhead: X: $x, width: $arrow_width\n";
    my $arrow_end;
    my $gdb = new GD::Image(1,2);
    $gdb->fill(0,0,$gdb->colorResolve(1,1,1));
    $gd->setBrush($gdb);
    unless ($self->no_three_D) {
      my $poly1 = GD::Polygon->new;
      my $poly2 = GD::Polygon->new;
      my @colors = $gd->rgb($self->get_color($self->color));
      my ($r, $g, $b) = ($colors[0],$colors[1],$colors[2]);
      if ($r -75 >= 0)
	     {$r -= 75;}
      if ($g - 75 >= 0)
	     {$g -= 75;}
      if ($b -75 >= 0)
	     {$b -= 75;}
      if ($self->strand =~ /-/i)
	{
	  $self->_make_3d(r=>$r, g=>$g, b=>$b, x1=>$x, x2=>$x+($arrow_width), y1=>0, y2=>$y, gray_lvl=>15);
	  $poly1->addPt($x, 0);
	  $poly1->addPt($x+($arrow_width), 0);
	  $poly1->addPt($x, $c-2);
	  $poly2->addPt($x, $c+2);
	  $poly2->addPt($x, $y);
	  $poly2->addPt($x+($arrow_width), $y);
	  $arrow_end = $x+($arrow_width);
	}
      else
	{
	  #$gd->filledRectangle($x+($w-$arrow_width),0, $x+$w, $y, $self->get_color(0,0,0));
	  $self->_make_3d(r=>$r, g=>$g, b=>$b, x1=>$x+($w-$arrow_width), x2=>$x+$w, y1=>0, y2=>$y, gray_lvl=>15);
	  #$poly->addPt($x+($w-$arrow_width), 0);
	  #$poly->addPt($x+($w-$arrow_width), $y);
	  #$poly->addPt($x+($w), $c);
	  $poly1->addPt($x+($w-$arrow_width), 0);
	  $poly1->addPt($x+$w, 0);
	  $poly1->addPt($x+($w), $c-2);
	  $poly2->addPt($x+($w), $c+2);
	  $poly2->addPt($x+($w), $y);
	  $poly2->addPt($x+($w-$arrow_width), $y);
	  $arrow_end=2;
	}
      #$gd->filledPolygon($poly, $self->get_color($self->color));
      $gd->filledPolygon($poly1, $self->get_color(255,255,255));
      $gd->filledPolygon($poly2, $self->get_color(255,255,255));
    }
    else {
      my $poly = GD::Polygon->new;
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
	  $gd->filledRectangle($x+($w-$arrow_width),0, $x+$w, $y, $self->get_color(0,0,0));
	  $poly->addPt($x+($w-$arrow_width), 0);
	  $poly->addPt($x+($w-$arrow_width), $y);
	  $poly->addPt($x+($w), $c);
	  $arrow_end=2;
	}
      
      my @tmp;
      foreach my $c (@{$self->color})
	{
	  $c-=100;
	  $c = 1 if $c < 1;
	  push @tmp, $c;
	}
      $gd->openPolygon($poly, $self->get_color([@tmp]));
    }
    return $arrow_end;
  }

sub _draw_join
  {
    my $self = shift;
    my $s = shift;
    my $e = shift;
  }
  
sub _make_3d
  {
    my $self = shift;
    my %opts = @_;
    my $gd = $self->gd;
    my $r = $opts{r};
    my $g = $opts{g};
    my $b = $opts{b};
    my ($x1, $x2, $y1, $y2, $gray_lvl) = ($opts{x1},$opts{x2},$opts{y1},$opts{y2},$opts{gray_lvl} || 15);
    my $draw_lines = $y1;
    my $color;
    while ($draw_lines <= $y2)
	{
	  #print STDERR "r: ", $r, ", g: ", $g, ", b: ", $b, "\n";
	  $color = $self->get_color($r,$g,$b);
	  $gd->line($x1, $draw_lines, $x2, $draw_lines, $color);
	  if ($draw_lines / ($y2 - $y1) < 0.34) {
	    if ($r + $gray_lvl <=255)
	     {$r += $gray_lvl;}
	    if ($g + $gray_lvl <=255)
	     {$g += $gray_lvl;}
	    if ($b + $gray_lvl <=255)
	     {$b += $gray_lvl;}
	    }
	  else {
	    if ($r - $gray_lvl >=0)
	     {$r -= $gray_lvl;}
	    if ($g - $gray_lvl >=0)
	     {$g -= $gray_lvl;}
	    if ($b - $gray_lvl >=0)
	     {$b -= $gray_lvl;}
	    }
	  $draw_lines++;
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

