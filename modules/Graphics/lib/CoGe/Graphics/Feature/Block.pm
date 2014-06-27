package CoGe::Graphics::Feature::Block;
use strict;
use base qw(CoGe::Graphics::Feature);

=head1 AUTHOR

	Eric Lyons
	elyons@nature.berkeley.edu

=head1 COPYRIGHT

Permission to use, copy, modify, and distribute this software and its documentation for educational, research, and not-for-profit purposes, without fee and without a signed licensing agreement, is hereby granted, provided that the above copyright notice, this paragraph and the following two paragraphs appear in all copies, modifications, and distributions. Contact The Office of Technology Licensing, UC Berkeley, 2150 Shattuck Avenue, Suite 510, Berkeley, CA 94720-1620, (510) 643-7201, for commercial licensing opportunities.

The full text of the license can be found in the
LICENSE file included with this module.

=cut

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
    unless ($self->segments)
      {
	$self->segments([]);
	push @{$self->segments}, [$self->start, $self->stop] if defined $self->start && defined $self->stop;
      }
    foreach my $seg (sort {$a->[0] <=> $b->[0]} @{$self->segments})
      {
	$s = $seg->[0] unless $s;
	$e = $seg->[1] unless $e;
	$s = $seg->[0] if $seg->[0] < $s;
	$e = $seg->[1] if $seg->[1] > $e;
      }
    my $w = $e-$s;
    $w =1 unless $w;
    $self->start($s);
    $self->stop($e);
    $self->image_width($w);
    $self->image_height($h);
    $self->merge_percent(100);
    $self->bgcolor([255,255,255]) unless $self->bgcolor;
    $self->color([255,100,255]) unless $self->bgcolor;
    $self->skip_overlap_search(0);
#    $self->font_size(1);
    $self->block_height($BLOCK_HEIGHT) unless $self->block_height;
    $self->print_label(0) unless defined $self->print_label();
    $self->type('block');
  }

sub _post_initialize
  {
    my $self = shift;
    my %opts = @_;
    $self->label($self->label." (".$self->type.")") if $self->add_type && $self->type;
#    my $label_loc = $self->strand =~ /-/ ? "bot" : "top";
#    $self->label_location($label_loc);
    my $gd = $self->gd;
    $gd->fill(0,0, $self->get_color($self->bgcolor));
#    $gd->transparent($self->get_color($self->bgcolor));
    my $s = $self->start;
    my $black = $self->get_color(0,0,0);
    my $color = $self->get_color($self->color);
    my $last;
    my $c = $self->ih()/2;
    my $bh = $self->image_height/2;
    my @sorted = sort {$a->[0] <=> $b->[0]} @{$self->segments};
    foreach my $seg (@sorted)
      {
	my $x1 = $seg->[0] - $s;
	my $x2 = $seg->[1] - $s;
	$x2 = $x1+1 if $x1 == $x2;
	my $y1 = $c-$bh;
	my $y2 = $c+$bh;
	$gd->filledRectangle($x1,$y1, $x2, $y2, $color);
	$gd->rectangle($x1,$y1, $x2, $y2, $black);
	$gd->setStyle($black, $black, $black, GD::gdTransparent, GD::gdTransparent);
	if ($last)
	  {
	    my $liney = $y1+($y2-$y1)/2;
	    my $mid = ($x1-$last)/2+$last;
	    $gd->line($last, $liney, $mid, 0, GD::gdStyled);
	    $gd->line($mid, 0, $x1, $liney, GD::gdStyled);
	  }
	$last = $x2;
      }
#    $self->_gd_string(y=>$c-$bh+2, x=>$x, text=>$self->label, size=>$self->block_height-4) if $self->print_label;
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
