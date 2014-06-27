package CoGe::Graphics::Feature::Link;
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
    use vars qw($VERSION $HEIGHT $EXTERNAL_IMAGE);
    $VERSION     = '0.1';
    $HEIGHT = 100;
    $EXTERNAL_IMAGE = '/opt/apache/CoGe/picts/bike_link.jpg';
    __PACKAGE__->mk_accessors(
"print_label", #flag for printing feature label in gene
);
}

sub add_segment
  {
    my $self = shift;
    my %opts = @_;
    my $start = $opts{start};
    my $stop = $opts{stop};
    $self->start($start) unless defined $self->start;
    $self->start($start) if $start < $self->start();
    $self->stop($stop) unless defined $self->stop;
    $self->stop($stop) if $stop > $self->stop();
  }

sub _initialize
  {
    my $self = shift;
    my %opts = @_;
    $self->image_width(abs($self->stop-$self->start+1)) unless $self->image_width();
    $self->image_height($HEIGHT) unless $self->image_height();
    $self->merge_percent(100);
    $self->bgcolor([255,255,255]) unless $self->bgcolor;
    $self->print_label(0) unless defined $self->print_label();
    $self->skip_duplicate_search(1);
  }

sub _post_initialize
  {
    my $self = shift;
    my %opts = @_;
    my $gd = $self->gd;
    $self->image_width(abs($self->stop-$self->start+1)) unless defined $self->image_width();
    my $bgcolor = $self->get_color($self->bgcolor);
    $gd->fill(0,0, $bgcolor);
#    $gd->transparent($self->get_color($self->bgcolor));
    my $s = $self->start;
    my $color = $self->get_color($self->color);
    my $black = $self->get_color(10,10,10);
    my $y = $self->ih()/2;
    my $x = $self->iw();
    $gd->filledArc($x/3, $y, $x/2-5, $y, 0,360,$color);
    $gd->filledArc($x/3*2, $y, $x/2-5, $y, 0,360,$color);
    $gd->setThickness(5);
    $gd->arc($x/3, $y, $x/2, $y, 0,360,$black);
    $gd->arc($x/3*2, $y, $x/2-5, $y-5, 0,360,$black);
#    $self->_gd_string(y=>$y, x=>10, text=>$self->label, size=>10) if $self->print_label;
#    if (-r $EXTERNAL_IMAGE)
#      {
#	my $ei= GD::Image->new($EXTERNAL_IMAGE);
#	$ei->transparent($ei->colorResolve(255,255,255));
#	my $ex_wid = $ei->width;
#	my $ex_hei = $ei->height;
	#2. copy, resize, and resample the feature onto the new image
#	$self->gd->copyResampled($ei, 0, 0, 0, 0, $self->gd->width, $self->gd->height, $ex_wid, $ex_hei);
#	$self->_gd($ei);
#      }
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
