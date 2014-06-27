package CoGe::Graphics::Feature::Circle;
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
    $HEIGHT = 20;
    $BLOCK_HEIGHT = 20;
    __PACKAGE__->mk_accessors(
"alignment",
"color_matches", #flag for whether to just fill the hsp with a color, or color each match individually
);
}

sub _initialize
  {
    my $self = shift;

    #my %opts = @_;
    my $h = $HEIGHT; #total image height
    my $s = $self->start;
    my $e = $self->stop;
    my $w = $e-$s;

    $self->image_width($w);
    $self->image_height($h);
    $self->merge_percent(100);
    $self->bgcolor([255,255,255]) unless $self->bgcolor;
    $self->skip_duplicate_search(1);
    $self->force_label(1);
    $self->type('Diamond') unless $self->type;
    $self->mag(1);
    $self->font_size(10);
    $self->color_matches(0) unless defined $self->color_matches();

  }

sub _post_initialize
  {
    my $self = shift;

    #my %opts = @_;

    my $w = $self->image_width;
    my $h = $self->image_height;
    my $s = $self->start;
    my $e = $self->stop;;
    my $color=$self->color;

    my @color=split(/\,/,$color);

    my $w2=$w-5;
    my $h2=$h-10;
    my $w3=$w-100;

    my $gd=$self->gd;

    #make a colored box around the hsp
    #$gd->fill(1,1, $self->get_color(@color));
    #$gd->rectangle(1,1,$w-1, $h-1, $self->get_color(@color));

    ##this works to draw a flag
    #$gd->filledRectangle(1,1,$w, $h, $self->get_color(@color));
    #$gd->filledRectangle(0,10,$w2, $h, $self->get_color([255,255,255]));
    ####

    ##This works for a diamond
    $gd->filledRectangle(1,1,$w, $h, $self->get_color([255,255,255]));
    #my $poly1 = GD::Polygon->new;
    #$poly1->addPt(0,10);
    #$poly1->addPt(100,0);
    #$poly1->addPt(200,10);
    #$gd->filledPolygon($poly1, $self->get_color(@color));

    #Trying to draw a circle
    my $poly3 = GD::Circle->new;
    $gd->Circle($gd,"100","10","50",[@color]);

    #my $poly2 = GD::Polygon->new;
    #$poly2->addPt(0,10);
    #$poly2->addPt(100,20);
    #$poly2->addPt(200,10);
    #$gd->filledPolygon($poly2, $self->get_color(@color));
  }

sub _rounded_edges
{
    my $self = shift;
    my %opts = @_;
    my $x1 = $opts{x1};
    my $y1 = $opts{y1};
    my $dist = $opts{dist} || 6;
    my $negx = $opts{negx} || 1;
    my $negy = $opts{negy} || 1;
    my $gd = $self->gd;
    foreach my $i (0..$dist-1)
    {
      #my $poly = GD::Polygon->new;
      #$poly->addPt($x1,$y1);
      #$poly->addPt($x1+($negx*($dist-$i)),$y1);
      #$poly->addPt($x1,$y1+($negy*($i+1)));
      print STDERR "($x1,$y1)\n";
      print STDERR "(",$x1+($negx*($dist-$i)),",$y1)\n";
      print STDERR "($x1,",$y1+($negy*($i+1)),")\n\n";
      #$gd->filledPolygon($poly, $self->get_color(255,255,255));
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
