package CoGe::Graphics::Feature;
use strict;
use base qw(Class::Accessor);
use Data::Dumper;
use GD;

BEGIN 
  {
    use vars qw($VERSION $DEFAULT_COLOR $FONTTT $FONT);
    $VERSION     = '0.1';
    $DEFAULT_COLOR = [0,0,0];
    $FONTTT = "/usr/lib/perl5/site_perl/CoGe/fonts/arial.ttf"; #path to true-type font
    $FONT = GD::Font->MediumBold; #default GD font
    __PACKAGE__->mk_accessors(
"DEBUG",
"_gd",
"image_height", "image_width", #generic image size for scaling later by Chromosome.pm
"label", #feature label
"label_location", #location to print label relative to image:  top, bottom, left, right, on.
"start", #chromosomal start position
"stop", #chromosomal stop position
"strand", #top strand ("1") or bottom strand ("-1")
"placement", #will feature be inside or outside the chromosome picture (may default to another place depending on options used and magnification of chromosome) ("in" or "out")
"fill", #should feature fill area (if possible?)
"type", #type of feature (e.g. gene)
"color", #color of feature e.g. [200,255,200]
"bgcolor", #background color of feature e.g.[255,255,255]
"order", #ordering number with which to display features (Features with the same order will be displayed at the same "level" or "track" on the image") 
"description", #feature description
);
  }

sub desc
  {
    my $self = shift;
    return $self->description(@_);
  }

sub ih 
  {
    my $self = shift;
    return $self->image_height(@_);
  }

sub iw
  {
    my $self = shift;
    return $self->image_width(@_);
  }

sub get_color
  {
    my $self = shift;
    my @colors;
    foreach (@_)
      {
	if (ref ($_) =~ /array/i)
	  {
	    push @colors, @$_;
	  }
	else
	  {
	    push @colors, $_;
	  }
      }
    return $self->get_color($DEFAULT_COLOR) unless (@colors >2 && @colors < 5);
    my $gd = $self->gd;
    if (@colors == 4)
      {
	return $gd->colorAllocateAlpha(@colors);
      }
    else
      {
	return $gd->colorResolve(@colors)
      }
  }

sub _gd_string
  {
    my $self = shift;
    my %opts = @_;
    my $text = $opts{text} || $opts{TEXT};
    return 0 unless $text;
    my $x = $opts{x} || $opts{X} || 0;
    my $y = $opts{'y'} || $opts{Y} || 0;
    my $color = $opts{color} || $opts{COLOR};
    my $size = $opts{size} || $opts{SIZE} || $self->padding;
    my $angle = $opts{angle} || $opts{ANGLE} || 0;
    $color = $self->get_color($color);
    my $gd = $self->gd;
    
    if (-r $FONTTT)
      {
	$gd->stringFT($color, $FONTTT, $size, $angle, $x, $y+$size, $text);
      }
    else
      {
	$gd->string($FONT, $x, $y, $text, $color);
      }
  }



sub name
  {
    my $self = shift;
    return $self->label(@_);
  }

sub begin 
  {
    my $self = shift;
    return $self->start(@_);
  }

sub end
  {
    my $self = shift;
    return $self->stop(@_);
  }

sub gd
  {
    my ($self) = shift;
    my %opts = @_;
    my $gd = $self->_gd();
    unless ($gd)
      {
	$self->_initialize(%opts);
	my ($wid, $hei) = ($self->image_width, $self->image_height);
	$gd = new GD::Image($wid, $hei,[1]);
	$self->_gd($gd);
	my $white = $self->get_color(255,255,255);
	$gd->transparent($white);
	$gd->interlaced('true');

	$self->_post_initialize(%opts);
      }
    return $gd;
  }
#this routine is meant to be overloaded by child classes
#this will be called before the GD object is created in order to set the width and height of the GD image
sub _initialize
  {
    my $self = shift;
    my %opts = @_;
  }

#this routine is meant to be overloaded by child classes
#this will be called after teh GD object is created in order to draw feature specific items
sub _post_initialize
  {
    my $self = shift;
    my %opts = @_;
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

CoGe::Graphics:Feature

=head1 SYNOPSIS

  use CoGe::Graphics::Feature;


=head1 DESCRIPTION

=head1 USAGE



=head1 BUGS



=head1 SUPPORT



=head1 AUTHOR

	Eric Lyons


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

