package CoGe::Graphics::Feature;
use strict;
use base qw(Class::Accessor);
use Data::Dumper;
use GD;
use Carp;
#################### main pod documentation begin ###################
##
##

=head1 NAME

CoGe::Graphics:Feature

=head1 SYNOPSIS

  use CoGe::Graphics::Feature;
  use CoGe::Graphics::Chromosome;

  #create a feature object for use in a chromosome object
  #in this case we are going to create a feature that will fill in the background
  #of the chromosome with a specific color to represent a single nucleotide "A" at
  #a specific chromosomal location
  my $a = CoGe::Graphics::Feature->new();

  #assign a start position for the feature in chromosomal coordinates (usually nucleotides)
  $a->start(50000);

  #you can assign an end if you want.  If you don't, the end will be assigned automatically
  #to equal the start.  In this case, we only want to create a single nucleotide feature
  #so we only need the start position
  #
  #Assign a relative height for the feature.  This will be scaled appropriately when the
  #feature is added to the chromosome object and the chromosome is drawn
  $a->ih(5);

  #assign a relative width for the feature
  $a->iw(5);

  #IMPORTANT: The actual feature image is drawn using a GD object.  Due to the nature
  #of GD, you must specify the height and width of the image before you can initialize
  #a GD object.

  #Features contain a GD object which we can use to draw our feature.  In this case, we
  #just want to fill it with some color.
  #$self->get_color is a method that calls the internal GD object to find the appropriate
  #GD color index for you specified color.  See get_color for more information
  $a->gd->fill(0,0,$a->get_color(200,200, 255));

  #set a label for the feature
  $a->label("A");

  #set the strand for the feature
  $a->strand("1");

  #Specify that this feature will be used to fill in the chromosome
  $a->fill(1);

  #Let's create another feature.  This one will be used to mark an HSP
  #since feature inherits from Class::Accessor, we can use its format to assign
  #values to accessor methods during object creation (See Class::Accessor for more
  #information)
  my $hsp = CoGe::Graphics::Feature->new({start=>49990, stop=>50010});

  #assign image's relative height and width
  $hsp->ih(10);
  $hsp->iw(10);

  #Fill the GD image with some color
  $hsp->gd->fill(0,0,$gd->get_color(250,20, 200));

  #assign a label to the feature
  $hsp->label("HSP1");

  #assign strand
  $hsp->strand("-1");

  #Now that we have a couple of features, we need to add them to a chromosome object
  #create chromosome object.  See CoGe::Graphics::Chromosome for more information regarding
  #this object
  my $c = CoGe::Graphics::Chromosome->new();

  #set the length of the chromosome in chromosomal units (usually nucleotides)
  $c->chr_length(100000)

  #add features to chromosome
  $c->add_feature($a, $hsp);

  #generate image
  $c->generate_png(file=>"./tmp/test.png");

=head1 DESCRIPTION

  In a genomic sense, features are "marked" regions of a genome that have some specific
  attribute defined or known.  For example, a gene is genomic feature that refers to
  a chromosomal region that can be transcribed by a cells transcriptional machinery
  to generate an RNA polymer.  Likewise, a single nucleotide is a feature that describes
  the base composition of a chromosome at a particular region.  The list of genomic
  features is limited to what information you want to attribute to chromosomal regions.

  The CoGe::Graphics::Chromsome object was designed to provide an easy-to-use (TM) and
  flexible way to paint information associated with chromosomes and generate images
  of the chromosome at various levels of magnification.  This object is thus designed
  to provide an easy-to-use an flexible way to create the features to be painted on
  a chromsome object.

  The core work-horse of both of these objects is GD as it provides many flexible routines
  to generating and manipulate images.  By taking advantage of GD's method copyResampled,
  CoGe::Grpahics::Chromosome can take any GD image and resize/resample it to any place.
  This methodology allows one to create any GD image and then use it to represent any
  feature on a chromosome.  In order to facility this process, CoGe::Graphics::Feature
  provides all the necessary routines, accessor methods, and default values, to make
  this process as easy as possible.

  Although this module can be used independently with the Chromosome object as shown in
  example able, its real power lies in a user's ability to use this as a base class
  for another module that generates the GD image in a specific way.  For example, you
  want to generate an image that represents an mRNA that has solid boxes to represent
  exons, an arrow head on one end to indicate directionality, and dotted lines
  connecting exons to indicate introns.  Instead of drawing this image from scratch
  every time, you could create a child class (e.g. CoGe::Graphics::Feature::Gene)
  that does all the nasty work for you each time you want to draw a mRNA.

  CoGe::Graphics::Feature contains two special metods to help with the derivation of
  child classes _initialize and _post_initialize.  _initialize is called BEFORE the gd
  object is created and is a good place to place the functions necessary to figure out
  dimensions of the feature in order to set the image's width and height.
  _post_initalize is called afterwards and is where you can place all the special
  drawing functions for GD.

  For examples of inheriting from CoGe::Graphics::Feature, please see
  CoGe::Graphics::Feature::Gene and CoGe::Graphics::Feature::NucTide.

=head1 USAGE

 use CoGe::Graphics::Feature'
 my $f = CoGe::Graphics::Feature->new()

=head1 BUGS

=head1 SUPPORT

 Please contact Eric Lyons with questions, comments, suggestions, etc.

=head1 AUTHOR

	Eric Lyons
	elyons@nature.berkeley.edu

=head1 COPYRIGHT

Permission to use, copy, modify, and distribute this software and its documentation for educational, research, and not-for-profit purposes, without fee and without a signed licensing agreement, is hereby granted, provided that the above copyright notice, this paragraph and the following two paragraphs appear in all copies, modifications, and distributions. Contact The Office of Technology Licensing, UC Berkeley, 2150 Shattuck Avenue, Suite 510, Berkeley, CA 94720-1620, (510) 643-7201, for commercial licensing opportunities.

The full text of the license can be found in the
LICENSE file included with this module.

=head1 SEE ALSO

perl(1).
CoGe::Graphics::Chromosome
GD

=cut

#################### main pod documentation end ###################
BEGIN
  {
    use vars qw($VERSION $DEFAULT_COLOR $FONTTT $FONT);
    $VERSION     = '0.1';
    $DEFAULT_COLOR = [0,0,0];
    $FONTTT = "/usr/local/fonts/arial.ttf"; #path to true-type font.  Needs to be fixed to get from a conf file
    $FONT = GD::Font->MediumBold; #default GD font
    __PACKAGE__->mk_accessors(
"DEBUG",
"_gd",
"image_height", "image_width", #generic image size for scaling later by Chromosome.pm
"label", #feature label
"label_location", #location to print label relative to image:  top, bottom, left, right, on.
"_start", #chromosomal start position
"_stop", #chromosomal stop position
"strand", #top strand ("1") or bottom strand ("-1")
#"placement", #will feature be inside or outside the chromosome picture (may default to another place depending on options used and magnification of chromosome) ("in" or "out")
"fill", #should feature fill area (if possible?)
"fill_height", #is there a limit to the fill area?
"type", #type of feature (e.g. gene)
"color", #color of feature e.g. [200,255,200]
"bgcolor", #background color of feature e.g.[255,255,255]
"order", #ordering number with which to display features (Features with the same order will be displayed at the same "level" or "track" on the image")
"description", #feature description
"link",  #a place to store an URI.
"_overlap", #place to store the number of regular features that overlap the same region at the same order
"_overlap_pos", #position for placement of overlapping region
"skip_overlap_search", #flag to skip overlap search on some objects (like fill objects)
"skip_duplicate_search", #flag to skip duplicate check on some objects
"merge_percent", #a number between 0-100 for how much to merge the feature onto the chromosome picture.  100 means to completely copy with no merging of colors.
"font_size", #scaling factor for font size.
"external_image", #external image for feature
"force_label", #flag to force the printing of the label
"magnification", #magnification level for increasing/decreasing the size of the feature
"overlay", #level of overlay for features such that one feature is drawn over another.
"use_external_image", #flag for using an external image when available
"options", #place to store additional options
"alt", #used for the alt tag in image maps
"image_coordinates", #used to store the coordinates (xmin, ymin, xmax, ymax) of the feature on the resulting image
"force_draw", #used to make sure that the feature is drawn with a width of 1px, even if it is smaller
);
  }

#################### subroutine header begin ####################

=head2 new

 Usage     : my $f = CoGe::Graphics::Feature->new();
 Purpose   : create a new Feature object for use with a CoGe::Graphics::Feature
           : object.
 Returns   : a new Feature object
 Argument  : a hash reference where keys are accessor function names and values
           : are the values to which they are set
 Throws    :
 Comment   : This new is inherited from Class::Accessor.  Please see for
           : more information

See Also   : Class::Accessor

=cut

#################### subroutine header end ####################

#################### subroutine header begin ####################

=head2 Accessor Methods

 These methods are provided by Class::Accessor and are used to get/set
 a variety of parameters used by the Feature object.  Each method is
 listed and described along with the default values.  The default values
 are often set by CoGe::Graphics::Chromosome->add_feature when the feature
 is added to the Chromosome object

 image_height => relative height of feature that will be used when gd object is created
 alias:  ih

 image_width  => relative width
 alias:  iw

 start        => start position in chromosomal units (usually nucleotides)
 alias:  begin

 stop         => step position in chromsomeal units (set to start if not defined)
 alias:  end

 strand       => (Default: 1) strand on chromosome (1 is top strand, -1 is bottome)

 order        => (auto increments from 1)  The order is the distance from the center
                 of the chromosome that a feature will be drawn.  For example, if set
		 to 1, it will be drawn as close to the center as possible (above
		 the center for the top strand, and below the center for the bottom
		 strand.)  If another feature has an order 2, it will be drawn one
		 "step" further away from the center of the chromosome.  The distance of
		 these steps is the height of the drawn feature plus a padding value.
		 These values are set and used by the Chromosome object.  If two
		 features have the same order value, they will both be drawn the same
		 distance from the center.
 alias:  track

 label        => the label of the feature (e.g. "Gene SuperFoo")
 alias:  name

 label_location => where the label of a feature is to be drawn:  top, bottom, left, right, on

 force_label  => Force the label to be printed (which sometimes is not if the feature is drawn too small)

 force_draw   => Force the drawing of features with a final width of < 1px to be 1px wide.

 description  => the description of the feature
 alias:  desc

 link         => a place to store a URI associated with the feature

 alt          => a place to store the alt tag for generating image maps from CoGe::Graphics::Chromosome

 image_coordinate => a place to store the coordinates (xmin, ymin, xmax, ymax) of the feature on the resulting
                     image generated by CoGe::Graphics::Chromosome

 type         => the type fo the feature (e.g. "gene")

 color        => a place to store a color of the feature.  This is usually important for
                 child classes

 bgcolor      => a place to store a background color of the feature.  Again, this is
                 usually used by child classes

 fill         => (Default 0)If set to true (1), this feature will be used to fill in
                 the background of the chomosome picture at the specified region
		 instead of drawn on top of the image at a particular distance from
		 the center of the chomosome as determined by the value of $self->order
 fill_height  => Percent height of a filled feature.  Valid values 0-1.  Nefault 1;

 skip_overlap_search => When set to true, no overlap search is performed by CoGe::Graphics::Chromosome->
                 add_feature.  This is often useful if there are so many features being searched
		 that the ordering algorithm slows down the performance of the object to intolerable
		 levels.  An example would be if you are adding nucleotide objects to your chromosome
		 and you KNOW that there will never be any overlap between those objects.

 DEBUG        => When true, debugging message will be printed.

 merge_percent=> A value set between 0 and 100 that dictates how the feature is merged with the
                 background image on the chromosome.  100 means to copy it completely over with
                 no merging.  0 means that it won't be seen.  50 means that half of the feature
                 will merge with 50 of the background image.  This is useful for highlighing and
                 blending features that will appear at the same place on the final image.
 alias:  transparency

 font_size    => Scaling factor for font size.  2 will multiple the font size by 2,
                 .5 will multiple the font size by .5 (divide it by 2)

 magnification=> Scaling factor for the size of the feature when used
                 by chromosome.pm for drawing the feature.  This
                 allows for increasing/decreasing the size of some
                 features relative to others
 alias:  mag

 overlay      => This allows features to be drawn over one another.  For example, if I want
                 to draw the CDS region of a gene on top of the mRNA region of a gene.  If
                 two features are drawn at the same order with the same overlay number,
                 they will be drawn next to one another.   However, if they have different
                 overlay numbers, the one with the lowest overlay number will be drawn first
                 followed by the others.
 alias:  layer

 external_image=> This is a place to store another GD object that is linked to an additional image to be
                 drawn in the feature.  An example of use in an image of a nucleotide that is drawn
                 on the background.  For exampe usage, see CoGe::Graphics::Feature::NucTide.pm

 use_external_image    => This flag determines whether or not an external image is used if one is available.

 skip_duplicate_search => this flag will ensure that this feature is not checked for duplication with
                          other features

 options      => place to store additional options for objects derived from this one

 _gd          => Internal place to store the GD object.
 _overlap     => Internal place to track the number of features that occure at the same position.
                 When CoGe::Graphics::Chromosome->add_feature is called, this is set to 1 unless
		 previously defined.  This value is modified by Chromosome->_check_overlap which
		 is called by add_feature.
 _overlap_pos => Internal place to store the relative placement of a feature.  This is initialized
                 by CoGe::Graphics::Chromosome->add_feature and set to one.  Each time an overlap is
		 detected, this value is incremented.  It is later used by Chromosome->_draw_features
		 to determine the relative placement of the feature in the final image.
=cut

#################### subroutine header end ####################

sub new
  {
    my $class = shift;
    my $self = {};
    bless $self , $class;
    my @opts = @_;
    my $opts = ref($opts[0])=~/hash/i ? $opts[0] : {@opts};
        my $start = $opts->{start};
        my $stop = $opts->{stop};
        delete $opts->{start};
        delete $opts->{stop};
        $opts->{_start}=$start if defined $start;
        $opts->{_stop}=$stop if defined $stop;
    while (my ($k, $v) = each %$opts)
     {
       $self->{$k} = $v;
     }
    return $self;
  }

sub start
  {
    my $self = shift;
    my $val = shift;
    return $self->_start($val) if $val;
    return $self->_start;
  }

sub stop
  {
    my $self = shift;
    my $val = shift;
    return $self->_stop($val) if $val;
    return $self->_stop();
  }

#################### subroutine header begin ####################

#=head2 transparency

# Purpose   : alias for $self->merge_percent

#=cut

#################### subroutine header end ####################

sub transparency
  {
    my $self = shift;
    return $self->merge_percent(@_);
  }

#################### subroutine header begin ####################

#################### subroutine header begin ####################

#=head2 layer

# Purpose   : alias for $self->overlay

#=cut

#################### subroutine header end ####################

sub layer
  {
    my $self = shift;
    return $self->overlay(@_);
  }

#################### subroutine header begin ####################

#=head2 track

# Purpose   : alias for $self->order
# Comment   : provides logistical compatibility with other genome visualization packages

#=cut

#################### subroutine header end ####################

sub track
  {
    my $self = shift;
    return $self->order(@_);
  }

#################### subroutine header begin ####################

#=head2 mag

# Purpose   : alias for $self->magnification

#=cut

#################### subroutine header end ####################
sub mag
  {
    my $self = shift;
    return $self->magnification(@_);
  }

#################### subroutine header begin ####################

#=head2 desc

# Purpose   : alias for $self->description

#=cut

#################### subroutine header end ####################

sub desc
  {
    my $self = shift;
    return $self->description(@_);
  }

#################### subroutine header begin ####################

#=head2 ih

# Purpose   : alias for $self->image_height

#=cut

#################### subroutine header end ####################
sub ih
  {
    my $self = shift;
    return $self->image_height(@_);
  }

#################### subroutine header begin ####################

#=head2 iw

# Purpose   : alias for $self->image_width

#=cut

#################### subroutine header end ####################
sub iw
  {
    my $self = shift;
    return $self->image_width(@_);
  }

#################### subroutine header begin ####################

=head2 get_color

 Usage     : my $color_index = $c->get_color([0,0,0]);
 Purpose   : Gets the color index from the internal GD object for your specified color
 Returns   : a GD color index (integer?)
 Argument  : an array or array ref of three to four integers between 0 and 255
 Throws    : this will return the index of the default color $DEFAULT_COLOR
           : if no color is specified or the wrong number of aruguments was passed
 Comment   : If three args are passed, GD->colorResolve is called
           : If four args are passed, the forth is assumed to by an alpha channel
	   : and GD->colorAllocateAlpha is called

See Also   : GD

=cut

#################### subroutine header end ####################
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

#################### subroutine header begin ####################

=head2 _gd_string

 Usage     : $self->_gd_string(text=>$text, x=>$x, y=>$y, $size=>$size);
 Purpose   : generate a string with gd for some text at some position specified
           : by x, y coordinates.
 Returns   : none
 Argument  : hash of key-value pairs where keys are:
             text    =>   text to be printed
             x       =>   x axis coordiate
             y       =>   y axis coordiate
             color   =>   (Optional) an array refof three integers between 1-255
                          This method calls $self->get_color to get the color from GD
                          and will return the default color if none was specified
             size    =>   For true type fonts, this will be the size of the font
             angle   =>   For true type fonts, this will be the angle offset for the font
  Throws    : 0 and a warning if X and Y are not defined
  Comment   : This will check to see if the file is readable as specified by $self->font.
            : If so, it will assume that file to be a true type font and use file in a call to
              GD->stringTF.  Otherwise, it will fallback on the global variable $FONT for the
              default GD font to use (GD::Font->MediumBold)

See Also   : GD

=cut

#################### subroutine header end ####################
sub _gd_string
  {
    my $self = shift;
    my %opts = @_;
    my $text = $opts{text} || $opts{TEXT};
    return 0 unless $text;
    my $x = $opts{x} || $opts{X} || 0;
    my $y = $opts{'y'} || $opts{Y} || 0;
    my $color = $opts{color} || $opts{COLOR};
    my $size = $opts{size} || $opts{SIZE} || 10;
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

#################### subroutine header begin ####################

#=head2 name

# Purpose   : alias for $self->label

#=cut

#################### subroutine header end ####################
sub name
  {
    my $self = shift;
    return $self->label(@_);
  }

#################### subroutine header begin ####################

#=head2 begin

# Purpose   : alias for $self->start

#See Also   :

#=cut

#################### subroutine header end ####################
sub begin
  {
    my $self = shift;
    return $self->start(@_);
  }

#################### subroutine header begin ####################

#=head2 end

# Purpose   : alias for $self->start

#=cut

#################### subroutine header end ####################

sub end
  {
    my $self = shift;
    return $self->stop(@_);
  }

#################### subroutine header begin ####################

=head2 gd

 Usage     : my $gd = $c->gd;
 Purpose   : initializes (if needed) and returns the GD object
 Returns   : GD object
 Argument  : none
 Throws    : none
 Comment   : This checks to see if a gd object has been previously created and stored
           : in $self->_gd.  If not, it creates the GD object using $self->image_width
           : and $self->image_height for dimensions.

 See Also   : GD (which is an excellent module to know if you need to generate images)

=cut

#################### subroutine header end ####################
sub gd
  {
    my ($self) = shift;
    my %opts = @_;
    my $gd = $self->_gd();
    unless ($gd)
      {
	$self->_initialize(%opts);
	my ($wid, $hei) = ($self->image_width, $self->image_height);
	unless ($wid && $hei)
	  {
	    carp "Warning:  image width ($wid) or height ($hei) is not defined\n";
	    ($wid, $hei) = (1,1);
	  }
	$gd = new GD::Image($wid, $hei,[1]);
	$self->_gd($gd);
#	print STDERR ref ($self->_gd),"\n";
	my $white = $self->get_color(255,255,255);
	$gd->transparent($white);
	$gd->interlaced('true');
	$self->_post_initialize(%opts);
      }
    return $gd;
  }
#################### subroutine header begin ####################

=head2 _initialize

 Usage     : $self->_initialize(%args)
 Purpose   : this method is called during the creation of the internal GD object
           : BEFORE the GD object is created.  It is intended to be overloaded by
	   : child classes in order to determine and set the image_height and
	   : image_width that will be used during the creation of the GD object.
 Comment   :
           :

See Also   : For an example of this, see CoGe::Graphics::Feature::Gene

=cut

#################### subroutine header end ####################
#this routine is meant to be overloaded by child classes
#this will be called before the GD object is created in order to set the width and height of the GD image
sub _initialize
  {
    my $self = shift;
    my %opts = @_;

  }

#################### subroutine header begin ####################

=head2 _post_initialize

 Usage     : $self->_post_initialize
 Purpose   : This internal method is called after the creation of the internal
           : GD object and is meant to be overloaded by child classes so custom
	   : GD drawing routines can be implemented

See Also   : See CoGe::Graphics::Feature::Gene for an example of this.

=cut

#################### subroutine header end ####################
#this routine is meant to be overloaded by child classes
#this will be called after teh GD object is created in order to draw feature specific items
sub _post_initialize
  {
    my $self = shift;
    my %opts = @_;
  }

sub feature_info
  {
    my $self = shift;
    return join ("\n", map {$_.": ".$self->{$_}} sort keys %$self);
  }

#################### subroutine header begin ####################

=head2

 Usage     :
 Purpose   :
 Returns   :
 Argument  :
 Throws    :
 Comment   :
           :

See Also   :

=cut

#################### subroutine header end ####################

1;
# The preceding line will help the module return a true value
