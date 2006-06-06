package CoGe::Graphics::Chromosome;
use strict;
use base qw(Class::Accessor);
use POSIX;
use Data::Dumper;
use GD;

#################### main pod documentation begin ###################
## Below is the stub of documentation for your module. 
## You better edit it!


=head1 NAME

CoGe::Graphics::Chromosome - Object for drawing chromosomes that provides functionality for painting chromosomes with location based features (such as genes), magnifying/zooming on particular regions, and printing pictures of chromosomes as pngs.

=head1 SYNOPSIS

  use CoGe::Graphics::Chromosome;
  use CoGe::Graphics::Feature::Gene;
  use CoGe::Graphics::Feature::NucTide;

  #create a chromosome object;
  my $c = CoGe::Graphics::Chromosome->new();
  #set the size of the chromosome in nucleotides
  $c->chr_length(250000000);
  #set the point (in nucleotides) where to center the view on the chromosome
  $c->set_point(10000)

  #create a gene feature that will be added to the chromosome
  #this feature object inherits from the base class CoGe::Graphics::Feature which 
  #provides the basic ties for generating features on the chromosome.  The base class
  #can be used alone and custom features designed ad hoc.  Addionally, you can create 
  #a new feature class by inheriting from the base class and designing you own custom
  #drawing routines.  Please see CoGe::Graphics::Feature for details and refer to 
  #CoGe::Graphics::Feature::Gene and others for examples.

  my $f = CoGe::Graphics::Feature::Gene->new();
  #set the strand of the feature
  $f->strand("-1");
  #add some segments for the gene (start and end locations are nucleotide positions)
  $f->add_segment(start=>8000, end=>9000);
  $f->add_segment(start=>9100, end=>9300);
  $f->add_segment(start=>9400, end=>9600);
  $f->add_segment(start=>9700, end=>9800);
  $f->add_segment(start=>10000, end=>10500);
  $f->add_segment(start=>11000, end=>12000);
  #give the feature a label
  $f->label("My special gene")
  #set the color of the feature (an array ref of RGB values where each is between 0 and 255)
  $f->color([255,0,0]); #RED!
  #set the display order of the feature on the chromosome.  "1" is closest to the center of the chromosome
  $f->order(1);

  #add the feature to the chromosome
  $c->add_feature($f);

  #next, let's add some nucleotide sequence data and use the Feature::NucTide object.
  #Remember, this object inherits from CoGe::Graphics::Feature, but has some special 
  #attributes and functionality to draw individual nucleotides in the background image
  #of the chromosome.  Refer to its documentation for more information
  
  #First, we'll need to get some DNA sequence covering the region of interest.  This is
  #NOT a subroutine of this object and is mearly provided to fill in code.
  my $seq = get_dna_sequence(start=>8000, end=>12000); #returns a string of DNA sequence (ATCGTC...)

  my %trans = (A=>'T',
               T=>'A',
       	       C=>'G',
	       G=>'C');

  my $ i = 0;
  foreach my $chr (split //, $seq)
   {
    $chr=uc($chr);
    my $rc_chr = $trans{$chr}
    my $f1 = CoGe::Graphics::Feature::NucTide->new({nt=>$chr, strand=>1, start =>$i+8000});
    my $f2 = CoGe::Graphics::Feature::NucTide->new({nt=>$rc_chr, strand=>-1, start =>$i+8000});
    $c->add_feature($f1, $f2);
    $i++;
   }

  #don't print labels of genes as the CoGe::Graphics::Feature::Gene object will take care of that
  $c->feature_labels(0);

  #turn on the flag for printing labels of the nucleotides (which are "fill" type features)
  $c->fill_labels(1);

  #One of the more powerful (aka "fun") aspects of this object is the idea of zooming in
  #and out on a chromosomal location.
  #By default, there are 10 steps of magnification.  1 is the lowest magnification and 10 is
  #the highest.  
  #You can set the number of steps with the method: $c->num_mag
  #You can set the max number of "units" (in this case nucleotides) seen at the highest 
  #magnification with $c->max_mag($number).  There is a default value if none is specified

  #Let's generate an image for each magnification step and save those to a file.
  foreach my $i (1..10)
   { 
     $c->mag($i);
     $c->generate_png(file=>"tmp/test$i.png");
   }

  #you are finished!


=head1 DESCRIPTION

The overall goal of this object is to create an easy-to-use thingie (TM) for generating images
of chromosomes on which enlightening "features" (genes, functional domains, expression data, 
whatever) are painted with the use of MINIMAL PROGRAMMING and dependencies, yet perserving as
much flexibility as possible so that the final image was fully customizable for advanced and
patient programmers.

Simply put, after working with existing genomics visualization tools and libraries written
in PERL, I felt that something new was needed.  Most of the existing tools were often hard to
work with, required vast knowledge of many other modules, and were rather inflexible towards
customization.  I wanted a module that would allow me to navigate a chromsome with as much ease
as Google Maps(tm) allowed me to navigate my local neighborhood, plop tags at specific
locales, and zoom in on things of interest. 

To this end, I hope this package helps others create views of genomes with the features 
and patterns they find interesting.

The specific aims of this package is to:

1. Create an object that represents a chromsome

2. Allow features (also objects) to easily be added to the chromosome

3. Generate a png of the chromosome at some level of magnification.

4. Gives the user power to customize many aspects of the final image if desired.


=head1 USAGE

 use CoGe::Graphics::Chromsome;
 my $c = CoGe::Graphics::Chromosome->new();

=head1 BUGS



=head1 SUPPORT

Please contact Eric Lyons with questions, comments, suggestions, and most importantly code 
improvements.

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

BEGIN {
    use vars qw($VERSION $DEFAULT_WIDTH $PADDING $DEFAULT_COLOR $MAX_MAG $MAG_SCALE_TYPE $CHR_MAG_HEIGHT $CHR_START_HEIGHT $CHR_INNER_COLOR $CHR_OUTER_COLOR $RULER_COLOR $TICK_COLOR $RULER_HEIGHT $FONT $FONTTT $NUM_MAG $FEATURE_START_HEIGHT $FEATURE_MAG_HEIGHT);
    $VERSION     = '0.1';
    $DEFAULT_WIDTH = 200;  #default image width pixels
    $PADDING = 15; #default value to pad image height
    $CHR_MAG_HEIGHT = 30; #the amount to increase the height of the chromosome picture for each magnification (features increase half this height)
    $CHR_START_HEIGHT = 30; #initial height of the chromosome image in pixels before magnification is added
    $MAX_MAG = 10;     #default number of "units" (e.g. base pairs) to show when maximally zoomed)
    $MAG_SCALE_TYPE = "log"; #default image scaling.  options "log", "linear"
    $DEFAULT_COLOR = [0,0,0];
    $CHR_INNER_COLOR = [220,255,220]; #inner color for chromosome
    $CHR_OUTER_COLOR = [0,0,0]; #border color for chromosome
    $RULER_COLOR = [0,0,255]; #color for measurement ruler
    $TICK_COLOR = [0,0,255]; #color for ticks on measurement ruler
    $RULER_HEIGHT = 20; #height (in pixels) of the ruler
    $FONTTT = "/usr/lib/perl5/site_perl/CoGe/fonts/arial.ttf"; #path to true-type font
    $FONT = GD::Font->MediumBold; #default GD font
    $FEATURE_START_HEIGHT = 5; #the starting height of a feature
    $FEATURE_MAG_HEIGHT = 5; #the magnification of height of a feature with zoom
    $NUM_MAG = 10; #number of magnification steps;
    __PACKAGE__->mk_accessors(
"DEBUG",
"chr_length",
"feature_start_height", #height of a feature
"feature_mag_height", #magnification of a feature with zoom
"draw_ruler", #flag for drawing ruler
"draw_chr_end", #flag for drawing "ends" of chromosome"
"ruler_color", "tick_color", #color for ruler and ticks on ruler respectively
"ruler_height", #height of ruler
"mag_scale_type", "max_mag", "chr_mag_height", "num_mag",
"image_width", "image_height", 
"padding",
"font",
"chr_start_height", #the starting height of the chromosome,
"feature_labels", "fill_labels", #flag to turn off the printing of labels, fill_lables are specifically for filled features;
"draw_chromosome", "chr_inner_color", "chr_outer_color",
"start_picture", #can be set to left, right or center and dictates at which point to 
"_region_start", "_region_stop", #image's start and stop (should be equal to or "larger" than the users
"_magnification", "_mag_scale", 
"_image_h_used", #storage for the amount of the image height used
"_gd", #store GD object
 "_chr_brush",
"_chr_center", "_chr_height", "_chr_h1", "_chr_h2", #interal storage of chromosome image height positions
"_features", #internal storage of features
"_fill_features", #internal storeage of fill features

#"start", "stop", #user defined start and stop.  Not sure if this is needed. . .
);
}

#################### subroutine header begin ####################

=head2 new

 Usage     : my $c = CoGe::Graphics::Chromosome->new()
 Purpose   : Creates a Chromsome object and set up the default parameters
 Returns   : a CoGe::Graphics::Chromosome object
 Argument  : Currently no paramters can be passed in and used to set the defaults. 
             However, you can use the objects Accessor functions to override the defaults
 Throws    : None
 Comment   : This is the mama-jama new.  If you don't know new, then you need to read up
           : on object oriented programming

See Also   : 

=cut

#################### subroutine header end ####################


sub new
{
    my ($class, %parameters) = @_;

    my $self = bless ({}, ref ($class) || $class);
    $self->mag_scale_type($MAG_SCALE_TYPE);
    $self->max_mag($MAX_MAG);
    $self->chr_mag_height($CHR_MAG_HEIGHT);
    $self->chr_start_height($CHR_START_HEIGHT);
    $self->image_width($DEFAULT_WIDTH);
    $self->padding ($PADDING);
    $self->chr_inner_color($CHR_INNER_COLOR);
    $self->chr_outer_color($CHR_OUTER_COLOR);
    $self->draw_chromosome(1);
    $self->draw_ruler(1);
    $self->draw_chr_end(1);
    $self->ruler_height($RULER_HEIGHT);
    $self->ruler_color($RULER_COLOR);
    $self->tick_color($TICK_COLOR);
    $self->num_mag($NUM_MAG);
    $self->feature_start_height($FEATURE_START_HEIGHT);
    $self->feature_mag_height($FEATURE_MAG_HEIGHT);
    $self->mag(5);
    $self->font($FONTTT);
    $self->_features([]);
    $self->_fill_features([]);
    return $self;
}


#################### subroutine header begin ####################

=head2 accessor methods

 These methods are provided by Class::Accessor and are used to get and set a variety of parameters
 used by the Chromosome object.  Each method is listed and described along with the default values
 set during when new is called.  Many of the defaults can be changed easily by looking at the
 BEGIN block of the module and finding the appropriate global variable.

 DEBUG            =>    (DEFAULT: 0) When set to 1, this will cause the object to print debugging 
                         messages
 chr_length       =>    This is used to set the length (usually in nucleotides) of the chromosome.
                        IMPORTANT:  You must set this, otherwise the object will complain.
 draw_chromosome  =>    (DEFAULT: 1) Flag (0 or 1) for whether or not the chromosome is
                        drawn on the final image
 draw_ruler       =>    (DEFAULT: 1) Flag (0 or 1) for whether or not the positional ruler
                        is drawn on the image
 draw_chr_end     =>    (DEFAULT: 1) Flag (0 or 1) for whether or not the rounded ends of the chromosome
                        are drawn where appropriate
 feature_start_height=> (DEFAULT: 5) This stores the starting height (in pixels) of a feature without
                        any increase due to zoom and feature_mag_height.  
 feature_mag_height=>   (DEFAULT: 5) This stores the feature's height in terms of how it 
                        is scaled as the magnification increases.  For example, if this is set
                        to 5 and the magnification is 5, then the resulting height of the feature
                        will be 25 pixels (5*5)

 ruler_color      =>    (DEFAULT: [0,0,255]) Defines the color of the positional ruler.
                        This is the an array reference of three values
                        that corresponds to RGB color values.  Each color ranges from 0-255
 tick_color       =>    (DEFAULT: [0,0,255]) Defines the color of ticks on the  positional ruler.
                        This is the an array reference of three values
                        that corresponds to RGB color values.  Each color ranges from 0-255
 ruler_height     =>    (DEFAULT: 20)  The heigth, in pixels of the positional ruler
 mag_scale_type   =>    (DEFAULT: log) The scaling that is used for the magnification steps.
                        The options are:  linear, log, constant_power

                        Linear scaling means that an equal 
                        numer of positional units (usually nucleotides) are add/removed for each
                        scaling level.  For example, with
                        a chromosome that is 10,000 nucleotides long, a linear scale of 
                        magnification for 10 steps could be:
                        1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000, 1000

                        Log scaling means that there is a logarithmic (or 
                        exponential) scaling between magnification steps.  
                        Log scale of magnification for 4 steps would be:
                        10, 100, 1000, 10000
                        Overall, the log scaling has a better "feel" than the linear scaling

                        Constant_power is designed such that for any chromosome of any length,
                        there is always the name number of chromosomal units at each magnification
                        level.  However, since the length of each chromosome can be different,
                        each chromosome object can have a different number of active magnification
                        levels.
 max_mag          =>    (DEFAULT: 10)  This is the limit of how many positional units (usually
                         nucleotides) are seen at the highest/most powerful magnification
 chr_mag_height   =>    (DEFAULT: 30)  This is the number, in pixels, by which the chromosome 
                        picture increases on the image for each level of magnification.  For 
                        example, if this is set to 30 and the magnification is 5, then the
                        chromosome image would be 150 pixels high.
 chr_start_height =>    (DEFAULT: 30)  This is the number, in pixels, of the starting height of the 
                        chromosome before magnification is applied
 num_mag          =>    (DEFAULT: 10) This is the number of magnification steps available.  If
                        this is set to 10, then there are 10 magnification steps (where the 1, 
                        the lowest magnification, shows a range equal to the lenght of the 
                        chromosome and 10, the highest magnification, shows a range equal to the 
                        value stored in max_mag.
 image_width      =>    (DEFAULT: 200) The width in pixels of the final image.  This value can be
                        modified by the user without undue (aka strange) effects
 image_height     =>    This holds the height of the image and is a value that is calculated
                        dynamically by the module (sub set_image_height) when the image is 
                        generated.  IMPORTANT:  THIS VALUE SHOULD NOT BE MODIFIED BY THE USER
                        DIRECTLY.  One thing to keep in mind is that, the height of the 
                        chromosomal images are dynamic.  This is due to the factors such as the
                        number and scaling aspects of features on the chromosome, how the 
                        scaling of the chromsome image changes with increased magnification, and
                        other factors such as the heigth of the positional ruler.  You may 
                        customize the final height of the image by specifying the scaling factors
                        and heights of the various image parts, but it is not recommended to 
                        change this value as strange(tm) things may happen.
 padding         =>     (DEFAULT: 15) This is the padding (in pixels) used between most items 
                        drawn on the final image.  
 font            =>     (DEFAULT: "/usr/lib/perl5/site_perl/CoGe/fonts/arial.ttf")
                        This is the path to a true-type font used for text labels on the image
 feature_labels  =>     (DEFAULT: 0) Flag used for whether or not to print feature labels.
                        Usually the feature object has already taken care of how to print a label
 fill_labels     =>     (DEFAULT: 1) Flag used for whether or not to print "fill" features labels.
                        A "fill feature" is one that is used to fill in a region on the chromosome
                        and is distinct from regular features.  An example of this would be 
                        the CoGe::Graphics::Feature::NucTide object with, by default, is a fill
                        feature.  This means that when one of these features is drawn, it fills
                        in the background area of the chromosome over the region is covers.  The
                        resulting image will then have individual regions of the chromsome colored
                        according the nucleotide composition and thus generates an easily viewed
                        image.
 chr_inner_color =>     (DEFAULT: [220,255,220]) Defines the interior color of the chromosome.
                        This is the an array reference of three values
                        that corresponds to RGB color values.  Each color ranges from 0-255
 chr_outer_color =>     (DEFAULT: [0,0,0]) Defines the border color of the chromosome.
                        This is the an array reference of three values
                        that corresponds to RGB color values.  Each color ranges from 0-255
 start_picture   =>     (Default: Center)  When setting the region on the chromosome for
                        viewing using $self->set_region($point), the valued contained in
                        this accessor function determines if that point is at the left
                        of the image, the right, or the center.  For example, if this is set
                        to "left" and the $point is equal to 1, the resulting image will
                        start at chomosomal location "1" and extend to the right however
                        many chromosomal units (usually nucleotides) it needs to for the
                        set magnification level.  Valid options for this are "left", "right"
                        and "center"

=cut

#################### subroutine header end ####################


#################### subroutine header begin ####################

=head2 set_region

 Usage     : $c->set_region(start=>$start, stop=>$stop);
 Purpose   : This routine sets the magnification to the appropriate level for viewing the
             selected region as well as the internal accessor functions to track the beginning
             and end of the viewable region.  
 Returns   : none
 Argument  : hash with at least one key-value pair for "start"
             accepts "start", "begin", "START", "BEGIN" to specify the beginning of the region
             accepts "stop", "end", "STOP", "END" to specify the end of the region
             Nominally, the values should be integers the correspond to a chromosomal location.
 Throws    : None
 Comment   : Since this object uses the concept of magnification to set the viewable range
           : on the chromosome, this routine will find the highest magnification that will 
           : encompase the requested region.  For example, if you request to see region
           : 300-400, you may actually see the region from 250-450.
           : 
           : If you only specify a start, this routine will act like set_point
See Also   : set_point

=cut

#################### subroutine header end ####################


sub set_region
  {
    #user sets either a range with a start and stop, or just a start point.
    #if only a start point is set, then we assume that will be the center of the view
    my $self = shift;
    my %opts = @_;
    my $start = $opts{start} || $opts{begin} || $opts{START} || $opts{BEGIN} 
      || $opts{center} || $opts{CENTER} 
	|| $opts{point} || $opts{POINT} || 0;
#	  ||$self->start || 0;
    my $end = $opts{stop} || $opts{end} || $opts{STOP} || $opts{END};# || $self->stop;
#    $self->start($start);
#    $self->stop($end);
    if (defined $start && $end)
      {
	$self->_set_region_start_stop($start, $end);
      }
    else 
      {
	$self->_set_region_for_point($start);
      }
  }

#################### subroutine header begin ####################

=head2 set_point

 Usage     : $c->set_point($position)
 Purpose   : this routine sets the chromosomal location on which to center the viewed region
 Returns   : 0 in case of error
 Argument  : integer that corresponds to a chromsomal location (usually in base_pairs)
 Throws    : Warning if no position is passed in
 Comment   : This routine along with set_region are the two ways in which to specify where
           : on the chromosome you wish to center your viewable image.

See Also   : set_region

=cut

#################### subroutine header end ####################

sub set_point
  {
    my $self = shift;
    my $point = shift;
    unless (defined $point)
      {
	warn "You must specify a chromosomal location!";
	return 0;
      }
    $self->set_region(start=>$point);
  }



#################### subroutine header begin ####################

=head2 add_feature

 Usage     : $c->add_feature($feat_obj);
 Purpose   : This method adds features to the chromosome.
           : Features represent anything that you want to draw on the chromosome and some
           : common examples are genes, mRNAs, nucleotides, etc.
 Returns   : a warning if a feature is not a Feature object (ref($feat) =~ /Feature/)
 Argument  : an array or array ref of CoGe::Graphics::Feature objects or child-class objects
 Throws    : warning
 Comment   : A few defaults will be set in the feature if they haven't been set:
           : strand         => 1
           : fill           => 0
           : stop           => start
           : merge_perecent => 100
           : magnification  => 1
           : overlay        => 1
           : Also, the feature's GD object will be initialized upon import.
	   : There is a check for whether the added feature overlaps other features.  
	   : If so, a counter, $feat->_overlap is incemented in the feature object.
	   : This is later used by the $self->_draw_feature algorithm to figure
	   : out how to best draw overlapping features.
See Also   : CoGe::Graphics::Feature

=cut

#################### subroutine header end ####################


sub add_feature
  {
    my $self = shift;
    my @feats;
    foreach (@_)
      {
	push @feats, ref($_) =~ /array/i ? @$_ : $_;
      }
    foreach my $feat (@feats)
      {
	$feat->strand(1) unless defined $feat->strand;
	$feat->fill(0) unless $feat->fill;
	$feat->merge_percent(100) unless defined $feat->merge_percent;
	$feat->stop($feat->start) unless defined $feat->stop;
	$feat->ih(0) unless defined $feat->ih;
	$feat->iw(0) unless defined $feat->iw;
	$feat->gd; #initialize feature;
	$feat->mag(1) unless defined $feat->mag;
	$feat->overlay(1) unless defined $feat->overlay();
	$feat->_overlap(1) unless $feat->_overlap;#detects overlapping feature on the same track 
	$feat->_overlap_pos(1) unless $feat->_overlap_pos; #placement for overlapping features
	if (ref($feat) =~ /Feature/i)
	  {
	    unless ($feat->order)
	      {
		my $last_feat = $self->get_feats(last_order=>1, strand=>$feat->strand, fill=>$feat->fill);
		my $order = $last_feat ? $last_feat->order()+1 : 1;
		$feat->order($order);
	      }
	    $self->_check_overlap($feat);
	    if ($feat->fill)
	      {
	        push @{$self->_fill_features}, $feat;
  	      }
	    else
	      {
	        push @{$self->_features},$feat;
	      }
 	  }
	else
	  {
	    warn "Feature ($feat) does not appear to be a feature object.  Skipping. . .\n";
	  }
  	
      }
  }


#################### subroutine header begin ####################

=head2 delete_features

 Usage     : $c->delete_features('fill');
 Purpose   : Deletes features from the object.  Either of a type (fill or regular), or all of them
 Returns   : none
 Argument  : string or none
               fill => 1 for deleting fill features
               regular => for deleting regular features
               all (or blank) => deletes all the features
 Throws    : 
 Comment   : 

=cut

#################### subroutine header end ####################


sub delete_features
  {
    my $self = shift;
    my $type = shift;
    my $fill = 1 if !$type || $type =~ /fill/i || $type =~ /all/;
    my $reg = 1 if !$type || $type =~ /regular/i || $type =~ /all/;
    if ($fill || $reg)
      {
	$self->_fill_features([]) if $fill;
	$self->_features([]) if $reg;
      }
    else
      {
	my @feats;
	foreach my $feat ($self->get_features(fill=>0))
	  {
	    push @feats, $feat unless ($feat->type && $feat->type =~ /$type/);
	  }
	$self->_features(\@feats);
      }
  }


#################### subroutine header begin ####################

=head2 _check_overlap

 Usage     : $self->_check_overlap($feature);
 Purpose   : This internal method is called by $self->add_feature in determine if the 
           : being added overlaps another feature on the same strand, order, overlay level, and fill
	   : type.  If so, it increments an internal counter in both features called
	   : _overlap. A positional counter called _overlap_pos is incremented in the feature 
	   : being searched.  This counter is later used by $self->_draw_feature to 
	   : determine the appropriate way to draw the overlapping features
 Returns   : none
 Argument  : a CoGe::Graphics::Feature object
 Throws    : none
 Comment   : this algorithm can get slow with lots of features and doing an overlap search.
           : The overlap search algorithm is a linear search through all previously entered features 
             for any that overlap the newly added feature.  This can probably go faster with a different 
             algo.

See Also   : $self->add_feature();

=cut

#################### subroutine header end ####################

sub _check_overlap
  {
    my $self = shift;
    my $feat = shift;
#    use Benchmark;
#    my $t0 = new Benchmark;    
    return if $feat->skip_overlap_search;
    my @feats = $self->get_feats(strand=>$feat->strand, fill=>$feat->fill, order=>$feat->order);
#    my $t1 = new Benchmark;    
    #@feats = sort {$a->stop <=> $b->stop} @feats;
#    my $t2 = new Benchmark;    
#    my @feats = sort {$a->stop <=> $b->stop} $self->get_feats(strand=>$feat->strand, fill=>$feat->fill, order=>$feat->order);
    return unless @feats;
    #since this can take a while, let's see if we can skip the search by checking if the feature to be checked is outside the bounds of existing features
#    return if ($feat->stop < $feats[0]->start || $feat->start > $feats[-1]->stop);
    foreach my $f (@feats)
    	{
	  next unless $feat->overlay() == $f->overlay();  #skip the check if the features are at different overlay levels.
	  unless ( ($feat->start > $f->stop) || ($feat->stop < $f->start) )
	    {
	      print STDERR "Overlap: ",$feat->name,"\t",$f->name,"\n" if $self->DEBUG; 
	      $feat->_overlap($feat->_overlap+1);
	      $f->_overlap($f->_overlap+1);
	      $feat->_overlap_pos($feat->_overlap_pos+1);
	    }
	}	
#    my $t3 = new Benchmark;
#    my $get_time = timestr(timediff($t1, $t0));
#    my $sort_time = timestr(timediff($t2, $t1));
#    my $find_time = timestr(timediff($t3, $t2));

#     print STDERR qq{
# get_feats:             $get_time
# sort_feats:            $sort_time
# find_time:             $find_time
# };
  }

sub _check_overlap_new
  {
    my $self = shift;
    my $feat = shift;
#    use Benchmark;
#    my $t0 = new Benchmark;    
    return if $feat->skip_overlap_search;
    my @feats = $self->get_feats(strand=>$feat->strand, fill=>$feat->fill, order=>$feat->order);
#    my $t1 = new Benchmark;
    @feats = sort {$a->stop <=> $b->stop} @feats;
    return unless @feats;
#    my $t2 = new Benchmark;
    #since this can take a while, let's see if we can skip the search by checking if the feature to be checked is outside the bounds of existing features
    return if ($feat->stop < $feats[0]->start || $feat->start > $feats[-1]->stop);
    my $start = 0;
    my $stop = $#feats;
    while ($stop-$start > 10)
      {
	my $check = int($stop/2+$start);
	last if ($check == $stop || $check == $start);
#	print STDERR "\tchecking $start - $check - $stop\n";
	if ($feats[$check]->start > $feat->stop)
	  {
	    $stop = $check;
	  }
	elsif ($feats[$check]->stop < $feat->start)
	  {
	    $start = $check;
	  }
	else
	  {
	    last;
	  }
      }
#    my $t3 = new Benchmark;
#    print STDERR "Overlap search $start - $stop\n";
    foreach my $i ($start..$stop)
    	{
	  my $f = $feats[$i];
	  next unless $feat->overlay() == $f->overlay();  #skip the check if the features are at different overlay levels.
	  unless ( ($feat->start > $f->stop) || ($feat->stop < $f->start) )
	    {
	      print STDERR "Overlap: ",$feat->name,"\t",$f->name,"\n" if $self->DEBUG; 
	      $feat->_overlap($feat->_overlap+1);
	      $f->_overlap($f->_overlap+1);
	      $feat->_overlap_pos($feat->_overlap_pos+1);
	    }
	}
#    my $t4 = new Benchmark;
#    my $get_time = timestr(timediff($t1, $t0));
#    my $sort_time = timestr(timediff($t2, $t1));
#    my $pos_time = timestr(timediff($t3, $t2));
#    my $find_time = timestr(timediff($t4, $t3));
#     print STDERR qq{
# get_feats:             $get_time
# sort_feats:            $sort_time
# pos_time:              $pos_time
# find_time:             $find_time
# };
  }

#################### subroutine header begin ####################

=head2 get_features

 Usage     : my @fill_feats = $c->get_features(fill=>1, strand=>1);
 Purpose   : find features that meet specific criteria such as their strand, type, order and fill
 Returns   : an array or array ref based on wantarray
 Argument  : optional hash with the following keys:
           : order => get features that are on that order.  Order is the order by which features
                      are drawn on the chromosome.  order=>1 is for features to be drawn closest 
                      to the center of the chromosome.  order=>2 is for the next layer of 
                      features, etc.
             type  => get features whose type match this value
             strand=> get features from that strand (1, -1, +, -) etc.  
                      this just searches for matching on "-"
             fill  => get "fill features".  Fill features are those feature that are drawn to 
                      "fill in" a region on a chromsome.  An example of this would be a 
                      nucleotide where you would want to color an entire region of the chromosome
                      for a specific nucleotide.
             start => get features that start at this position
             stop  => get features that stop at this position
	     last_order => flag for retrieving only the feature with the highest order
             overlay    => get features at a particular overlay level
 Throws    : none
 Comment   : This is mostly used internally, but is provided in case you want to retrieve a 
           : feature that was previously added

See Also   : CoGe::Graphics::Feature

=cut

#################### subroutine header end ####################


sub get_features
  {
    my $self = shift;
    my %opts = @_;
    my $order = $opts{order} || $opts{ORDER};
    my $type = $opts{type} || $opts{TYPE};
    my $last = $opts{last_order} || $opts{LAST_ORDER}; #flag to get the highest order feature for a location
    my $strand = $opts{strand} || $opts{STRAND};
    my $fill = $opts{fill};
    $fill = $opts{FILL} unless defined $fill; #get filled features?
    my $start = $opts{start};
    my $stop  = $opts{stop};
    my $overlay  = $opts{overlay};
    my @rfeats;
    my @feat_refs;
    push @feat_refs, $self->_fill_features if $fill || !(defined $fill);
    push @feat_refs, $self->_features if (defined $fill && $fill == 0) || !(defined $fill);
    foreach my $ref (@feat_refs)
     {
#       foreach my $feat (sort {$a->overlay <=> $b->overlay} @$ref)
       foreach my $feat (@$ref)
       {
	if ($strand)
	  {
	    if ($strand =~ /-/)
	      {
		next unless $feat->strand =~ /-/;
	      }
	    else
	      {
		next if $feat->strand =~ /-/;
	      }
	  }
	if ($order)
	  {
	    next unless $feat->order eq $order;
	  }
	if ($type)
	  {
	    next unless $feat->type && $feat->type eq $type;
	  }
	if ($start) {next unless $feat->start eq $start;}
	if ($stop) {next unless $feat->stop eq $stop;}
	if ($overlay) {next unless $feat->overlay == $overlay;}
	push @rfeats, $feat
      }
     }
#    @rfeats = sort {$a->order <=> $b->order} @rfeats;

    if ($last)
      {
	return $rfeats[-1];
      }
    return wantarray ? @rfeats : \@rfeats;
  }


#################### subroutine header begin ####################

=head2 get_feature

 Purpose   : alias for get_features

=cut

#################### subroutine header end ####################


sub get_feature
  {
    my $self = shift;
    my %opts = @_;
    return $self->get_features(%opts);
  }

#################### subroutine header begin ####################

=head2 get_feats

 Purpose   : alias for get_features

=cut

#################### subroutine header end ####################


sub get_feats
  {
    my $self = shift;
    my %opts = @_;
    return $self->get_features(%opts);
  }


#################### subroutine header begin ####################

=head2 find_magnification_for_size

 Usage     : my $mag = $c->find_magnification_for_size(10000)
 Purpose   : Find the highest magnification for the size to be viewed
 Returns   : an integer
 Argument  : an integer
 Throws    : Will warn you if the chromosome length has not been specified first.
           : This is needed in order to calculate the magnification scale
 Comment   : When this routine is set, the internal magnification is automatically set
           : with a call to $self->magnification($mag);

See Also   : sub magnification

=cut

#################### subroutine header end ####################


sub find_magnification_for_size
  {
    my $self = shift;
    my $len = shift;
    my $chr_len = $self->chr_length;
    unless ($chr_len)
      {
	warn 'You must set the length of the chromosome before setting the region.\n$self->chr_length($length)';
	return 0;
      }
    my $mag_scale = $self->mag_scale;
    #we want the highest magnification that contains the entire region
    my $mag = $self->magnification;
    foreach my $magt (sort {$a <=> $b} keys %$mag_scale)
      {
	$mag = $magt if $mag_scale->{$magt} >= $len;
      }
#    return $mag;
    $self->magnification($mag);
    return $self->magnification();
  }

#################### subroutine header begin ####################

=head2 find_size_for_magnification

 Usage     : my $size = $c->find_size_for_magnification(5);
 Purpose   : This find the size of the window in chromosomal units (usually nucleotides)
           : 
 Returns   : an integer
 Argument  : an integer or will call $self->magnification to get the the currently set 
           : magnification
 Throws    : none, but will return 0 if the internal magnification scale cannot be defined
           : This requires that the chromsomal length has been previously set (chr_length)
 Comment   : 
           : 

See Also   : 

=cut

#################### subroutine header end ####################


sub find_size_for_magnification
  {
    my $self = shift;
    my $mag = shift || $self->mag();
    my $scale = $self->mag_scale;
    return 0 unless $scale && $mag;
    return $scale->{$mag};
  }

#################### subroutine header begin ####################

=head2 magnification

 Usage     : my $mag = $c->magnification(); ##getting magnification
           : ## OR ##
           : $c->magnification(5); #change magnification
 Purpose   : get, set, and change the magnification of the object
 Returns   : an integer
 Argument  : an integer (optional)
 Throws    : none
 Comment   : When the magnification is changed, the viewable range on the chromosome is
           : changed by a call to the method set_point.

See Also   : 

=cut

#################### subroutine header end ####################


sub magnification
  {
    my $self = shift;
    my $mag = shift;
    #are we changing magnification?  if so, we need to set the region start and end points

    $mag = $self->num_mag if $mag && $mag > $self->num_mag;
    $mag = 1 if defined $mag && $mag < 1;
    $self->_magnification($mag) if $mag;
    if ($self->_region_start && $mag)
      {
	my $point = $self->_region_start;
	$point += ceil (($self->_region_stop - $self->_region_start)/2) if $self->_region_stop;
	$self->set_point($point);
      }
    
    return $self->_magnification();
  }

#################### subroutine header begin ####################

=head2 mag

 Purpose   : alias for $c->magnification

=cut

#################### subroutine header end ####################


sub mag
  {
    my $self = shift;
    return $self->magnification(@_);
  }

#################### subroutine header begin ####################

=head2 generate_png

 Usage     : $c->generate_png(file=>$file_name); #generates THE png by name $file_name
             $c->generate_png();  #generates THE png to STDOUT
 Purpose   : This generates the picture of the chromosome!
 Returns   : none
 Argument  : optional hash where:
             file => is the path for the output png
 Throws    : none
 Comment   : This routine calls the method generate_region to render the image in GD
           : and then calls GD->png to generate the png.  You may wish to generate the 
           : picture in another format or do additional modifications on the GD object.
           : If so, you can call generate_region and then access the gd object direcly.
           : When this routine is finished, the gd object is cleared (set to undef) so
           : that the same object may be used again to generate another image 
See Also   : 

=cut

#################### subroutine header end ####################


sub generate_png
  {
    my $self = shift;
    my %opts = @_;
    my $file_name = $opts{file_name} || $opts{file} || $opts{filename};
    $self->generate_region();
    if ($file_name)
      {
	open (OUT, ">$file_name") || die "Can't open $file_name for writing: $!";
	binmode OUT;
	print OUT $self->gd->png;
	close OUT;
      }
    else
      {
	binmode STDOUT;
	print $self->gd->png;
      }
    $self->_gd(undef);
  }


#################### subroutine header begin ####################

=head2 generate_region

 Usage     : $c->generate_region
 Purpose   : Initializes the GD object and creates the image in GD
 Returns   : none
 Argument  : none
 Throws    : none
 Comment   : This is the necessary method to call to actually initialize the GD object and create
           : the image.  It calls (in order):
             set_image_height (to calculate the height of the image)
             _draw_ruler      (to generate the ruler at the top of the image)
             _draw_chromosome (to generate the chromosome image)
             _draw_features   (add features to the image)
See Also   : 

=cut

#################### subroutine header end ####################


sub generate_region
  {
    my $self = shift;
    my %opts = @_;
    print STDERR "\n", join "\t", "mag: ".$self->mag, "region start: ".$self->_region_start,"region end: ".$self->_region_stop,"\n" if $self->DEBUG;
#    print Dumper $self if $self->DEBUG;
    $self->set_image_height();
    $self->gd->fill(0,0,$self->get_color([255,255,255]));
    $self->_draw_ruler;
    $self->_draw_chromosome;
    $self->_draw_features;
#    $self->gd->fill($self->iw/2,$self->_image_h_used +( $self->ih - $self->_image_h_used)/2+1, $self->get_color($self->chr_inner_color));
  }

#################### subroutine header begin ####################

=head2 generate_imagemap

 Usage     : $c->generate_imagemap
 Purpose   : Generates an image for the features on the chromosome
 Returns   : html image map string
 Argument  : mapname => name of the image map to use the HTML tag: <map name="mapname">             
 Throws    : none
 Comment   : This is also designed to use a javascript function called "change" to display
           : information stored in the feature object description accessor function 
             ($feat->description)in a textarea box.  For example:
           : <script type="text/javascript">
             function change( info ) {	document.info.info.value = info }
             </script>
             <FORM NAME="info"><TEXTAREA NAME="info" cols=50 rows=12></TEXTAREA></FORM>
             Will cause the textarea to change to the features description as the mouse 
             is moved over the feature on the image.
             Also, this usll generate a link in the imagemap using the URI stored in the 
             feature's link accessor function ($feat->link)
             The alt field is the label of the feature.

             Example line from the imagemap:
             <area coords="562,220,663,227" href="annotation_lookup.pl?id=At1g79180.1" onMouseOver="change('Locus: At1g79180.1')"  alt="ID:At1g79180.1">

           : NOTE: Currently skips filled features
See Also   : CoGe::Graphics::Feature

=cut

#################### subroutine header end ####################


sub generate_imagemap
  {
    my $self = shift;
    my %opts = @_;
    my $name = $opts{name}|| $opts{mapname};
    my $html = "<map name=\"$name\">\n";
    my $c = $self->_image_h_used+($self->ih - $self->_image_h_used)/2;
    foreach my $feat ( $self->get_feature(fill=>1), $self->get_features(fill=>0))
      {
	next if $feat->fill;
	#skip drawing features that are outside (by two times the range being viewed) the view
	if ($feat->start)
	  {
	    next if $feat->start > $self->_region_end+2*($self->_region_stop - $self->_region_start );
	  }
	if ($feat->stop > 0)
	  {
	    next if $feat->stop < $self->_region_start-2*($self->_region_stop - $self->_region_start );
	  }
	my $feat_height = ($self->feature_start_height+$self->feature_mag_height*$self->mag);
	my $feat_h = $feat_height/$feat->_overlap;
	my $offset = ($feat->order-1)*($feat_height+$self->padding/1.5)+$self->padding;
	$offset = 0 if $feat->fill;
	$feat_h = ($self->_chr_height-$self->mag-1)/2 if $feat->fill;
	my $y = $feat->strand =~ /-/ ? $c+ $offset+1+($feat_h)*($feat->_overlap_pos-1): $c - $offset-$feat_h*$feat->_overlap_pos;
	my $rb = $self->_region_start;
	my $re = $self->_region_stop;
	my $range = $re-$rb;
	my $w = $self->iw;
	$feat->gd;
	$feat->stop($feat->start) unless defined $feat->stop;
	my $feat_range = $feat->stop-$feat->start;
	my $unit = $self->_calc_unit_size;
	my $fs = $unit*($feat->start-$rb);
	my $fe = $unit*($feat->end-$rb+1);
	my $fw = sprintf("%.1f",$fe - $fs)+1; #calculate the width of the feature;
	
	next if $fw < 1; #skip drawing if less than one pix wide
	my $link = $feat->link;
	my $anno = $feat->description;
	my $alt = $feat->label;
	my $x2 = ceil($fs+$fw);
	my $y2 = ceil($y+$feat_h);
	$y = floor $y;
	$fs = floor($fs);
	$html .= qq{
<area coords="$fs, $y, $x2, $y2"};
	$html .= qq{ href="$link" } if $link;
	$html .= qq{ onMouseOver="change('$anno')"} if $anno;
	$html .= qq{ alt="$alt"} if $alt;
	$html .= ">\n";
      }
    $html .= "</map>\n";
    return $html;
  }

#################### subroutine header begin ####################

=head2 ih

 Purpose   : alias for Accessor method $self->image_height

=cut

#################### subroutine header end ####################


sub ih 
  {
    my $self = shift;
    return $self->image_height(@_);
  }

#################### subroutine header begin ####################

=head2 iw

 Purpose   : alias for Accessor method $self->image_width

=cut

#################### subroutine header end ####################


sub iw
  {
    my $self = shift;
    return $self->image_width(@_);
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
    my ($self) = @_;
    my $gd = $self->_gd();
    unless ($gd)
      {
	my ($wid, $hei) = ($self->image_width, $self->image_height);
	$gd = new GD::Image($wid, $hei,[1]);
	my $white = $gd->colorAllocate(255,255,255);
#	$gd->transparent($white);
	$gd->interlaced('true');
	$self->_gd($gd);
      }
    return $gd;
  }

#################### subroutine header begin ####################

=head2 get_color

 Usage     : my $color_index = $c->get_color([0,0,0]);
 Purpose   : gets the color index from the GD object for your specified color
 Returns   : a GD color index (integer?)
 Argument  : an array or array ref of three to four integers between 0 and 255
 Throws    : this will return the index of the default color $DEFAULT_COLOR if no color
           : is specified or it was passed the wrong number of arguments
 Comment   : If three arguments are passed in, GD->colorResolve is called.
           : If four arguments are passed in, the forth is assumed to be an alpha channel
           : and GD->colorAllocateAlpha is called.
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
    return $gd->colorResolve(@colors)
  }

#################### subroutine header begin ####################

=head2 set_image_height

 Usage     : $c->set_image_height
 Purpose   : This routine figures out how tall the final image will be and sets 
             $self->image_height with that value.  The height of the image depends on
             a number of factors including magnification level, feature height, the number 
             type and placement of features, the height of the positional ruler, the padding
             between picture elements, etc.
 Returns   : none
 Argument  : none
 Throws    : none
 Comment   : This is called internally when generate_region is called.  This must be called
           : before the GD object is created internally by $self->gd in order for the correct
             image height to be passed to the GD object during creation

See Also   : $self->gd

=cut

#################### subroutine header end ####################


sub set_image_height
  {
    #this sub calculates how much height out picture needs based on options and features to be drawn
    my $self = shift;
    my $feat_height = $self->feature_start_height+$self->feature_mag_height*$self->mag;
#    $feat_height = $feat_height*($self->mag);
    my $h=0;# = $self->padding; #give use some padding
    $h += $self->ruler_height;#+$self->padding;
    my $ch = ($self->chr_start_height+$self->mag * $self->chr_mag_height); #chromosome image height
    $self->_chr_height($ch);
    my $chrh = $ch+$self->padding;# if $self->draw_chromosome; #chromosome image height
    my $top_feat = $self->get_feats(last_order=>1, strand=>1, fill=>0);
    my $tfh = $top_feat->order * ($feat_height+$self->padding)if $top_feat;
    $tfh = 0 unless $tfh;
    my $bot_feat = $self->get_feats(last_order=>1, strand=>-1, fill=>0);
    my $bfh = $bot_feat->order * ($feat_height+$self->padding)if $bot_feat;
    $bfh = 0 unless $bfh;
    $h += $tfh > $chrh/2 ? $tfh : $chrh/2;
    $self->_chr_center($h);
    $h += $bfh > $chrh/2 ? $bfh : $chrh/2;
    $h += $self->padding;
#    print STDERR join("\t", $tfh, $bfh, $chrh/2),"\n";
    print STDERR "Image Height: $h\n" if $self->DEBUG;
    $self->ih($h);
    $self->_image_h_used($self->padding);

  }


#################### subroutine header begin ####################

=head2 chr_brush

 Usage     : my $chr_brush = $c->chr_brush
 Purpose   : returns a GD object that is used to generate the border of the chromsome
 Returns   : a GD object
 Argument  : none, but uses $self->chr_outer_color and $self->chr_inner_color to figure out the
             colors needed for the GD image
 Throws    : none
 Comment   : This routine generates a GD object that is used to paint the border of the 
           : chromosome.  It makes a smooth blend from the interior color to the exterior color
           : of the chromosome.  The actual obejct is stored in $self->_chr_brush.  If that
           : exists, that object is returned, otherwise a GD object is create, the image 
           : generated, the object store in $self->_chr_brush, and then returned.
See Also   : Accessor functions $self->chr_outer_color, $self->chr_inner_color

=cut

#################### subroutine header end ####################



sub chr_brush
  {
    my $self = shift;
#    return $self->_chr_brush if $self->_chr_brush;
    my $mag = $self->mag;
    my @bc = @{$self->chr_outer_color()};
    my @ic = @{$self->chr_inner_color()};
    my $edge_brush = new GD::Image(1,$mag+1);
    $edge_brush->colorResolve(255,255,255);
    $edge_brush->line(0,0,0,0,$edge_brush->colorResolve(@bc));
    $edge_brush->line(0,$mag+1,0,$mag+1,$edge_brush->colorResolve(@ic));
    for my $i (1..$mag)
      {
	my @c;
	for my $j (0..2)
	  {
	    push @c, $bc[$j]+($ic[$j]-$bc[$j])/($mag+1)*$i;
	  }
	$edge_brush->line(0,$i,1,$i,$edge_brush->colorResolve(@c));
      }
    $self->_chr_brush($edge_brush);
    return $self->_chr_brush;
  }

#################### subroutine header begin ####################

=head2 mag_scale

 Usage     : my $mag_scale = $c->mag_scale
 Purpose   : The magnification of the object is based on several "steps" of magnification
           : where each step shows a range (or window) of the chromosome.  This routine
           : generates the magnification scale based on the length of the chromosome, 
           : the number of step of magnification as defined by $self->num_mag, and the type
           : of scaling to use (linear, log, or constant_power) as defined by $self->mag_scale_type
           : for generating the magnification scale.  This magnification scale is stored in
           : an Accessor function called $self->_mag_scale and is a hash ref of values where
           : the keys is the magnification step number and the value is the number of 
           : chromosomal units (usually nucleotides) that are visable at the magnification step.
           :
           : If constant_power is used to calculate the scaling effects, $self->num_mag is
           : then set to the number of magnification steps used to cover the full chromsome.
 Returns   : a hash ref where the key is the magnification step number and the value is 
           : the number of chromosomal units (usually nucleotides) that are visable at 
           : the magnification step
 Argument  : none
 Throws    : none
 Comment   : The formulas for determining the magnification steps are
           : linear:  $step = ceil(($self->chr_length-$self->max_mag)/($self->num_mag-1));
           : log   :  $step = 10**(log10($rang)/($self->num_mag));
                      where $rang = $self->chr_length-$self->max_mag;
           : constant_power: $range = $self->max_mag * (2**$step)
           :         This is repeated until $range is larger than $eslf->chr_length
See Also   : 

=cut

#################### subroutine header end ####################


sub mag_scale
  {
    #returns a hash where each key refers to a magnification and the value is the range of viewing (in base pairs)
    my $self = shift;
    my %opts = @_;
    unless ($self->_mag_scale)
      {
	my %scale;
	if ($self->mag_scale_type =~ /lin/i)
	  {
	    my $step = ceil(($self->chr_length-$self->max_mag)/($self->num_mag-1));
	    for my $i (1..$self->num_mag)
	      {
		my $rang = $self->chr_length - $i * $step + $step;
		$rang = $self->max_mag if $rang < $self->max_mag;
		$rang = $self->chr_length + $step if $rang >= $self->chr_length;
		$scale{$i} = $rang;
	      }
	  }
	elsif ($self->mag_scale_type =~ /log/i)
	  {
	    my $rang = $self->chr_length-$self->max_mag;;
	    my $step = 10**(log10($rang)/($self->num_mag));
	    print STDERR "Log Magnification Step: $step, Num Mags: ", $self->num_mag,"\n" if $self->DEBUG;
	    $scale{1} = $self->chr_length;
	    for my $i (2..$self->num_mag-1)
	      {
		$rang = ceil ($rang)/$step;
		#$rang = $self->max_mag if $rang < $self->max_mag;
		$scale{$i} = $rang+$self->max_mag;
		$scale{$i} = $self->max_mag if $scale{$i} < $self->max_mag;
	      }
	    $scale{$self->num_mag} = $self->max_mag;
	  }
	elsif ($self->mag_scale_type =~ /constant_power/i)
	  {
	    my $i = 1;
	    my $range = $self->max_mag;
	    my %tmp;
	    while ($range < $self->chr_length)
	      {
		$tmp{$i}=$range;
		$range = $self->max_mag*(2**$i);
		$i++;
	      }
	    $tmp{$i} = $range;
	    $self->num_mag($i);
	    $self->mag($i) if $i < $self->mag;
	    foreach my $key (keys %tmp)
	      {
		$scale{$i-$key+1} = $tmp{$key};
	      }
	  }
	$self->_mag_scale(\%scale);
      }
    return $self->_mag_scale;
  }



#################### subroutine header begin ####################

=head2 _draw_chromosome

 Usage     : $c->_draw_chromosome
 Purpose   : this internal routine draws the chromosome picture if $self->draw_chromosome is true
 Returns   : none
 Argument  : none
 Throws    : none
 Comment   : This routine will generate the chromsome background picture as well as calculating
           : the center and height of the chromosome picture.  These latter values are important
           : for the drawing of features.  This method is called internally by 
           : $self->generate_region
See Also   : 

=cut

#################### subroutine header end ####################


sub _draw_chromosome
  {
    my $self = shift;
    my $gd = $self->gd;
    my $ch = $self->_chr_height;
    $ch = $ch/2; #want half the height for the rest of the calcs
    my $hc = $self->_image_h_used+($self->ih-$self->_image_h_used)/2; #height of center of chromsome image
    my $w = $self->iw; #width of image
    my $xs = $self->_region_start < 1 ? $w*abs(1-$self->_region_start)/($self->_region_stop - $self->_region_start): 0;
    my $xe = $self->_region_stop > $self->chr_length ? $w-$w*($self->_region_stop - $self->chr_length-1)/($self->_region_stop - $self->_region_start) : $w;
    print STDERR "\nChromosome image: Height/2: $ch, Height Center: $hc, xs: $xs, xe:\n" if $self->DEBUG;
    
    
    return unless $self->draw_chromosome; #do we draw the chromsome picture?
    my $brush = $self->chr_brush;
    $gd->setBrush($brush);
    $gd->line($xs, $hc-$ch, $xe, $hc-$ch, gdBrushed) unless $xe < 0;
    $brush->flipVertical();
    $gd->setBrush($brush);
    $gd->line($xs, $hc+$ch, $xe, $hc+$ch, gdBrushed) unless $xe < 0;
    my $color = $self->get_color($self->chr_outer_color);
    $gd->setStyle($color, $color, $color, GD::gdTransparent, GD::gdTransparent);
    $gd->line($xs, $hc, $xe, $hc, gdStyled) unless $xe < 0;
    $self->_draw_chr_end (x=>$xs, dir=>"left", 'y'=>$hc) if ($xs > 0) && $self->draw_chr_end;
    $self->_draw_chr_end (x=>$xe, dir=>"right",'y'=> $hc) if ($xe < $w) && $self->draw_chr_end;
    
  }

#################### subroutine header begin ####################

=head2 _draw_chr_end

 Usage     : $c->_draw_chr_end (x=>$x_pos, dir=>"left", y=>$y_pos) 
 Purpose   : this internal method draws a semi-circle end to the chromosome picture (if needed)
 Returns   : none
 Argument  : hash of key-value pairs where:
             x   => is the x coordinate and where the open half of the semi-circle should lie
             y   => is the y coordinate and the center of the chromosome
             dir => ('left' or 'right') for which end of the chromosome this will lie
 Throws    : none
 Comment   : this is called internall by $self->_draw_chromosome
           : 

See Also   : 

=cut

#################### subroutine header end ####################


sub _draw_chr_end
  {
    my $self = shift;
    my %opts = @_;
    my $x = $opts{x};
    my $dir = $opts{dir} || "left";
    my $y = $opts{'y'}; #$ch/2+$self->_image_h_used+$self->mag/2; #height of center of chromsome image
    my $gd = $self->gd;
    my $ch = $self->_chr_height;
    
    my @arc1 = $dir =~ /left/i ? (90, 180) : (0, 90);
    my @arc2 = $dir =~ /left/i ? (180, 270) : (270, 0);
    my @arc3 = $dir =~ /left/i ? (90, 270) : (270,90);
    $gd->filledArc($x, $y, $ch, $ch, @arc3, $self->get_color($self->chr_inner_color));
    $self->chr_brush->flipVertical;
    $gd->arc($x, $y, $ch, $ch, @arc1, gdBrushed);
    $self->chr_brush->flipVertical;
    $gd->setBrush($self->chr_brush);
    $gd->arc($x, $y, $ch, $ch, @arc2, gdBrushed);
  }


#################### subroutine header begin ####################

=head2 _draw_features

 Usage     : $self->_draw_features
 Purpose   : this routine parses all the feature objects store internally, determines all the
           : necessary positional information for where they are to be drawn, and sents the
           : information to $self->_draw_feature for rendering
 Returns   : none
 Argument  : none
 Throws    : none
 Comment   : Has an internal method to skip rendering features that are not within the visable
           : window

See Also   : $self>_draw_feature for individual feature rendering
           : $self->add_feature for adding features
           : $self->get_feature for retrieving features

=cut

#################### subroutine header end ####################


sub _draw_features
  {
    my $self = shift;
    my $c = $self->_image_h_used+($self->ih - $self->_image_h_used)/2;
    print STDERR "Image used: ".$self->_image_h_used."  Image Height: ".$self->ih."  Center: $c\n" if $self->DEBUG;
    foreach my $feat ( $self->get_feature(fill=>1), sort {$a->overlay <=> $b->overlay} $self->get_features(fill=>0))
      {
	#skip drawing features that are outside (by two times the range being viewed) the view
	if ($feat->start)
	  {
	    next if $feat->start > $self->_region_end+2*($self->_region_stop - $self->_region_start );
	  }
	if ($feat->stop > 0)
	  {
	    next if $feat->stop < $self->_region_start-2*($self->_region_stop - $self->_region_start );
	  }
	my $feature_height = ($self->feature_start_height+$self->feature_mag_height*$self->mag);
	my $feat_h = $feature_height/$feat->_overlap;#*$feat->mag;
	my $offset = ($feat->order-1)*($feature_height+$self->padding/1.5)+$self->padding;
	$offset = 0 if $feat->fill;
	$feat_h = ($self->_chr_height-$self->mag-1)/2 if $feat->fill;
	my $feat_mag_offset = ($feat_h*$feat->mag - $feat_h)/2;
	my $y = $feat->strand =~ /-/ ? 
	  $c + $offset + $feat_h  * ($feat->_overlap_pos-1) - $feat_mag_offset+ 1  :
	  $c - $offset - $feat_h  *  $feat->_overlap_pos    - $feat_mag_offset;
	$feat_h *= $feat->mag;
        my $sy;
	if ($feat->fill)
	  {
	    $sy = $feat->strand =~ /-/ ? $c+2 : $c-$self->padding;
	  }
	elsif ($feat->label_location && $feat->label_location =~ /bot/)
	  {
	    $sy = $y+$feature_height+1;
	  }
	elsif ($feat->label_location && $feat->label_location =~ /top/)
	  {
	    $sy = $y;
	  }
	$self->_draw_feature_fast(feat=>$feat, 'y'=>$y, ih=>$feat_h, 'sy'=>$sy);
      }
  }

#################### subroutine header begin ####################

=head2 _draw_feature

 Usage     : $self->_draw_feature(feat=>$feat, 'y'=>$y, ih=>$feat_h, 'sy'=>$sy);
 Purpose   : draws a feature at specific y axis position with a particular height
 Returns   : 0 if a valid feature object was not specified
 Argument  : hash of key value pairs where the keys are:
              feat         => (or "FEAT", "f") a CoGe::Graphics::Feature object
              image_height => (or "ih", "IH", or takes height from feature object) the height
                              at which the feature will be drawn
              y            => (or "Y") the y axis position from which the feature will be drawn
              string_y     => (or "sy") the y axis position from which the feature label will
                              be drawn if the chromosome object permits the drawing of labels.
                              (As determined from $self->feature_labels and $self->fill_labels)
 Throws    : 0 if a valid feature object was not specified
 Comment   : This uses GD->copyResampled to resample the gd image from the feature object onto
           : the chromosome gd objects.  The feature height is determined by either a specified
             parameter or by the feature object.  The width of the feature is calculated based 
             on the chromosomal location of the feature (usually in nucleotides).  Together
             this easily allows for the generation of a feature image on the chromosome image
             that scales smoothly at the requested magnification.  This routine is called by
             $self->_draw_features

See Also   : $self->_draw_features

=cut

#################### subroutine header end ####################


sub _draw_feature_slow
  {
    my $self = shift;
    my %opts = @_;
    my $feat = $opts{feat} || $opts{FEAT} || $opts{f};
    return 0 unless ref($feat) =~ /Feature/i;
    my $y = $opts{'y'} || $opts{Y};
    my $ih = $opts{'image_height'} || $opts{'ih'} || $opts{'IH'} || $feat->ih;
    my $sy = $opts{'string_y'} || $opts{'sy'};#label y axis
    my $rb = $self->_region_start;
    my $re = $self->_region_stop;
    my $range = $re-$rb;
    my $w = $self->iw;
    $feat->gd;
    $feat->stop($feat->start) unless defined $feat->stop;
    my $feat_range = $feat->stop-$feat->start;
    my $unit = $self->_calc_unit_size;
    my $fs = $unit*($feat->start-$rb);
    my $fe = $unit*($feat->end-$rb+1);
    my $fw = sprintf("%.1f",$fe - $fs)+1; #calculate the width of the feature;
    return if $fw < 1; #skip drawing if less than one pix wide
    print STDERR "Drawing feature ".$feat->label." Order: ".$feat->order." Overlap: ".$feat->_overlap." : ", $feat->start, "-", $feat->end," Dimentions:",$fw,"x",$ih, " at position: $fs,$y"."\n" if $self->DEBUG;
    if ($feat->fill)
      {
	$self->gd->copyResampled($feat->gd, $fs, $y,0,0, $fw, $ih, $feat->iw, $feat->ih);
	if ($feat->external_image && $fw > 10) #if we have an external image and the feature width is greater than 10. . .
	  {
	    my $ei = $feat->external_image;
#	    $ei->transparent($ei->colorResolve(255,255,255));
	    my $ex_wid = $ei->width;
	    my $ex_hei = $ei->height;
	    my $scale = $fw/$ex_wid; #scaling factor for external image
	    my $hei = $feat->strand =~ /-/ ? $y : $y+$ih-($ex_hei*$scale);
	    #okay, we are going to need to do some fancy stuff in order to smoothly resize and paste the feature onto the main image
	    #1. create a blank image of the appropriate size
	    my $newgd = GD::Image->new ($fw, $ex_hei*$scale,[1]);
	    #2. copy, resize, and resample the feature onto the new image
	    $newgd->copyResampled($ei, 0, 0, 0, 0, $newgd->width, $newgd->height, $ex_wid, $ex_hei);  
	    #3. find any colors that are close to white and set them to white.
	    my $max = 200;
	    for my $x (0..$fw)
	      {
		for my $y (0..$ih)
		  {
		    my ($r, $g, $b) = $newgd->rgb($newgd->getPixel($x, $y));
		    if ($r > $max && $g > $max && $b > $max)
		      {
			$newgd->setPixel($x, $y, $newgd->colorResolve(255,255,255));
		      }
		  }
	      }
	    #4. make white transparent
	    $newgd->transparent($newgd->colorResolve(255,255,255));
	    #5. copy new image into appropriate place on image.
#	    $self->gd->copyMerge($newgd, $fs, $hei, 0, 0, $fw, $ih, $feat->merge_percent);
 	    $self->gd->copy($newgd, $fs, $hei, 0, 0, $newgd->width, $newgd->height);
	  }
      }
    else
      {
	#okay, we are going to need to do some fancy stuff in order to smoothly resize and paste the feature onto the main image
	#1. create a blank image of the appropriate size
	my $newgd = GD::Image->new ($fw, $ih,[1]);
	#2. copy, resize, and resample the feature onto the new image
	$newgd->copyResampled($feat->gd, 0, 0, 0, 0, $fw, $ih, $feat->iw, $feat->ih);  
	#3. find any colors that are close to white and set them to white.
	my $max = 240;
	for my $x (0..$fw)
	  {
	    for my $y (0..$ih)
	      {
		my ($r, $g, $b) = $newgd->rgb($newgd->getPixel($x, $y));
		if ($r > $max && $g > $max && $b > $max)
		  {
		    $newgd->setPixel($x, $y, $newgd->colorResolve(255,255,255));
		  }
 	      }
	  }
	#4. make white transparent
	$newgd->transparent($newgd->colorResolve(255,255,255));
	#5. copy new image into appropriate place on image.
	$self->gd->copyMerge($newgd, $fs, $y, 0, 0, $fw, $ih, $feat->merge_percent);
      }
    
    my $size;
    if ($self->fill_labels && $feat->fill) {$size=$fw >= 15 ? 15 : $fw;}
    elsif ($self->feature_labels && defined $feat->label) 
      {
        $size = $ih > 13 ? 13 : $ih; 
	$size=$size/2 if $fw <$size * (length $feat->label)/1.5;
	#print STDERR $feat->label,": $fw, $size\n";
        $sy=$y+$ih/2-$size/2 unless $sy;
	my $adjust = 0;
	$adjust = $fw/10;
	$fs+=$adjust;
      }
    $size = $size*$feat->font_size if $size && $feat->font_size;

    $self->_gd_string(y=>$sy, x=>$fs, text=>$feat->label, size=>$size) if ( ($self->feature_labels || $self->fill_labels)&& ($fw>5 || $feat->force_label)); #don't make the string unless the feature is at least 5 pixels wide
  }
sub _draw_feature_fast
  {
    my $self = shift;
    my %opts = @_;
    my $feat = $opts{feat} || $opts{FEAT} || $opts{f};
    return 0 unless ref($feat) =~ /Feature/i;
    my $y = $opts{'y'} || $opts{Y};
    my $ih = $opts{'image_height'} || $opts{'ih'} || $opts{'IH'} || $feat->ih;
    my $sy = $opts{'string_y'} || $opts{'sy'};#label y axis
    my $rb = $self->_region_start;
    my $re = $self->_region_stop;
    my $range = $re-$rb;
    my $w = $self->iw;
    $feat->gd;
    $feat->stop($feat->start) unless defined $feat->stop;
    my $feat_range = $feat->stop-$feat->start;
    my $unit = $self->_calc_unit_size;
    my $fs = $unit*($feat->start-$rb);
    my $fe = $unit*($feat->end-$rb+1);
    my $fw = sprintf("%.1f",$fe - $fs)+1; #calculate the width of the feature;
    return if $fw < 1; #skip drawing if less than one pix wide
    print STDERR "Drawing feature ".$feat->label." Order: ".$feat->order." Overlap: ".$feat->_overlap." : ", $feat->start, "-", $feat->end," Dimentions:",$fw,"x",$ih, " at position: $fs,$y"."\n" if $self->DEBUG;
    if ($feat->fill)
      {
	$self->gd->copyResampled($feat->gd, $fs, $y,0,0, $fw, $ih, $feat->iw, $feat->ih);
	if ($feat->external_image && $fw > 10) #if we have an external image and the feature width is greater than 10. . .
	  {
	    my $ei = $feat->external_image;
#	    $ei->transparent($ei->colorResolve(255,255,255));
	    my $ex_wid = $ei->width;
	    my $ex_hei = $ei->height;
	    my $scale = $fw/$ex_wid; #scaling factor for external image
	    my $hei = $feat->strand =~ /-/ ? $y : $y+$ih-($ex_hei*$scale);
	    #okay, we are going to need to do some fancy stuff in order to smoothly resize and paste the feature onto the main image
	    #1. create a blank image of the appropriate size
	    my $newgd = GD::Image->new ($fw, $ex_hei*$scale,[1]);
	    #2. copy, resize, and resample the feature onto the new image
	    $newgd->copyResampled($ei, 0, 0, 0, 0, $newgd->width, $newgd->height, $ex_wid, $ex_hei);  
	    #3. find any colors that are close to white and set them to white.
	    my $max = 200;
	    for my $x (0..$fw)
	      {
		for my $y (0..$ih)
		  {
		    my ($r, $g, $b) = $newgd->rgb($newgd->getPixel($x, $y));
		    if ($r > $max && $g > $max && $b > $max)
		      {
			$newgd->setPixel($x, $y, $newgd->colorResolve(255,255,255));
		      }
		  }
	      }
	    #4. make white transparent
	    $newgd->transparent($newgd->colorResolve(255,255,255));
	    #5. copy new image into appropriate place on image.
#	    $self->gd->copyMerge($newgd, $fs, $hei, 0, 0, $fw, $ih, $feat->merge_percent);
 	    $self->gd->copy($newgd, $fs, $hei, 0, 0, $newgd->width, $newgd->height);
	  }
      }
    else
      {
# 	#okay, we are going to need to do some fancy stuff in order to smoothly resize and paste the feature onto the main image
# 	#1. create a blank image of the appropriate size
 	my $newgd = GD::Image->new ($fw, $ih,[1]);
# 	#2. copy, resize, and resample the feature onto the new image
 	$newgd->copyResized($feat->gd, 0, 0, 0, 0, $fw, $ih, $feat->iw, $feat->ih);  
# 	#3. find any colors that are close to white and set them to white.
# 	my $max = 240;
# 	for my $x (0..$fw)
# 	  {
# 	    for my $y (0..$ih)
# 	      {
# 		my ($r, $g, $b) = $newgd->rgb($newgd->getPixel($x, $y));
# 		if ($r > $max && $g > $max && $b > $max)
# 		  {
# 		    $newgd->setPixel($x, $y, $newgd->colorResolve(255,255,255));
# 		  }
#  	      }
# 	  }
# 	#4. make white transparent
 	$newgd->transparent($newgd->colorResolve(0,0,0));
# 	#5. copy new image into appropriate place on image.
 	$self->gd->copyMerge($newgd, $fs, $y, 0, 0, $fw, $ih, $feat->merge_percent);
#	$self->gd->copyResized($feat->gd, $fs, $y, 0, 0, $fw, $ih, $feat->iw, $feat->ih);
      }
    
    my $size;
    if ($self->fill_labels && $feat->fill) {$size=$fw >= 15 ? 15 : $fw;}
    elsif ($self->feature_labels && defined $feat->label) 
      {
        $size = $ih > 13 ? 13 : $ih; 
	$size=$size/2 if $fw <$size * (length $feat->label)/1.5;
	#print STDERR $feat->label,": $fw, $size\n";
        $sy=$y+$ih/2-$size/2 unless $sy;
	my $adjust = 0;
	$adjust = $fw/10;
	$fs+=$adjust;
      }
    $size = $size*$feat->font_size if $size && $feat->font_size;

    $self->_gd_string(y=>$sy, x=>$fs, text=>$feat->label, size=>$size) if ( ($self->feature_labels || $self->fill_labels)&& ($fw>5 || $feat->force_label)); #don't make the string unless the feature is at least 5 pixels wide
  }

#################### subroutine header begin ####################

=head2 _calc_unit_size

 Usage     : my $unit = $self->_calc_unit_size();
 Purpose   : returns the width, in pixels, of one chromsomal unit (usually nucleotides)
 Returns   : an number
 Argument  : none
 Throws    : none
 Comment   : formula is image_width/visable_region_size
           : 

See Also   : 

=cut

#################### subroutine header end ####################


sub _calc_unit_size
  {
    my $self = shift;
    return ($self->iw/($self->_region_stop-$self->_region_start));
    return sprintf("%.5f",($self->iw/($self->_region_stop-$self->_region_start)));
  }

#################### subroutine header begin ####################

=head2 _draw_ruler

 Usage     : $self->_draw_ruler;
 Purpose   : generates the positional ruler at the top of the image
 Returns   : none
 Argument  : none
 Throws    : none
 Comment   : 
           : called by $self->generate_region

See Also   : 

=cut

#################### subroutine header end ####################


sub _draw_ruler
  {
    my $self = shift;
    my $gd = $self->gd;
    my $c = $self->ruler_height/2+$self->_image_h_used; #center of ruler
    $self->_image_h_used($self->_image_h_used + $self->ruler_height+$self->padding/2);
    return unless $self->draw_ruler;
    my $w = $self->iw; #width of image
    my $xb = $self->_region_start < 1 ? $w*abs(1-$self->_region_start)/($self->_region_stop - $self->_region_start): 0; #x position begin
    my $xe = $self->_region_stop > $self->chr_length ? $w-$w*($self->_region_stop - $self->chr_length-1)/($self->_region_stop - $self->_region_start) : $w; #x position end
    return unless ($xe>0);
    $gd->line($xb,$c,$xe,$c,$self->get_color($self->ruler_color));
    my $mtyb = $c-$self->ruler_height/2; #major tick y begin
    my $mtye = $c+$self->ruler_height/2; #major tick y end
    my $styb = $c-$self->ruler_height/4; #minor tick y begin
    my $stye = $c+$self->ruler_height/4; #minor tick y end
    my $rb = $self->_region_start;
    $rb = 1 if $rb < 1;
    my $re = $self->_region_end;
    $re = $self->chr_length if $re > $self->chr_length;
    my $range = $re-$rb+1; #chromosomal positional range
    $rb = $rb-($range/10); #back up a bit
    my $div = "1"."0"x int log10($self->find_size_for_magnification); #determine range scale (10, 100, 1000, etc)
    print STDERR "\nRULER: Center: $c, Start $xb, Stop: $xe, Ticks: $div, \n" if $self->DEBUG;
    $self->_make_ticks(scale=>$div*10, y1=>$mtyb, y2=>$mtye, range_begin=>$rb, range_end=>$re,text_loc=>1);
#    $div /= 10;
    $self->_make_ticks(scale=>$div, y1=>$styb, y2=>$stye, range_begin=>$rb, range_end=>$re, text_loc=>-1);

  }

#################### subroutine header begin ####################

=head2 _make_ticks

 Usage     : $self->_make_ticks(scale=>1000, y1=>$y1, y2=>$y2, range_begin=>1, range_end=>10000)
 Purpose   : generate tick marks of height ($y2-$y1) using the range positions and the scale
           : to calcualte where the ticks should be.
 Returns   : none
 Argument  : hash of key-value pairs where keys are:
              y1   => starting y axis position of tick
              y2   => ending y axis position of tick
              range_begin => number representing the starting point of the ruler
              range_end   => number representing the ending poitn of the ruler
              scale       => points along the range at which to generate a tick mark
 Throws    : 0
 Comment   : This method will convert numbers drawn at the tick marks to use 
              "K" if the number ends in 000
              "M" if the nubmer ends in 000000
              "G" if the number ends in 000000000
           : This method also generates ticks that are the size of one chromsomal unit 
             (usually nucleotides) if the magnification is high enough.
           : This method is called by $self->_generate_ruler
See Also   : 

=cut

#################### subroutine header end ####################


sub _make_ticks
  {
    my $self = shift;
    my %opts = @_;
    my $div = $opts{scale};
    my $y1 = $opts{y1}; #top of tick
    my $y2 = $opts{y2}; #bottom of tick
    my $rb = $opts{range_begin};
    my $re = $opts{range_end};
    my $text_loc  = $opts{text_loc}; #location of text -- flag (0 off, 1 above line, -1 below line, default 1)
    $text_loc = 1 unless defined $text_loc;
    my $tick = $div;
    my $w = $self->iw;
    my $gd = $self->gd;
    #find first tick regional posistion
    while ($tick < $rb){$tick += $div;}
    #make ticks
    my $unit = $self->_calc_unit_size;
    print STDERR "UNIT Size: $unit\n" if $self->DEBUG;
    if ($rb <= 1 && $re >= 1)
      {
	my $x = $w *(1 - $self->_region_start)/($self->_region_stop - $self->_region_start);
	$gd->filledRectangle($x, $y1, $x+$unit, $y2, $self->get_color($self->tick_color));
	if ($text_loc)
	  {
	    my $h = $text_loc =~ /-/ ? $y2-1: $y1-$self->padding/2;
	    $self->_gd_string(text=>'1',x=>$x+$unit+2, 'y'=>$h, size => ($y2-$y1)/1.25);
	  }
      }
    while ($tick <= $re)
      {
	my $x = $w *($tick - $self->_region_start)/($self->_region_stop - $self->_region_start);
	print STDERR "Generating tick at $tick ($x)\n" if $self->DEBUG;
	$gd->filledRectangle($x, $y1, $x+$unit, $y2, $self->get_color($self->tick_color));
	if ($text_loc)
	  {
	    my ($key) = $tick =~ /(0+)$/;
	    $key = defined $key ? "1".$key : "1";
	    my %end = (
		       1=>"",
		       10=>"0",
		       100=>"00",
		       1000=>"K",
		       10000=>"0K",
		       100000=>"00K",
		       1000000=>"M",
		       10000000=>"0M",
		       100000000=>"00M",
		       1000000000=>"G",
		       10000000000=>"0G",
		       100000000000=>"00G",
		      );
	    my $t = $tick/$key . $end{$key};
	    my $h = $text_loc =~ /-/ ? $y2-1: $y1-$self->padding/2;
	    $self->_gd_string(text=>$t,x=>$x+$unit+2, 'y'=>$h, size => ($y2-$y1)/1.25);#/1.5  );
	  }
	$tick+= $div;
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
    my $x = $opts{x};
    $x = $opts{X} unless defined $x;
    my $y = $opts{'y'};
    $y = $opts{Y} unless defined $y;
    unless (defined $x && defined $y)
      {
	warn ("X: $x or Y: $y is not defined.  Can't generate string without coordinates");
	return 0;
      }
    my $color = $opts{color} || $opts{COLOR};
    my $size = $opts{size} || $opts{SIZE} || $self->padding;
    my $angle = $opts{angle} || $opts{ANGLE} || 0;
    $color = $self->get_color($color);
    my $gd = $self->gd;
    
    if (-r $self->font)
      {
	$gd->stringFT($color, $self->font, $size, $angle, $x, $y+$size, $text);
      }
    else
      {
	$gd->string($FONT, $x, $y, $text, $color);
      }
  }

#################### subroutine header begin ####################

=head2 _set_region_for_point

 Usage     : $self->_set_region_for_point($chromosomal_location);
 Purpose   : internal method to set the ranges of the viewable region of the chromosome
           : it checks the accessor function $self->start_picture in order to determine
           : if this point is at the left, right or center of the image.
 Returns   : none
 Argument  : integer which cooresponds to the chromosomal location on which the view will
           : will be centered
 Throws    : 0 and a warning if no point was specified
 Comment   : called internally by $self->set_region
           : 

See Also   : 

=cut

#################### subroutine header end ####################


sub _set_region_for_point
  {
    my $self = shift;
    my $point = shift;
    unless (defined $point)
      {
	warn 'You must specify a point!\n';
	return (0);
      }
    my ($start, $stop);
    if ($self->start_picture && $self->start_picture =~ /^l/)
      {
	$start = $point;
	$stop = $start + $self->find_size_for_magnification();
      }
    elsif ($self->start_picture && $self->start_picture =~ /^r/)
      {
	$stop = $point;
	$start = $stop - $self->find_size_for_magnification();
      }
    else
      {
	my $range = ceil ($self->find_size_for_magnification()/2);
	$start = $point - $range;
	$stop = $point + $range;
      }
    $self->_region_start($start);
    $self->_region_stop($stop);
  }


#################### subroutine header begin ####################

=head2 _set_region_start_stop

 Usage     : $self->_set_region_start_stop($start, $end);
 Purpose   : uses the $start and $end chromosomal positions to determin the highest magnification
           : that will show the entire specified range
 Returns   : 0 unless $start and $end are defined
 Argument  : two integers that correspond to the start and end chomosomal positions to be viewed
 Throws    : 0 unless $start and $end are defined
 Comment   : called internally by $self->set_region
           : 

See Also   : 

=cut

#################### subroutine header end ####################


sub _set_region_start_stop
  {
    my $self = shift;
    my ($start, $stop) = @_;
    return 0 unless (defined $start && defined $stop);
    if ($stop < $start)
      {
	my $tmp = $start;
	$start = $stop;
	$stop = $tmp;
      }
    $self->_set_region_for_point($start) if $start == $stop;
    my $len = $stop - $start+1;
    my $mag = $self->find_magnification_for_size($len);
    $self->mag($mag);
#    $self->start($start);
#    $self->stop($stop) if $stop;
    my $size = $self->mag_scale->{$mag};
    my $diff = ceil( ($size-$len)/2);
    my $rstart = $start-$diff;
    my $rstop = $stop + $diff;
    $self->_region_start($rstart);
    $self->_region_stop($rstop);
#    print STDERR Dumper $self->mag_scale;
#     print STDERR qq{
# IN Chomosome.pm sub _set_region_start_stop
#    requested ($start - $stop)
#    mag : $mag
#    len : $len
#    size: $size
#    diff: $diff
#    set ($rstart - $rstop)
# };
  }


#################### subroutine header begin ####################

=head2 _region_end

 Purpose   : alias for accessor function $self->_region_stop

=cut

#################### subroutine header end ####################


sub _region_end
  {
    my $self = shift;
    return $self->_region_stop(@_);
  }

#################### subroutine header begin ####################

=head2 _region_begin

 Purpose   : alias for accessor function $self->_region_start

=cut

#################### subroutine header end ####################


sub _region_begin
  {
    my $self = shift;
    return $self->_region_start(@_);
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

