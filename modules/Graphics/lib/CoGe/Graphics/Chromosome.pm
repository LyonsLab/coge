package CoGe::Graphics::Chromosome;
use strict;
use base qw(Class::Accessor);
use POSIX;
use Data::Dumper;
use Benchmark;
use GD;

#################### main pod documentation begin ###################
##
##

=head1 NAME

CoGe::Graphics::Chromosome - Object for drawing chromosomes that provides functionality for painting chromosomes with location based features (such as genes), and printing pictures of chromosomes as pngs.

=head1 SYNOPSIS

  #!/usr/bin/perl -w

  use strict;

  use CoGe::Graphics::Chromosome;
  use CoGe::Graphics::Feature::Gene;
  use CoGe::Graphics::Feature::NucTide;

  #create a chromosome object;
  my $c = CoGe::Graphics::Chromosome->new();
  #set the size of the chromosome in nucleotides
  $c->set_region(start=>1, stop=>100000);

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
  $f->add_segment(start=>11000, end=>12000); #feature goes off region, but no worries!
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
  my $seq = get_dna_sequence(); #returns a string of DNA sequence (ATCGTC...) for the region

  my $seq_len = length $seq;
  my $chrs = int (($c->region_length)/$c->iw); #number of characters to use per pixel
  $chrs = 1 if $chrs < 1;
  my $pos = 0; #position in sequence string
  while ($pos < $seq_len)my $ i = 0;
   {
     my $subseq = substr ($seq, $pos, $chrs);
     my $rcseq = substr ($seq, $pos, $chrs);
     $rcseq =~ tr/ATCG/TAGC/;
     next unless $subseq && $rcseq;
     my $f1 = CoGe::Graphics::Feature::NucTide->new({nt=>$subseq, strand=>1, start =>$pos+1, options=>"gc"});
     my $f2 = CoGe::Graphics::Feature::NucTide->new({nt=>$rcseq, strand=>-1, start =>$pos+1, options=>"gc"});
     $f1->show_label(1);
     $f2->show_label(1);
     $c->add_feature($f1) if $f1;
     $c->add_feature($f2) if $f2;
     $pos+=$chrs;
   }

  #set to "1" to print labels of features;
  $c->feature_labels(0);

  #turn on the flag for printing labels of the nucleotides (which are "fill" type features)
  $c->fill_labels(1);

   $c->generate_png(file=>"tmp/test$i.png");

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

2. Allow features (also objects) to easily be added/painted on to the chromosome

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

Permission to use, copy, modify, and distribute this software and its documentation for educational, research, and not-for-profit purposes, without fee and without a signed licensing agreement, is hereby granted, provided that the above copyright notice, this paragraph and the following two paragraphs appear in all copies, modifications, and distributions. Contact The Office of Technology Licensing, UC Berkeley, 2150 Shattuck Avenue, Suite 510, Berkeley, CA 94720-1620, (510) 643-7201, for commercial licensing opportunities.

The full text of the license can be found in the
LICENSE file included with this module.

=head1 SEE ALSO

perl(1).

=cut

#################### main pod documentation end ###################

BEGIN {
    use vars qw($VERSION $DEFAULT_WIDTH $PADDING $DEFAULT_COLOR $CHR_HEIGHT $CHR_INNER_COLOR $CHR_OUTER_COLOR $RULER_COLOR $TICK_COLOR $RULER_HEIGHT $FONT $FONTTT $FEATURE_HEIGHT);
    $VERSION     = '0.1';
    $DEFAULT_WIDTH = 200;  #default image width pixels
    $PADDING = 10; #default value to pad image height
    $CHR_HEIGHT = 30; #initial height of the chromosome image in pixels before magnification is added
    $DEFAULT_COLOR = [0,0,0];
    $CHR_INNER_COLOR = [255,255,255]; #inner color for chromosome
    $CHR_OUTER_COLOR = [0,0,0]; #border color for chromosome
    $RULER_COLOR = [200,200,255]; #color for measurement ruler
    $TICK_COLOR = [0,0,255]; #color for ticks on measurement ruler
    $RULER_HEIGHT = 20; #height (in pixels) of the ruler
    $FONTTT = "/usr/local/fonts/arial.ttf";  #needs to be fixed to get from conf file
    $FONT = GD::Font->MediumBold; #default GD font
    $FEATURE_HEIGHT = 20; #the height of a feature in pixels
    __PACKAGE__->mk_accessors(
"DEBUG",
"feature_height", #height of a feature in pixels
"draw_ruler", #flag for drawing ruler
"draw_chr_end", #flag for drawing "ends" of chromosome"
"ruler_color", "tick_color", #color for ruler and ticks on ruler respectively
"ruler_height", #height of ruler
"image_width", "image_height",
"padding",
"font",
"chr_height", #the starting height of the chromosome,
"feature_labels", "fill_labels", #flag to turn off the printing of labels, fill_lables are specifically for filled features;
"draw_chromosome", "chr_inner_color", "chr_outer_color",
"overlap_adjustment", #flag to turn on/off the overlap adjustment for overlapping features
"skip_duplicate_features", #flag to turn on/off skipping duplicate feature
"major_tick_labels", #whether to and where to draw major tick lables (1 above, -1 below, 0 none)
"minor_tick_labels", #whether to and where to draw minor tick lables (1 above, -1 below, 0 none)
"region_start", "region_stop", #image's start and stop (should be equal to or "larger" than the users
"_image_h_used", #storage for the amount of the image height used
"_gd", #store GD object
 "_chr_brush",
"_chr_center", "_chr_height", "_chr_h1", "_chr_h2", #interal storage of chromosome image height positions
"_features", #internal storage of features
"_fill_features", #internal storeage of fill features
"chr_end",#end position of chromosome -- this lets you stop the chromosome in the middle of the image
"_max_track", #place to store the maximum number of tracks on which to draw genomic features
"region_generated",#place to store a flag for whether or not the region has been previously generated;
"benchmark", #stores a flag for printing benchmark information for image generation
"invert_chromosome", #flag to draw the image in reverse so that the chromosome has been "flipped" 180 degrees
"draw_hi_qual", #flag to draw high quality features on chromosome (but slower)
"top_padding", #how many pixels to pad the top of the image with white space
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
    my ($class, %opts) = @_;
    #print STDERR "CHR IS BEING CALLED\n";
    my $self = bless ({}, ref ($class) || $class);

    $self->chr_height($CHR_HEIGHT);
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
    $self->feature_height($FEATURE_HEIGHT);
    $self->overlap_adjustment(0);
    $self->skip_duplicate_features(0);
    $self->font($FONTTT);
    $self->major_tick_labels(-1);
    $self->minor_tick_labels(0);
    $self->_features({});
    $self->_fill_features([]);
    $self->_image_h_used(0);
    $self->draw_hi_qual(0);
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

 benchmark        =>    (DEFAULT: 0) Output benchmarking on image generation

 region_start     =>    starting position of the chromosome (USER SPECIFIED)
 alias: start

 region_stop      =>    stopping position of the chromosome (USER SPECIFIED)
 alias: stop

 draw_chromosome  =>    (DEFAULT: 1) Flag (0 or 1) for whether or not the chromosome is
                        drawn on the final image
 draw_ruler       =>    (DEFAULT: 1) Flag (0 or 1) for whether or not the positional ruler
                        is drawn on the image
 draw_chr_end     =>    (DEFAULT: 1) Flag (0 or 1) for whether or not the rounded ends of the chromosome
                        are drawn where appropriate
 feature_height   =>    Height of a feature in pixels.  (Default: 20).  This is used if automatic zoom
                        is not used

 ruler_color      =>    (DEFAULT: [0,0,255]) Defines the color of the positional ruler.
                        This is the an array reference of three values
                        that corresponds to RGB color values.  Each color ranges from 0-255
 tick_color       =>    (DEFAULT: [0,0,255]) Defines the color of ticks on the  positional ruler.
                        This is the an array reference of three values
                        that corresponds to RGB color values.  Each color ranges from 0-255
 ruler_height     =>    (DEFAULT: 20)  The heigth, in pixels of the positional ruler

 major_tick_labels=>    Options for drawing major tick lables.  1 draws them above the tick, -1 draws them below the tick,
                        0 for not drawing tick labels.  (DEFAULT: -1)

 minor_tick_labels=>    Options for drawing minor tick lables.  1 draws them above the tick, -1 draws them below the tick,
                        0 for not drawing tick labels.  (DEFAULT: 0)

 chr_height       =>    (DEFAULT: 30)  This is the number, in pixels, of the starting height of the
                        chromosome before adjustments for featurs are made

 image_width      =>    (DEFAULT: 200) The width in pixels of the final image.
 alias:  iw

 image_height     =>    This holds the height of the image and is a value that is calculated
                        dynamically by the module (sub set_image_height) when the image is
                        generated.  IMPORTANT:  THIS VALUE SHOULD NOT BE MODIFIED BY THE USER
                        DIRECTLY.  One thing to keep in mind is that, the height of the
                        chromosomal images are dynamic.  This is due to the factors such as the
                        number and scaling aspects of features on the chromosome,and
                        customize the final height of the image by specifying the scaling factors
                        and heights of the various image parts, but it is not recommended to
                        change this value as strange(tm) things may happen.
 alias:  ih

 padding         =>     (DEFAULT: 15) This is the padding (in pixels) used between most items
                        drawn on the final image.

 font            =>     (DEFAULT: "/usr/lib/perl5/site_perl/CoGe/fonts/arial.ttf")
                        This is the path to a true-type font used for text labels on the image

 feature_labels  =>     (DEFAULT: 0) Flag used for whether or not to print feature labels.

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

 invert_chromosome  =>  Draw the chromosome such that it has been inverted 5' => 3'

 overlap_adjustment =>  flag for whether overlapping features are rescaled and position such that
                        they don't overlap when the image is generated.  Default: 1

 skip_duplicate_features => flag for whether to skip two featrues if they are identical.  Default: 0

 draw_hi_qual       =>  (DEFAULT: 0)This flag determines if the high quality mapping function for drawing features
                        on the chromosome is used or the low quality mapping.  The cost, of course, is
                        speed (roughly twice as long for high quality).  Overall, there is only minor difference
                        between hi qual and low qual images

 top_padding        =>  (DEFAULT: 0) Amount of whitespace padding added to the top of the final image.

=cut

#################### subroutine header end ####################

#################### subroutine header begin ####################

#=head2 start

# Usage     : $c->start(100);
#             my $start = $c->start;
# Purpose   : alias for region_start

#=cut

#################### subroutine header end ####################

sub start
  {
    return shift->region_start(@_);
  }

#################### subroutine header begin ####################

#=head2 stop

# Usage     : $c->stop(100);
#             my $stop = $c->stop;
# Purpose   : alias for region_stop

#=cut

#################### subroutine header end ####################

sub stop
  {
    return shift->region_stop(@_);
  }

sub _region_start
  {
    return shift->region_start(@_);
  }

sub _region_stop
  {
    return shift->region_stop(@_);
  }

#################### subroutine header begin ####################

=head2 set_region

 Usage     : $c->set_region(start=>$start, stop=>$stop);
 Purpose   : This routine sets the region by define region_start and region_stop
 Returns   : none
 Argument  : hash with at least one key-value pair for "start"
             accepts "start", "begin", "START", "BEGIN" to specify the beginning of the region
             accepts "stop", "end", "STOP", "END" to specify the end of the region
             Nominally, the values should be integers the correspond to a chromosomal location.
 Throws    : None

=cut

#################### subroutine header end ####################

sub set_region
  {
    my $self = shift;
    my %opts = @_;
    my $start = $opts{start} || $opts{begin} || $opts{START} || $opts{BEGIN};
    my $end = $opts{stop} || $opts{end} || $opts{STOP} || $opts{END};
    $self->region_start($start);
    $self->region_stop($end);
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
           : fill_height    => 1
           : stop           => start
           : merge_perecent => 100
           : magnification  => 1
           : overlay        => 1
           : mag            => 1
           : layer          => 1
           : type           => "unknown"
           : Also, the feature's GD object will be initialized upon import.
	   : There is a check for whether the added feature overlaps other features.
	   : If so, a counter, $feat->_overlap is incemented in the feature object.
	   : This is later used by the $self->_draw_feature algorithm to figure
	   : out how to best draw overlapping features.  The overlap check is skipped
           : unless $self->overlap_adjustment is true.
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
	unless (ref ($feat) =~ /Feature/i)
	  {
	    warn "Feature ($feat) does not appear to be a feature object.  Skipping. . .\n";
	    next;
	  }
	if ($feat->start || $feat->stop)
	  {
	    my $start = $feat->start;
	    my $stop = $feat->stop;
	    $stop = $start unless $stop;
	    $start = $stop unless $start;
	    ($start, $stop) = ($stop, $start) if $start> $stop;
	    $feat->start($start);
	    $feat->stop($stop);
	  }
	$feat->strand(1) unless defined $feat->strand;
	$feat->fill(0) unless $feat->fill;
	$feat->fill_height(1) unless defined $feat->fill_height && $feat->fill_height =~ /\d/ && $feat->fill_height <=1 && $feat->fill_height >= 0;

	$feat->merge_percent(100) unless defined $feat->merge_percent;
	$feat->stop($feat->start) unless defined $feat->stop;
	$feat->ih(0) unless defined $feat->ih;
	$feat->iw(0) unless defined $feat->iw;
	$feat->gd; #initialize feature;
	$feat->mag(1) unless defined $feat->mag;
	$feat->overlay(1) unless defined $feat->overlay();
	$feat->_overlap(1) unless defined $feat->_overlap;#detects overlapping feature on the same track
	$feat->_overlap_pos(1) unless $feat->_overlap_pos; #placement for overlapping features
	$feat->layer(1) unless $feat->layer;
	$feat->type("unknown") unless $feat->type;
	unless (defined $feat->order)
	  {
	    my $last_feat = $self->get_feats(last_order=>1, strand=>$feat->strand, fill=>$feat->fill);
	    my $order = $last_feat ? $last_feat->order()+1 : 1;
	    $feat->order($order);
	  }
#	$self->_check_overlap($feat) if $self->overlap_adjustment;
	if ($self->skip_duplicate_features)
	  {
	    next if $self->_check_duplicate($feat);   #should implement this
	  }
	my $feats = $self->_features;
	my $strand = $feat->strand =~ /-/ ? "-1" : "1";
	$feats->{$feat->type}{$strand}{$feat->order}{$feat->layer}{$feat->fill}{$feat->start}{$feat}=$feat;
      }
  }

#################### subroutine header begin ####################

=head2 delete_feature

 Usage     : $c->delete_featuer($feat);
 Purpose   : deletes a feature from the chromosome graphics object
 Returns   : nothing
 Argument  : a CoGe::Graphics::Feature object or derivative object
 Throws    :
 Comment   : Features are stored in a complex hash for quick and speedy retrieval

See Also   :

=cut

#################### subroutine header end ####################

sub delete_feature
  {
    my $self = shift;
    my $f = shift;
    delete $self->_features->{$f->type}{$f->strand}{$f->order}{$f->layer}{$f->fill}{$f->start}{$f};
  }

#################### subroutine header begin ####################

=head2 delete_features

 Usage     : $c->delete_features('nt');
 Purpose   : Deletes features from the object by the type of feature
 Returns   : none
 Argument  : string or none
               all (or blank) => deletes all the features
               <name of feature type> => e.g. "gene", "tRNA", "aa", "nt", etc.  depends on what feature derivatives used
 Throws    :
 Comment   :

=cut

#################### subroutine header end ####################

sub delete_features
  {
    my $self = shift;
    my $type = shift;
    if ($type eq "all" || !defined $type)
      {
	$self->_features({});
	return;
      }
    my $feats = $self->_features;
    delete $feats->{$type};
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
    my $order = $opts{order} || $opts{ORDER} || $opts{track};
    my $type = $opts{type} || $opts{TYPE};
    my $last = $opts{last_order} || $opts{LAST_ORDER}; #flag to get the highest order feature for a location
    my $strand = $opts{strand} || $opts{STRAND};
    if ($strand)
      {
	$strand = $strand =~ /-/ ? "-1" : "1"; #make this consistent with add_feature
      }
    my $fill = $opts{fill};
    $fill = $opts{FILL} unless defined $fill; #get filled features?
    my $start = $opts{start};
    my $stop  = $opts{stop};
    my $layer  = $opts{overlay} || $opts{layer};

    my $feats = $self->_features;
    my @return_feats;

    foreach my $t (keys %$feats)
      {
	if ($type)
	  {
	    next unless $t eq $type;
	  }
	foreach my $s (keys %{$feats->{$t}})
	  {
	    if ($strand)
	      {
		next unless $strand eq $s;
	      }
	    foreach my $o (keys %{$feats->{$t}{$s}})
	      {
		if ($order)
		  {
		    next unless $o eq $order;
		  }
		foreach my $l (keys %{$feats->{$t}{$s}{$o}})
		  {
		    if (defined $layer)
		      {
			next unless $l eq $layer;
		      }
		    foreach my $f (keys %{$feats->{$t}{$s}{$o}{$l}})
		      {
			if (defined $fill)
			  {
			    next unless $f eq $fill;
			  }
			foreach my $sp (keys %{$feats->{$t}{$s}{$o}{$l}{$f}})
			  {
			    if (defined $start)
			      {
				$stop = $start unless $stop;
				if ($start > $stop)
				  {
				    my $tmp = $stop;
				    $stop = $start;
				    $start = $tmp;
				  }
				next unless $sp >= $start && $sp <=$stop;
			      }
			    push @return_feats, values %{$feats->{$t}{$s}{$o}{$l}{$f}{$sp}};
			  }
		      }
		  }
	      }
	  }
      }
    if ($last)
      {
	my @return_feats = sort {abs($b->order) <=> abs($a->order)} @return_feats;
	return $return_feats[0];
      }
#    print "Found feats: ", scalar @return_feats,"\n";
#    print Dumper \@return_feats;
    return wantarray ? @return_feats : \@return_feats;
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
    $self->generate_region() unless $self->region_generated;
#    $self->gd->transparent($self->gd->colorResolve(255,255,255));
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
    print STDERR "\n", join "\t",  "region start: ".$self->region_start,"region end: ".$self->region_stop,"\n" if $self->DEBUG;
#    print Dumper $self if $self->DEBUG;
    $self->chr_end($self->stop) unless $self->chr_end;
    my $t0 = new Benchmark;
    $self->set_image_height();
    my $t1 = new Benchmark;
    $self->gd->fill(0,0,$self->get_color([255,255,255]));
    my $t2 = new Benchmark;
    $self->_draw_ruler;
    my $t3 = new Benchmark;
    $self->_draw_chromosome;
    my $t4 = new Benchmark;
    $self->_draw_features;
    my $t5 = new Benchmark;
    $self->region_generated(1);
    my $ih_time = timestr(timediff($t1, $t0));
    my $fl_time = timestr(timediff($t2, $t1));
    my $rl_time = timestr(timediff($t3, $t2));
    my $ch_time = timestr(timediff($t4, $t3));
    my $ft_time = timestr(timediff($t5, $t4));
    print STDERR qq {
         ---------------------------------------
         Time to set image height:   $ih_time
         Time to full region bg:     $fl_time
         Time to draw ruler:         $rl_time
         Time to draw chromosome:    $ch_time
         Time to draw features:      $ft_time
} if $self->benchmark;
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
    foreach my $feat ( $self->get_features(fill=>0))
      {
	next if $feat->fill;
	#skip drawing features that are outside (by two times the range being viewed) the view
	if ($feat->start)
	  {
	    next if $feat->start > $self->region_stop+2*($self->region_length );
	  }
	if ($feat->stop > 0)
	  {
	    next if $feat->stop < $self->region_start-2*($self->region_length );
	  }
	my $anno = $feat->description;
	next unless $anno;
	#clean this for javascript
	$anno =~ s/\n/<br>/g;
	$anno =~ s/\t/&nbsp;&nbsp;/g;
	$anno =~ s/"//g;
	$anno =~ s/'//g;
	my $feat_height = $self->feature_height;
	my $feat_h = $feat_height/$feat->_overlap;
	my $offset = ($feat->order-1)*($feat_height+$self->padding/1.5)+$self->padding/2;
	$offset = 0 if $feat->fill;
	$feat_h = ($self->_chr_height)/2 * $feat->fill_height if $feat->fill;
	my $y = $feat->strand =~ /-/ ? $c+ $offset+1+($feat_h)*($feat->_overlap_pos-1): $c - $offset-$feat_h*$feat->_overlap_pos;
	my $rb = $self->region_start;
	my $re = $self->region_stop;
	my $range = $re-$rb;
	my $w = $self->iw;
	$feat->gd;
	$feat->stop($feat->start) unless defined $feat->stop;
	my $feat_range = $feat->stop-$feat->start;
	my $unit = $self->_calc_unit_size;
	my $fs = $unit*($feat->start-$rb);
	my $fe = $unit*($feat->end-$rb+1);
	#skip stuff that is outside of the image
	next if $fe < 1;
	next if $fs > $self->iw;
	my $fw = sprintf("%.1f",$fe - $fs)+1; #calculate the width of the feature;

	next if $fw < 1; #skip drawing if less than one pix wide
	my $link = $feat->link;
	my $alt = $feat->alt || $feat->label;
	my $x2 = ceil($fs+$fw);
	my $y2 = ceil($y+$feat_h);
	$y = floor $y;
	$fs = floor($fs);
	$html .= qq{
<area coords="$fs, $y, $x2, $y2"};
	$html .= qq{ href="$link" target=_new } if $link;
	$html .= qq{ onMouseOver="change('$anno')"} if $anno;
	$html .= qq{ alt="$alt"} if $alt;
	$html .= ">\n";
      }
    $html .= "</map>\n";
    return $html;
  }

#################### subroutine header begin ####################

#=head2 ih

# Purpose   : alias for Accessor method $self->image_height

#=cut

#################### subroutine header end ####################

sub ih
  {
    my $self = shift;
    return $self->image_height(@_);
  }

#################### subroutine header begin ####################

#=head2 iw

# Purpose   : alias for Accessor method $self->image_width

#=cut

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
	$gd->transparent($white);
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
             a number of factors including feature height, the number
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
    unless ($self->_max_track)
      {
	my $top_feat = $self->get_feats(last_order=>1, strand=>1, fill=>0);
	my $bot_feat = $self->get_feats(last_order=>1, strand=>-1, fill=>0);
	my $max = $top_feat->order if $top_feat;
	$max = $bot_feat->order unless !$bot_feat || $max;
	$max = $bot_feat->order if $bot_feat && $max && $bot_feat->order > $max;
	$max = 1 unless $max;
	$self->_max_track($max);
      }
    my $feat_space = $self->_max_track * ($self->feature_height+$self->padding);
    my $ch = $feat_space*2 > $self->chr_height ? $feat_space*2 : $self->chr_height; #chromosome height

    my $ruler_h = $self->ruler_height + $self->padding;
    my $h=$ruler_h;
    $h+=$ch;
    if ($self->top_padding)#add space to top
      {
	$h+= $self->top_padding;
	$self->_image_h_used($self->top_padding);
      }
    print STDERR "Image Height: $h\n" if $self->DEBUG;
    if ($self->ih)
      {
	my $scale = $self->ih/$h;
	$self->padding($self->padding*$scale);
	$self->feature_height($self->feature_height*$scale);
	$self->_chr_height($ch*$scale);
	$self->_chr_center(($h-$ch/2)*$scale);
      }
    else
      {
	$self->ih($h);
	$self->_chr_height($ch);
	$self->_chr_center($h-$ch/2);
      }
    $self->ih($self->ih+$self->chr_brush->height*3);# if $self->draw_chromosome;;
#    $self->_image_h_used(0);

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
    return $self->_chr_brush if $self->_chr_brush;
    my @bc = @{$self->chr_outer_color()};
    my @ic = @{$self->chr_inner_color()};
    my $size = 5;
    my $edge_brush = new GD::Image(1,$size);
    $edge_brush->colorResolve(255,255,255);
    $edge_brush->line(0,0,0,0,$edge_brush->colorResolve(@bc));
    $edge_brush->line(0,$size,0,$size,$edge_brush->colorResolve(@ic));
    for my $i (1..$size-1)
      {
	my @c;
	for my $j (0..2)
	  {
	    push @c, $bc[$j]+($ic[$j]-$bc[$j])/($size)*$i;
	  }
	$edge_brush->line(0,$i,1,$i,$edge_brush->colorResolve(@c));
      }
    $self->_chr_brush($edge_brush);
    return $self->_chr_brush;
  }

#################### subroutine header begin ####################

=head2 region_length

 Usage     : my $length = $self->region_length()
 Purpose   : returns the length of the chromosomal region
 Returns   : int
 Argument  : none
 Throws    : none
 Comment   : return the value of $self->region_stop - $self->region_start + 1;
           :

See Also   :

=cut

#################### subroutine header end ####################

sub region_length
  {
    my $self = shift;
#    my $stop = $self->chr_end ? $self->chr_end : $self->stop;
    return ($self->stop - $self->start+1);
  }

#################### subroutine header begin ####################

=head2 chr_length

 Usage     : my $length = $c->chr_length
 Purpose   : alias for region_length

=cut

#################### subroutine header end ####################

sub chr_length
  {
    shift->region_length(@_);
  }

#################### subroutine header begin ####################

=head2 _invert_chromosome

 Usage     : $c->_invert_chromosome
 Purpose   : makes up->down, down->up, etc.
 Returns   :
 Argument  :
 Throws    :
 Comment   : Should not be called directly.  Set invert_chromosome flag to 1 and this will be called by _draw_features
           :

See Also   :

=cut

#################### subroutine header end ####################

sub _invert_chromosome
  {
    my $self = shift;
    #make up->down, down->up, left->right and right->left . . . up, up, down, down, left, right, left, right, B, A start!
    foreach my $feat ($self->get_features())
      {
	my $strand = $feat->strand =~ /-/ ? 1 : "-1";
	my $start = $self->region_stop - $feat->stop;
	my $stop = $self->region_stop - $feat->start;
	$feat->strand($strand);
	$feat->start($start);
	$feat->stop($stop);
	$feat->gd->flipHorizontal();
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
    return if $feat->skip_overlap_search;
    my @feats = $self->get_feats(strand=>$feat->strand, fill=>$feat->fill, order=>$feat->order);
    return unless @feats;
    my @overlapped_feats;
#    print STDERR "Checking ",$feat->start,"-",$feat->stop,"\n";
    foreach my $f (@feats)
    	{
	  next if $feat eq $f;
	  next unless $feat->overlay() == $f->overlay();  #skip the check if the features are at different overlay levels.
	  push @overlapped_feats, $f unless ( ($feat->start > $f->stop) || ($feat->stop < $f->start) );
	}
    #let's figure out the number of overlaps for each sequence
    my @overlap_tracker;
    my $i =0;
#    foreach my $f1 (sort {abs($b->stop-$b->start) <=> abs($a->stop-$a->start)} @overlapped_feats, $feat)
    foreach my $f1 (sort {$a->start <=> $b->start} @overlapped_feats, $feat)
      {
	my $count=1;
	my $prevfeat;
	foreach my $f2 (sort {$a->start <=> $b->start} @overlapped_feats, $feat)
#	foreach my $f2 (sort {abs($b->stop-$b->start) <=> abs($a->stop-$a->start)} @overlapped_feats, $feat)
	  {
	    my $s1 = scalar $f1;
	    my $s2 = scalar $f2;
	    next if $s1 eq $s2;
	    next if ( ($f1->start > $f2->stop) || ($f1->stop < $f2->start) );# skip if there ain't no overlap
	    if ($prevfeat)
	      {
		if ($prevfeat->stop < $f2->start)
		  {
		    $prevfeat = $f2;
		    next;
		  }
	      }
	    $count++;
	    $prevfeat = $f2;
	  }
	$i++;
	push @overlap_tracker,[$f1,$count, $i];
      }
    my ($item1, $item2) = sort {$b->[1] <=> $a->[1]} @overlap_tracker;
    $item1->[1] = $item2->[1] if $item2->[1];
    foreach my $item (@overlap_tracker)
      {
	my ($f,$count,$pos) = @$item;
	$f->_overlap($count) unless $f->_overlap() > $count;
	$f->_overlap_pos($pos);
#	print STDERR $f->start,"-",$f->stop," ",$f->_overlap," ",$f->_overlap_pos,"\n";
      }
#    print STDERR "\n";
  }

sub _check_duplicate
  {
    my $self = shift;
    my $feat = shift;
    return 0 if $feat->skip_duplicate_search;
    my @feats = $self->get_feats(strand=>$feat->strand, fill=>$feat->fill, order=>$feat->order, type=>$feat->type, start=>$feat->start);
    return 0 unless @feats;
    my $check = 0;
    foreach my $f (@feats)
      {
	$check = 1 if $f->stop eq $feat->stop && $f->layer eq $feat->layer && $f->order eq $feat->order;
      }
    return $check;
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
    return unless $self->draw_chromosome; #do we draw the chromsome picture?

    my $gd = $self->gd;
    my $ch = $self->_chr_height;
    $ch = $ch/2; #want half the height for the rest of the calcs
    my $hc = $self->_image_h_used+($self->ih-$self->_image_h_used)/2; #height of center of chromsome image
    my $w = $self->iw; #width of image
    my $xs = $self->region_start < 1 ? $w*abs(1-$self->region_start)/($self->region_length): 0;
    my $xe = $self->region_stop > $self->chr_end ? $w-$w*($self->region_stop-$self->chr_end -1)/($self->region_stop-$self->start-1) : $w;
#    print STDERR join ("\t",$self->start, $self->stop, $self->chr_end, $self->chr_length,($self->chr_end-$self->region_stop - 1)/($self->chr_end-$self->start-1)),"\n";
    print STDERR "Chromosome image: Height/2: $ch, Height Center: $hc, xs: $xs, xe: $xe  draw_chromsome: ".$self->draw_chromosome."\n" if $self->DEBUG;

    my $brush = $self->chr_brush;
    $gd->setBrush($brush);
    $gd->line($xs, $hc-$ch-$brush->height, $xe, $hc-$ch-$brush->height, gdBrushed) unless $xe < 0;
    $brush->flipVertical();
    $gd->setBrush($brush);
    $gd->line($xs, $hc+$ch+$brush->height, $xe, $hc+$ch+$brush->height, gdBrushed) unless $xe < 0;
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
    my $ch = $self->_chr_height+$self->chr_brush->height*2.5;

    my @arc1 = $dir =~ /left/i ? (90, 180) : (0, 90);
    my @arc2 = $dir =~ /left/i ? (180, 270) : (270, 360);
    my @arc3 = $dir =~ /left/i ? (90, 270) : (270,90);
    $gd->filledArc($x, $y, $ch, $ch, @arc3, $self->get_color($self->chr_inner_color));
#    $self->chr_brush->flipVertical;
    $gd->setBrush($self->chr_brush);
    $gd->arc($x, $y, $ch, $ch, @arc1, gdBrushed);
    $self->chr_brush->flipVertical;
    $gd->setBrush($self->chr_brush);
    $gd->arc($x, $y, $ch, $ch, @arc2, gdBrushed);
    $self->chr_brush->flipVertical;
    $gd->setBrush($self->chr_brush);

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
    $self->_invert_chromosome if $self->invert_chromosome;
    my $c = $self->_image_h_used+($self->ih - $self->_image_h_used)/2;
    print STDERR "Image used: ".$self->_image_h_used."  Image Height: ".$self->ih."  Center: $c\n" if $self->DEBUG;
    foreach my $feat(sort {$a->start <=> $b->start || $a->stop <=> $b->stop}$self->get_features(fill=>0))
      {
	$self->_check_overlap($feat) if $self->overlap_adjustment;
      }
    foreach my $feat ( (sort {abs($b->stop-$b->start) <=> abs($a->stop-$a->start) || $b->type cmp $a->type} $self->get_feature(fill=>1)), sort {$a->overlay <=> $b->overlay || $a->start <=> $b->start || abs($a->start-$a->stop) <=> abs($b->start-$b->stop)} $self->get_features(fill=>0))
      {
	my $feature_height = $self->feature_height;
	my $feat_h = $feature_height/$feat->_overlap;
	$feat_h = 1 if $feat_h < 1;
	my $offset = ($feat->order-1)*($feature_height+$self->padding/1.5)+$self->padding/2;
	$offset = 0 if $feat->fill;
	$feat_h = ($self->_chr_height)/2 * $feat->fill_height if $feat->fill;
	my $feat_mag_offset = ($feat_h*$feat->mag - $feat_h)/2;
	my $y = $feat->strand =~ /-/ ?
	  $c + $offset + $feat_h  * ($feat->_overlap_pos-1) - $feat_mag_offset+ 1  :
	  $c - $offset - $feat_h  *  $feat->_overlap_pos - $feat_mag_offset;
#	print STDERR "overlap: ", $feat->_overlap," ",$feat->_overlap_pos," ",$feat->type," ",$feat->start," ",$y," ",$feat_h,"\n" if $feat->_overlap > 1;
	$y = sprintf("%.0f",$y);
#	print STDERR "\t",$y,"\n"if $feat->_overlap > 1;
	$feat_h *= $feat->mag;
        my $sy;
	if ($feat->fill)
	  {
	    $sy = $feat->strand =~ /-/ ? $c+2 : $c-$self->padding-$self->feature_height;
	  }
	elsif ($feat->label_location && $feat->label_location =~ /bot/)
	  {
	    $sy = $y+$feat_h*.75;
	  }
	elsif ($feat->label_location && $feat->label_location =~ /top/)
	  {
	    $sy = $y-$feat_h*.25;
	  }
	$self->_draw_feature(feat=>$feat, 'y'=>$y, ih=>$feat_h, 'sy'=>$sy, highqual=>$self->draw_hi_qual);
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

sub _draw_feature
  {
    my $self = shift;
    my %opts = @_;
    my $feat = $opts{feat} || $opts{FEAT} || $opts{f};
    return 0 unless ref($feat) =~ /Feature/i;
    my $y = $opts{'y'} || $opts{Y};
    my $ih = $opts{'image_height'} || $opts{'ih'} || $opts{'IH'} || $feat->ih;
    my $sy = $opts{'string_y'} || $opts{'sy'};#label y axis
    my $highqual = $opts{highqual};
    my $rb = $self->region_start;
    my $re = $self->region_stop;
    my $w = $self->iw;
    $feat->gd;
    $feat->stop($feat->start) unless defined $feat->stop;
    my $unit = $self->_calc_unit_size;
    my $fs = $unit*($feat->start-$rb);
    $fs = 0 if $fs < 0 && !$feat->force_draw;
    my $fe = $unit*($feat->end-$rb+1);
    $fe = $w if $fe > $w && !$feat->force_draw;
    my $fw = sprintf("%.1f",$fe - $fs)+1; #calculate the width of the feature;
    $fw = 1 if $feat->force_draw() && $fw < 1;
    unless ($feat->force_draw)
      {
	return if $fe < 0;
	return if $fw < 1; #skip drawing if less than one pix wide
      }
#    return if $fw > 1000000;
    my ($xmin, $xmax, $ymin , $ymax);

    print STDERR "Drawing feature ".$feat->label." Order: ".$feat->order." Overlap: ".$feat->_overlap." : ", $feat->start, "-", $feat->end," Dimentions:",$fw,"x",$ih, " at position: $fs,$y"."\n" if $self->DEBUG;
    if ($feat->fill)
      {
	$self->gd->copyResized($feat->gd, $fs, $y,0,0, $fw, $ih, $feat->iw, $feat->ih);
	if ($feat->external_image && $fw > 10) #if we have an external image and the feature width is greater than 10. . .
	  {
	    my $ei = $feat->external_image;
	    my $ex_wid = $ei->width;
	    my $ex_hei = $ei->height;
	    $ei->fill(0,0,$ei->colorResolve(@{$feat->color}));
	    my $scale = $fw/$ex_wid; #scaling factor for external image
	    my $hei = $feat->strand =~ /-/ ? $y : $y+$ih-($ex_hei*$scale);
	    #okay, we are going to need to do some fancy stuff in order to smoothly resize and paste the feature onto the main image
	    #1. create a blank image of the appropriate size
	    my $newgd = GD::Image->new ($fw, $ex_hei*$scale,[1]);
	    #2. copy, resize, and resample the feature onto the new image
	    $newgd->copyResampled($ei, 0, 0, 0, 0, $newgd->width, $newgd->height, $ex_wid, $ex_hei);
	    if ($highqual)
 	      {
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
 	      }#4. make white transparent
 	    $newgd->transparent($newgd->colorResolve(255,255,255));
# 	    #5. copy new image into appropriate place on image.
 	    $self->gd->copyMerge($newgd, $fs, $hei, 0, 0, $newgd->width, $newgd->height,100);
	  }
      }
    elsif ($fe > 0)
      {
# 	#okay, we are going to need to do some fancy stuff in order to smoothly resize and paste the feature onto the main image
# 	#1. create a blank image of the appropriate size
 	my $newgd = GD::Image->new ($fw, $ih,[1]);
	$newgd->fill(0,0,$newgd->colorResolve(255,255,255));
# 	#2. copy, resize, and resample the feature onto the new image
	my $srcx = $feat->start < $rb ? $rb-$feat->start : 0;
	$rb = 0 if $rb < 1;
	my $srcw = $feat->stop-$feat->start+1 - $srcx-1;
	$srcw -= ($feat->stop-$re) if $re< $feat->stop;
	# source image   destx dest y  srcx srcy  destw desth  srcw  srch
 	$newgd->copyResized($feat->gd, 0, 0, $srcx, 0, $fw, $ih, $srcw, $feat->ih);
	if ($highqual)
	      {
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
	      }#4. make white transparent# 	#4. make white transparent
 	$newgd->transparent($newgd->colorResolve(255,255,255));
# 	#5. copy new image into appropriate place on image.
	#srcimg, destx, desty srcx, srcy, wid, hei, mergepercent
 	$self->gd->copyMerge($newgd, $fs, $y, 0, 0, $fw, $ih, $feat->merge_percent);
      }
    $xmin = sprintf("%.0f",$fs);
    $ymin = sprintf("%.0f",$y);
    $xmax = sprintf("%.0f",$fe);
    $ymax = sprintf("%.0f",$y+$ih-1);

    my $size = $feat->font_size ? $feat->font_size : $ih;
    if ($self->fill_labels && $feat->fill) {$size=$fw >= 15 ? 15 : $fw;}
    $sy=$y+$ih/2-$size/2 unless $sy;
    my $adjust = 2;
	$fs+=$adjust unless $fs+$adjust > $self->iw;
    $self->_gd_string(y=>$sy, x=>$fs, text=>$feat->label, size=>$size) if $feat->label && ( ($self->feature_labels || $self->fill_labels)&& ($fw>5 || $feat->force_label)); #don't make the string unless the feature is at least 5 pixels wide
#    print STDERR $feat->type," ","$xmin, $ymin, $xmax, $ymax\n" if $feat->{anchor};
    $feat->image_coordinates("$xmin, $ymin, $xmax, $ymax") if defined $xmin && defined $ymin && defined $xmax && defined $ymax;
  }

#################### subroutine header begin ####################

=head2 _calc_unit_size

 Usage     : my $unit = $self->_calc_unit_size();
 Purpose   : returns the width, in pixels, of one chromsomal unit (usually nucleotides)
 Returns   : an number
 Argument  : none
 Throws    : none
 Comment   : formula is image_width/visable_region_size (nt)
           :

See Also   :

=cut

#################### subroutine header end ####################

sub _calc_unit_size
  {
    my $self = shift;
    return ($self->iw/($self->region_stop-$self->region_start+1));
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
    my %opts = @_;
    my $major_tick_labels = $opts{major_tick_labels};
    $major_tick_labels = $self->major_tick_labels unless defined $major_tick_labels;
    my $minor_tick_labels = $opts{minor_tick_labels};
    $minor_tick_labels = $self->minor_tick_labels unless defined $minor_tick_labels;
    my $gd = $self->gd;
    my $c = $self->ruler_height/2+$self->_image_h_used; #center of ruler
    $self->_image_h_used($self->_image_h_used + $self->ruler_height+$self->padding/2);
    return unless $self->draw_ruler;
    my $w = $self->iw; #width of image
    my $xb = $self->region_start < 1 ? $w*abs(1-$self->region_start)/($self->region_length): 0; #x position begin
    my $xe = $self->region_stop > $self->chr_end ? $w-$w*($self->region_stop-$self->chr_end -1)/($self->region_stop-$self->start-1) : $w;
#    print STDERR join ("\t", $self->start, $self->stop, $self->chr_end, $self->region_length),"\n";
    #$self->region_stop > $self->chr_length ? $w-$w*($self->region_stop - $self->chr_length)/($self->region_length) : $w; #x position end
    return unless ($xe>0);
    $gd->line($xb,$c,$xe,$c,$self->get_color($self->ruler_color));
    my $mtyb = $c-$self->ruler_height/2; #major tick y begin
    my $mtye = $c+$self->ruler_height/2; #major tick y end
    my $styb = $c-$self->ruler_height/4; #minor tick y begin
    my $stye = $c+$self->ruler_height/4; #minor tick y end
    my $region_width = $self->region_length;
    my $rb = $self->region_start-$region_width/4;
    $rb = 1 if $rb < 1;
    my $re = $self->region_stop;
    #$re = $self->chr_length if $re > $self->chr_length;
    my $range = $re-$rb+1; #chromosomal positional range
    $rb = $rb-($range/10); #back up a bit
    my $div = "1"."0"x int (log10($self->region_length)+.5); #determine range scale (10, 100, 1000, etc)
    print STDERR "\nRULER: Center: $c, Start $xb, Stop: $xe, range: $range, Ticks: $div, \n" if $self->DEBUG;
    $self->_make_ticks(scale=>$div, y1=>$mtyb, y2=>$mtye, range_begin=>$rb, range_end=>$re,text_loc=>$major_tick_labels);
    $self->_make_ticks(scale=>$div/10, y1=>$styb, y2=>$stye, range_begin=>$rb, range_end=>$re, text_loc=>$minor_tick_labels, label_size=>($mtye-$mtyb)/2.5);
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
    my $label = $opts{label}; #write labels on ticks?
    my $label_size = $opts{label_size};
    $label_size = ($y2-$y1)/2.5 unless defined $label_size;
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
    my @text;
    if ($rb <= 1 && $re >= 1)
      {
	my $x = $w *(1 - $self->region_start)/($self->stop-$self->start);
	print STDERR "Generating tick '1' at $tick ($x)\n" if $self->DEBUG;
	$gd->filledRectangle($x, $y1, $x+1, $y2, $self->get_color($self->tick_color));
	my $h = $text_loc =~ /-/ ? $y2-1: $y1;#$y1-$self->padding/2;
	push @text, {text=>'1',x=>$x+2, 'y'=>$h, size => $label_size};
	$self->_gd_string(text=>'1',x=>$x+1, 'y'=>$h, size => $label_size) if $text_loc;
      }
    while ($tick <= $re)
      {
	last if $tick > $self->chr_end;
	my $x = $w *($tick - $self->region_start)/($self->stop-$self->start);
	print STDERR "Generating tick at $tick ($x)\n" if $self->DEBUG;
	$gd->filledRectangle($x, $y1, $x+1, $y2, $self->get_color($self->tick_color)) unless $x<0;
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
	my $h = $text_loc =~ /-/ ? ($y2-$y1+1): $y1;#-$self->padding/2;
	push @text, {text=>$t,x=>$x+2, 'y'=>$h, size => $label_size};
	$self->_gd_string(text=>$t,x=>$x+2, 'y'=>$h, size => $label_size) if ($text_loc);
	$tick+= $div;
      }
    return \@text;
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
#	warn ("X: $x or Y: $y is not defined.  Can't generate string without coordinates");
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

sub round
  {
    my $self = shift;
    my $num = shift;
    return $num unless $num;
    return floor($num+0.5);
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
