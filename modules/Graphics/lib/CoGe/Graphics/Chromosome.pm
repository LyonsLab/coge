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

CoGe::Graphics::Chromosome - Object for drawing chromosomes that provides functionality for painting chromosomes with location based features (such as genes), magnification/zooming on particular regions, and printing pictures of chromosomes as pngs.

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
  $c->labels(0);

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
    use vars qw($VERSION $DEFAULT_WIDTH $PADDING $DEFAULT_COLOR $MAX_MAG $MAG_SCALE_TYPE $MAG_STEP_HEIGHT $CHR_INNER_COLOR $CHR_OUTER_COLOR $SCALE_COLOR $TICK_COLOR $SCALE_HEIGHT $FONT $FONTTT $NUM_MAG $FEATURE_HEIGHT);
    $VERSION     = '0.1';
    $DEFAULT_WIDTH = 200;  #default image width pixels
    $PADDING = 15; #default value to pad image height
    $MAG_STEP_HEIGHT = 30; #the amount to increase the height of the chromosome picture for each magnification (features increase half this height)
    $MAX_MAG = 10;     #default number of "units" (e.g. base pairs) to show when maximally zoomed)
    $MAG_SCALE_TYPE = "log"; #default image scaling.  options "log", "linear"
    $DEFAULT_COLOR = [0,0,0];
    $CHR_INNER_COLOR = [220,255,220]; #inner color for chromosome
    $CHR_OUTER_COLOR = [0,0,0]; #border color for chromosome
    $SCALE_COLOR = [0,0,255]; #color for measurement scale
    $TICK_COLOR = [0,0,255]; #color for ticks on measurement scale
    $SCALE_HEIGHT = 20; #height (in pixels) of the scale
    $FONTTT = "/usr/lib/perl5/site_perl/CoGe/fonts/arial.ttf"; #path to true-type font
    $FONT = GD::Font->MediumBold; #default GD font
    $FEATURE_HEIGHT = 4; #the heigth of a feature is determined by this number * $self->magnification.  if the magnification is 5 and the feature_height is set to 4, then the resulting feature will be 20 pixels high.
    $NUM_MAG = 10; #number of magnification steps;
    __PACKAGE__->mk_accessors(
"DEBUG",
"chr_length",
"feature_height", #height of a feature
"draw_scale", #flag for drawing scale
"scale_color", "tick_color", #color for scale and ticks on scale respectively
"scale_height", #height of scale
"mag_scale_type", "max_mag", "mag_step_height", "num_mag",
"image_width", "image_height", 
"padding",
"font",
"labels", "fill_labels", #flag to turn off the printing of labels, fill_lables are specifically for filled features;

"_region_start", "_region_stop", #image's start and stop (should be equal to or "larger" than the users
"_magnification", "_mag_scale", 
"_image_h_used", #storage for the amount of the image height used
"_gd", #store GD object
"draw_chromosome", "chr_inner_color", "chr_outer_color",
 "_chr_brush",
"_chr_center", "_chr_height", "_chr_h1", "_chr_h2", #interal storage of chromosome image height positions
"_features", #internal storage of features


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
    $self->mag_step_height($MAG_STEP_HEIGHT);
    $self->image_width($DEFAULT_WIDTH);
    $self->padding ($PADDING);
    $self->chr_inner_color($CHR_INNER_COLOR);
    $self->chr_outer_color($CHR_OUTER_COLOR);
    $self->draw_chromosome(1);
    $self->draw_scale(1);
    $self->scale_height($SCALE_HEIGHT);
    $self->scale_color($SCALE_COLOR);
    $self->tick_color($TICK_COLOR);
    $self->num_mag($NUM_MAG);
    $self->feature_height($FEATURE_HEIGHT);
    $self->mag(5);
    $self->font($FONTTT);
    $self->_features([]);
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


=cut

#################### subroutine header end ####################


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
	print $self->gd->png;
      }
    $self->_gd(undef);
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
    return $gd->colorResolve(@colors)
  }

sub set_image_height
  {
    #this sub calculates how much height out picture needs based on options and features to be drawn
    my $self = shift;
    my $feat_height = $self->feature_height || $FEATURE_HEIGHT;
    $feat_height = $feat_height*($self->mag);
    my $h = $self->padding; #give use some padding
    $h += $self->scale_height+$self->padding;# if $self->draw_scale; #scale image height
    my $chrh = $self->mag * $self->mag_step_height+$self->mag/2+$self->padding;# if $self->draw_chromosome; #chromosome image height
    my $top_feat = $self->get_feats(last=>1, strand=>1, fill=>0);
    my $tfh = $top_feat->order * ($feat_height+$self->padding)+2*$self->padding if $top_feat;
    my $bot_feat = $self->get_feats(last=>1, strand=>-1, fill=>0);
    my $bfh = $bot_feat->order * ($feat_height+$self->padding)+2*$self->padding if $bot_feat;
    $h += $tfh > $chrh/2 ? $tfh : $chrh/2;
    $self->_chr_center($h);
    $h += $bfh > $chrh/2 ? $bfh : $chrh/2;
    print STDERR "Image Height: $h\n" if $self->DEBUG;
    $self->ih($h);
    $self->_image_h_used($self->padding);

  }


sub generate_region
  {
    my $self = shift;
    my %opts = @_;
    print STDERR "\n", join "\t", "mag: ".$self->mag, "region start: ".$self->_region_start,"region end: ".$self->_region_stop,"\n" if $self->DEBUG;
#    print Dumper $self if $self->DEBUG;
    $self->set_image_height();
    $self->gd->fill(0,0,$self->get_color([255,255,255]));
    $self->_draw_scale;# if $self->draw_scale;
    $self->_draw_chromosome;# if $self->draw_chromosome;
    $self->_draw_features;
#    $self->gd->fill($self->iw/2,$self->_image_h_used +( $self->ih - $self->_image_h_used)/2+1, $self->get_color($self->chr_inner_color));
  }

sub _draw_features
  {
    my $self = shift;
    my $c = $self->_image_h_used+($self->ih - $self->_image_h_used)/2;
    print STDERR "Image used: ".$self->_image_h_used."  Image Height: ".$self->ih."  Center: $c\n" if $self->DEBUG;
    foreach my $feat (sort {$b->fill <=> $a->fill} $self->get_features)
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
	my $feat_h = $self->feature_height*$self->mag;
	my $offset = ($feat->order-1)*($feat_h+$self->padding)+$self->padding;
	$offset = 0 if $feat->fill;
	$feat_h = ($self->_chr_height-$self->mag)/2 if $feat->fill;
	my $y = $feat->strand =~ /-/ ? $c+ $offset+1: $c - $offset-$feat_h-1;
	
#	print STDERR "Feature offset: $y, Order: ", $feat->order,"\n";
        my $sy;
	if ($feat->fill)
	  {
	    $sy = $feat->strand =~ /-/ ? $c+2 : $c-$self->padding;
	  }
	$self->_draw_feature(feat=>$feat, 'y'=>$y, ih=>$feat_h, 'sy'=>$sy);
	#may need to fill in color with $bgcolor. . .
	if ($y > ($self->_chr_center-$self->_chr_height/2) && $y < ($self->_chr_center+$self->_chr_height/2))
	  {
	    #$self->gd->fill($self->iw/2,$y-1, $self->get_color($self->chr_inner_color));
	  }
      }
  }
sub _draw_feature
  {
    my $self = shift;
    my %opts = @_;
    my $feat = $opts{feat} || $opts{FEAT} || $opts{f};
    return 0 unless ref($feat) =~ /Feature/i;
    my $y = $opts{'y'} || $opts{Y};
    my $ih = $opts{'image_height'} || $opts{'ih'} || $opts{'IH'} || $feat->ih;
    my $sy = $opts{'string_y'} || $opts{'sy'} || $y; #label y axis
    my $rb = $self->_region_start;
    my $re = $self->_region_stop;
    my $range = $re-$rb;
    my $w = $self->iw;
    $feat->gd;
    $feat->stop($feat->start) unless defined $feat->stop;
    my $feat_range = $feat->stop-$feat->start;

    my $unit = $w/$range;
    my $fs = int($w* ($feat->start-$rb)/$range);#+($unit/2);
    my $fe = int($w* ($feat->end-$rb)/$range+$unit); 
    my $fw = int($fe - $fs); #calculate the width of the feature;
    return if $fw < 1; #skip drawing if less than one pix wide
    print STDERR "Drawing feature ".$feat->label.": ", $feat->start, "-", $feat->end,"Dimentions:",$fw,"x",$ih, " at position: $fs,$y"."\n" if $self->DEBUG;
    $self->gd->copyResampled($feat->gd, $fs, $y,0,0, $fw, $ih, $feat->iw, $feat->ih);
  #  $self->_gd_string(y=>$y+$feat->ih*.9, x=>$fs, text=>$feat->label, size=>10) if ($self->labels && $fw>5); #don't make the string unless the feature is at least 5 pixels wide
    $self->_gd_string(y=>$y, x=>$fs, text=>$feat->label, size=>15) if ($self->labels && $fw>5); #don't make the string unless the feature is at least 5 pixels wide
    $self->_gd_string(y=>$sy, x=>$fs, text=>$feat->label, size=>15) if ($self->fill_labels && $feat->fill && $fw>5); #don't make the string unless the feature is at least 5 pixels wide
    return $y+$feat->ih+$self->padding;
  }

sub _calc_unit_size
  {
    my $self = shift;
    return (($self->iw/($self->_region_stop-$self->_region_start)));
  }

sub add_feature
  {
    my $self = shift;
    my @feats;
    my $sort = 0;
    foreach (@_)
      {
	push @feats, ref($_) =~ /array/i ? @$_ : $_;
      }
    foreach my $feat (@feats)
      {
	$feat->strand(1) unless defined $feat->strand;
	$feat->fill(0) unless $feat->fill;
	$feat->stop($feat->start) unless defined $feat->stop;
	$feat->gd; #initialize feature;
	if (ref($feat) =~ /Feature/i)
	  {
	    unless ($feat->order)
	      {
		my $last_feat = $self->get_feats(last=>1, strand=>$feat->strand, fill=>$feat->fill);
		my $order = $last_feat ? $last_feat->order()+1 : 1;
		$feat->order($order);
	      }
	    if ($sort)
	      {
	        my @feats = sort {$a->order <=> $b->order} @{$self->_features}, $feat;
	        $self->_features(\@feats);
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

sub get_features
  {
    my $self = shift;
    my %opts = @_;
    my $order = $opts{order} || $opts{ORDER};
    my $type = $opts{type} || $opts{TYPE};
    my $last = $opts{last} || $opts{LAST}; #flag to get the highest order feature for a location
    my $strand = $opts{strand} || $opts{STRAND};
    my $fill = $opts{fill};
    $fill = $opts{FILL} unless defined $fill; #get filled features?
    return unless $self->_features;
    my @rfeats;
    foreach my $feat (@{$self->_features})
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
	    next unless $feat->type eq $type;
	  }
	if (defined $fill)
	 {
	   next unless $feat->fill eq $fill;
	 }
	push @rfeats, $feat
      }
    @rfeats = sort {$a->order <=> $b->order} @rfeats;

    if ($last)
      {
	return $rfeats[-1];
      }
    return wantarray ? @rfeats : \@rfeats;
  }


sub get_feature
  {
    my $self = shift;
    my %opts = @_;
    return $self->get_features(%opts);
  }

sub get_feats
  {
    my $self = shift;
    my %opts = @_;
    return $self->get_features(%opts);
  }

sub _draw_scale
  {
    my $self = shift;
    my $gd = $self->gd;
    my $c = $self->scale_height/2+$self->_image_h_used; #center of scale
    my $w = $self->iw; #width of image
    my $xb = $self->_region_start < 1 ? $w*abs($self->_region_start)/($self->_region_stop - $self->_region_start): 0; #x position begin
    my $xe = $self->_region_stop > $self->chr_length ? $w-$w*($self->_region_stop - $self->chr_length)/($self->_region_stop - $self->_region_start) : $w; #x position end
    $gd->line($xb,$c,$xe,$c,$self->get_color($self->scale_color));
    my $mtyb = $c-$self->scale_height/2; #major tick y begin
    my $mtye = $c+$self->scale_height/2; #major tick y end
    my $styb = $c-$self->scale_height/4; #minor tick y begin
    my $stye = $c+$self->scale_height/4; #minor tick y end
    my $rb = $self->_region_start;
    $rb = 1 if $rb < 1;
    my $re = $self->_region_end;
    $re = $self->chr_length if $re > $self->chr_length;
    my $range = $re-$rb; #chromosomal positional range
    my $div = "1"."0"x int log10($range); #determine range scale (10, 100, 1000, etc)
    print STDERR "\nSCALE: Center: $c, Start $xb, Stop: $xe, Ticks: $div, \n" if $self->DEBUG;
    $self->_make_ticks(scale=>$div*10, y1=>$mtyb, y2=>$mtye, range_begin=>$rb, range_end=>$re,text_loc=>1);
#    $div /= 10;
    $self->_make_ticks(scale=>$div, y1=>$styb, y2=>$stye, range_begin=>$rb, range_end=>$re, text_loc=>-1);
    $self->_image_h_used($self->_image_h_used + $self->scale_height+$self->padding/2);
  }

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
    if ($rb == 1)
      {
	my $x = $w *($rb - $self->_region_start)/($self->_region_stop - $self->_region_start);
	$gd->filledRectangle($x, $y1, $x+$unit, $y2, $self->get_color($self->tick_color));
	if ($text_loc)
	  {
	    my $h = $text_loc =~ /-/ ? $y2-1: $y1-$self->padding/2;
	    $self->_gd_string(text=>$rb,x=>$x+$unit+2, 'y'=>$h, size => ($y2-$y1)/1.5);
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
	    $key = "1".$key;
	    my %end = (10=>"0",
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
	    $self->_gd_string(text=>$t,x=>$x+$unit+2, 'y'=>$h, size => ($y2-$y1));#/1.5  );
	  }
	$tick+= $div;
      }

  }

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

sub _draw_chromosome
  {
    my $self = shift;
    my $gd = $self->gd;
    my $ch = $self->mag * $self->mag_step_height/2; #half chromosome image height
    my $hc = $self->_image_h_used+($self->ih-$self->_image_h_used)/2; #height of center of chromsome image
    my $w = $self->iw; #width of image
    return unless $self->draw_chromosome;
    print STDERR "\nChromosome image: Height/2: $ch, Height Center: $hc\n" if $self->DEBUG;
    my $xs = $self->_region_start < 1 ? $w*abs($self->_region_start)/($self->_region_stop - $self->_region_start): 0;
    my $xe = $self->_region_stop > $self->chr_length ? $w-$w*($self->_region_stop - $self->chr_length)/($self->_region_stop - $self->_region_start) : $w;
#    $gd->filledRectangle($xs,$hc-$ch,$xe, $hc+$ch,$self->get_color(@{$self->chr_inner_color}));
    $self->_chr_height($ch*2);
    $gd->setBrush($self->chr_brush);
    $gd->line($xs, $hc-$ch, $xe, $hc-$ch, gdBrushed);
    $self->chr_brush->flipVertical();
    $gd->setBrush($self->chr_brush);
    $gd->line($xs, $hc+$ch, $xe, $hc+$ch, gdBrushed);
    my $color = $self->get_color($self->chr_outer_color);
    $gd->setStyle($color, $color, $color, GD::gdTransparent, GD::gdTransparent);
    $gd->line($xs, $hc, $xe, $hc, gdStyled);
    $self->draw_chr_end ($xs, "left", $hc) if ($xs > 0);
    $self->draw_chr_end ($xe, "right", $hc) if ($xe < $w);
    
  }

sub chr_brush
  {
    my $self = shift;
    return $self->_chr_brush if $self->_chr_brush;
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

sub draw_chr_end
  {
    my $self = shift;
    my $x = shift;
    my $dir = shift || "left";
    my $hc = shift; #$ch/2+$self->_image_h_used+$self->mag/2; #height of center of chromsome image
    my $gd = $self->gd;
    my $ch = $self->mag * $self->mag_step_height; #chromosome image height
    
    my @arc1 = $dir =~ /left/i ? (90, 180) : (0, 90);
    my @arc2 = $dir =~ /left/i ? (180, 270) : (270, 0);
    my @arc3 = $dir =~ /left/i ? (90, 270) : (270,90);
    $gd->filledArc($x, $hc, $ch, $ch, @arc3, $self->get_color($self->chr_inner_color));
    $gd->arc($x, $hc, $ch, $ch, @arc1, gdBrushed);
    $self->chr_brush->flipVertical;
    $gd->setBrush($self->chr_brush);
    $gd->arc($x, $hc, $ch, $ch, @arc2, gdBrushed);
  }

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
	else #log
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
	$self->_mag_scale(\%scale);
      }
    return $self->_mag_scale;
  }

sub set_region
  {
    #user sets either a range with a start and stop, or just a start point.
    #if only a start point is set, then we assume that will be the center of the view
    my $self = shift;
    my %opts = @_;
    my $start = $opts{start} || $opts{begin} || $opts{START} || $opts{BEGIN} 
      || $opts{center} || $opts{CENTER} 
	|| $opts{point} || $opts{POINT} || 1;
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

sub _set_region_for_point
  {
    my $self = shift;
    my $point = shift;
    unless (defined $point)
      {
	warn 'You must specify a point!\n';
	return (0);
      }
    my $range = ceil ($self->find_size_for_magnification()/2);
    my $start = $point - $range;
    $self->_region_start($start);
    my $stop = $point + $range;
    $self->_region_stop($stop);
  }


sub _set_region_start_stop
  {
    my $self = shift;
    my ($start, $end) = @_;
    return 0 unless (defined $start && defined $end);
    if ($end < $start)
      {
	my $tmp = $start;
	$start = $end;
	$end = $tmp;
      }
    my $len = $end - $start;
    my $mag = $self->find_magnification_for_size($len);
    $self->mag($mag);
#    $self->start($start);
#    $self->stop($end) if $end;
    my $size = $self->mag_scale->{$mag};
    my $diff = ceil( ($size-$len)/2);
    my $rstart = $start-$diff;
    my $rend = $end + $diff;
    $self->_region_start($rstart);
    $self->_region_stop($rend);
  }


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
	$mag = $magt if $mag_scale->{$magt} > $len;
      }
#    return $mag;
    $self->magnification($mag);
    return $self->magnification();
  }

sub find_size_for_magnification
  {
    my $self = shift;
    my $mag = shift || $self->mag();
    my $scale = $self->mag_scale;
    return 0 unless $scale && $mag;
    return $scale->{$mag};
  }

sub magnification
  {
    my $self = shift;
    my $mag = shift;
    #are we changing magnification?  if so, we need to set the region start and end points

    $mag = $self->num_mag if $mag && $mag > $self->num_mag;
    $self->_magnification($mag) if $mag;
    if ($self->_region_start && $mag)
      {
	my $point = $self->_region_start;
	$point += ceil (($self->_region_stop - $self->_region_start)/2) if $self->_region_stop;
	$self->set_point($point);
      }
    
    return $self->_magnification();
  }

sub mag
  {
    my $self = shift;
    return $self->magnification(@_);
  }

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


sub _region_end
  {
    my $self = shift;
    return $self->_region_stop(@_);
  }

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

