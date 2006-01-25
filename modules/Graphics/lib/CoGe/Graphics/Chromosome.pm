package CoGe::Graphics::Chromosome;
use strict;
use base qw(Class::Accessor);
use POSIX;
use Data::Dumper;
use GD;

BEGIN {
    use vars qw($VERSION $DEFAULT_WIDTH $PADDING $DEFAULT_COLOR $MAX_MAG $MAG_SCALE_TYPE $MAG_STEP_HEIGHT $CHR_INNER_COLOR $CHR_OUTER_COLOR $SCALE_COLOR $TICK_COLOR $SCALE_HEIGHT $FONT $FONTTT $NUM_MAG $FEATURE_HEIGHT);
    $VERSION     = '0.1';
    $DEFAULT_WIDTH = 200;  #default image width pixels
    $PADDING = 10; #default value to pad image height
    $MAG_STEP_HEIGHT = 20; #the amount to increase the height of the chromosome picture for each magnification
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
    $FEATURE_HEIGHT = 10;
    $NUM_MAG = 20; #number of magnification steps;
    __PACKAGE__->mk_accessors(
"DEBUG",
"chr_length",
"start", "stop", #user defined start and stop
"_region_start", "_region_stop", #image's start and stop (should be equal to or "larger" than the users
"_magnification", "_mag_scale", "mag_scale_type", "max_mag", "mag_step_height", "num_mag",
"features", "feature_height",
"image_width", "image_height", "padding",
"_image_h_used", #storage for the amount of the image height used
"_gd", #store GD object
"draw_chromosome", "chr_inner_color", "chr_outer_color", "_chr_brush",
"_chr_h1", "_chr_h2", #interal storage of chromosome image height positions
"draw_scale", "scale_color", "tick_color", "scale_height",
"_features",
);
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
    return $gd->colorResolve(@colors)
  }

sub set_image_height
  {
    #this sub calculates how much height out picture needs based on options and features to be drawn
    my $self = shift;
    my $h = 3*$self->padding; #give use some padding
    $h += $self->mag * $self->mag_step_height+$self->mag/2+$self->padding;# if $self->draw_chromosome; #chromosome image height
    $h += $self->scale_height+$self->padding;# if $self->draw_scale; #scale image height
    my $top_feat = $self->get_feats(last=>1, strand=>1, loc=>"out");
    $h += $top_feat->order * ($self->feature_height+$self->padding)+$self->padding*2;
    my $bot_feat = $self->get_feats(last=>1, strand=>-1, loc=>"out");
    $h += $bot_feat->order * ($self->feature_height+$self->padding)+$self->padding*2;

#    my %max;
#    foreach my $feat ($self->get_features(loc=>'out'))
#      {
#	$seen{$feat->order} = $feat->ih unless $seen{$feat->order};
#	$seen{$feat->order} = $feat->ih if $seen{$feat->order} < $feat->ih;
#      }
#    foreach my $ih (values %seen)
#      {
#	$h += $ih+$self->padding;
#      }
    print STDERR "\nImage Height: $h\n" if $self->DEBUG;
    $self->ih($h);
    $self->_image_h_used($self->padding);

  }


sub generate_region
  {
    my $self = shift;
    my %opts = @_;
    print STDERR "\n", join "\t", "mag: ".$self->mag, "region start: ".$self->_region_start,"region end: ".$self->_region_stop,"\n" if $self->DEBUG;
    print STDERR Dumper $self if $self->DEBUG;
    $self->set_image_height();
    $self->gd->fill(0,0,$self->get_color([255,255,255]));
    $self->_draw_scale;# if $self->draw_scale;
    $self->_draw_outer_top_features();
    $self->_draw_chromosome;# if $self->draw_chromosome;
    $self->_draw_inner_features();
    $self->_draw_outer_bottom_features();
  }

sub _draw_inner_features
  {
    my $self = shift;
    my $ch1 = $self->_chr_h1; #chromosome start height
    my $ch2 = $self->_chr_h2; #chromosome end height
    my $ih = ($ch2-$ch1) / 2; #height of chromosome picture
    my $ty = $ch1+$self->padding; #y position for top half
    my $by = $ih+$ch1+$self->padding; #y position for bottom half 
    foreach my $feat (sort {$b->fill <=> $a->fill} @{$self->get_features(loc=>'in')})
      {
	my $strand = $feat->strand =~ /-/ ? 0 : 1;
	my $h = $feat->fill ? $ih :$feat->ih;
	my $y = $strand ? $ty: $by;
	if ($feat->fill)
	  {
	    $y = $strand ? $ch1: $ch1+$ih;
	  }
	$self->_draw_feature(feat=>$feat, 'y'=>$y, ih=>$h);
	$ty += $feat->ih unless $feat->fill;
	$by += $feat->ih unless $feat->fill;
      }
  }

sub _draw_outer_top_features
  {
    my $self = shift;
    my $h = $self->_image_h_used;
    foreach my $feat (sort {$a->order <=> $b->order} @{$self->get_features(loc=>'out', strand=>"1")})
      {
	$h = $self->_image_h_used + ($self->feature_height+$self->padding) * $feat->order;
	print STDERR "Top feature height: $h\n";
	$self->_draw_feature(feat=>$feat, 'y'=>$h)
      }
    $h += $self->feature_height+$self->padding*3;
    print STDERR "End of top feature height: $h\n";
    $self->_image_h_used($h);
  }

sub _draw_outer_bottom_features
  {
    my $self = shift;
    my $h = $self->_image_h_used;
    foreach my $feat (@{$self->get_features(loc=>'out')})
      {
	if ($feat->strand =~ /-/)
	  {
	    $h = $self->_draw_feature(feat=>$feat, 'y'=>$h)
	  }
      }
    $self->_image_h_used($h);
  }

sub _draw_feature
  {
    my $self = shift;
    my %opts = @_;
    my $feat = $opts{feat} || $opts{FEAT} || $opts{f};
    return 0 unless ref($feat) =~ /Feature/i;
    my $y = $opts{'y'} || $opts{Y};
    my $ih = $opts{'image_height'} || $opts{'ih'} || $opts{'IH'} || $feat->ih;
    my $rb = $self->_region_start;
    my $re = $self->_region_stop;
    my $range = $re-$rb;
    my $w = $self->iw;
    $feat->stop($feat->start) unless defined $feat->stop;
    my $feat_range = $feat->stop-$feat->start;

    my $unit = $w/$range;
    my $fs = $w* ($feat->start-$rb)/$range;#+($unit/2);
    my $fe = $w* ($feat->end-$rb)/$range+$unit; 
    my $fw = $fe - $fs; #calculate the width of the feature;
    print STDERR "Drawing feature ".$feat->label.": ", $feat->start, "-", $feat->end,"Dimentions:",$fw,"x",$ih, "\n" if $self->DEBUG;
    $self->gd->copyResized($feat->gd, $fs, $y,0,0, $fw, $ih, $feat->iw, $feat->ih);
    $self->_gd_string(y=>$y+$feat->ih*.9, x=>$fs, text=>$feat->label, size=>10) if $fw>5; #don't make the string unless the feature is at least 5 pixels wide
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
    $self->_features({'in'=>[],
		     'out'=>[],
		    }) unless $self->_features(); #initialize if needed
    my @feats;
    foreach (@_)
      {
	push @feats, ref($_) =~ /array/i ? @$_ : $_;
      }
    foreach my $feat (@feats)
      {
	$feat->strand(1) unless defined $feat->strand;
	$feat->fill(0) unless $feat->fill;
	$feat->gd; #initialize gd object;
	if (ref($feat) =~ /Feature/i)
	  {
	    unless ($feat->order)
	      {
		my $last_feat = $self->get_feats(loc=>$feat->placement, last=>1, strand=>$feat->strand);
		my $order = $last_feat ? $last_feat->order()+1 : 1;
		$feat->order($order);
	      }
	    if ($feat->placement =~ /in/i)
	      {
		push @{$self->_features->{in}},$feat;
	      }
	    else
	      {
		push @{$self->_features->{out}},$feat;
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
    my $loc = $opts{location} || $opts {loc} || $opts{LOCATION} || $opts{LOC};
    my $last = $opts{last} || $opts{LAST}; #flag to get the highest order feature for a location
    my $strand = $opts{strand} || $opts{STRAND};
    my @feats = $loc ? @{$self->_features->{$loc}} : @{$self->_features->{in}}, @{$self->_features->{out}};
    return wantarray ? @feats : \@feats unless ($order || $type || $strand || $last);
    my @rfeats;
    foreach my $feat (@feats)
      {
	if ($strand)
	  {
	    next unless $feat->strand eq $strand;
	  }
	if ($order)
	  {
	    next unless $feat->order eq $order;
	  }
	if ($type)
	  {
	    next unless $feat->type eq $type;
	  }
	push @rfeats, $feat
      }
    if ($last)
      {
	($last) = sort {$b->order <=> $a->order} @rfeats;
	return $last;
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
    
    $self->_make_ticks($div*10, $mtyb, $mtye, $rb, $re,1);
#    $div /= 10;
    $self->_make_ticks($div, $styb, $stye, $rb, $re, -11);
    $self->_image_h_used($self->_image_h_used + $self->scale_height+$self->padding/2);
  }

sub _make_ticks
  {
    my $self = shift;
    my ($div, $y1, $y2, $rb, $re, $text) = @_;
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
	if ($text)
	  {
	    my $h = $text =~ /-/ ? $y2-1: $y1-1;
	    $self->_gd_string(text=>$rb,x=>$x+$unit+2, 'y'=>$h, size => ($y2-$y1)/2  )
	  }
      }
    while ($tick <= $re)
      {
	my $x = $w *($tick - $self->_region_start)/($self->_region_stop - $self->_region_start);
	print STDERR "\nGenerating tick at $tick ($x)\n" if $self->DEBUG;
	$gd->filledRectangle($x, $y1, $x+$unit, $y2, $self->get_color($self->tick_color));
	if ($text)
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
	    my $h = $text =~ /-/ ? $y2-1: $y1-1;
	    $self->_gd_string(text=>$t,x=>$x+$unit+2, 'y'=>$h, size => ($y2-$y1)/2  );
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
    my $x = $opts{x} || $opts{X};
    my $y = $opts{'y'} || $opts{Y};
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

sub _draw_chromosome
  {
    my $self = shift;
    my $gd = $self->gd;
    my $ch = $self->mag * $self->mag_step_height/2; #half chromosome image height
    my $hc = $ch+$self->_image_h_used+$self->mag/2; #height of center of chromsome image
    my $w = $self->iw; #width of image
    $self->_image_h_used($self->_image_h_used + $ch*2 + $self->mag + $self->padding/2);
    return unless $self->draw_chromosome;
    print STDERR "\nChromosome image: Height/2: $ch, Height Center: $hc\n" if $self->DEBUG;
    my $xs = $self->_region_start < 1 ? $w*abs($self->_region_start)/($self->_region_stop - $self->_region_start): 0;
    my $xe = $self->_region_stop > $self->chr_length ? $w-$w*($self->_region_stop - $self->chr_length)/($self->_region_stop - $self->_region_start) : $w;
    $gd->filledRectangle($xs,$hc-$ch,$xe, $hc+$ch,$self->get_color(@{$self->chr_inner_color}));
    $self->_chr_h1($hc-$ch+$self->mag/2+1);
    $self->_chr_h2($hc+$ch-$self->mag/2);
    $gd->setBrush($self->chr_brush);
    $gd->line($xs, $hc-$ch, $xe, $hc-$ch, gdBrushed);
    $self->chr_brush->flipVertical();
    $gd->setBrush($self->chr_brush);
    $gd->line($xs, $hc+$ch, $xe, $hc+$ch, gdBrushed);
    $self->draw_chr_end ($xs, "left") if ($xs > 0);
    $self->draw_chr_end ($xe, "right") if ($xe < $w);

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
    my $gd = $self->gd;
    my $ch = $self->mag * $self->mag_step_height; #chromosome image height
    my $hc = $ch/2+$self->_image_h_used+$self->mag/2; #height of center of chromsome image
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
	    my $rang = $self->chr_length;
	    for my $i (1..$self->num_mag)
	      {
		$rang = ceil ($rang)/3;
		$rang = $self->max_mag if $rang < $self->max_mag;
		$scale{$i} = $rang;
	      }
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
	|| $opts{point} || $opts{POINT} 
	  ||$self->start || 0;
    my $end = $opts{stop} || $opts{end} || $opts{STOP} || $opts{END} || $self->stop;
    $self->start($start);
    $self->stop($end);
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
#    print STDERR "\nsetting point: $point\n";
#    print STDERR "\n", $self->_region_start,"\t",$self->_region_stop,"\n";
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
#    $start = 1 if $start < 1;
    $self->_region_start($start);
    my $stop = $point + $range;
#    $stop = $self->chr_length if $stop > $self->chr_length;
    $self->_region_stop($stop);
  }


sub _set_region_start_stop
  {
    my $self = shift;
    my ($start, $end) = @_;
    return 0 unless ($start && $end);
    if ($end < $start)
      {
	my $tmp = $start;
	$start = $end;
	$end = $tmp;
      }
    my $len = $end - $start;
    my $mag = $self->find_magnification_for_size($len);
    $self->mag($mag);
    $self->start($start);
    $self->stop($end) if $end;
    my $size = $self->mag_scale->{$mag};
    my $diff = ceil( ($size-$len)/2);
    my $rstart = $start-$diff;
#    $rstart = 1 if $rstart < 1;
    my $rend = $end + $diff;
#    $rend = $self->chr_length if $rend > $self->chr_length;
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
    return $self;
}


#################### main pod documentation begin ###################
## Below is the stub of documentation for your module. 
## You better edit it!


=head1 NAME

CoGe::Graphics::Chromosome - CoGe::Graphics::Chromosome

=head1 SYNOPSIS

  use CoGe::Graphics::Chromosome;
  blah blah blah


=head1 DESCRIPTION

Stub documentation for this module was created by ExtUtils::ModuleMaker.
It looks like the author of the extension was negligent enough
to leave the stub unedited.

Blah blah blah.


=head1 USAGE



=head1 BUGS



=head1 SUPPORT



=head1 AUTHOR

	Eric Lyons
	CPAN ID: MODAUTHOR
	XYZ Corp.
	elyons@nature.berkeley.edu
	http://a.galaxy.far.far.away/modules

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

