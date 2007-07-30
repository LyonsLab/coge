package CoGe::Graphics::GenomeView;

use base qw(Class::Accessor);
use strict;
use Data::Dumper;
use GD;


#for best performance, create all the chromosomes before generating the features.

__PACKAGE__->mk_accessors(qw(organism chromosomes features image_width image_height legend_height _default_feature_color _gd _color_set color_band_flag legend chromosome_height));

my $DEFAULT_COLOR = [255,0,0];
my $FONT = GD::Font->MediumBold;
my $FONTTT="/usr/lib/perl5/site_perl/CoGe/fonts/arial.ttf";


sub generate_imagemap
  {
    my ($self) = shift @_;
    my %opts = @_;
    my $map_name = $opts{map_name};
    $map_name = $opts{mapname} unless $map_name;
    $map_name = $opts{name} unless $map_name;
    $map_name = shift unless $map_name;
    my $map = qq!<map name="$map_name">\n!;
    my $max;
    my $chrs = $self->chr;
    my $vert_spacer = ($self->image_height()-$self->legend_height())/(keys (%$chrs) +1);
    my $horz_spacer = $self->image_width()/10;
    my $count = 1;
    foreach my $name (sort {$chrs->{$b}->{end} <=> $chrs->{$a}->{end}} keys %$chrs)
      {
	my $height = $vert_spacer-$vert_spacer/1.5;
	$max = $chrs->{$name}->{end} unless $max;
	my $width = ($horz_spacer*8)*($chrs->{$name}->{end}/$max);
	my $locs = $self->imagemap_features($chrs->{$name}, $horz_spacer, $count*$vert_spacer, $width, $height);
	$map .= $locs if $locs;
	$count++;
      }
    $map .= qq!</map>\n!;
    return $map;
  }

sub imagemap_features
  {
    my ($self, $chr, $x, $y, $width, $height) = @_;
    my $feats = $self->features();
    my $up = 1;
    my $color_band_flag = $self->color_band_flag;
    my $js_1 = qq!onMouseOver="change('!;
    my $js_2 = qq!')" !;
    my $map;
    return unless $feats->{$chr->{name}};
    foreach my $feat (sort {$a->{end} <=> $b->{end} }@{$feats->{$chr->{name}}})
      {
	my $x1 = ($x+$feat->{end}/$chr->{end}*$width)-$height/5;
	my $x2 = ($x+$feat->{end}/$chr->{end}*$width)+$height/5;
	$up = $feat->{up} if defined $feat->{up};
	if ($up)
	  {
	    my $y1 = $y;
	    my $y2 = $y-$height/1.3;
	    $map .= qq!<area coords="$x1, $y1, $x2, $y2" !;
	    $up = 0;
	  }
	else
	  {
	    my $y1 = $y+$height;
	    my $y2 = $y+$height+$height/1.3;
	    $map .= qq!<area coords="$x1, $y1, $x2, $y2" !;
	    $up = 1;
	  }
	$map .= "\n";
	$map .= qq!href="$feat->{link}" ! if $feat->{link};
	$map .= $js_1;
	$map .= $feat->{name} ."\\n";
	$map .= $feat->{desc} ."\\n" if $feat->{desc};
	$map .= "Start:  ".$feat->{start}."\\n";
	$map .= "End:    ".$feat->{end}."\\n";
	$map .= "Length: ".($feat->{end}-$feat->{start})."\\n";
	$map .= $js_2."\n";
	$map .= qq!alt="$feat->{name}">\n!;
	if ($color_band_flag)
	  {
	    my $xt = ($x+$feat->{end}/$chr->{end}*$width)+2;
	    my $yt = $y+$height;
	    $map .= qq!<area coords="$xt, $y, $xt, $yt" !;
	    $map .= "\n";
	    $map .= qq!href="$feat->{link}" ! if $feat->{link};
	    $map .= $js_1;
	    $map .= $feat->{name} ."\\n";
	    $map .= $feat->{desc} ."\\n" if $feat->{desc};
	    $map .= "Start:  ".$feat->{start}."\\n";
	    $map .= "End:    ".$feat->{end}."\\n";
	    $map .= "Length: ".($feat->{end}-$feat->{start})."\\n";
	    $map .= $js_2."\n";
	    $map .= qq!alt="$feat->{name}">\n!;
	  }
      }
    return $map;
  }

sub generate_png
  {
    my ($self) = shift@_;
    my %opts = @_;
    my $file_name = $opts{file_name};
    $file_name = $opts{filename} unless $file_name;
    $file_name = $opts{name} unless $file_name;
    $file_name = shift unless $file_name;
    my $color_band_flag = $opts{band_flag};
    $color_band_flag = $opts{bandflag} unless defined $color_band_flag;
    $color_band_flag = shift unless defined $color_band_flag;
    $color_band_flag = $self->color_band_flag unless defined $color_band_flag;
    $color_band_flag = 0 unless defined $color_band_flag;
    $self->color_band_flag($color_band_flag);
    open (OUT, ">$file_name") || die "Can't open $file_name for writing: $!";
    binmode OUT;
#    my $gd = $self->gd;
    $self->generate_chromosomes();
    $self->generate_legend() if $self->legend;
    print OUT $self->gd->png;
    close OUT;
  }

sub get_color
  {
    my ($self, $color) = @_;
    $self->_color_set({}) unless $self->_color_set(); #initialize
    my $color_set = $self->_color_set();
    $color_set->{$color->[0]}{$color->[1]}{$color->[2]} = $self->gd->colorAllocate(@$color) unless ($color_set->{$color->[0]}{$color->[1]}{$color->[2]});
    return $color_set->{$color->[0]}{$color->[1]}{$color->[2]};
  }

sub generate_chromosomes
  {
    my $self = shift;
    my $gd = $self->gd;
    my $max;
    my $chrs = $self->chr;
    my $vert_spacer = ($self->image_height()-$self->legend_height)/(keys (%$chrs) +1);
    my $horz_spacer = $self->image_width()/10;
    my $count = 1;
    my $black = $self->get_color([0,0,0]);
    my $white = $self->get_color([255,255,255]);
    my $TITLE_FONT = gdGiantFont;
    $gd->string($TITLE_FONT, $self->image_width/2-(length ($self->organism)/2*$TITLE_FONT->width), $TITLE_FONT->height, $self->organism, $black) if $self->organism;
    

    foreach my $name (sort {$chrs->{$b}->{end} <=> $chrs->{$a}->{end}} keys %$chrs)
      {
	my $height = $vert_spacer-$vert_spacer/1.5;
	$max = $chrs->{$name}->{end} unless $max;
	my $width = ($horz_spacer*8)*($chrs->{$name}->{end}/$max);
	my ($cstart, $cend) = ($chrs->{$name}->{centromere_start}, $chrs->{$name}->{centromere_end});
	$gd->rectangle($horz_spacer, $count*$vert_spacer, $horz_spacer+$width, $count*$vert_spacer+$height, $black);
	$gd->arc($horz_spacer, $count*$vert_spacer+$height/2, $height, $height, 90, 270, $black);
	$gd->arc($horz_spacer+$width, $count*$vert_spacer+$height/2, $height, $height, 270, 90, $black);
	$gd->line($horz_spacer, $count*$vert_spacer+1, $horz_spacer, $count*$vert_spacer+$height-1, $white);
	$gd->line($horz_spacer+$width, $count*$vert_spacer+1, $horz_spacer+$width, $count*$vert_spacer+$height-1, $white);
        $gd->stringFT($black, $FONTTT, $vert_spacer/5, 0, 5, $count*$vert_spacer, $name);
#	$gd->string($FONT, 5, $count*$vert_spacer-$FONT->height, $name, $black);
#	$gd->string($FONT, $horz_spacer-$FONT->width*(length($chrs->{$name}{start})+1), $count*$vert_spacer-$FONT->height, $chrs->{$name}->{start}, $black);
	$gd->stringFT($black, $FONTTT, $vert_spacer/5, 0, $horz_spacer+$width, $count*$vert_spacer, $chrs->{$name}->{end});
	#$gd->string($FONT, $horz_spacer+$width+10, $count*$vert_spacer-$FONT->height, $chrs->{$name}->{end}, $black);
	$self->draw_centromere($horz_spacer, $count*$vert_spacer, $width, $height, $cstart/$chrs->{$name}->{end}, $cend/$chrs->{$name}->{end}, $black) if $cstart && $cend;
	$self->draw_features($chrs->{$name}, $horz_spacer, $count*$vert_spacer, $width, $height);

	$count++;
      }
    $self->_gd($gd);
  }

sub draw_centromere
  {
    my ($self, $x, $y, $width, $height, $cstart, $cend, $color) =@_;
    my $gd = $self->gd;
    my $cwidth = $cend-$cstart;
    $gd->filledEllipse($x+($width*(($cend-$cstart)/2+$cstart)), $y+$height/2, $width*($cend-$cstart), $height, $color);
  }

sub draw_features
  {
    my ($self, $chr, $x, $y, $width, $height) = @_;
    my $feats = $self->features();
    my $gd = $self->gd;
    my $up = 1;
    my $color_band_flag = $self->color_band_flag;
    return unless $feats->{$chr->{name}};
    my $black = $gd->colorAllocate(0,0,0);
    foreach my $feat (sort {$a->{end} <=> $b->{end} }@{$feats->{$chr->{name}}})
      {
	my $color = $feat->{color};
	$color = $self->get_color($color) if ref ($color) =~ /array/i;
	$color = $self->default_feature_color unless $color;

	my $x1 = $x+$feat->{end}/$chr->{end}*$width;
	my $poly = new GD::Polygon;
	$up = $feat->{up} if defined $feat->{up};
	if ($up)
	  {
	    $poly->addPt($x1, $y);
	    $poly->addPt($x1-$height/5, $y-$height/1.3);
	    $poly->addPt($x1+$height/5, $y-$height/1.3);
	    $up = 0;
	  }
	else
	  {
	    $poly->addPt($x1, $y+$height);
	    $poly->addPt($x1-$height/5, $y+$height+$height/1.3);
	    $poly->addPt($x1+$height/5, $y+$height+$height/1.3);
	    $up = 1;
	  }
	$gd->filledPolygon($poly, $color);
	$gd->openPolygon($poly, $black);
	if ($color_band_flag)
	  {
	    $gd->line($x1, $y, $x1, $y+$height, $color);
	  }
      }
  }

sub generate_legend
  {
    my $self = shift;
    my $h = 40 * (scalar (keys %{$self->chr})*2);
    my $legend = $self->legend;
    my $black = $self->get_color([0,0,0]);
    foreach my $str (keys %$legend)
      {
	my $color = $self->get_color($legend->{$str});
	my $poly = new GD::Polygon;
	$poly->addPt(5, $h);
	$poly->addPt(5, $h+$FONT->height);
	$poly->addPt(5+$FONT->height*2, $h+$FONT->height/2);
	$self->gd->filledPolygon($poly, $color);
	$self->gd->openPolygon($poly, $black);
	$self->gd->string($FONT, 5+$FONT->height*2, $h, $str, $black);
	$h += $FONT->height;
      }
  }
sub org
  {
    my ($self, $id) = @_;
    $self->organism($id) if $id;
    return ($self->organism );
  }

sub species
  {
    my ($self, $id) = @_;
    $self->organism($id) if $id;
    return ($self->organism );
  }

sub chr
  {
    my ($self, $id) = @_;
    $self->chromosomes($id) if $id;
    return ($self->chromosomes );
  }

sub gd
  {
    my ($self) = @_;
    my $gd = $self->_gd();
    unless ($gd)
      {
	$self->image_width(800) unless $self->image_width;
	if ($self->legend)
	  {
	    my $legend_hei = scalar keys (%{$self->legend}) * $FONT->height if keys %{$self->legend};
	    $self->legend_height($legend_hei);
	    $self->image_height($self->image_height() + $self->legend_height);
	  }
	else
	  {
	    $self->legend_height(0);
	  }
	$self->image_height(40 * (scalar (keys %{$self->chr})*2)) unless $self->image_height;
        $self->image_height(((keys (%{$self->chr}) +1)*$self->chromosome_height+$self->legend_height)*1.5) if $self->chromosome_height;

	my ($wid, $hei) = ($self->image_width, $self->image_height);
	$gd = new GD::Image($wid, $hei);
	$gd->colorAllocate(255,255,255);
	$self->_gd($gd);
	
      }

    return $gd;
  }

sub default_feature_color
  {
    my ($self, $color) = @_;
    $self->_default_feature_color($color) if $color && ref ($color) =~ /GD/;
    unless ($self->_default_feature_color)
      {
	$self->_default_feature_color($self->get_color($DEFAULT_COLOR));
      }
    return $self->_default_feature_color;
  }

sub add_feature
  {
    my ($self) = shift;
    my (%opts) = @_;
    my $name = $opts{'name'};
    my $start = $opts{'start'};
    $start = $opts{'begin'} unless $start;
    $start = 1 unless $start;
    my $end = $opts{'end'};
    $end = $opts{'stop'} unless $end;
    my $chr = $opts{'chr'};
    $chr = $opts{'chrom'} unless $chr;
    $chr = $opts{'chromosome'} unless $chr;
    $chr = $self->get_chromosome($chr);
    my $link = $opts{'link'};
    my $color = $opts{'color'};
#    $color = $self->get_color($color) if ref ($color) =~ /array/i;
#    $color = $self->default_feature_color unless $color;
    my $desc = $opts{'desc'};
    my $up = $opts{'up'};
    unless ($start && $end && $chr)
      {
	warn "add_feature call failed -- must have valid start, end, and chromosome";
	return 0;
      }
    $self->features({}) unless $self->features;
    my $feats = $self->features();
    push @{$feats->{$chr->{name}}}, {
				     name=>$name,
				     start=>$start,
				     end=>$end,
				     chr=>$chr,
				     color=>$color,
				     link=>$link,
				     desc=>$desc,
				     up=>$up,
				 };
    return 1;
  }

sub add_chromosome
  {
    my ($self) = shift;
    my (%opts) = @_;
    my $name = $opts{'name'};
    my $start = $opts{'start'};
    $start = $opts{'begin'} unless $start;
    $start = 1 unless $start;
    my $end = $opts{'end'};
    $end = $opts{'stop'} unless $end;
    my $cstart = $opts{'centromere_start'};
    $cstart = $opts{'cen_start'} unless $cstart;
    $cstart = $opts{'cstart'} unless $cstart;
    my $cend = $opts{'centromere_end'};
    $cend = $opts{'cen_end'} unless $cend;
    $cend = $opts{'cend'} unless $cend;
     #are the start/end positions based on morgans or nucleotide (nt) positions
#    my $metric = $opts{'metric'};
#    $metric = shift unless $metric;
#    $metric = "nt" unless $metric;
    unless ($name && $start && $end)
      {
	warn "add_chromosome call failed -- must have valid name, start and end\n";
	return 0;
      }
    $self->chr({}) unless $self->chr;
    my $chr = $self->chr();
    $chr->{$name} = {
		     start=>$start, 
		     end=>$end,
		     name=>$name,
		     centromere_start=>$cstart,
		     centromere_end=>$cend,
		    };
    return 1;
  }

sub get_chromosome
  {
    my ($self, $chr) = @_;
    my $chrs = $self->chr;
    return $chrs->{$chr};
  }

1;
