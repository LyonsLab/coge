package CoGe::Graphics::GenomeView;

use base qw(Class::Accessor);
use strict;
use Data::Dumper;
use GD;

#################### main pod documentation begin ###################
##
##

=head1 NAME  GenomeView

 CoGe::Graphics::GenomeView

=head1 SYNOPSIS

 my $gv = new CoGe::Graphics::GenomeView({color_band_flag=>1, image_width=>800, chromosome_height=>20}) ;
 $gv->add_chromosome(name=>"Chr: 1",
		     end=>30000000,
		     );
 $gv->add_chromosome(name=>"Chr: 10",
		     end=>10000000,
		     );
 #add a series of tick marks at various positions
 my $up = 0; #is the tick on top or bottom of chromosome?
 foreach (my $i=1;$i<=30000000,$i+=100000)
  {
    $gv->add_feature(name=>"Mark at position $i,
  		     start=>$i,
  		     stop=>$i,
  		     chr=>"Chr: 1", #must match one of the specified chromosomes
  		     imagemap=>qq/"Mark at position $i" onclick="your_js_function(1, $i);"/,
  		     up=>$up,
  		     color=>[255,0,0],
  		    );
    $up = $up == 1 ? 0 : 1;
  }

 foreach (my $i=1;$i<=10000000,$i+=1000000)
  {
    $gv->add_feature(name=>"Mark at position $i,
  		     start=>$i,
  		     stop=>$i,
  		     chr=>"Chr: 10",
  		     imagemap=>qq/"Mark at position $i" onclick="your_js_function(1, $i);"/,
  		     up=>$up,
  		     color=>[0,255,0],
  		    );
    $up = $up == 1 ? 0 : 1;
  }

 my $file = $gv->generate_png(filename=>"image.png");
 my $map = $gv->generate_imagemap(mapname=>"genomeview_imagemap");

 print qq{Content-Type: text/html

<html>
<head></head>
<body>
<img src=$file ismap usemap='#genomeview_imagemap' border=0>$map
</body>
</html>
};

=head1 DESCRIPTION

 GenomeView creates images that give a graphical overview of an entire chromosome that can be painted
 with tick marks.  These tick marks can be used to show where blast hits fall on a chromosome, where
 genes reside, etc.  Each tick mark can also be used as a link in an imagemap to other web-pages or
 javascript functions.

=head1 USAGE

 use CoGe::Graphics::GenomeView;
 my $gv = CoGe::Graphics::GenomeView->new();

=head1 BUGS
 Documentation for this module is lacking

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
 GD
=cut

#for best performance, create all the chromosomes before generating the features.

__PACKAGE__->mk_accessors(qw(organism chromosomes features image_width image_height legend_height _default_feature_color _gd _color_set color_band_flag legend chromosome_height show_count draw_ends max_count));

my $DEFAULT_COLOR = [255,100,100];
my $FONT = GD::Font->MediumBold;
my $FONTTT="/usr/local/fonts/arial.ttf"; #needs to be fixed to get from conf file

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
	my $locs = $self->imagemap_features(chr=>$chrs->{$name}, x=>$horz_spacer, 'y'=>$count*$vert_spacer, width=>$width, height=>$height);
	$map .= $locs if $locs;
	$count++;
      }
    $map .= qq!</map>\n!;
    return $map;
  }

sub imagemap_features
  {
    my ($self) = shift;
    my %opts = @_;
    my $chr =$opts{chr};
    my $x = sprintf("%.0f",$opts{x});
    my $y = sprintf("%.0f",$opts{'y'});
    my $width = $opts{width};
    my $height = $opts{height};
    my $onchange = $opts{onchange} || 0;
    my $feats = $self->features();
    my $up = 1;
    my $color_band_flag = $self->color_band_flag;
    my $js_1 = qq!onMouseOver="change('!;
    my $js_2 = qq!')" !;
    my $map;
    return unless $feats->{$chr->{name}};
    my %points;
    foreach my $feat (sort {$a->{end} <=> $b->{end} }@{$feats->{$chr->{name}}})
      {
	my $x = sprintf("%.0f",($x+$feat->{end}/$chr->{end}*$width));
	$points{$x}{0}{count}=0 unless defined $points{$x}{0}{count};
	$points{$x}{1}{count}=0 unless defined $points{$x}{1}{count};
	$points{$x}{$feat->{up}}{count}++;
	$points{$x}{$feat->{up}}{feat} = $feat;

      }
    foreach my $x (sort {$points{$a}{0}{count}+$points{$a}{1}{count} <=> $points{$b}{0}{count}+$points{$b}{1}{count}} keys %points)
      {
	my $x1 = sprintf("%.0f",$x-$height/5);
	my $x2 = sprintf("%.0f",$x+$height/5);
	foreach my $up (keys %{$points{$x}})
	 {
	   next unless $points{$x}{$up}{count};

	   my $feat = $points{$x}{$up}{feat};
	   if ($up)
	     {
	       my $y1 = sprintf("%.0f",$y);
	       my $y2 = sprintf("%.0f",$y-$height/1.3);
	       $map .= qq!<area shape="rect" coords="$x1,$y1, $x2,$y2" !;
	       $up = 0;
	     }
	   else
	     {
	       my $y1 = sprintf("%.0f",$y+$height);
	       my $y2 = sprintf("%.0f",$y+$height+$height/1.3);
	       $map .= qq!<area shape="rect" coords="$x1,$y1, $x2,$y2" !;
	       $up = 1;
	     }
	   $map .= "\n";
	   $map .= qq!href="$feat->{link}" ! if $feat->{link};
	   if ($onchange)
	     {
	       $map .= $js_1;
	       $map .= $feat->{name} ."\\n";
	       $map .= $feat->{desc} ."\\n" if $feat->{desc};
	       $map .= "Start:  ".$feat->{start}."\\n";
	       $map .= "End:    ".$feat->{end}."\\n";
	       $map .= "Length: ".($feat->{end}-$feat->{start})."\\n";
	       $map .= $js_2."\n";
	     }
	   $map .= " ".$feat->{imagemap}." " if $feat->{imagemap};
	   $map .= qq!alt="$feat->{name}">\n!;

	   if ($color_band_flag)
	     {
	       my $xt = sprintf("%.0f",($x+$feat->{end}/$chr->{end}*$width));
	       my $yt = sprintf("%.0f",($y+$height));
	       $map .= qq!<area shape="rect" coords="$xt,$y, $xt,$yt" !;
	       $map .= "\n";
	       $map .= qq!href="$feat->{link}" ! if $feat->{link};
	       if ($onchange)
		 {
		   $map .= $js_1;
		   $map .= $feat->{name} ."\\n";
		   $map .= $feat->{desc} ."\\n" if $feat->{desc};
		   $map .= "Start:  ".$feat->{start}."\\n";
		   $map .= "End:    ".$feat->{end}."\\n";
		   $map .= "Length: ".($feat->{end}-$feat->{start})."\\n";
		   $map .= $js_2."\n";
		 }
	       $map .= " ".$feat->{imagemap}." " if $feat->{imagemap};
	       $map .= qq!alt="$feat->{name}">\n!;
	     }
	  }
      }
    return $map;
  }

sub generate_png
  {
    my ($self) = shift@_;
    my %opts = @_;
    my $gif = $opts{gif} || 0;
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
    $self->_gd(undef);
    $self->_color_set(undef);
    $self->generate_chromosomes();
    $self->generate_legend() if $self->legend;
    if ($gif)
      {
	print OUT $self->gd->gif;
      }
    else
      {
	print OUT $self->gd->png;
      }
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
	my $real_width = $horz_spacer+$width;
	my $pos_word_length = (length $chrs->{$name}->{end})*($vert_spacer/5);
	my $offset_width = ($real_width+$pos_word_length) - $self->image_width;
	#print STDERR "offset width: ",$offset_width,"\n";
	$real_width -= ($offset_width - 15) if $offset_width > 0;
	my ($cstart, $cend) = ($chrs->{$name}->{centromere_start}, $chrs->{$name}->{centromere_end});
	$gd->rectangle($horz_spacer, $count*$vert_spacer, $horz_spacer+$width, $count*$vert_spacer+$height, $black);
	$gd->arc($horz_spacer, $count*$vert_spacer+$height/2, $height, $height, 90, 270, $black) unless defined $self->draw_ends && $self->draw_ends == 0;
	$gd->arc($horz_spacer+$width, $count*$vert_spacer+$height/2, $height, $height, 270, 90, $black) unless defined $self->draw_ends && $self->draw_ends == 0;
	$gd->line($horz_spacer, $count*$vert_spacer, $horz_spacer, $count*$vert_spacer+$height, $white);
	$gd->line($horz_spacer+$width, $count*$vert_spacer, $horz_spacer+$width, $count*$vert_spacer+$height, $white);
        $gd->stringFT($black, $FONTTT, $vert_spacer/5, 0, 5, $count*$vert_spacer, $name);
	$gd->stringFT($black, $FONTTT, $vert_spacer/5, 0, $real_width, $count*$vert_spacer, $chrs->{$name}->{end});
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
    my $black = $self->get_color([0,0,0]);
    my $white = $self->get_color([255,255,255]);
    my %points;
    my $max_count =1;
    foreach my $feat (@{$feats->{$chr->{name}}})
      {
	my $x1 = sprintf("%.0f",$x+$feat->{end}/$chr->{end}*$width);
	$points{$x1}{0}{count}=0 unless defined $points{$x1}{0}{count};
	$points{$x1}{1}{count}=0 unless defined $points{$x1}{1}{count};
	$points{$x1}{$feat->{up}}{count}++;
	$max_count = $points{$x1}{$feat->{up}}{count} if $points{$x1}{$feat->{up}}{count} > $max_count;
	$points{$x1}{$feat->{up}}{color} = $feat->{color};
	$points{$x1}{$feat->{up}}{heatmap} = $feat->{heatmap};

      }
    $self->max_count($max_count);
    foreach my $x1 (sort {$points{$a}{0}{count}+$points{$a}{1}{count} <=> $points{$b}{0}{count}+$points{$b}{1}{count}} keys %points)
     {
       foreach my $up (sort keys %{$points{$x1}})
	 {
	   my $count = $points{$x1}{$up}{count};
	   next unless $count;

	   my $color;
	   if (defined $points{$x1}{$up}{heatmap})
	     {
	       my $c1 = sprintf("%.0f",100 + 155*$count/$max_count);
	       $color = [0,0,0];
	       $color->[$points{$x1}{$up}{heatmap}]=$c1;
	     }
	   else
	     {
	       $color = $points{$x1}{$up}{color};
	     }
	   $color = $self->get_color($color) if ref ($color) =~ /array/i;
	   $color = $self->default_feature_color unless $color;

	   $y = sprintf("%.0f", $y);
	   my $poly = new GD::Polygon;
	   if ($up)
	     {
	       $poly->addPt($x1, $y);
	       $poly->addPt(sprintf("%.0f",$x1-$height/5), sprintf("%.0f",$y-$height/1.3));
	       $poly->addPt(sprintf("%.0f",$x1+$height/5), sprintf("%.0f",$y-$height/1.3));
	       $up = 0;
	     }
	   else
	     {
	       $poly->addPt($x1, sprintf("%.0f",$y+$height));
	       $poly->addPt(sprintf("%.0f",$x1-$height/5), sprintf("%.0f",$y+$height+$height/1.3));
	       $poly->addPt(sprintf("%.0f",$x1+$height/5), sprintf("%.0f",$y+$height+$height/1.3));
	       $up = 1;
	     }
	   $gd->filledPolygon($poly, $color);
	   $gd->openPolygon($poly, $black);
	   my $font_size = $self->chromosome_height/5;
	   my $h = $up? sprintf("%.0f",$height+$height/1.3) : sprintf("%.0f",0-$height/1.3+$font_size);
	   if ($self->show_count)
	     {
	       $gd->stringFT($white, $FONTTT, $font_size+6, 0, sprintf("%.0f",$x1-$height/5), $y+$h, $count) unless $count == 1;
	       $gd->stringFT($black, $FONTTT, $font_size+3, 0, sprintf("%.0f",$x1-$height/5+1), $y+$h-1, $count) unless $count == 1;
	     }
	   if ($color_band_flag)
	     {
	       $gd->line($x1, $y, $x1, $y+$height, $color);
	     }
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
    $self->chromosomes({}) unless defined $self->chromosomes;
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
    my $desc = $opts{'desc'};
    my $up = $opts{'up'} || 0;
    my $imagemap = $opts{imagemap};
    my $heatmap = $opts{heatmap};
    $self->features({}) unless $self->features;

    unless ($start && $end && $chr)
      {
	warn "add_feature call failed -- must have valid start, end, and chromosome";
	return 0;
      }
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
				     imagemap=>$imagemap,
				     heatmap=>$heatmap,
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
