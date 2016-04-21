#!/usr/bin/perl -w

use strict;
use GD;
use CoGe::Accessory::histogram;
use Data::Dumper;
use DBI;
use Getopt::Long;
use POSIX;

use vars qw($help $outfile $width $height $log $max $min $file $x_title $hist_type $calc_mean);

GetOptions(
	   "out|o=s"=>\$outfile,
	   "help|h"=>\$help,
	   "width|w=i"=>\$width,
	   "height|hi=i"=>\$height,
	   "log"=>\$log,
	   "max=s"=>\$max,
	   "min=s"=>\$min,
	   "file|f=s"=>\$file,
	   "title|t=s"=>\$x_title,
	   "hist_type|ht=s"=>\$hist_type,
	   "mean"=>\$calc_mean,
	   );
$hist_type = "count" unless $hist_type && $hist_type eq "percentage";
usage() if $help;
usage("Unable to read data file\n") unless -r $file;
$width = 500 unless $width;
$height = 500 unless $height;

my $data = get_data(file=>$file);
my $mean;
if ($calc_mean)
  {
    $mean = 0;
    map {$mean+=$_} @$data;
    $mean = $mean/scalar(@$data);
  }
if ($log)
  {
    my ($max, $min, $non_zero_min) = range( $data);
#    print STDERR join ("\t", $min, $max),"\n";
    my @data;
    foreach my $val (@$data)
      {
	$val = $non_zero_min if ($val == 0);
	push @data, log10($val);
      }
    $x_title = "log10() ".$x_title;
    $data = \@data;
  }
#check data if necessary
if (defined $min || defined $max)
  {
    my @data;
    foreach my $val (@$data)
      {
	next if defined $min && $val < $min;
	next if defined $max && $val > $max;
	push @data, $val;
      }
    $data = \@data;
  }

my $hist = new CoGe::Accessory::histogram($width, $height);
my $bins = CoGe::Accessory::histogram::_histogram_bins($data, 100);
my $colors = gen_color_list($bins);
my $count=0;
my @color_names;
foreach my $color (@$colors)
  {
    my $name = "c$count";
    $hist->add_colour($name=>$color);
    push @color_names, $name;
    $count++;
  }
my $title= "Mean: $mean" if defined $mean;

$hist->set(histogram_bins=>100,
	   x_label_skip=>5,
	   y_label =>$hist_type,
	   x_label=>$x_title,
	   cycle_clrs=>1,
	   y_long_ticks=>1,
	   y_tick_number=>10,
	   histogram_type=>$hist_type,
	   title=>$title,
	  );
$hist->set(dclrs=>[@color_names]) if scalar @color_names;
$hist->set_x_label_font("/usr/lib/perl5/site_perl/CoGe/fonts/arial.ttf",8);
$hist->set_y_label_font("/usr/lib/perl5/site_perl/CoGe/fonts/arial.ttf",8);
$hist->set_x_axis_font("/usr/lib/perl5/site_perl/CoGe/fonts/arial.ttf",8);
$hist->set_y_axis_font("/usr/lib/perl5/site_perl/CoGe/fonts/arial.ttf",8);
$hist->set_values_font("/usr/lib/perl5/site_perl/CoGe/fonts/arial.ttf",8);
my $gd = $hist->plot(data=>$data, min=>$min, max=>$max);
if ($outfile)
  {
    $outfile .= ".png" unless $outfile =~ /png$/;
    open (OUT, ">$outfile");
    binmode OUT;
    print OUT $gd->png;
    close OUT;
  }
else
  {
    print $gd->png;
  }

sub gen_color_list
  {
    my $bins = shift;
    my $range = $bins->[-1][0]-$bins->[0][0];
    return unless $range;
    my @colors;
    foreach my $item (@$bins)
      {
	my $val = sprintf("%.4f", ($item->[0]-$bins->[0][0])/$range);
	push @colors, get_color(val=>$val);
      }
    return \@colors;
  }

sub get_data
  {
    my %opts = @_;
    my $file = $opts{file};
    my @data;
    return \@data unless -r $file;
    open (IN, $file);
    while (<IN>)
      {
	chomp;
	next unless $_;
	next if /^#/;
	push @data, $_;
      }
    return \@data;
  }

sub range
  {
    my $data = shift;
    my $min =undef;
    my $max =undef;
    my $non_zero_min=undef;
    foreach my $val (@$data)
      {
	$val = 0 if $val == 0;
	$min = $val unless defined $min;
	$non_zero_min = $val unless defined $non_zero_min || $val == 0;
	$max = $val unless defined $max;
	$min = $val if $min > $val;
	$max = $val if $max < $val;
	$non_zero_min = $val if $val > 0 && $non_zero_min > $val;
      }
    return ($max, $min, $non_zero_min);
  }

sub get_color
  {
    my %opts = @_;
    my $val = $opts{val};
    return [0,0,0] unless defined $val;
    my @colors = (
		  [255,0,0], #red
		  [255,126,0], #orange
		  [255,255,0], #yellow
		  [0,255,0], # green
		  [0,255,255], # cyan
		  [0,0,255], # blue
#		  [255,0,255], #magenta
		  [126,0,126], #purple
		 );
    @colors = reverse @colors;
    my ($index1, $index2) = ((floor((scalar(@colors)-1)*$val)), ceil((scalar(@colors)-1)*$val));

    my $color=[];
    my $step = 1/(scalar (@colors)-1);
    my $scale = $index1*$step;
    my $dist = ($val-$scale)/($step);
    for (my $i=0; $i<=2; $i++)
      {
	my $diff = ($colors[$index1][$i]-$colors[$index2][$i])*$dist;
	push @$color, sprintf("%.0f", $colors[$index1][$i]-$diff);
      }
    return $color;
  }

sub get_pairs
  {
    my %opts = @_;
    my $file = $opts{file};
    my $chr1 = $opts{chr1};
    my $chr2 = $opts{chr2};
    my %data;
    open (IN, $file) || die $!;
    while (<IN>)
      {
	chomp;
	next if /^#/;
	my @line = split/\t/;
	my @item1 = split/\|\|/, $line[1];
	my @item2 = split/\|\|/, $line[5];
	if ($chr1)
	  {
	    next unless $item1[0] eq $chr1 || $item2[0] eq $chr1;
	  }
	if ($chr2)
	  {
	    next unless $item1[0] eq $chr2 || $item2[0] eq $chr2;
	  }
	$data{$item1[6]}{$item2[6]}=1;
	$data{$item2[6]}{$item1[6]}=1;
      }
    close IN;
    return \%data;
  }
sub usage
  {
    print qq{
welcome to $0!

This program generates a histogram of numerical data stored in a file.

Options:

-file     |    f       data file.  Each data line contains a number.  Lines beginning with '#'
                       are skipped

-width    |    w       graph width in pixels (DEFAULT 500)

-height   |    hi      graph height in pixels  (DEFAULT 500)

-out      |    o       output file.  STDOUT by default

-log                   log transfor output (NO OPTIONS-- just flag to set)

-hist_type |   ht      type for histogram:  "count" or "percentage"  DEFAULT count

-title    |    t       Title for histogram

-min                   Minimum value for data

-max                   Maximum value for data

-mean                  flag to calculate the mean of the values

-hist_type |  ht       Histogram of data counts or distribution of percentages.  Values: count percentage.  Default: count

-help     |    h       print this message

};
    exit;
  }
