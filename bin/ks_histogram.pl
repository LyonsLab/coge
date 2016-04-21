#!/usr/bin/perl -w

use strict;
use GD;
use CoGe::Accessory::histogram;
use Data::Dumper;
use DBI;
use Getopt::Long;
use POSIX;
use vars qw($help $ks_db $ks_type $outfile $width $height $log $pair_file $chr1 $chr2 $max $min $color_scheme $fontsize);

GetOptions(
	   "ks_db|db=s"=>\$ks_db,
	   "ks_type|kst=s"=>\$ks_type,
	   "out|o=s"=>\$outfile,
	   "help|h"=>\$help,
	   "width|w=i"=>\$width,
	   "height|hi=i"=>\$height,
	   "log"=>\$log,
	   "pair_file|pf=s"=>\$pair_file,
	   "chr1|c1=s"=>\$chr1,
	   "chr2|c2=s"=>\$chr2,
	   "max=s"=>\$max,
	   "min=s"=>\$min,
	   "color_scheme=s"=>\$color_scheme,
	   "fontsize|fs=i" => \$fontsize,
	   );

usage() if $help;
usage("Unable to read database\n") unless -r $ks_db;
$ks_type = "kS" unless $ks_type;
$ks_type = uc($ks_type) if $ks_type;
$width = 500 unless $width;
$height = 500 unless $height;
$fontsize = 8 unless $fontsize;
my $pairs = get_pairs(file=>$pair_file, chr1=>$chr1, chr2=>$chr2) if $pair_file && -r $pair_file;
my $x_title = "substitution per site for $ks_type";
$x_title = "log10() ".$x_title if $log;

my ($data) = get_ksdata(ks_db=>$ks_db, type=>$ks_type, pairs=>$pairs);
unless (@$data) {
    warn("No data was found in the database");
    system("touch $outfile");
    return;
}

my @data;
my ($val_max, $val_min, $non_zero_min) = range( $data);
foreach my $val (@$data)
  {
    my $tmp = $val;
    if ($log)
      {
	$tmp = $non_zero_min if ($val == 0);
	$tmp = log10($tmp);
      }
    next if defined $min && $tmp < $min;
    next if defined $max && $tmp > $max;
    push @data, $tmp;
  }

my $mean =0;
map {$mean+=$_} @data;
$mean = sprintf("%.4f",$mean/scalar(@data));
@data = sort {$a<=>$b} @data;

my $median = sprintf("%.4f",$data[floor(scalar(@data/2))]);

my $hist = new CoGe::Accessory::histogram($width, $height);
my $bins = CoGe::Accessory::histogram::_histogram_bins($data, 100);
my $colors = gen_color_list(bins=>$bins, color_scheme=>$color_scheme);
my $count=0;
my @color_names;
foreach my $color (@$colors)
  {
    my $name = "c$count";

    $hist->add_colour($name=>$color);
    push @color_names, $name;
    $count++;
  }
my $title = "Mean: $mean" if defined $mean;
$title .= " Median: $median" if defined $median;
$hist->set(histogram_bins=>100,
	   x_label_skip=>5,
	   y_label =>'counts',
	   x_label=>$x_title,
	   cycle_clrs=>1,
	   y_long_ticks=>1,
	   y_tick_number=>10,
	   title=>$title,
	  );
$hist->set(dclrs=>[@color_names]) if scalar @color_names;
$hist->set_title_font("/usr/lib/perl5/CoGe/fonts/arial.ttf",$fontsize);
$hist->set_x_label_font("/usr/lib/perl5/CoGe/fonts/arial.ttf",$fontsize);
$hist->set_y_label_font("/usr/lib/perl5/CoGe/fonts/arial.ttf",$fontsize);
$hist->set_x_axis_font("/usr/lib/perl5/CoGe/fonts/arial.ttf",$fontsize);
$hist->set_y_axis_font("/usr/lib/perl5/CoGe/fonts/arial.ttf",$fontsize);
$hist->set_values_font("/usr/lib/perl5/CoGe/fonts/arial.ttf",$fontsize);
my $gd = $hist->plot(data=>\@data);
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
    my %opts = @_;
    my $bins = $opts{bins};
    my $color_scheme = $opts{color_scheme};
    my $range = $bins->[-1][0]-$bins->[0][0];
    return unless $range;
    my @colors;
    foreach my $item (@$bins)
      {
	my $val = sprintf("%.4f", ($item->[0]-$bins->[0][0])/$range);
	push @colors, get_color(val=>$val, color_scheme=>$color_scheme);
      }
    return \@colors;
  }

sub get_ksdata
  {
    my %opts = @_;
    my $ks_db = $opts{ks_db};
    my $type= $opts{type};
    my $pairs = $opts{pairs};
    return [] unless -r $ks_db;
    my $select = "select * from ks_data";
    my $dbh = DBI->connect("dbi:SQLite:dbname=$ks_db","","");
    my $sth = $dbh->prepare($select);
    $sth->execute();
    my @data;
    while (my $data = $sth->fetchrow_arrayref)
      {
	if ($pairs)
	  {
	    next unless $pairs->{$data->[1]}{$data->[2]};
	  }
#	print STDERR $data->[3],"\n" unless $data->[3];
	next unless $data->[3] =~ /\d/;
	my %data = (
		    KS=>$data->[3],
		    KN=>$data->[4],
		    'KN_KS'=>$data->[5],
		   );
	push @data, $data{$type} if $data{$type};
      }
    $sth->finish();
    undef $sth;
    $dbh->disconnect();
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
    my $color_scheme = $opts{color_scheme};
    $color_scheme=1 unless defined $color_scheme;
    return [0,0,0] unless defined $val;
   my $schemes = [
		   [
		    [255,255,0], #yellow
		    [200,200,0], # green
		    [0,200,0], # green
		    [0,100,100], # green
		    [0,200,200], # cyan
		    [0,0,200], # blue
		    [100,0,100], #purple
		    [200,0,200], #magenta
		    [200,0,0], #red
		    [100,0,0], #red
		    [200,100,0], #orange
		    [255,126,0], #orange
		    ],
		   [
		    [255,255,0], #yellow
		    [255,0,0], #red
		    [0,255,0], # green
		    [0,255,255], # cyan
		    [220,0,220], #magenta
		    [0,0,255], # blue
		   ],
		   [
		    [0,0,150], # blue
		    [220,220,20], #yellow
		    [255,0,0], #red
		   ],
		   [
		    [0,200,0], # green
		    [0,0,200], # blue
		    [220,220,20], #yellow
		    [255,0,0], #red
		   ],
		   [
		    [255,0,0], # red
		    [0,0,0], #black
		   ],
		   [
		    [255,0,0], #red
		    [255,255,0], #yellow
		    [0,255,0], # green
		    [0,255,255], # cyan
		    [220,0,220], #magenta
		    [0,0,255], # blue
		    [255,0,0], #red
		    [255,255,0], #yellow
		    [0,255,0], # green
		    [0,255,255], # cyan
		    [220,0,220], #magenta
		    [0,0,255], # blue
		   ],
		   [
		    [255,0,0], #red
		    [255,255,0], #yellow
		    [0,255,0], # green
		    [0,255,255], # cyan
		    [220,0,220], #magenta
		    [0,0,255], # blue
		    [255,0,0], #red
		    [255,255,0], #yellow
		    [0,255,0], # green
		    [0,255,255], # cyan
		    [220,0,220], #magenta
		    [0,0,255], # blue
		    [255,0,0], #red
		    [255,255,0], #yellow
		    [0,255,0], # green
		    [0,255,255], # cyan
		    [220,0,220], #magenta
		    [0,0,255], # blue
		   ],
		   [
		    [255,0,0], #red
		    [255,155,0], #yellow
		    [0,255,0], # green
		    [0,255,255], # cyan
		    [220,0,220], #magenta
		    [0,0,255], # blue
		    [255,0,0], #red
		    [255,155,0], #yellow
		    [0,255,0], # green
		    [0,255,255], # cyan
		    [220,0,220], #magenta
		    [0,0,255], # blue
		   ],
		   [
		    [255,0,0], #red
		    [255,155,0], #yellow
		    [0,255,0], # green
		    [0,255,255], # cyan
		    [220,0,220], #magenta
		    [0,0,255], # blue
		    [255,0,0], #red
		    [255,155,0], #yellow
		    [0,255,0], # green
		    [0,255,255], # cyan
		    [220,0,220], #magenta
		    [0,0,255], # blue
		    [255,0,0], #red
		    [255,155,0], #yellow
		    [0,255,0], # green
		    [0,255,255], # cyan
		    [220,0,220], #magenta
		    [0,0,255], # blue
		   ],
		  [
		   [0,255,255], # blue
		   [255,0,0], #orange
		   [0,255,255], # blue
		   [255,0,0], #orange
		   [0,255,255], # blue
		   [255,0,0], #orange
		  ],
		  [
		   [0,0,255], # blue
		   [255,99,33], #orange
		   [0,0,255], # blue
		   [255,99,33], #orange
		   [0,0,255], # blue
		   [255,99,33], #orange
		  ],
		 ];
    my @colors = @{$schemes->[$color_scheme]};
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
	next unless $_;
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

This program generates a histogram of kS data stored in an sqlite database generated by CoGe's SynMap program.

Options:

-ks_db    |    db      the sqlite database

-ks_type  |    kst     what data to plot (default:  dS)  Options dS, dN, dS/dN

-out      |    o       output file.  STDOUT by default

-help     |    h       print this message
};
    exit;
  }
