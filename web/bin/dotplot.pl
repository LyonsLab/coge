#!/usr/bin/perl -w

use strict;
use GD;
use Getopt::Long;
use CoGeX;
use Data::Dumper;

use vars qw($dagfile $alignfile $width $link $min_chr_size $orgid1 $orgid2 $help $coge $gd $CHR1 $CHR2 $basename $link_type);


GetOptions(
	   "dagfile|d=s"=>\$dagfile,
	   "alignfile|a=s"=>\$alignfile,
	   "width|w=i"=>\$width,
	   "link|l=s"=>\$link,
	   "link_type|lt=i"=>\$link_type,
	   "min_chr_size|mcs=i"=>\$min_chr_size,
	   "orgid1|o1=i"=>\$orgid1,
	   "orgid2|o2=i"=>\$orgid2,
	   "help|h" =>\$help,
	   "chr1|c1=s"=>\$CHR1,
	   "chr2|c2=s"=>\$CHR2,
	   "basename|b=s"=>\$basename,
	   
	   );
usage() if $help;
usage() unless -r $dagfile;

$basename = "test" unless $basename;
$width = 1024 unless $width;

$coge = CoGeX->dbconnect();


my $org1info = get_org_info(oid=>$orgid1, chr=>$CHR1, minsize=>$min_chr_size);
my $org1length =0;
map {$org1length+=$_->{length}} values %$org1info;
my $org2info = get_org_info(oid=>$orgid2, chr=>$CHR2, minsize=>$min_chr_size);
my $org2length =0;
map {$org2length+=$_->{length}} values %$org2info;

my $height = sprintf("%.0f", $width*$org2length/$org1length);
$height = $width if ($height > 5*$width);
my $x_bp_per_pix = sprintf("%.0f", $org1length/$width);
my $x_pix_per_bp = 1/$x_bp_per_pix;
my $y_bp_per_pix = sprintf("%.0f", $org2length/$height);
my $y_pix_per_bp = 1/$y_bp_per_pix;
#print STDERR join ("\n", $width."x". $height, $org1length."x".$org2length, $x_bp_per_pix, $y_bp_per_pix, $x_pix_per_bp, $y_pix_per_bp);

my $gd = new GD::Image($width, $height);
my $white = $gd->colorResolve(255,255,255);
$gd->fill(1,1,$white);

draw_chromosome_grid(gd=>$gd, org1=>$org1info, org2=>$org2info, x_pix_per_bp=>$x_pix_per_bp, y_pix_per_bp=>$y_pix_per_bp, link=>$link);
draw_dots(gd=>$gd, file=>$dagfile, org1=>$org1info, org2=>$org2info, x_pix_per_bp=>$x_pix_per_bp, y_pix_per_bp=>$y_pix_per_bp);
my $add = 1 if $orgid1 eq $orgid2;
draw_dots(gd=>$gd, file=>$alignfile, org1=>$org1info, org2=>$org2info, x_pix_per_bp=>$x_pix_per_bp, y_pix_per_bp=>$y_pix_per_bp, color=>$gd->colorResolve(0,150,0), size=>2, add_inverse=>$add);

open (OUT, ">".$basename.".png") || die "$!";
binmode OUT;
print OUT $gd->png;
close OUT;
sub draw_dots
  {
    my %opts = @_;
    my $file = $opts{file};
    my $gd = $opts{gd};
    my $org1 = $opts{org1};
    my $org2 = $opts{org2};
    my $x_pix_per_bp=$opts{x_pix_per_bp};
    my $y_pix_per_bp=$opts{y_pix_per_bp};
    my $size = $opts{size} || 1;
    my $color = $opts{color};
    my $add_inverse = $opts{add_inverse};
    $color = $gd->colorResolve(150,150,150) unless $color;
    open (IN, $file) || die "$!";
    while (<IN>)
      {
	chomp;
	next if /^#/;
	next unless $_;
	my @line = split /\t/;
	my ($a, $chr1) = split /_/,$line[0],2;
	my ($b, $chr2) = split /_/,$line[4],2;
	next if $CHR1 && $chr1 ne $CHR1;
	next if $CHR2 && $chr2 ne $CHR2;
	my $org1 = $opts{org1};
	my $org2 = $opts{org2};
	my ($xmin) = sort ($line[2], $line[3]);
	my $x = sprintf("%.0f",($org1->{$chr1}{start}+$xmin+abs($line[3]-$line[2])/2)*$x_pix_per_bp);
	my ($ymin) = sort ($line[6], $line[7]);
	my $y = sprintf("%.0f",($org2->{$chr2}{start}+$ymin+abs($line[7]-$line[6])/2)*$y_pix_per_bp);
#	print STDERR $x,"x", $y,"\n";
	$gd->arc($x, $gd->height-$y, $size, $size, 0, 360, $color);
	$gd->arc($y, $gd->height-$x, $size, $size, 0, 360, $color) if $add_inverse && $x ne $y;
      }
    close IN;
  }

sub draw_chromosome_grid
  {
    my %opts = @_;
    my $gd = $opts{gd};
    my $org1 = $opts{org1};
    my $org2 = $opts{org2};
    my $x_pix_per_bp=$opts{x_pix_per_bp};
    my $y_pix_per_bp=$opts{y_pix_per_bp};
    my $link = $opts{link};
    my $black = $gd->colorResolve(0,0,0);
    $gd->line(0,0, $gd->width, 0, $black);
    $gd->line(0, $gd->height-1, $gd->width, $gd->height-1, $black);
    $gd->line(0,0, 0, $gd->height, $black);
    $gd->line($gd->width-1,0, $gd->width-1, $gd->height, $black);
    #make x-axis chromosome deliniators
    my %data;
    my $pchr;
    my $pv;
    foreach my $chr (sort {$org1->{$a}{start}<=>$org1->{$b}{start} } keys %$org1)
      {
	my $x = sprintf("%.0f",$org1->{$chr}{start}*$x_pix_per_bp);
	$gd->line($x, 0, $x, $gd->height, $black);
	$gd->string(gdSmallFont, $x+2, $gd->height-15, $chr, $black);
	$data{x}{$pchr}=[$pv,$x] if $pchr;
	$pchr=$chr;
	$pv = $x+1;
      }
    $data{x}{$pchr}=[$pv,$gd->width];
    $pchr = undef;
    $pv=undef;
    foreach my $chr (sort {$org2->{$a}{start}<=>$org2->{$b}{start} } keys %$org2)
      {
	my $y = $gd->height-sprintf("%.0f",$org2->{$chr}{start}*$y_pix_per_bp);
	
	$gd->line(0, $y, $gd->width, $y, $black);
	$gd->string(gdSmallFont, 2, $y-15, $chr, $black);
	$data{y}{$pchr}=[$y, $pv] if $pchr;
	$pchr = $chr;
	$pv =$y+1;
      }
    $data{y}{$pchr}=[0,$pv-1];
    if ($link)
      {
	open (OUT, ">".$basename.".html") || die "$!";
	print OUT "<html><head></head><body>\n";
	my ($img) = $basename =~ /([^\/]*$)/;
	print OUT qq{
<IMG SRC="$img.png" usemap="#points" border="0">
<map name="points">
};
	foreach my $xchr (keys %{$data{x}})
	  {
	    my ($x1, $x2) = @{$data{x}{$xchr}};
	    foreach my $ychr (keys %{$data{y}})
	      {
		my $tmp = $link;
		$tmp =~ s/XCHR/$xchr/;
		$tmp =~ s/YCHR/$ychr/;
		my ($y1, $y2) = @{$data{y}{$ychr}};
		print OUT qq{
<area href='$tmp' alt="$xchr vs $ychr" title="$xchr vs $ychr" shap=rect coords="$x1, $y1, $x2, $y2">
};
	      }
	  }
	
	print OUT qq{
</map>
</body></html>
};
	close OUT;
      }
  }

sub get_org_info
  {
    my %opts = @_;
    my $oid = $opts{oid};
    my $chr = $opts{chr};
    my $minsize = $opts{minsize};
    my $org = $coge->resultset('Organism')->find($oid);
    unless ($org)
      {
	warn "No organism found with dbid $oid\n";
	return;
      }
    my %data;
    foreach my $ds ($org->current_datasets)
      {
	foreach my $chrtmp ($ds->chromosomes)
	  {
	    next if $chr && $chr ne $chrtmp;
	    my $last = $ds->last_chromosome_position($chrtmp);
	    next if $minsize && $minsize > $last;
	    if ($data{$chrtmp})
	      {
		warn "Duplicate chromosome: $chrtmp\n";
	      }
	    $data{$chrtmp}{length}=$last;
	  }
      }
    my $pos = 1;
    foreach my $item (sort {$data{$b}{length} <=> $data{$a}{length} } keys %data )
      {
	$data{$item}{start} = $pos;
	$pos += $data{$item}{length};
      }
		      
    return \%data;
  }

sub usage
  {
    print qq{
Welcome to $0

dagfile      | d       path to dag file containing all the hits

alignfile    | a       path to .aligncoords file generated by 
                       dagchainer containing just the diags

width        | w       width of image in pixels (1024)

min_chr_size | mcs     minimim size of chromosome to be drawn (0)

orgid1       | oid1    database id of organism on x-axis

orgid2       | oid2    database id of organism on y-axis

chr1         | c1      only show data for this chromosome on the x axis

chr2         | c2      only show data for this chromosome on the y axis

basename     | b       base path and name for output files (test)

link         | l       link to be used in image map

link_type    | lt      are image map links for chromosome blocks or points:
                       1  ::   blocks  (default if link is specified)  (Use "XCHR","YCHR" which will get the appropriate chr substituted in) 
                       2  ::   points

help         | h       print this message

};
    exit;
  }
