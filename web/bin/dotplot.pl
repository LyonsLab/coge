#!/usr/bin/perl -w

use strict;
use GD;
use Getopt::Long;

use CoGeX;
use CoGe::Accessory::Web qw(url_for);
use CoGe::Accessory::SynMap_report;
use CoGe::Accessory::Utils qw( commify );
use CoGe::Core::Chromosomes;
use DBI;
use Data::Dumper;
use DBI;
use POSIX;
use Sort::Versions;
use HTML::Template;

use vars
  qw($P $dagfile $alignfile $width $link $min_chr_size $dsgid1 $dsgid2 $help $coge $graphics_context $CHR1 $CHR2 $basename $link_type $flip $grid $ks_db $ks_type $log $MAX $MIN $assemble $axis_metric $color_type $box_diags $fid1 $fid2 $selfself $labels $color_scheme $chr_sort_order $font $GZIP $GUNZIP $conffile $skip_random $force_box $chr_order $dotsize $linesize $STANDALONE);

GetOptions(
    "dagfile|d=s"            => \$dagfile,        #all dots
    "alignfile|a=s"          => \$alignfile,      #syntenic dots
    "width|w=i"              => \$width,
    "link|l=s"               => \$link,
    "link_type|lt=i"         => \$link_type,
    "min_chr_size|mcs=i"     => \$min_chr_size,
    "dsgid1|dsg1=i"          => \$dsgid1,
    "dsgid2|dsg2=i"          => \$dsgid2,
    "help|h"                 => \$help,
    "chr1|c1=s"              => \$CHR1,
    "chr2|c2=s"              => \$CHR2,
    "basename|b=s"           => \$basename,
    "flip|f=i"               => \$flip,
    "grid|g=i"               => \$grid,
    "ksdb|ks_db=s"           => \$ks_db,
    "ks_type|kst=s"          => \$ks_type,
    "color_diags_type|cdt=s" => \$color_type,     #how to color diagonals
      #NOTE this is trumped by ksdb and the diag dots will be colored according to ks_type
    "max=s"            => \$MAX,
    "min=s"            => \$MIN,
    "log=s"            => \$log,
    "assemble=s"       => \$assemble,      #syntenic path assembly option
    "axis_metrix|am=s" => \$axis_metric,
    "box_diags|bd=i"   => \$box_diags,
    "fid1|f1=i"        => \$fid1,
    "fid2|f2=i"        => \$fid2,
    "selfself"         => \$selfself,      #draw diag for self self comparison
    "labels=s" => \$labels,    #do we draw labels for chromosomes, on my default
    "color_scheme=s"       => \$color_scheme,
    "font=s"               => \$font,
    "chr_sort_order|cso=s" => \$chr_sort_order
    ,    #display sort order for chromosomes, "Name|N" || "Size|S";
    "config_file|cf=s" => \$conffile,
    "skip_random|sr=i" => \$skip_random
    ,    #flag or skipping chromosome containing the word 'random' in their name
    "force_box|fb" => \$force_box
    , #flag to make dotplot dimensions a box instead of relative on genomic content
    "chr_order|co=s" => \$chr_order
    , #string of ":" delimted chromosome name order to display for the genome with fewer chromosomes
    "dotsize|ds=i"  => \$dotsize,     #size of dots in dotplot
    "linesize|ls=i" => \$linesize,    #size of lines in dotplot
    "standalone|sa=i"  => \$STANDALONE
);
$selfself       = 1      unless defined $selfself;
$labels         = 1      unless defined $labels;
$chr_sort_order = "size" unless defined $chr_sort_order;

$P      = CoGe::Accessory::Web::get_defaults($conffile);
$font   = $P->{FONT} unless $font && -r $font;
$GZIP   = $P->{GZIP};
$GUNZIP = $P->{GUNZIP};

usage() if $help;
unless ( ( defined $dagfile && -r $dagfile )
    || -r $alignfile
    || -r "$alignfile.gz" )
{
  print qq{
Need to define an input file for dots.  Either the syntenic pairs file or the alignment file.
};
    usage();
}

if ( defined $dagfile and !( -r $dagfile || -r $dagfile . ".gz" ) ) {
    warn "dagfile specified but not present or readable: $!";
}

$dagfile = CoGe::Accessory::Web::gunzip($dagfile)
  if $dagfile;    # && $dagfile =~ /\.gz$/;
$alignfile = CoGe::Accessory::Web::gunzip($alignfile)
  if $alignfile;    # && $alignfile =~ /\.gz$/;

if ( $alignfile && -r $alignfile && $alignfile =~ /\.gz$/ ) {
    print "Problem decompressing $alignfile.\n";
    exit;
}

if ( $dagfile && -r $dagfile && $dagfile =~ /\.gz$/ ) {
    print "Problem decompressing $dagfile.\n";
    exit;
}

#set a default for this, make sure it is uppercase
$ks_type = "kS" unless $ks_type;
$ks_type = uc($ks_type) if $ks_type;

my $ks_hist = $P->{KS_HISTOGRAM};

$basename = "test" unless $basename;
$width    = 1024   unless $width;
$coge = CoGeX->dbconnect($P);

my $synmap_report = new CoGe::Accessory::SynMap_report;

my ($dsg1) = $coge->resultset('Genome')->find($dsgid1);
my ($dsg2) = $coge->resultset('Genome')->find($dsgid2);
unless ($dsg1) {
    warn "No genome found with id $dsgid1\n";
    return;
}
unless ($dsg2) {
    warn "No genome found with id $dsgid2\n";
    return;
}

#get display order of chromosomes, get genome information
my ( $org1_order, $org1info ) = get_dsg_order(
    dsg            => $dsg1,
    chr            => $CHR1,
    minsize        => $min_chr_size,
    chr_sort_order => $chr_sort_order,
    skip_random    => $skip_random
);
my ( $org2_order, $org2info ) = get_dsg_order(
    dsg            => $dsg2,
    chr            => $CHR2,
    minsize        => $min_chr_size,
    chr_sort_order => $chr_sort_order,
    skip_random    => $skip_random
);

add_user_order(
    order1     => $org1_order,
    order2     => $org2_order,
    user_order => $chr_order
) if $chr_order;

if ( $axis_metric && $axis_metric =~ /gene/ )    #add in gene information
{
    get_gene_info( dsgid => $dsgid1, info => $org1info );
    get_gene_info( dsgid => $dsgid2, info => $org2info );
}

# Reorder whichever genome has more chromosomes/contigs
my $spa_info_file = $basename . ".spa_info.txt";
if ($assemble) {
    my $skip_non_aligned_contigs = $assemble && $assemble =~ /2/ ? 1 : 0;
    print STDERR 'alignfile=', $alignfile, "\n";
    my ( $org1_association, $org2_association );
    ( $org1_association, $org2_association ) =
        $synmap_report->parse_syn_blocks( file => $alignfile ) if $assemble;
    #print STDERR Dumper $org1_association, $org2_association;
    
    my $output;
    if (   (@$org1_association > @$org2_association && $assemble > 0)
        || (@$org1_association < @$org2_association && $assemble < 0) )
    {
        $output = reord(
            reorder     => $org1_order,
            order       => $org2_order,
            assoc       => $org2_association,
            skip        => $skip_non_aligned_contigs,
            info        => $org1_association,
            skip_random => $skip_random
        );
    }
    else {
        $output = reord(
            reorder     => $org2_order,
            order       => $org1_order,
            assoc       => $org1_association,
            skip        => $skip_non_aligned_contigs,
            info        => $org2_association,
            skip_random => $skip_random
        );
    }

    open OUT, ">$spa_info_file";
    print OUT $output;
    close OUT;
    
    add_rev_info( info => $org1info, assoc => $org1_association );
    add_rev_info( info => $org2info, assoc => $org2_association );
}

calc_abs_start_pos( order => $org1_order, info => $org1info );
calc_abs_start_pos( order => $org2_order, info => $org2info );

#my $org1info = get_dsg_info(dsgid=>$dsgid1, chr=>$CHR1, minsize=>$min_chr_size, order=>$org1_order, metric=>$axis_metric, chr_sort_order=>$chr_sort_order, skip_random=>$skip_random);
my $org1length = 0;
map { $org1length += $_->{length} } values %$org1info;
#my $org2info = get_dsg_info(dsgid=>$dsgid2, chr=>$CHR2, minsize=>$min_chr_size, order=>$org2_order, metric=>$axis_metric, chr_sort_order=>$chr_sort_order, skip_random=>$skip_random);
my $org2length = 0;
map { $org2length += $_->{length} } values %$org2info;

unless ( $org1length && $org2length ) {
    print STDERR "Error: one or both of the genomes has no effective length: Org1: $org1length Org2: $org2length\n";
    if ( $axis_metric && $axis_metric =~ /gene/ ) {
        print STDERR qq{Possible source of error is using genes as axis metric when one or both genomes has no CDSs.};
    }
    exit;
}

( $org1info, $org1length, $dsgid1, $org2info, $org2length, $dsgid2 ) =
  ( $org2info, $org2length, $dsgid2, $org1info, $org1length, $dsgid1 )
  if $flip;
( $CHR1, $CHR2 ) = ( $CHR2, $CHR1 ) if $flip && ( $CHR1 || $CHR2 );

my $height = sprintf( "%.0f", $width * $org2length / $org1length );
$height = $width
  if ( $height > 20 * $width ) || ( $height < $width / 20 ) || $force_box;
my $x_bp_per_pix = $org1length / $width;   #sprintf("%.0f", $org1length/$width);
#$x_bp_per_pix = 1 if $x_bp_per_pix < 1;
my $x_pix_per_bp = 1 / $x_bp_per_pix;
my $y_bp_per_pix = $org2length / $height; #sprintf("%.0f", $org2length/$height);

#$y_bp_per_pix = 1 if $y_bp_per_pix < 1;
my $y_pix_per_bp = 1 / $y_bp_per_pix;

#Generate new graphics context and fill the background
my $graphics_context = new GD::Image( $width, $height );
my $white       = $graphics_context->colorResolve( 255, 255, 255 );
my $black       = $graphics_context->colorResolve( 0,   0,   0 );
my $grey        = $graphics_context->colorResolve( 200, 200, 200 );
my $red         = $graphics_context->colorResolve( 255, 0,   0 );
my $green       = $graphics_context->colorResolve( 0,   100, 0 );
my $alert_color = $graphics_context->colorResolve( 255, 0,   0 );

$graphics_context->fill( 1, 1, $white );

if ($ks_db) {
    my $cmd = $ks_hist;
    $cmd .= " -db $ks_db";
    $cmd .= " -ks_type $ks_type";
    $cmd .= " -log" if $log;
    $cmd .= " -pf $alignfile";
    $cmd .= " -chr1 $CHR1" if defined $CHR1;
    $cmd .= " -chr2 $CHR2" if defined $CHR2;
    $cmd .= " -min $MIN" if defined $MIN;
    $cmd .= " -max $MAX" if defined $MAX;
    $cmd .= " -color_scheme $color_scheme" if defined $color_scheme;
    $cmd .= " -o $basename";
    #    $cmd .= ".log" if $log;
    #    $cmd .= ".$MIN" if defined $MIN;
    #    $cmd .= ".$MAX" if defined $MAX;
    $cmd .= ".hist.png";
    #print STDERR "HIST: ",$cmd,"\n";
    `$cmd &`;
}

#$basename .= ".$MIN" if defined $MIN;
#$basename .= ".$MAX" if defined $MAX;

my $pairs = get_pairs( file => $alignfile, chr1 => $CHR1, chr2 => $CHR2 )
  if $alignfile && $ks_db && -r $alignfile && -r $ks_db;

#Magic happens here.
#Link type seems to indicate the type of tile; i.e. a 'master' (a large, all chromosome) or a blow up of a two chromosome intersection
#draw_chromosome_grid draws either the black chomosome lines, or the light green tile lines, so its always called in addition to the draw_dots function.
my $coords = draw_chromosome_grid(
    gd           => $graphics_context,
    org1         => $org1info,
    org2         => $org2info,
    x_pix_per_bp => $x_pix_per_bp,
    y_pix_per_bp => $y_pix_per_bp,
    link         => $link,
    link_type    => $link_type,
    flip         => $flip,
    grid         => $grid,
    color        => $grey,
    linesize     => $linesize
);

my $x_labels_gd = draw_x_labels( coords => $coords, axis => 'x' ) if $labels;
my $y_labels_gd = draw_y_labels( coords => $coords, axis => 'y' ) if $labels;

#get syntenic gene pairs for ks_data (if needed)
my $ksdata = get_ksdata(
    ks_db   => $ks_db,
    ks_type => $ks_type,
    chr1    => $CHR1,
    chr2    => $CHR2,
    pairs   => $pairs
) if $ks_db && -r $ks_db;

my $size = 1;
$size = int($x_pix_per_bp) if $x_pix_per_bp > 1;
$size = int($y_pix_per_bp) if $y_pix_per_bp > 1;
$size = 4                  if $x_pix_per_bp > 4;
$size = 4                  if $y_pix_per_bp > 4;
$size = $dotsize           if $dotsize;

#draw dots for all matches
draw_dots(
    gd           => $graphics_context,
    file         => $dagfile,
    org1         => $org1info,
    org2         => $org2info,
    x_pix_per_bp => $x_pix_per_bp,
    y_pix_per_bp => $y_pix_per_bp,
    link_type    => $link_type,
    dsgid1       => $dsgid1,
    dsgid2       => $dsgid2,
    flip         => $flip,
    metric       => $axis_metric,
    fid1         => $fid1,
    fid2         => $fid2,
    size         => $size
) if defined $dagfile && -r $dagfile;

#color_scheme
my @colors;
if ( $color_type && $color_type eq "inv" ) {
    push @colors, $graphics_context->colorResolve( 0, 150, 0 );   #forward green
    push @colors, $graphics_context->colorResolve( 0, 0,   150 ); #reverse blue
}
elsif ( $color_type && $color_type eq "diag" ) {
    push @colors, $graphics_context->colorResolve( 150, 0,   0 );
    push @colors, $graphics_context->colorResolve( 0,   150, 0 );
    push @colors, $graphics_context->colorResolve( 0,   0,   150 );
    push @colors, $graphics_context->colorResolve( 150, 150, 0 );
    push @colors, $graphics_context->colorResolve( 150, 0,   150 );
    push @colors, $graphics_context->colorResolve( 0,   150, 150 );
}
elsif ( $color_type && $color_type eq "chr" ) {
    my $color_num =
        scalar keys %$org1info > scalar keys %$org2info
      ? scalar keys %$org2info
      : scalar keys %$org1info;
    push @colors, $graphics_context->colorResolve( 200, 0,   0 );
    push @colors, $graphics_context->colorResolve( 200, 200, 0 );
    push @colors, $graphics_context->colorResolve( 0,   200, 0 );
    push @colors, $graphics_context->colorResolve( 0,   0,   200 );
    push @colors, $graphics_context->colorResolve( 0,   200, 200 );
    push @colors, $graphics_context->colorResolve( 200, 0,   200 );
    push @colors, $graphics_context->colorResolve( 100, 0,   0 );
    push @colors, $graphics_context->colorResolve( 100, 100, 0 );
    push @colors, $graphics_context->colorResolve( 0,   100, 0 );
    push @colors, $graphics_context->colorResolve( 0,   0,   100 );
    push @colors, $graphics_context->colorResolve( 0,   100, 100 );
    push @colors, $graphics_context->colorResolve( 100, 0,   100 );
    push @colors, $graphics_context->colorResolve( 200, 100, 0 );
    push @colors, $graphics_context->colorResolve( 0,   200, 100 );
    push @colors, $graphics_context->colorResolve( 200, 0,   100 );
    push @colors, $graphics_context->colorResolve( 100, 200, 0 );
    push @colors, $graphics_context->colorResolve( 0,   100, 200 );
    push @colors, $graphics_context->colorResolve( 100, 0,   200 );
    my @tmp_colors;
    my $index = 0;

    for ( my $i = 0 ; $i < $color_num ; $i++ ) {
        push @tmp_colors, $colors[$index];
        $index++;
        $index = 0 if $index > ( scalar @colors - 1 );
    }
    @colors = @tmp_colors;
}
else {
    push @colors, $graphics_context->colorResolve( 0, 150, 0 );
}

#add syntenic gene pairs
my $add;    # = 1 if $dsgid1 eq $dsgid2;

$size = 2;
$size = int($x_pix_per_bp) + 2 if $x_pix_per_bp > 1;
$size = int($y_pix_per_bp) + 2 if $y_pix_per_bp > 1;
$size = 5 if $x_pix_per_bp > 5;
$size = 5 if $y_pix_per_bp > 5;
$size = $dotsize if $dotsize;

#print STDERR $alignfile,"\n";
my $box_coords = draw_dots(
    gd           => $graphics_context,
    file         => $alignfile,
    org1         => $org1info,
    org2         => $org2info,
    x_pix_per_bp => $x_pix_per_bp,
    y_pix_per_bp => $y_pix_per_bp,
    size         => $size,
    add_inverse  => $add,
    flip         => $flip,
    ksdata       => $ksdata,
    ks_type      => $ks_type,
    log          => $log,
    metric       => $axis_metric,
    colors       => \@colors,
    color_type   => $color_type,
    color_scheme => $color_scheme,
    fid1         => $fid1,
    fid2         => $fid2
);

draw_boxes( gd => $graphics_context, boxes => $box_coords )
  if $box_diags && $box_coords && @$box_coords;

#draw self-self line?
my $draw_selfself = 0;

if ( $selfself && ( $dsgid1 == $dsgid2 ) ) {
    if ( $CHR1 && $CHR2 ) {
        $draw_selfself = 1 if $CHR1 eq $CHR2;
    }
    elsif ( !$assemble ) {
        $draw_selfself = 1;
    }
}
if ($draw_selfself)
  {
    $graphics_context->line( 0, $height, $width, 0, $green );
    if ($linesize)
      {
	my $size = ceil($linesize/2);
	for (my $i=1; $i< $size; $i++)
	  {
	        $graphics_context->line( 0+$i, $height+$i, $width+$i, 0+$i, $green );
	        $graphics_context->line( 0-$i, $height-$i, $width-$i, 0-$i, $green );
	  }
      }
  }

#Write out graphics context - the generated dot plot - to a .png file
open( OUT, ">" . $basename . ".png" ) || die "$!";
binmode OUT;
print OUT $graphics_context->png;
close OUT;

if ($labels) {
    #Write out graphics context - the generated x axis - to a .png file
    open( OUT, ">" . $basename . ".x.png" ) || die "$!";
    binmode OUT;
    print OUT $x_labels_gd->png;
    close OUT;

    #Write out graphics context - the generated y axis - to a .png file
    open( OUT, ">" . $basename . ".y.png" ) || die "$!";
    binmode OUT;
    print OUT $y_labels_gd->png;
    close OUT;
}

#CoGe::Accessory::Web::gzip($dagfile) if $dagfile && -r $dagfile;
#CoGe::Accessory::Web::gzip($alignfile) if $alignfile && -r $alignfile;
#generate_historgram of ks values if necessary

#This function appears to parse dagchainer output, generated in SynMap.pl, and draw the results to the GD graphics context.
sub draw_dots {
    my %opts             = @_;
    my $file             = $opts{file};
    my $graphics_context = $opts{gd};
    my $org1             = $opts{org1};
    my $org2             = $opts{org2};
    my $x_pix_per_bp     = $opts{x_pix_per_bp};
    my $y_pix_per_bp     = $opts{y_pix_per_bp};
    my $size             = $opts{size} || 1;
    my $colors           = $opts{colors};
    my $add_inverse      = $opts{add_inverse};
    my $link_type        = $opts{link_type} || 0;
    my $dsgid1           = $opts{dsgid1};
    my $dsgid2           = $opts{dsgid2};
    my $flip             = $opts{flip};
    my $ksdata           = $opts{ksdata};
    my $log        = $opts{log};          #log normalize ksdata for display?
    my $metric     = $opts{metric};
    my $color_type = $opts{color_type};
    my $ks_type    = $opts{ks_type};
    my $fid1       = $opts{fid1};
    my $fid2       = $opts{fid2};
    my $order1     = $opts{order1};       #chromosome display order for genome 1
    my $order2     = $opts{order2};       #chromosome display order for genome 2
    my $color_scheme =
      $opts{color_scheme}; #color pattern for Ks/Kn caluclations
                           #easier lookup, can scale to more pairs in the future
    my %fids;
    $fids{$fid1} = 1 if $fid1;
    $fids{$fid2} = 1 if $fid2;
    my $has_ksdata = keys %$ksdata ? 1 : 0;

    #min and max will be log normalized if log flag is set
    my ( $max, $min ) =
      get_range( data => $ksdata, min => $MIN, max => $MAX, log => $log )
      if $has_ksdata;
    my $range = $max - $min if $has_ksdata;
    my $default_color = $graphics_context->colorResolve( 150, 150, 150 );
    $colors = [$default_color] unless ref($colors) =~ /array/i && @$colors;
    open( IN, $file )
      || die "Can't open $file: $!";    #this is where the problem lies!
    my $use_color;
    $use_color = $colors->[0] unless $color_type;
    my $color_index = 0;

    my @feats;
    my %points;
    my @points;
    my @boxes;
    my ( $min_x, $min_y, $max_x, $max_y );
    my $count = 0;
    while (<IN>) {
        my $tuse_color = $use_color;
        my $tsize      = $size
          ; #might want to dynmaically change the size of the dot.  Reset to default after each line
        chomp;
        if ( /^#/ && $color_type && $color_type eq "diag" ) {
            $tuse_color = $colors->[$color_index];
            $color_index++;
            $color_index = 0 if $color_index >= @$colors;
            $use_color = $tuse_color;
        }
        if (/^#/) {
            push @boxes, [ $min_x - 1, $min_y - 1, $max_x + 1, $max_y + 1 ]
              if defined $min_x
                  && defined $min_y
                  && defined $max_x
                  && defined $max_y;
            $min_x = undef;
            $min_y = undef;
            $max_x = undef;
            $max_y = undef;
        }

        next if /^#/;
        next unless $_;

        my @line = split /\t/;
        my $val;
        my @item1 = split(/\|\|/, $line[1]);
        my @item2 = split(/\|\|/, $line[5]);
        my $fid1  = $item1[6];
        my $fid2  = $item2[6];
        if ( $color_type && $color_type eq "inv" && $item1[4] && $item2[4] ) {
            $tuse_color = $item1[4] eq $item2[4] ? $colors->[0] : $colors->[1];
        }
        if ( $fid1 && $fid2 && $fids{$fid1} && $fids{$fid2} ) {
            $tsize      = 10;
            $tuse_color = $alert_color;
        }

        if ($has_ksdata) {
            $val = $ksdata->{$fid1}{$fid2};
            if ( defined $val && $val =~ /\d/ ) {
                my $orig_val = $val;
                $val = $min if $log && $val == 0;
                if ($log) {
                    if ( $val <= 0 ) {
                        $val = 0;    #set to minimum color value
                    }
                    else {
                        $val = ( log10($val) - $min ) / $range;
                    }
                }
                else {
                    $val = ( $val - $min ) / $range;
                }
                $val = sprintf( "%.4f", $val );
                $tuse_color =
                  get_color( val => $val, color_scheme => $color_scheme )
                  ;                  #val is 0<=x<=1
                $tuse_color = $graphics_context->colorResolve(@$tuse_color);
            }
            else {
            #		print STDERR "Skipping due to no ks data for dot: $fid1 $fid2\n";
            #don't have ks data -- skip drawing this dot!
                next;
                #		print Dumper $ksdata->{$item1[6]}{$item2[6]};
                #		$tuse_color = $graphics_context->colorResolve(0,0,0);
            }
        }
        if ($flip) {
            my @tmp = @line[ 0 .. 3 ];
            @line[ 0 .. 3 ] = @line[ 4 .. 7 ];
            @line[ 4 .. 7 ] = @tmp;
            ( $fid1, $fid2 ) = ( $fid2, $fid1 );
        }
        my ( $a, $chr1 ) = split /_/, $line[0], 2;
        my ( $b, $chr2 ) = split /_/, $line[4], 2;
        my $special = 0
          ; #stupid variable name so that if we are viewing a single chr to single chr comparison within the same organism, this will make collinear matches appear on the inverse section

        if ( $CHR1 && $CHR2 ) {
            if ( $add_inverse && $CHR1 eq $chr2 && $CHR2 eq $chr1 ) {
                $special = 1;
                ( $chr1, $chr2 ) = ( $chr2, $chr1 );
            }
            else {
                next if $chr1 ne $CHR1;
                next if $chr2 ne $CHR2;
            }
        }
        next
          unless $org1->{$chr1}
              && $org2->{ $chr2
              }; #sometimes there will be data that is skipped, e.g. where chromosome="random";

        my ( $xmin, $ymin );
        my ( $midx, $midy );

        if ( $metric && $metric =~ /gene/i ) {
            next unless $fid1          && $fid2;
            next unless $org1->{$chr1} && $org1->{$chr1}{gene};
            next unless $org2->{$chr2} && $org2->{$chr2}{gene};
            $xmin = $org1->{$chr1}{gene}{$fid1};
            $ymin = $org2->{$chr2}{gene}{$fid2};
            next unless $xmin && $ymin;
            $midx =
              $org1->{$chr1}{rev}
              ? sprintf( "%.0f",
                $org1->{$chr1}{start} + $org1->{$chr1}{length} - ($xmin) )
              : sprintf( "%.0f", $org1->{$chr1}{start} + $xmin );
            $midy =
              $org2->{$chr2}{rev}
              ? sprintf( "%.0f",
                $org2->{$chr2}{start} + $org2->{$chr2}{length} - ($ymin) )
              : sprintf( "%.0f", $org2->{$chr2}{start} + $ymin );

        }
        else    #using nucleotides
        {
            ($xmin) = sort ( $line[2], $line[3] );
            ($ymin) = sort ( $line[6], $line[7] );
            $midx =
              $org1->{$chr1}{rev}
              ? sprintf( "%.0f",
                $org1->{$chr1}{start} +
                  $org1->{$chr1}{length} -
                  ( $xmin + abs( $line[3] - $line[2] ) / 2 ) )
              : sprintf( "%.0f",
                $org1->{$chr1}{start} +
                  $xmin +
                  abs( $line[3] - $line[2] ) / 2 );
            $midy =
              $org2->{$chr2}{rev}
              ? sprintf( "%.0f",
                $org2->{$chr2}{start} +
                  $org2->{$chr2}{length} -
                  ( $ymin + abs( $line[7] - $line[6] ) / 2 ) )
              : sprintf( "%.0f",
                $org2->{$chr2}{start} +
                  $ymin +
                  abs( $line[7] - $line[6] ) / 2 );

        }

        my $x = sprintf( "%.0f", $midx * $x_pix_per_bp );
        my $y = sprintf( "%.0f", $midy * $y_pix_per_bp );

        ( $x, $y ) = ( $y, $x ) if $special;
        $val = 0 unless $val;    #give it some value for later sorting
        $x = $width - ceil( $tsize / 2 ) if $x >= $width;
        my $y_real = $graphics_context->height - $y;
        push @points,
          [ $x, $y_real, $tsize, $tsize, 0, 360, $tuse_color, $val ];
        $min_x = $x unless $min_x;
        $min_x = $x if $x < $min_x;
        $min_y = $y_real unless $min_y;
        $min_y = $y_real if $y_real < $min_y;
        $max_x = $x unless $max_x;
        $max_x = $x if $x > $max_x;
        $max_y = $y_real unless $max_y;
        $max_y = $y_real if $y_real > $max_y;

        $y_real = $graphics_context->height - $x;
        $tuse_color = $colors->[0] unless $tuse_color;    #default val just in case
        push @points, [ $y, $y_real, $tsize, $tsize, 0, 360, $tuse_color, $val ]
          if ( $add_inverse && !$CHR1 && $x ne $y );
        push @points, [ $y, $y_real, $tsize, $tsize, 0, 360, $tuse_color, $val ]
          if ( $add_inverse && $chr1 eq $chr2 && $x ne $y );

        my %URL_PARAMS = (
            drup1 => 50000,
            drdown1 => 50000,
            drup2 => 50000,
            drdown2 => 50000
        );

        if ( $link_type == 1 ) {
#working here.  Need to build a GEvo link using datasets/chr/position if dealing with genomic data.
            if ($fid1) {
                $URL_PARAMS{fid1} = $fid1;
            }
            else {
                #just added the coge->ds object to $org
                $item1[6] = "Chr: "
                  . $item1[0] . " "
                  . commify( $item1[1] ) . " - "
                  . commify( $item1[2] );
                my $chr = $item1[0];
                @URL_PARAMS{qw(dsgid1 x1 chr1)} = ($dsgid1, $midx, $chr);
                #$link .= qq{;dsgid1=$dsgid1;x1=$midx;chr1=$chr};
            }
            if ($fid2) {
                #$link .= qq{;fid2=$fid2};
                $URL_PARAMS{fid2} = $fid2;
            }
            else {
                #just added the coge->ds object to $org
                $item2[6] = "Chr: "
                  . $item2[0] . " "
                  . commify( $item2[1] ) . " - "
                  . commify( $item2[2] );
                my $chr = $item2[0];
                @URL_PARAMS{qw(dsgid2 x2 chr2)} = ($dsgid2, $midy, $chr);
            }
            unless ( $points{$x}{$y} )   #cuts down on the size of the image map
            {
                push @feats,
                  [
                    $x, $graphics_context->height - $y,
                    $item1[6], $item2[6], url_for("GEvo.pl", %URL_PARAMS)
                  ];
                $points{$x}{$y} = 1;
            }
        }
        $count++;
    }
    close IN;
    push @boxes, [ $min_x - 1, $min_y - 1, $max_x + 1, $max_y + 1 ]
      if defined $min_x && defined $min_y && defined $max_x && defined $max_y;
#@data = $ks_type && $ks_type eq "KN" ? sort {$b<=>$a} @data : sort {$a<=>$b} @data;
    if ($has_ksdata) {
        if ( $ks_type && $ks_type eq "KN" ) {
            @points = sort { $a->[-1] <=> $b->[-1] } @points;
        }
        else {
            @points = sort { $b->[-1] <=> $a->[-1] } @points;
        }
    }
    foreach my $point (@points) {
        my $val   = pop @$point;
        my $tsize = $point->[2];
        if ( $tsize > 3 ) {
            $graphics_context->filledArc(@$point);
        }
        else {
            $graphics_context->arc(@$point);
        }
    }
    if ( $link_type == 1 ) {
        #Okay, now generate the HTML document that contains the click map for the image
        open( OUT, ">" . $basename . ".html" ) || die "$basename.html\n$!";
        my ( $org1name, $org2name );
        if ($dsgid1) {
            my $dsg = $coge->resultset('Genome')->find($dsgid1);
            $org1name = $dsg->organism->name . " (v" . $dsg->version . ")";
        }
        if ($dsgid2) {
            my $dsg = $coge->resultset('Genome')->find($dsgid2);
            $org2name = $dsg->organism->name . " (v" . $dsg->version . ")";
        }

        my $org1length = commify( $org1->{$CHR1}->{length} ) . " bp";
        my $org2length = commify( $org2->{$CHR2}->{length} ) . " bp";

        my ($img) = $basename =~ /([^\/]*$)/;

        my $template = HTML::Template->new(
            filename => $P->{TMPLDIR} . '/widgets/Dotplot.tmpl'
        );

        #Loop through the features list and print out the click map info
        my ($index, $clickmap) = (0, "");

        foreach my $item (@feats) {
            my ( $x, $y, $f1, $f2, $link ) = @$item;
#<area shape='circle' coords='$x, $y, 2' href='$link' onMouseOver="get_pair_info(['args__$f1','args__$f2'],['pair_info']);" target='_blank' >
            $clickmap .= qq{location_list[$index]}
            . qq{ = ['circle', [$x, $y, 2],}
            . qq {'$link', 'get_pair_info(["args__$f1","args__$f2"]}
            . qq {,["pair_info"])'];\n};

            $index++;
        }

        #Print out the canvas tag, onLoad call, and other stuff.
        my $y = $y_labels_gd->width . "px";

        #Print out the nametag
        my $pos = int($graphics_context->height + 45) . "px";

        # Include histogram if the file exists
        my $include_histogram = (-e $basename . ".hist.png") ? 1 : 0;

        $template->param(
            url         => url_for(""),
            img         => $img,
            chr1        => $CHR1,
            chr2        => $CHR2,
            org1name    => $org1name,
            org2name    => $org2name,
            length1     => $org1length,
            length2     => $org2length,
            standalone  => $STANDALONE,
            histogram   => $include_histogram,
            click_map   => $clickmap,
        );

        print OUT $template->output;
        close OUT;
    }

    return \@boxes;
}

#This method draws the grid lines indicating cromosome locations, and write out the click map to the HTML file.
#On 'zoom in' tile, this method draws the light green seperator lines.
sub draw_chromosome_grid {
    my %opts             = @_;
    my $graphics_context = $opts{gd};
    my $org1             = $opts{org1};
    my $org2             = $opts{org2};
    my $x_pix_per_bp     = $opts{x_pix_per_bp};
    my $y_pix_per_bp     = $opts{y_pix_per_bp};
    my $link             = $opts{link};
    my $link_type        = $opts{link_type};
    my $grid             = $opts{grid};
    my $color            = $opts{color};
    my $linesize         = $opts{linesize};

    $color = $black unless $color;
    my $height     = $graphics_context->height;
    my $width      = $graphics_context->width;
    my $span_color = $graphics_context->colorResolve( 200, 255, 200 );
    $graphics_context->line( 0, 0, $graphics_context->width, 0, $black );
    $graphics_context->line( 0, $graphics_context->height - 1,
        $graphics_context->width, $graphics_context->height - 1, $black );
    $graphics_context->line( 0, 0, 0, $graphics_context->height, $black );
    $graphics_context->line(
        $graphics_context->width - 1,
        0, $graphics_context->width - 1,
        $graphics_context->height, $black
    );
    #make x-axis chromosome deliniators
    my %data;
    my $pchr;
    my $pv;
    my @x_label_coords;
    my @y_label_coords;
    my $str_color;
    foreach
      my $chr ( sort { $org1->{$a}{start} <=> $org1->{$b}{start} } keys %$org1 )
    {
        my $x = $org1->{$chr}{start} * $x_pix_per_bp;
        next if $pv && $x == $pv - 1;
        if ($grid) {
            my $tmp  = $x;
            my $span = $org1->{$chr}{length} / 10 * $x_pix_per_bp;
            for ( 1 .. 9 ) {
                $tmp += $span;
                $graphics_context->line( $tmp, 0, $tmp,
                    $graphics_context->height, $span_color ); #Draw green lines?
            }
        }
        if ($linesize) {
            $graphics_context->filledRectangle(
                $x - ( $linesize / 2 ),
                0, $x + ( $linesize / 2 ),
                $graphics_context->height, $color
            );                                                #Draw black line
        }
        else {
            $graphics_context->line( $x, 0, $x, $graphics_context->height,
                $color );                                     #Draw black line
        }
        $str_color = $org1->{$chr}{rev} ? $red : $black;
#	if ($labels)
# {#Draws name of chromosome
#   if (-r $font)
#     {
# 	$graphics_context->stringFT($str_color, $font, 12, 0, $x+2, $graphics_context->height-15,$chr);
#     }
#   else
#     {
# 	$graphics_context->string(gdSmallFont, $x+2, $graphics_context->height-15, $chr, $str_color);
#     }
# }
        $data{x}{$pchr} = [ $pv, $x, $str_color ] if defined $pchr;
        $pchr           = $chr;
        $pv             = $x + 1;
    }
    $data{x}{$pchr} = [ $pv, $graphics_context->width, $str_color ];
    $pchr           = undef;
    $pv             = undef;
    foreach
      my $chr ( sort { $org2->{$a}{start} <=> $org2->{$b}{start} } keys %$org2 )
    {
        my $y =
          $graphics_context->height - $org2->{$chr}{start} * $y_pix_per_bp;
        next if $pv && $y == $pv - 1;
        if ($grid) {
            my $tmp  = $y;
            my $span = $org2->{$chr}{length} / 10 * $y_pix_per_bp;
            for ( 1 .. 9 ) {
                $tmp -= $span;
                $graphics_context->line( 0, $tmp, $graphics_context->width,
                    $tmp, $span_color );
            }
        }
        if ($linesize) {
            $graphics_context->filledRectangle( 0, $y - ( $linesize / 2 ),
                $graphics_context->width, $y + ( $linesize / 2 ), $color );
        }
        else {
            $graphics_context->line( 0, $y, $graphics_context->width, $y,
                $color );
        }
        $str_color = $org2->{$chr}{rev} ? $red : $black;
#	$graphics_context->string(gdSmallFont, 2, $y-15, $chr, $str_color) if $labels;
        $data{y}{$pchr} = [ $y, $pv, $str_color ] if defined $pchr;
        $pchr           = $chr;
        $pv             = $y + 1;
    }
    $data{y}{$pchr} = [ 0, $pv - 1, $str_color ];

    if ( $link_type == 2 ) {
        my $template = HTML::Template->new(
            filename => $P->{TMPLDIR} . '/widgets/Dotplot.tmpl'
        );

        #Loop through our list of chromosomes and print out the corrosponding click map tag
        my ($index, $map) = (0, "");

        foreach my $xchr ( sort keys %{ $data{x} } ) {
            my ( $x1, $x2 ) = @{ $data{x}{$xchr} };
            next if abs( $x1 - $x2 ) < 3;
            foreach my $ychr ( sort keys %{ $data{y} } ) {
                my $tmp = $link;
                $tmp =~ s/XCHR/$xchr/;
                $tmp =~ s/YCHR/$ychr/;
                my ( $y1, $y2 ) = @{ $data{y}{$ychr} };
                next if abs( $y1 - $y2 ) < 3;

                $map .= qq{location_list[$index] = ['rect'}
                     .  qq{, [$x1, $y1, $x2, $y2], '$tmp', ''];\n};

                $index++;
            }
        }

        my ($img) = $basename =~ /([^\/]*$)/;

        $template->param(
            url        => url_for(""),
            img        => $img,
            histogram  => 0,
            click_map  => $map,
        );

        open( OUT, ">" . $basename . ".html" ) || die "$!";
        print OUT $template->output;
        close OUT;
    }
    return \%data;
}

sub draw_x_labels {
    my %opts   = @_;
    my $coords = $opts{coords};
    my $axis   = $opts{axis};
    my $gd     = $opts{gd};
    unless ($gd) {
        #calculate size of figure by legends sizes
        my ( $x, $y ) = ( 0, 0 );
        foreach my $chr ( keys %{ $coords->{$axis} } ) {
            my ( $start, $stop, $color ) = @{ $coords->{$axis}{$chr} };
            my $tx = ( $stop - $start ) / 2;
            my $fsize = $tx > 12 ? 12 : $tx;
            $x     = $stop if $stop > $x;
            $fsize = 4     if $fsize < 4;
            my @bounds =
              GD::Image->stringFT( $color, $font, $fsize, 45, $x, 0, $chr );
            my $ty = ( $bounds[1] - $bounds[5] );
            $y = $ty if $ty > $y;
        }
        $gd = new GD::Image( $x, $y + 5 );
        return unless $gd;

        my $white       = $gd->colorResolve( 255, 255, 255 );
        my $black       = $gd->colorResolve( 0,   0,   0 );
        my $red         = $gd->colorResolve( 255, 0,   0 );
        my $green       = $gd->colorResolve( 0,   150, 0 );
        my $alert_color = $gd->colorResolve( 255, 0,   0 );
    }
    foreach my $chr ( keys %{ $coords->{$axis} } ) {
        my ( $start, $stop, $color ) = @{ $coords->{$axis}{$chr} };
        my $x = ( $stop - $start ) / 2;
        my $fsize = $x > 12 ? 12 : $x;
        $fsize = 4 if $fsize < 4;
        $x += $start;
        if ( -r $font ) {
            $gd->stringFT( $color, $font, $fsize, 45, $x, $gd->height - 5,
                $chr );
        }
        else {
            $gd->string( gdSmallFont, $start, $gd->height - 15, $chr, $color );
        }
    }
    return $gd;
}

sub draw_y_labels {
    my %opts   = @_;
    my $coords = $opts{coords};
    my $axis   = $opts{axis};
    my $gd     = $opts{gd};
    unless ($gd) {
        #calculate size of figure by legends sizes
        my ( $x, $y ) = ( 0, 0 );
        foreach my $chr ( keys %{ $coords->{$axis} } ) {
            my ( $start, $stop, $color ) = @{ $coords->{$axis}{$chr} };
            my $ty = ( $stop - $start ) / 2;
            my $fsize = $ty > 12 ? 12 : $ty;
            $y     = $stop if $stop > $y;
            $fsize = 4     if $fsize < 4;
            my @bounds =
              GD::Image->stringFT( $color, $font, $fsize, 45, 0, 0, $chr );
            my $tx = ( $bounds[2] - $bounds[6] );
            $x = $tx if $tx > $x;
        }
        $gd = new GD::Image( $x, $y + 5 );
        return unless $gd;

        my $white       = $gd->colorResolve( 255, 255, 255 );
        my $black       = $gd->colorResolve( 0,   0,   0 );
        my $red         = $gd->colorResolve( 255, 0,   0 );
        my $green       = $gd->colorResolve( 0,   150, 0 );
        my $alert_color = $gd->colorResolve( 255, 0,   0 );
    }
    foreach my $chr ( keys %{ $coords->{$axis} } ) {
        my ( $start, $stop, $color ) = @{ $coords->{$axis}{$chr} };
        my $y = ( $stop - $start ) / 2;
        my $fsize = $y > 12 ? 12 : $y;
        $fsize = 4 if $fsize < 4;
        $y += $start;
        if ( -r $font ) {
            $gd->stringFT( $color, $font, $fsize, 45, 0 + 10, $y, $chr );
        }
        else {
            $gd->string( gdSmallFont, 0, $y, $chr, $color );
        }
    }
    return $gd;
}

sub add_user_order {
    my %opts       = @_;
    my $order1     = $opts{order1};
    my $order2     = $opts{order2};
    my $user_order = $opts{user_order};
    my $to_order =
      @$order1 > @$order2 ? $order2 : $order1;    #get the smaller list
    my %chrs = map { $_, 1 } @$to_order;
    my %seen;
    my @new_order;

    foreach my $item ( split ":", $user_order ) {
        if ( $chrs{$item} ) {
            push @new_order, $item;
        }
        else {
            warn "$item is not in chromosome list!";
        }
        $seen{$item} = 1;
    }
    foreach my $item (@$to_order) {
        next if $seen{$item};
        push @new_order, $item;
        $seen{$item} = 1;
    }
    @$to_order = @new_order;
}

sub reord {
    my %opts        = @_;
    my $order       = $opts{order};
    my $reorder     = $opts{reorder};
    my $association = $opts{assoc};
    my $skip        = $opts{skip};
    my $skip_random = $opts{skip_random};
    my $info        = $opts{info};          #for determining orientation
    print STDERR 'reord: info=', scalar(@$info), ' order=', scalar(@$order), ' reorder=', scalar(@$reorder), ' assoc=', scalar(@$association), ' skip=', $skip, ' skip_random=', $skip_random, "\n";

    #create mapping hash of which contigs are in the reverse orientation
    my %rev_info;
    foreach my $item (@$info) {
        my $rev = $item->{rev} ? -1 : 1;
        my $chr = $item->{chr};
        $rev_info{$chr} = $rev;
    }

    my %mapped_association;
    map { $mapped_association{ $_->{chr} } = $_->{matching_chr} } @$association;
    my @new_order;

    my $output = join( "\t", ( "#CHR1", "CHR2", "ORIENTATION" ) ) . "\n";
    foreach my $chr (@$order) {
        next if check_random($chr) && $skip_random;
        if ( $mapped_association{$chr} ) {
            foreach my $item ( @{ $mapped_association{$chr} } ) {
                next if check_random($item) && $skip_random;
                push @new_order, $item;
            }
        }
        $output .= join( "\n",
            map { join( "\t", $chr, $_, $rev_info{$_} ) }
              @{ $mapped_association{$chr} } )
          . "\n";
    }
    unless ($skip) {
        my %seen = map { $_ => 1 } @new_order;
        foreach my $item (@$reorder) {
            $output .= join( "\t", "unmapped", $item ) . "\n"
              unless $seen{$item};
            push @new_order, $item unless $seen{$item};
        }
    }
    @$reorder = @new_order;
    return $output;
}

sub calc_abs_start_pos {
    my %opts  = @_;
    my $order = $opts{order};
    my $info  = $opts{info};
    my $pos   = 1;
    my %seen;
    foreach my $chr (@$order) {
        next if $seen{$chr};
        $seen{$chr} = 1;
        $info->{$chr}{start} = $pos - 1;
        $pos += $info->{$chr}{length};
    }
    #delete items from info that are not in the ordered list
    foreach my $chr ( keys %$info ) {
        delete $info->{$chr} unless $seen{$chr};
    }
}

sub get_dsg_info {
    my %opts             = @_;
    my $dsgid            = $opts{dsgid};
    my $chr              = $opts{chr};
    my $minsize          = $opts{minsize};
    my $order            = $opts{order};
    my $skip_non_ordered = $opts{skip_non_ordered};
    my $metric           = $opts{metric};
    my $chr_sort_order   = $opts{chr_sort_order};
    my $skip_random =
      $opts{skip_random}; #skip "random" chromosome where sequences are added ad hoc
    my $dsg = $coge->resultset('Genome')->find($dsgid);

    unless ($dsg) {
        warn "No genome found with dbid $dsgid\n";
        return;
    }
    my %rev;              #store chromosomes to be reversed in display
    my @ordered;
    my %data;
    if ($order)           #chromsomes have a prespecified order
    {
        my %seen;
        foreach my $chr (@$order) {
            #	    my $chr = $item->{chr};
            next unless $data{$chr};
            push @ordered, $chr;
            #	    $rev{$chr}=1 if $item->{rev};
            $seen{$chr} = 1;
        }
        my @chr;
        #get any that were not in @$order
        #how to sort chromosomes for diplay?
        if ( $chr_sort_order =~ /^n/i )    #sorting by name
        {
            @chr = sort { $a cmp $b } keys %data;
        }
        elsif ( $chr_sort_order =~ /^s/i )    #sorting by size
        {
            @chr =
              sort { $data{$b}{chr_length} <=> $data{$a}{chr_length} }
              keys %data;
        }
        foreach my $chr (@chr) {
            next if $seen{$chr};
            if ($skip_non_ordered) {
                delete $data{$chr};
            }
            else {
                push @ordered, $chr;
            }
        }
    }
    else    #no predefined order
    {
        #how to sort chromosomes for diplay?
        if ( $chr_sort_order =~ /^n/i )    #sorting by name
        {
            my @numbered;
            my @lettered;
            foreach my $chr ( keys %data ) {
                if ( $chr =~ /\d+/ ) {
                    push @numbered, $chr;
                }
                else {
                    push @lettered, $chr;
                }
            }
            @ordered = (
                ( sort { chr_sort($a) <=> chr_sort($b) } @numbered ),
                ( sort { $a cmp $b } @lettered )
            );
        }
        elsif ( $chr_sort_order =~ /^s/i )    #sorting by size
        {
            @ordered =
              sort { $data{$b}{chr_length} <=> $data{$a}{chr_length} }
              keys %data;
        }
    }
    my $pos = 1;
    foreach my $item (@ordered) {
        $data{$item}{start} = $pos - 1;
        $pos += $data{$item}{length};
        $data{$item}{rev} = 1 if $rev{$item};
    }
    return \%data;
}

sub chr_sort {
    my $item = shift;
    if ( $item =~ /(\d+)(.*)/ ) {
        return $1 . "." . ord($2);
    }
    return $item;
}

sub get_ksdata {
    my %opts  = @_;
    my $ks_db = $opts{ks_db};
    my $pairs = $opts{pairs};
    my $type  = $opts{ks_type};
    my $min   = $opts{min};
    my $max   = $opts{max};
    my $log   = $opts{log};
    my %data;
    return \%data unless -r $ks_db;
    my $select = "select * from ks_data";
    my $dbh    = DBI->connect( "dbi:SQLite:dbname=$ks_db", "", "" );
    my $sth    = $dbh->prepare($select);
    $sth->execute();

    while ( my $data = $sth->fetchrow_arrayref ) {
        if ($pairs) {
            unless ( $pairs->{ $data->[1] }{ $data->[2] } ) {
                #		print STDERR "ks pair is not in pairs list\n";
                next;
            }
        }
        my %item = (
            KS      => $data->[3],
            KN      => $data->[4],
            'KN_KS' => $data->[5],
        );
        my $val = $item{$type};
        next unless defined $val && $val =~ /\d/;
        #	if (defined $min || defined $max)
        #	  {
        #	  }
        $data{ $data->[1] }{ $data->[2] } = $val;

    }
    $sth->finish();
    undef $sth;
    $dbh->disconnect();
    return \%data;
}

sub get_pairs {
    my %opts = @_;
    my $file = $opts{file};
    my $chr1 = $opts{chr1};
    my $chr2 = $opts{chr2};
    my %data;
    open( IN, $file ) || die $!;
    while (<IN>) {
        chomp;
        next if /^#/;
        next unless $_;
        my @line  = split /\t/;
        my @item1 = split(/\|\|/, $line[1]);
        my @item2 = split(/\|\|/, $line[5]);
        next unless $item1[6] && $item2[6];
        if ($chr1) {
            next unless $item1[0] eq $chr1 || $item2[0] eq $chr1;
        }
        if ($chr2) {
            next unless $item1[0] eq $chr2 || $item2[0] eq $chr2;
        }
        $data{ $item1[6] }{ $item2[6] } = 1;
        $data{ $item2[6] }{ $item1[6] } = 1;
    }
    close IN;
    return \%data;
}

sub get_dsg_order {
    my %opts           = @_;
    my $dsg            = $opts{dsg};
    my $chr            = $opts{chr};
    my $minsize        = $opts{minsize};
    my $chr_sort_order = $opts{chr_sort_order};
    my $skip_random    = $opts{skip_random}; #skip "random" chromosome where sequences are added ad hoc

    my %data;
#    foreach my $gs ( $dsg->genomic_sequences ) {
#        next if check_random( $gs->chromosome ) && $skip_random;
#        next if defined $chr && $chr ne $gs->chromosome;
#        my $len = $gs->sequence_length;
#        next if $minsize && $minsize > $len;
#        if ( $data{ $gs->chromosome } ) {
#            warn "Duplicate chromosome:" . $gs->chromosome . "\n";
#        }
#        $data{ $gs->chromosome }{chr_length} = $len;
#        $data{ $gs->chromosome }{length}     = $len;
#    }
	my $c = CoGe::Core::Chromosomes->new($dsg->id);
	while ($c->next) {
        next if check_random( $c->name ) && $skip_random;
        next if defined $chr && $chr ne $c->name;
        my $len = $c->length;
        next if $minsize && $minsize > $len;
        if ( $data{ $c->name } ) {
            warn "Duplicate chromosome:" . $c->name . "\n";
        }
        $data{ $c->name }{chr_length} = $len;
        $data{ $c->name }{length}     = $len;
    }
    #how to sort chromosomes for display?
    my @ordered;
    if ( $chr_sort_order =~ /^n/i )    #sorting by name
    {
        @ordered = sort { versioncmp( $a, $b ) } keys %data;
        #old below
        #	my @numbered;
        #	my @lettered;
        #	foreach my $chr (keys %data) {
        #	    if ($chr =~ /\d+/) {
        #		    push @numbered, $chr;
        #	    }
        #	    else {
        #		    push @lettered, $chr;
        #	    }
        #	}
        #	@ordered = ( (sort {chr_sort($a) <=> chr_sort($b) } @numbered), (sort { $a cmp $b } @lettered));
    }
    elsif ( $chr_sort_order =~ /^s/i ) #sorting by size
    {
        @ordered =
          sort { $data{$b}{chr_length} <=> $data{$a}{chr_length} } keys %data;
    }
    return \@ordered, \%data;
}

sub check_random {
    my $chr = shift;
    return 1 if ( $chr =~ /random/i || $chr =~ /unknown/i || $chr =~ /^un$/i );
    return 0;
}

sub get_gene_info {
    my %opts  = @_;
    my $dsgid = $opts{dsgid};
    my $info  = $opts{info};
    my $dbh   = $coge->storage->dbh;
    foreach my $tmp_chr ( keys %$info ) {
        my $query = qq{
SELECT count(distinct(feature_id))
  FROM feature
  JOIN dataset_connector dc using (dataset_id)
 WHERE genome_id = $dsgid
   AND feature_type_id = 3
   AND feature.chromosome = '$tmp_chr'

};
        my ($res) = $dbh->selectrow_array($query);

        $info->{$tmp_chr}{gene_length} = $res;
        next unless defined $res;
        $info->{$tmp_chr}{length} = $res;
        #get gene order
        $query = qq{
SELECT feature_id
  FROM feature
  JOIN dataset_connector dc using (dataset_id)
 WHERE genome_id = $dsgid
   AND feature_type_id = 3
   AND feature.chromosome = '$tmp_chr'
 ORDER BY feature.start

};
        my $sth = $dbh->prepare($query);
        $sth->execute();
        my $i = 1;
        while ( my $row = $sth->fetchrow_arrayref ) {
            $info->{$tmp_chr}{gene}{ $row->[0] } = $i;
            $i++;
        }
    }
}

sub add_rev_info {
    my %opts  = @_;
    my $info  = $opts{info};
    my $assoc = $opts{assoc};
    foreach my $item (@$assoc) {
        $info->{ $item->{chr} }{rev} = $item->{rev};
    }
}

sub get_range {
    my %opts = @_;
    my $data = $opts{data};
    my $min  = $opts{min};
    my $max  = $opts{max};
    my $log  = $opts{log};

    my ( $set_max, $set_min, $non_zero_min );
    foreach my $k1 ( keys %$data ) {
        foreach my $k2 ( keys %{ $data->{$k1} } ) {
            next unless defined $data->{$k1}{$k2};
            $non_zero_min = $data->{$k1}{$k2}
              unless defined $non_zero_min || $data->{$k1}{$k2} == 0;
            $set_max = $data->{$k1}{$k2} unless defined $set_max;
            $set_min = $data->{$k1}{$k2} unless defined $set_min;
            $set_max = $data->{$k1}{$k2} if $data->{$k1}{$k2} > $set_max;
            $set_min = $data->{$k1}{$k2} if $data->{$k1}{$k2} < $set_min;
            $non_zero_min = $data->{$k1}{$k2}
              if $non_zero_min
                  && $data->{$k1}{$k2} < $non_zero_min
                  && $data->{$k1}{$k2} > 0;
        }
    }
    #let's minimize the data, if needed
    if ( defined $min || defined $max ) {
        foreach my $k1 ( keys %$data ) {
            foreach my $k2 ( keys %{ $data->{$k1} } ) {
                my $val = $data->{$k1}{$k2};
                next unless defined $val;
                if ($log) {
                    $val = $non_zero_min if $val == 0;
                    $val = log10($val);
                }
                if (   ( defined $min && $val < $min )
                    || ( defined $max && $val > $max ) )
                {
                    $data->{$k1}{$k2} = undef;
                    delete $data->{$k1}{$k2};
                    delete $data->{$k1} unless keys %{ $data->{$k1} };
                }
            }
        }
    }
    if ($log) {
        $set_max = log10($set_max)      if $set_max > 0;
        $set_min = log10($non_zero_min) if $non_zero_min;
    }
    $set_max = $max if defined $max && $max < $set_max;
    $set_min = $min if defined $min && $min > $set_min;
    return ( $set_max, $set_min );
}

sub get_color {
    my %opts         = @_;
    my $val          = $opts{val};
    my $color_scheme = $opts{color_scheme};
    $color_scheme = 1 unless defined $color_scheme;
    unless ( defined $val && $val >= 0 && $val <= 1 ) {
        print STDERR "in sub get_color val is not [0,1]: $val\n";
        return [ 0, 0, 0 ];
    }
    my $schemes = [
        [
            [ 255, 255, 0 ],      #yellow
            [ 200, 200, 0 ],      # green
            [ 0,   200, 0 ],      # green
            [ 0,   100, 100 ],    # green
            [ 0,   200, 200 ],    # cyan
            [ 0,   0,   200 ],    # blue
            [ 100, 0,   100 ],    #purple
            [ 200, 0,   200 ],    #magenta
            [ 200, 0,   0 ],      #red
            [ 100, 0,   0 ],      #red
            [ 200, 100, 0 ],      #orange
            [ 255, 126, 0 ],      #orange
        ],
        [
            [ 255, 255, 0 ],      #yellow
            [ 255, 0,   0 ],      #red
            [ 0,   255, 0 ],      # green
            [ 0,   255, 255 ],    # cyan
            [ 220, 0,   220 ],    #magenta
            [ 0,   0,   255 ],    # blue
        ],
        [
            [ 0,   0,   150 ],    # blue
            [ 220, 220, 20 ],     #yellow
            [ 255, 0,   0 ],      #red
        ],
        [
            [ 0,   200, 0 ],      # green
            [ 0,   0,   200 ],    # blue
            [ 220, 220, 20 ],     #yellow
            [ 255, 0,   0 ],      #red
        ],
        [
            [ 255, 0, 0 ],        # red
            [ 0,   0, 0 ],        #black
        ],
        [
            [ 255, 0,   0 ],      #red
            [ 255, 255, 0 ],      #yellow
            [ 0,   255, 0 ],      # green
            [ 0,   255, 255 ],    # cyan
            [ 220, 0,   220 ],    #magenta
            [ 0,   0,   255 ],    # blue
            [ 255, 0,   0 ],      #red
            [ 255, 255, 0 ],      #yellow
            [ 0,   255, 0 ],      # green
            [ 0,   255, 255 ],    # cyan
            [ 220, 0,   220 ],    #magenta
            [ 0,   0,   255 ],    # blue
        ],
        [
            [ 255, 0,   0 ],      #red
            [ 255, 255, 0 ],      #yellow
            [ 0,   255, 0 ],      # green
            [ 0,   255, 255 ],    # cyan
            [ 220, 0,   220 ],    #magenta
            [ 0,   0,   255 ],    # blue
            [ 255, 0,   0 ],      #red
            [ 255, 255, 0 ],      #yellow
            [ 0,   255, 0 ],      # green
            [ 0,   255, 255 ],    # cyan
            [ 220, 0,   220 ],    #magenta
            [ 0,   0,   255 ],    # blue
            [ 255, 0,   0 ],      #red
            [ 255, 255, 0 ],      #yellow
            [ 0,   255, 0 ],      # green
            [ 0,   255, 255 ],    # cyan
            [ 220, 0,   220 ],    #magenta
            [ 0,   0,   255 ],    # blue
        ],
        [
            [ 255, 0,   0 ],      #red
            [ 255, 155, 0 ],      #yellow
            [ 0,   255, 0 ],      # green
            [ 0,   255, 255 ],    # cyan
            [ 220, 0,   220 ],    #magenta
            [ 0,   0,   255 ],    # blue
            [ 255, 0,   0 ],      #red
            [ 255, 155, 0 ],      #yellow
            [ 0,   255, 0 ],      # green
            [ 0,   255, 255 ],    # cyan
            [ 220, 0,   220 ],    #magenta
            [ 0,   0,   255 ],    # blue
        ],
        [
            [ 255, 0,   0 ],      #red
            [ 255, 155, 0 ],      #yellow
            [ 0,   255, 0 ],      # green
            [ 0,   255, 255 ],    # cyan
            [ 220, 0,   220 ],    #magenta
            [ 0,   0,   255 ],    # blue
            [ 255, 0,   0 ],      #red
            [ 255, 155, 0 ],      #yellow
            [ 0,   255, 0 ],      # green
            [ 0,   255, 255 ],    # cyan
            [ 220, 0,   220 ],    #magenta
            [ 0,   0,   255 ],    # blue
            [ 255, 0,   0 ],      #red
            [ 255, 155, 0 ],      #yellow
            [ 0,   255, 0 ],      # green
            [ 0,   255, 255 ],    # cyan
            [ 220, 0,   220 ],    #magenta
            [ 0,   0,   255 ],    # blue
        ],
        [
            [ 0,   255, 255 ],    # blue
            [ 255, 0,   0 ],      #orange
            [ 0,   255, 255 ],    # blue
            [ 255, 0,   0 ],      #orange
            [ 0,   255, 255 ],    # blue
            [ 255, 0,   0 ],      #orange
        ],
        [
            [ 0,   0,  255 ],     # blue
            [ 255, 99, 33 ],      #orange
            [ 0,   0,  255 ],     # blue
            [ 255, 99, 33 ],      #orange
            [ 0,   0,  255 ],     # blue
            [ 255, 99, 33 ],      #orange
        ],

    ];
    my @colors = @{ $schemes->[$color_scheme] };
    @colors = reverse @colors;
    my ( $index1, $index2 ) = (
        ( floor( ( scalar(@colors) - 1 ) * $val ) ),
        ceil( ( scalar(@colors) - 1 ) * $val )
    );

    my $color = [];
    my $step  = 1 / ( scalar(@colors) - 1 );
    my $scale = $index1 * $step;
    my $dist  = ( $val - $scale ) / ($step);
    #    print join ("\t", $val, $step, $scale, $dist),"\n\t";
    for ( my $i = 0 ; $i <= 2 ; $i++ ) {
        my $diff = ( $colors[$index1][$i] - $colors[$index2][$i] ) * $dist;
        push @$color, sprintf( "%.0f", $colors[$index1][$i] - $diff );
    }
    return $color;
}

sub draw_boxes {
    my %opts  = @_;
    my $gd    = $opts{gd};
    my $boxes = $opts{boxes};
    #    my $color = $gd->colorAllocate(0,0,0);
    foreach my $box (@$boxes) {
        $gd->rectangle( @$box, $black );
    }

}

#Print out info on script usage
sub usage {
    print qq{
Welcome to $0

dagfile      | d       path to dag file containing all the hits

alignfile    | a       path to .aligncoords file generated by
                       dagchainer containing just the diags

width        | w       width of image in pixels (1024)

min_chr_size | mcs     minimim size of chromosome to be drawn (0)

dsgid1       | dsg1 | gid1    database id of genome on x-axis

dsgid2       | dsg2 | gid1    database id of genome on y-axis

chr1         | c1      only show data for this chromosome on the x axis

chr2         | c2      only show data for this chromosome on the y axis

basename     | b       base path and name for output files (test)

link         | l       link to be used in image map

link_type    | lt      are image map links for chromosome blocks or points:
                       1  ::   blocks  (Use "XCHR","YCHR" which will get the appropriate chr substituted in)
                       2  ::   points

flip         | f       flip axes (1 flip, 0 don't flip [DEFAULT])

grid         | g       add a positional grid to dotplot

ks_db        | ksdb    specify a sqlite database with synonymous/nonsynonymous data
                       to color syntenic points

ks_type| kst           specify the synonymous data to use (kS, kN, kn_ks) for coloring syntenic points

color_diags_type|cdt   specify how the diagonals are to be colored.
                       inv == colored by inversions
                       diag == colored by diagonals
                       --none-- == all diags are colored the same
                       NOTE:  this option is superceded by ks_db if that is specified

log                    log10 transform ks datadata (val = 0 set to minimum non-zero val)

max                    max ks val cutoff

min                    min ks val cutoff

assemble               Syntenic Path Assembly (SPA)
                       Will use the genome with FEWER pieces as the reference genome.
                       If set to 1, output will try to be assembled based on syntenic path
                       If set to 2, will not add any pieces that are not syntenic

                       REVERSE SPA: (genome with MORE pieces is reference genome)
                       If set to -1, output will try to be assembled based on syntenic path
                       If set to -2, will not add any pieces that are not syntenic

axis_metric  | am      metic (units) for dotplot axes:  nucleotides or genes.  Default: nucleotide distances

skip_random | sr       flag for skipping chromomses with the word "random" in their name

help         | h       print this message

};
    exit;
}
