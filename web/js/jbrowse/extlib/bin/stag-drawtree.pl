#!/usr/bin/perl -w

eval 'exec /usr/bin/perl -w -S $0 ${1+"$@"}'
    if 0; # not running under some shell

# POD docs at end

use strict;

use Data::Stag qw(:all);
use Getopt::Long;

use FileHandle;
use GD;

my $e = "";
my $fmt = 'png';
my $out;
my $parser;
my $fontchoice;
GetOptions("element|e=s"=>\$e,
	   "out|o=s"=>\$out,
           "parser|format|p=s" => \$parser,
	   "font=s"=>\$fontchoice,
	   "help|h"=>sub { system("perldoc $0"); exit },
	  );

my $fn = shift @ARGV;
my $fh;
if ($fn eq '-') {
    $fh = \*STDIN;
}
else {
    $fh = FileHandle->new($fn) || die $fn;
}
my $tab = 4;
my $stag = Data::Stag->parse(-file=>$fn, -format=>$parser);
my $box = box($stag, 0, 0);
my $FONT = gdLargeFont;
if ($fontchoice) {
    if ($fontchoice eq 'gdMediumBoldFont') {
	$FONT = gdMediumBoldFont;
    }
}
my $CW = $FONT->width;
my $CH = $FONT->height;
my $width = $CW * $box->{x2};
my $height = $CH * $box->{y2};
my $im = new GD::Image($width, $height);

my $bgcol = $im->colorAllocate(255, 255, 255);

my $fgcol = $im->colorAllocate(0, 0, 0);
my $linecol = $im->colorAllocate(0, 0, 0);
my $namecol = $im->colorAllocate(0, 0, 0);
my $datacol = $im->colorAllocate(255, 0, 0);
#$im->transparent($bgcol);
$im->interlaced('true');


drawbox($im, $box);
if ($out) {
    open(F, ">$out") || die($out);
    print F $im->$fmt();
    close(F);
}
else {
    print $im->$fmt();
}

sub box {
    my $stag = shift;
    my ($x, $y)= @_;

    my $y1 = $y;
    my $x2 = $x;
    my @kids = $stag->kids;
    my @subboxes = ();
    my $data;
    if ($stag->isterminal) {
	$data = $stag->data;
	$x2 += length($data) + length($stag->element) + 1;
	@kids = ();
    }
    $y++;
    foreach my $kid (@kids) {
	my $subbox = box($kid, $x + $tab, $y);
	push(@subboxes, $subbox);
	$y = $subbox->{y2};
	if ($subbox->{x2} > $x2) {
	    $x2 = $subbox->{x2};
	}
    }
    return 
      {boxes=>\@subboxes,
       name=>$stag->element,
       data=>$data,
       x1=>$x,
       x2=>$x2,
       y1=>$y1,
       y2=>$y,
      };
}

sub drawbox {
    my $im = shift;
    my $box = shift;
    my $sx1 = $box->{x1} * $CW;
    my $sy1 = $box->{y1} * $CH;
    my $name = $box->{name};
    my $data = $box->{data};
    
    $im->string($FONT,
		$sx1,
		$sy1,
		$name,
		$namecol);
    if ($data) {
	$im->string($FONT,
		    $sx1 + length($name) * $CW + 1 * $CW,
		    $sy1,
		    $data,
		    $datacol);
    }
    my $subboxes = $box->{boxes};
    my $XOFF = int($CW / 2);
    my $YOFF = int($CH / 2);
    foreach (@$subboxes) {
	drawbox($im, $_);
	my $ix = $sx1 + $XOFF;
	my $iy = $_->{y1} * $CH + $YOFF;
	my $jx = $_->{x1} * $CW;
	$im->line($ix, $sy1 + $CH, $ix, $iy, $linecol);
	$im->line($ix, $iy, $jx, $iy, $linecol);
    }
}

exit 0;

__END__

=head1 NAME 

stag-drawtree.pl - draws a stag file (xml, itext, sxpr) as a PNG diagram

=head1 SYNOPSIS

  stag-drawtree.pl -o my.png myfile.xml

  stag-drawtree.pl -p My::MyFormatParser -o my.png myfile.myfmt

=head1 DESCRIPTION

requires GD library and GD perl module

=head1 ARGUMENTS

=over

=item -p|parser FORMAT

FORMAT is one of xml, sxpr or itext, or the name of a perl module

xml assumed as default

=back


=head1 SEE ALSO

L<Data::Stag>

=cut

