package CoGe::Graphics::Feature::AminoAcid;
use strict;
use base qw(CoGe::Graphics::Feature);

=head1 AUTHOR

	Eric Lyons
	elyons@nature.berkeley.edu

=head1 COPYRIGHT

Permission to use, copy, modify, and distribute this software and its documentation for educational, research, and not-for-profit purposes, without fee and without a signed licensing agreement, is hereby granted, provided that the above copyright notice, this paragraph and the following two paragraphs appear in all copies, modifications, and distributions. Contact The Office of Technology Licensing, UC Berkeley, 2150 Shattuck Avenue, Suite 510, Berkeley, CA 94720-1620, (510) 643-7201, for commercial licensing opportunities.

The full text of the license can be found in the
LICENSE file included with this module.

=cut

BEGIN {
    use vars qw($VERSION $HEIGHT $WIDTH $FOB $PRO $GLY $BAS $ACD $RSH $ROH $ALA $ARO $AMN $HIS $TRP %AA $STOP);
    $VERSION     = '0.1';
    $HEIGHT = 25;
    $WIDTH = 25;
    $GLY = [235, 235, 235];
    $FOB = [15, 130, 15];
    $PRO = [220, 150, 130];
    $BAS = [20, 90, 255];
    $ACD = [230, 10, 10]; #pink
    $RSH = [230, 230, 0]; #yellow
    $STOP= [255, 255, 255];
    $ARO = [50, 50, 170];
    $ALA = [200, 200, 200];
    $HIS = [130, 130, 210];
    $ROH = [250,150, 0];
    $AMN = [0, 220, 220];
    $TRP = [180, 90, 180];

    %AA = (
    	   A =>'ala',
	   R =>'bas',
	   N =>'amn',
	   D =>'acd',
	   C =>'rsh',
	   E =>'acd',
	   Q =>'amn',
	   G =>'gly',
	   H =>'his',
	   I =>'fob',
	   L =>'fob',
	   K =>'bas',
	   M =>'rsh',
	   F =>'aro',
	   P =>'pro',
	   S =>'roh',
	   T =>'roh',
	   W =>'trp',
	   Y =>'aro',
	   V =>'fob',
	   '*' =>'stop',
          );
    __PACKAGE__->mk_accessors(
"aa",
"lastaa",
"no_three_D" #switch between flat images and "3D" images, if given value, use flat images
);
}

sub _initialize
  {
    #print STDERR "AA IS BEING CALLED\n";
    my $self = shift;
    my %opts = @_;
    my $h = $HEIGHT; #total image height
    my $w = $WIDTH;
    $self->image_width($w);
    $self->image_height($h);
    $self->merge_percent(73);
    $self->bgcolor([255,255,255]) unless $self->bgcolor;
    $self->label($self->aa) if $self->aa;
    $self->label(reverse($self->label)) if $self->strand =~ /-/;
    $self->stop($self->start + length($self->aa)*3-1);
    $self->type('aa');
    $self->skip_overlap_search(0) unless $self->skip_overlap_search;
    my ($sum, $fob, $pro, $gly, $bas, $acd, $rsh, $roh, $aro, $trp, $his, $amn, $ala, $stop)= (0,0,0,0,0,0,0,0,0,0,0,0,0,0);
    my @color;
    foreach my $aa (split //, $self->aa)
      {
        $aa = uc $aa;
	next unless $AA{$aa};
	if ($AA{$aa} eq 'fob')
	  {
	    $fob++;
	  }
	elsif ($AA{$aa} eq'gly')
	  {
	    $gly++;
	  }
	elsif ($AA{$aa} eq'pro')
	  {
	    $pro++;
	  }
	elsif ($AA{$aa} eq'bas')
	  {
	    $bas++;
	  }
	elsif ($AA{$aa} eq'acd')
	  {
	    $acd++;
	  }
	elsif ($AA{$aa} eq'rsh')
	  {
	    $rsh++;
	  }
	elsif ($AA{$aa} eq'roh')
	  {
	    $roh++;
	  }
	elsif ($AA{$aa} eq'ala')
	  {
	    $ala++;
	  }
	elsif ($AA{$aa} eq'aro')
	  {
	    $aro++;
	  }
	elsif ($AA{$aa} eq'trp')
	  {
	    $trp++;
	  }
	elsif ($AA{$aa} eq'his')
	  {
	    $his++;
	  }
	elsif ($AA{$aa} eq'amn')
	  {
	    $amn++;
	  }
	elsif ($AA{$aa} eq'stop')
	  {
	    $stop++;
	  }
        $sum++;
      }
    return 0 unless $sum;
    for my $i (0..2)
      {
	push @color, $FOB->[$i]*$fob/$sum+$GLY->[$i]*$gly/$sum+$PRO->[$i]*$pro/$sum+$BAS->[$i]*$bas/$sum
	             +$ACD->[$i]*$acd/$sum+$RSH->[$i]*$rsh/$sum+$ROH->[$i]*$roh/$sum+$ALA->[$i]*$ala/$sum+$ARO->[$i]*$aro/$sum+$TRP->[$i]*$trp/$sum+$HIS->[$i]*$his/$sum+$AMN->[$i]*$amn/$sum+$STOP->[$i]*$stop/$sum;
      }
    print "Amino acid color for ",join (":", $self->aa, @color),"\n" if $self->DEBUG;# if $self->aa =~ /\*/;
    @color = (255,255,255) if ($color[0] == 0 && $color[1] == 0 && $color[2] == 0);
    $self->color(\@color);
  }

sub _post_initialize
  {
    my $self = shift;
    my %opts = @_;
    my $gd = $self->gd;
    my ($x1, $x2, $y1, $y2) = (0,16,5,19);
    #print STDERR "la check ", $self->lastaa, "\n";
    my $gray_lvl = 21 || $opts{gray_lvl};
    $self->get_color(255,255,255);
    $self->gd->fill(12,12,$self->get_color(255,255,255));
    unless ($self->lastaa)
     {$gd->line(0,12,23,12,$self->get_color(0,0,0));}
    unless($self->no_three_D)
     {
       my @colors = $gd->rgb($self->get_color($self->color));
       my ($r, $g, $b) = ($colors[0],$colors[1],$colors[2]);
       if ($r + 80 <= 255)
	     {$r += 80;}
       if ($g + 80 <= 255)
	     {$g += 80;}
       if ($b + 80 <= 255)
	     {$b += 80;}
       $self->_make_3d(r=>$r, g=>$g, b=>$b, x1=>$x1, x2=>$x2, y1=>$y1, y2=>$y2, gray_lvl=>$gray_lvl);
     }
    else {
    $gd->filledRectangle($x1,$y1,$x2,$y2,$self->get_color($self->color));}
    $self->_rounded_edges(x1=>$x1, y1=>$y1);
    $self->_rounded_edges(x1=>$x2, y1=>$y1, negx=>-1);
    $self->_rounded_edges(x1=>$x2, y1=>$y2, negx=>-1, negy=>-1);
    $self->_rounded_edges(x1=>$x1, y1=>$y2, negy=>-1);

    #$gd->rectangle(3,14,19,30,$self->get_color(60,60,60));
    #$gd->rectangle(12,12,25,25,$self->get_color(0,0,0));

  }

sub _make_3d
  {
    my $self = shift;
    my %opts = @_;
    my $gd = $self->gd;
    my $r = $opts{r};
    my $g = $opts{g};
    my $b = $opts{b};
    my ($x1, $x2, $y1, $y2, $gray_lvl) = ($opts{x1},$opts{x2},$opts{y1},$opts{y2},$opts{gray_lvl} || 15);
    my $draw_lines = $y1;
   # my ($tmp1, $tmp2, $tmp3, $tmp4) = (0,0,0,0);
    my $color;
     while ($draw_lines <= $y2)
	{
	  $color = $self->get_color($r,$g,$b);
          $gd->line($x1, $draw_lines, $x2, $draw_lines, $color);
	  if ($draw_lines / ($y2 - $y1) < 0.38) {
	    if ($r + $gray_lvl <=255)
	     {$r += $gray_lvl;}
	    if ($g + $gray_lvl <=255)
	     {$g += $gray_lvl;}
	    if ($b + $gray_lvl <=255)
	     {$b += $gray_lvl;}
	    }
	  else {
	    if ($r - $gray_lvl >=0)
	     {$r -= $gray_lvl;}
	    if ($g - $gray_lvl >=0)
	     {$g -= $gray_lvl;}
	    if ($b - $gray_lvl >=0)
	     {$b -= $gray_lvl;}
	    }
 	  $draw_lines++;
	}
}

sub _rounded_edges
{
	my $self = shift;
	my %opts = @_;
	my $x1 = $opts{x1};
	my $y1 = $opts{y1};
	my $neg_x = $opts{negx} || 1;
	my $neg_y = $opts{negy} || 1;
	my $gd = $self->gd;
	my $poly1 = GD::Polygon->new;
        my $poly2 = GD::Polygon->new;
        my $poly3 = GD::Polygon->new;
        $poly1->addPt($x1,$y1);
        $poly1->addPt($x1+($neg_x*3),$y1);
        $poly1->addPt($x1,$y1+($neg_y*1));
        $poly2->addPt($x1,$y1);
        $poly2->addPt($x1+($neg_x*2),$y1);
        $poly2->addPt($x1,$y1+($neg_y*2));
        $poly3->addPt($x1,$y1);
        $poly3->addPt($x1+($neg_x*1),$y1);
        $poly3->addPt($x1,$y1+($neg_y*3));
        $gd->filledPolygon($poly1, $self->get_color(255,255,255));
        $gd->filledPolygon($poly2, $self->get_color(255,255,255));
        $gd->filledPolygon($poly3, $self->get_color(255,255,255));
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

#################### main pod documentation begin ###################
## Below is the stub of documentation for your module.
## You better edit it!

=head1 NAME

CoGe::Graphics::Feature::Base

=head1 SYNOPSIS

  use CoGe::Graphics::Feature::Base

=head1 DESCRIPTION

=head1 USAGE

=head1 BUGS

=head1 SUPPORT

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

1;
# The preceding line will help the module return a true value
