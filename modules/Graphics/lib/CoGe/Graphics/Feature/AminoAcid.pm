package CoGe::Graphics::Feature::AminoAcid;
use strict;
use base qw(CoGe::Graphics::Feature);


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
    $self->merge_percent(72);
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
    #print STDERR "la check ", $self->lastaa, "\n";
    my $gray_lvl = 5 || $opts{gray_lvl};;
    $self->get_color(255,255,255);
    $self->gd->fill(12,12,$self->get_color(255,255,255));
    unless ($self->lastaa)
     {$gd->line(0,12,23,12,$self->get_color(0,0,0));}
    unless($self->no_three_D)
     {
       my @colors = $gd->rgb($self->get_color($self->color));
       my ($r, $g, $b) = ($colors[0],$colors[1],$colors[2]);
       $self->_make_3d(r=>$r, g=>$g, b=>$b, x1=>0, x2=>16, y1=>5, y2=>19, gray_lvl=>$gray_lvl);
     }
    else {
     $gd->filledRectangle(0,5,16,19,$self->get_color($self->color));}

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
    my $color;
    while ($draw_lines <= $y2)
	{
	  #print STDERR "r: ", $r, ", g: ", $g, ", b: ", $b, "\n";
	  $color = $self->get_color($r,$g,$b);
	  $gd->line($x1, $draw_lines, $x2, $draw_lines, $color);
	  if ($draw_lines / ($y2 - $y1) < 0.42) {
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

