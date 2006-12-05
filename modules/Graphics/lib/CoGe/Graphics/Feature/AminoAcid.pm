package CoGe::Graphics::Feature::AminoAcid;
use strict;
use base qw(CoGe::Graphics::Feature);


BEGIN {
    use vars qw($VERSION $HEIGHT $WIDTH $FOB $PRO $FIL $BAS $ACD $CYS %AA);
    $VERSION     = '0.1';
    $HEIGHT = 25;
    $WIDTH = 25;
    $FIL = [100, 100, 100];
    $FOB = [200, 200, 200];
    $PRO = [100, 255, 100];
    $BAS = [255, 100, 100];
    $ACD = [100, 100, 255]; #pink
    $CYS = [255, 255, 100]; #yellow
    %AA = (
    	   A =>'fob',
	   R =>'bas',
	   N =>'fil',
	   D =>'acd',
	   C =>'cys',
	   E =>'acd',
	   Q =>'fil',
	   G =>'fob',
	   H =>'bas',
	   I =>'fob',
	   L =>'fob',
	   K =>'bas',
	   M =>'fob',
	   F =>'fob',
	   P =>'pro',
	   S =>'fil',
	   T =>'fil',
	   W =>'fob',
	   Y =>'fob',
	   V =>'fob',
          );	
    __PACKAGE__->mk_accessors(
"aa",
);
}

sub _initialize
  {
    my $self = shift;
    my %opts = @_;
    my $h = $HEIGHT; #total image height 
    my $w = $WIDTH;
    $self->image_width($w);
    $self->image_height($h);
    $self->bgcolor([255,255,255]) unless $self->bgcolor;
    $self->label($self->aa) if $self->aa;
    $self->stop($self->start + length($self->aa)*3-1);
    $self->type('aa');
    $self->skip_overlap_search(0) unless $self->skip_overlap_search;
    my ($sum, $fob, $pro, $fil, $bas, $acd, $cys)= (0,0,0,0,0,0,0);
    my @color;
    foreach my $aa (split //, $self->aa)
      {
        $aa = uc $aa;
	next unless $AA{$aa};
	if ($AA{$aa} eq 'fob')
	  {
	    $fob++;
	  }
	elsif ($AA{$aa} eq'fil')
	  {
	    $fil++;
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
	elsif ($AA{$aa} eq'cys')
	  {
	    $cys++;
	  }
        $sum++;
      }
    return 0 unless $sum;
    for my $i (0..2)
      {
	push @color, $FOB->[$i]*$fob/$sum+ $FIL->[$i]*$fil/$sum+$PRO->[$i]*$pro/$sum+$BAS->[$i]*$bas/$sum
	             +$ACD->[$i]*$acd/$sum+$CYS->[$i]*$cys/$sum;
      }
    @color= (255,255,255) if ($color[0] == 0 && $color[1] == 0 && $color[2] == 0);
    $self->color(\@color);
  }

sub _post_initialize
  {
    my $self = shift;
    my %opts = @_;
    my $gd = $self->gd;
    $self->get_color(255,255,255);
    $self->gd->fill(12,12,$self->get_color(255,255,255));
    $gd->filledArc(12,12, 25,25,0,360,$self->get_color($self->color));
    #$gd->arc(5,5, 9,9,0,360,$self->get_color(0,0,0));

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

