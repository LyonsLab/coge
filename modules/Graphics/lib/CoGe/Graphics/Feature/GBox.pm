package CoGe::Graphics::Feature::GBox;
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
    use vars qw($VERSION $HEIGHT $WIDTH %EXTERNAL_IMAGES);
    $VERSION     = '0.1';
    $HEIGHT = 5;
    $WIDTH = 5;
    %EXTERNAL_IMAGES = (
			A=>'/opt/apache/CoGe/picts/A.png',
			T=>'/opt/apache/CoGe/picts/T.png',
			C=>'/opt/apache/CoGe/picts/C.png',
			G=>'/opt/apache/CoGe/picts/G.png',
		       );
    __PACKAGE__->mk_accessors(
"nt",
"show_label",
"extra",
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
    $self->fill(1);
    $self->order(1);
    $self->stop($self->start + length $self->nt-1) unless $self->stop;
    $self->skip_overlap_search(1); #make sure to skip searching for overlap for these guys.  Search can be slow
    my $gbox = 0;
    my $sum = 0;
    my $seq = $self->nt;
#    print STDERR length ($seq),"\n";
    if (length ($seq) > 5)
      {
	while ($seq=~ /(CACGTG)/ig)
	  {
	    $sum += length $1;
	    $gbox++;
	  }
	while ($seq=~ /(GTGCAC)/ig)
	  {
	    $sum += length $1;
	    $gbox++;
	  }
      }
    my $p = ($sum)/(length ($seq));
    my @color;
    my $red = 255;
    $red -= ($p*200);
    my $blue = 255;
    my $green = 255;

    @color = ($red, $red, $blue);

    $self->color(\@color);
    $self->label("GBOX\n".$self->nt) if $self->nt && $self->show_label && $sum;
    $self->type('gbox');
  }

sub _post_initialize
  {
    my $self = shift;
    my %opts = @_;
    my $gd = $self->gd;
    $gd->fill(0,0, $self->get_color($self->color));
    if ($self->label && length ($self->label) == 1 && $EXTERNAL_IMAGES{uc($self->label)} && -r $EXTERNAL_IMAGES{uc($self->label)})
      {
	my $ei= GD::Image->new($EXTERNAL_IMAGES{uc($self->label)});
#	print STDERR Dumper ($ei);
	if ($self->strand =~ /-/ || $self->strand =~ /bot/i)
	  {
	    $ei->rotate180();
	  }
	$self->external_image($ei) if $self->use_external_image;
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
