##########################################################################################################
#	GD::Graph::histogram
#
#	Copyright 2007, Snehanshu Shah
#modified by eric lyons
#
##########################################################################################################
package CoGe::Accessory::histogram;
use strict;

BEGIN {
	use Exporter ();
	use vars qw ($VERSION @ISA @EXPORT @EXPORT_OK %EXPORT_TAGS);
	$VERSION     = 1.10;
	@ISA         = qw (Exporter);
	#Give a hoot don't pollute, do not export more than needed by default
	@EXPORT      = qw ();
	@EXPORT_OK   = qw ();
	%EXPORT_TAGS = ();
}

use GD::Graph::bars;
use GD::Graph::Data;
use Data::Dumper;

@CoGe::Accessory::histogram::ISA = qw(GD::Graph::bars);

my %Defaults = (
	histogram_bins => undef,
	histogram_type => 'count'
);

sub plot
{
    my $self = shift;
    my %opts = @_;
    my $dataArrRef = $opts{data};
    my $min = $opts{min};
    my $max = $opts{max};

    my $histogram_bins = $self->get('histogram_bins');

    my $cp = _histogram_bins( $dataArrRef, $histogram_bins, $min, $max );

    my $binArrRef = _histogram_frequency( $dataArrRef, $cp );

    if (defined $self->get('histogram_type') && $self->get('histogram_type') eq 'percentage')
    {
	    my $total = 0;
	    grep($total += $_ , @$binArrRef);

	    if ($total > 0)
	    {
	    	for (my $i = 0; $i < scalar(@$binArrRef); $i++)
	    	{
		    	$binArrRef->[$i] = 100 * $binArrRef->[$i] / $total;
	    	}
	    }
    }

    my @labelArr;
    for my $bin (@$cp)
    {
	    push(@labelArr, _numformat( $bin->[0] + ($bin->[1] - $bin->[0])/2 ) );
    }

    # Display the labels veritcally for histogram
    $self->set( x_labels_vertical => 1 );

    my $data = GD::Graph::Data->new([ \@labelArr, $binArrRef ]) or die GD::Graph::Data->error;
    return $self->SUPER::plot($data);

}

sub add_colour
  {
    my $self = shift;
    GD::Graph::bars::add_colour(@_);
  }

###################################################################
# _histogram_bins - calculates the bins usings Scott's algorithm
#
# Arguements:
#
# $data  (Vector)  - Data values
#
# $nbins (Integer) - Number of bins to create. If $nbins is undef
#                    the number of bins is calculated using Scott's
#                    algorithm
#
###################################################################
sub _histogram_bins {

	my ( $data, $nbins, $min_bin, $max_bin ) = @_;
	if( !defined $data ) { return; }

	my $calcBins = ( defined $nbins )? 0 : 1;
	my $cnt = 0;
	my $mean= 0;
	my $max = my $min = $data->[0];
	foreach (@$data) {
		$mean += $_;
		$min = ( $_ < $min )? $_ : $min;
		$max = ( $_ > $max )? $_ : $max;
		$cnt++;
	}
	$min = $min_bin if defined $min_bin && $min_bin < $min;
	$max = $max_bin if defined $max_bin && $max_bin > $max;
	$mean /= $cnt if( $cnt > 1 );

	my $sumsq = 0;
	$nbins = 1 if( $calcBins );
	my $s = 0;
	if( $cnt > 1 ) {
		foreach (@$data) {
			$sumsq += ( $_ - $mean )**2;
		}
		$s = sqrt( $sumsq / ($cnt - 1));
		$nbins = 3.49 * $s / $cnt**0.33 if( $s > 0 && $calcBins );
	}

	my $binwidth = ( $max - $min ) / $nbins;

	my $lower = $min;
	my $upper = $lower;

	my $bins;
	my @cutPoints;
	my $cntr = 0;
	while ( $upper <= $max && $cntr < $nbins) {
		$upper = $lower + $binwidth;
		push( @cutPoints, [$lower, $upper] );
		$lower = $upper;
		$cntr++;
	}

	return \@cutPoints;
}

###################################################################
# _histogram_frequency - bins the data
#
#     Lower Boundry <= data value < Upper Boundry
#
# Arguements:
#
# $data  (Vector)  - Data values
#
# $nbins (Integer) - Vector containing the cutpoints to bin the data
#
###################################################################
sub _histogram_frequency {
	my ( $data, $cutPoints ) = @_;

	if( !defined $data || !defined $cutPoints ) { return; }

	my @freqs;
	foreach (@$cutPoints) {
		push( @freqs, 0 );
	}

	foreach (@$data)
	{
		for( my $i = 0; $i < scalar( @$cutPoints ); $i++ )
		{
		if( ($_ >= $cutPoints->[$i]->[0] && $_ < $cutPoints->[$i]->[1])
			||
			($i == (scalar (@$cutPoints) - 1) && $_ >= $cutPoints->[$i]->[1]) )
			{

				$freqs[$i]++;
			}
		}
	}
	return \@freqs;
}

sub _numformat {
	my ($v, $f1, $f2) = @_;

	unless(defined $v) { return undef; }

	unless(defined $f1) { $f1 = "%.4e"; }

	unless(defined $f2) {
		if ($v < 1) {
			$f2 = "%.5f";
		} else {
			$f2 = "%.3f";
		}
	}

    	## To display no for eg. 22.50 as 22.5
    	if ($v =~ /^([+-]?\d+)\.(\d+)$/) {
    		my $no = $1;
		my $fraction = $2;
		$fraction =~ s/0+$//;
		$v = (length($fraction) == 0) ? $no : "$no.$fraction";
	}

	if ($v =~ /\./){
		if ($v == 0) {
			$v = 0;
		} elsif (($v > -0.001) and ($v < 0.001)) {
			$v = sprintf($f1, $v);
		} else {
			$v = sprintf($f2, $v);
		}
	}

	return $v;
}

sub _has_default {
    my $self = shift;
    my $attr = shift || return;
    exists $Defaults{$attr} || $self->SUPER::_has_default($attr);
}

########################################### main pod documentation begin ##
# Below is the stub of documentation for your module. You better edit it!

=head1 NAME

GD::Graph::histogram - Histogram plotting module for Perl5

=head1 SYNOPSIS

  use GD::Graph::histogram;

=head1 DESCRIPTION

GD::Graph::histogram extends the GD::Graph module to create histograms.
The module allow creation of count or percentage histograms.

=head1 USAGE

Fill an array with all the data values that are to be plotted. Note that
GD::Graph::histogram unlike the other GD::Graph modules can only plot one
data set at a time.

	$data = [1,5,7,8,9,10,11,3,3,5,5,5,7,2,2];

Create the graph

	my $graph = new GD::Graph::histogram(400,600);

Set graph options

	$graph->set(
		x_label         => 'X Label',
		y_label         => 'Count',
		title           => 'A Simple Count Histogram Chart',
		x_labels_vertical => 1,
		bar_spacing     => 0,
		shadow_depth    => 1,
		shadowclr       => 'dred',
		transparent     => 0,
	)
	or warn $graph->error;

plot the graph

	my $gd = $graph->plot($data) or die $graph->error;

save the graph to a file

	open(IMG, '>histogram.png') or die $!;
	binmode IMG;
	print IMG $gd->png;

=head1 METHODS

GD::Graph::histogram supports all the methods support by GD::Graph.
Please refer to the GD::Graph documentation for more information.

The only method that behaves differently is I<plot>

The I<plot> method provided by GD::Graph::histogram expects a
reference to an array of numbers.

Based on the input data, GD::Graph::histogram will generate the
appropriate labels for the X axis. The X axis label represent the center
point of the range of each histogram bin.

=head1 OPTIONS

GD::Graph::histogram supports all the options supported by GD::Graph::bars.
Please refer to the GD::Graph documentation for more information.

The two additional options that are specific to GD::Graph::histogram are:

	histogram_bins
		Specify the number of histogram bins to bucket the data into.
		The default is for the module to automatically computed the
		histogram bins based on the data.

	histogram_type
		Can be set to either 'percentage' or 'count'. By default the module
		will create 'count' histograms.

=head1 NOTES

As with all Modules for Perl: Please stick to using the interface. If
you try to fiddle too much with knowledge of the internals of this
module, you could get burned. I may change them at any time.

=head1 AUTHOR

	Snehanshu Shah
	perl@whizdog.com
	http://www.whizdog.com

=head1 ACKNOWLEDGEMENTS

Thanks for all the feedback, bug reports and bug fixes

Martin Corley
Jonathan Barber
William Miller

=head1 COPYRIGHT

This program is free software; you can redistribute
it and/or modify it under the same terms as Perl itself.

The full text of the license can be found in the
LICENSE file included with this module.

=head1 SEE ALSO

perl(1), GD::Graph

=cut

__END__
