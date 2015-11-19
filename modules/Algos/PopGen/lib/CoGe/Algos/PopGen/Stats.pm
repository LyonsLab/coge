#-------------------------------------------------------------------------------
# Purpose:	Generic stats functions
# Author:	Matt Bomhoff
# Created:	8/10/11, imported into CoGe 9/21/15
#-------------------------------------------------------------------------------
package CoGe::Algos::PopGen::Stats;

use warnings;
use strict;
use POSIX qw(ceil floor);
use base 'Exporter';

our @EXPORT = qw(min max mean stats);

sub min {
	my $x = shift;
	my $y = shift;
	
	return $x if (defined $x and $x <= $y);
	return $y;
}

sub max {
	my $x = shift;
	my $y = shift;
	
	return $x if (defined $x and $x >= $y);
	return $y;
}

sub bound {
	my $x = shift;
	my $min = shift;
	my $max = shift;
	
	return $min if ($x <= $min);
	return $max if ($x >= $max);
	return $x;
}

sub mean {
	my $pArray = shift; # ref to array of scalars
	my $n = scalar @$pArray;
	die if (not defined $pArray or ref($pArray) ne 'ARRAY');
	return if ($n == 0);
	
	my $total = 0;
	foreach my $x (@$pArray) {
		die if (not defined $x);
		$total += $x;
	}
	
	return ($total / $n);
}

sub stats {
	my $pArray = shift; # ref to array of scalars
	my $n = scalar @$pArray;
	die if (not defined $pArray or ref($pArray) ne 'ARRAY');
	return if ($n == 0);
	
	# Compute average
	my ($min, $max);
	my $total = 0;
	foreach my $x (@$pArray) {
		die if (not defined $x);
		$total += $x;
		$min = $x if (not defined $min or $x < $min);
		$max = $x if (not defined $max or $x > $max);
	}
	my $avg = $total / $n;
	
	# Compute variance
	my $var = 0;
	foreach my $x (@$pArray) {
		my $diff = $x - $avg;
		$var += $diff * $diff / $n;
	}
	
	# Compute median
	my @a = sort {$a<=>$b} @$pArray;
	my $middle = (@a+1)/2 - 1;
	my $a1 = $a[floor($middle)];
	my $a2 = $a[ceil($middle)];
	my $med = ($a1+$a2)/2;
	
	# Find 95% confidence interval
	my $ci_lower = percentile(0.025, \@a);
	my $ci_upper = percentile(0.975, \@a);
	
	return ($avg, $med, $var, $min, $max, $ci_lower, $ci_upper);
}

# Linear interpolation method, source: http://en.wikipedia.org/wiki/Percentile
sub percentile {
	my $p = shift; # percentile
	my $a = shift; # ref to sorted (ascending) array of scalars
	my $n = scalar @$a;
	
	my $i1 = floor(($n+1) * $p - 1);
	$i1 = bound($i1, 0, $n-1);
	my $i2 = ceil(($n+1) * $p - 1);
	$i2 = bound($i2, 0, $n-1);
	return $a->[$i1] if ($i1 == $i2);
	my $p1 = (100/$n)*(($i1+1)-0.5);
	my $v = $a->[$i1] + $n * ($p*100-$p1) * ($a->[$i2]-$a->[$i1]) / 100;
	
	return $v; # the value of percentile $p
}

1;