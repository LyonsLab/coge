#-------------------------------------------------------------------------------
# Purpose:	Functions for dealing with results files
# Author:	Sean Davey
# Created:	11/17/2015
#-------------------------------------------------------------------------------
package CoGe::Algos::PopGen::Results;

use warnings;
use strict;
use Data::Dumper;
use File::Spec::Functions qw( catfile );
use JSON::XS;
use base 'Exporter';

our @EXPORT = qw( export get_data get_plot_data );

sub export {
	my $file = shift;
	my $type = shift;
	my @columns = map {substr $_, 3} split ',', shift;

	print "Content-Disposition: attachment; filename=PopGen.txt\n\n";
	my $in_type = 0;
	open my $fh, '<', $file;
	while (my $row = <$fh>) {
		chomp $row;
		if (substr($row, 0, 1) eq '#') {
			return if $in_type;
			my @tokens = split '\t', $row;
			$in_type = 1 if (substr($tokens[0], 1) eq $type);
			$row = <$fh>;
			chomp $row;
		}
		if ($in_type) {
			my @tokens = split '\t', $row;
			my $line;
			foreach (@columns) {
				$line .= "\t" if $line;
				$line .= $tokens[$_] if ($tokens[$_]);
			}
			print $line . "\n";
		}
	}
}

sub get_data {
	my $file = shift;
	my $columns = {};
	my $data;
	my $type = '';
	my @types;
	my $set;
	open my $fh, '<', $file;
	while (my $row = <$fh>) {
		chomp $row;
		my @tokens = split '\t', $row;
		if (substr($tokens[0], 0, 1) eq '#') {
			if ($type) {
				$data .= ',' if $data;
				$data .= '\'';
				$data .= $type;
				$data .= '\':[';
				$data .= $set;
				$data .= ']';
				$set = '';
			}
			$type = substr $tokens[0], 1;
			push @types, $type;
			shift @tokens;
			$columns->{$type} = [@tokens];
			$row = <$fh>;
			chomp $row;
			@tokens = split '\t', $row;
		}
		$set .= ',' if $set;
		$set .= encode_json(\@tokens);
	}
	$data .= ',' if $data;
	$data .= '\'';
	$data .= $type;
	$data .= '\':[';
	$data .= $set;
	$data .= ']';
	my $options;
	for (@types) {
		$options .= '<option>';
		$options .= $_;
		$options .= '</option>';
	}
	return encode_json($columns), '{' . $data . '}', $options;
}

sub get_plot_data {
    my $file = shift;
    my $type = shift;
    my $column = shift;
    my $in_type = 0;
    my @p;
    open my $fh, '<', $file;
    while (my $row = <$fh>) {
        chomp $row;
        if (substr($row, 0, 1) eq '#') {
            return if $in_type;
            my @tokens = split '\t', $row;
            $in_type = 1 if (substr($tokens[0], 1) eq $type);
            $row = <$fh>;
            chomp $row;
        }
        if ($in_type) {
            my @tokens = split '\t', $row;
            push @p, [($tokens[1] + $tokens[2]) / 2, $tokens[$column]];
        }
    }
    @p = sort { $a->[0] <=> $b->[0] } @p;
    my @x;
    my @y;
    for (@p) {
        push @x, shift @$_;
        push @y, shift @$_;
    }
    return \@x, \@y;
}

1;