package FastBit;

=head1 NAME

FastBit

=head1 SYNOPSIS

Functions to interact with FastBit indexes

=head1 DESCRIPTION

=head1 AUTHOR

Sean Davey

=head1 COPYRIGHT

The full text of the license can be found in the
LICENSE file included with this module.

=head1 SEE ALSO

=cut

use strict;
use warnings;

use CoGe::Accessory::Web qw(get_defaults get_command_path);
use CoGe::Core::Storage qw(get_experiment_path);
use Data::Dumper;

BEGIN {
	use Exporter 'import';
	our @EXPORT_OK =
	  qw( max min query );
}

sub max {
	my $eid = shift;
	my $chr = shift;
	my $where = '0.0=0.0';
	$where .= " and chr='$chr'" if $chr;
	my $max = query('select max(value1) where ' . $where, $eid);
	return 0 + $max->[0];
}

sub min {
	my $eid = shift;
	my $chr = shift;
	my $where = '0.0=0.0';
	$where .= " and chr='$chr'" if $chr;
	my $min = query('select min(value1) where ' . $where, $eid);
	return 0 + $min->[0];
}

sub query {
	my $query = shift;
	my $eid = shift;

	my $cmdpath = get_command_path('FASTBIT_QUERY', 'ibis');
	my $storage_path = get_experiment_path($eid);
	my $cmd = "$cmdpath -v 1 -d $storage_path -q \"$query\" 2>&1";
#	warn $cmd;
	my @cmdOut = qx{$cmd};

	my $cmdStatus = $?;
	if ( $? != 0 ) {
		warn "FastBit::query: error $? executing command: $cmd";
	}
	my @lines;
	foreach (@cmdOut) {
		chomp;
		if (/^(\"|\d)/) {
			push @lines, $_;
		}
	}
	return \@lines;
}

1;
