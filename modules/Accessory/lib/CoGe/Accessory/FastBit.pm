package CoGe::Accessory::FastBit;

=head1 NAME

CoGe::Accessory::FastBit

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

use CoGe::Accessory::Web qw(get_defaults);
use CoGe::Core::Storage qw(get_experiment_path);
use Data::Dumper;

BEGIN {
	use Exporter 'import';
	our @EXPORT_OK =
	  qw( query );
}

sub debug {
	my $data = shift;
	my $new_file = shift;
	my $OUTFILE;
	open $OUTFILE, ($new_file ? ">/tmp/sean" : ">>/tmp/sean");
	print {$OUTFILE} Dumper $data;
	print {$OUTFILE} "\n";
	close $OUTFILE;
}

sub query {
	my $query = shift;
	my $eid = shift;

    my $cmdpath = get_defaults()->{FASTBIT_QUERY};
    my $storage_path = get_experiment_path($eid);
    my $cmd = "$cmdpath -v 1 -d $storage_path -q \"$query\" 2>&1";
debug $cmd;
    my @cmdOut = qx{$cmd};

    my $cmdStatus = $?;
    if ( $? != 0 ) {
        warn "CoGe::Accessory::FastBit::query: error $? executing command: $cmd";
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
