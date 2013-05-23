# $Id: ITextParser.pm,v 1.19 2008/06/03 17:31:15 cmungall Exp $
#
# Copyright (C) 2002 Chris Mungall <cjm@fruitfly.org>
#
# See also - http://stag.sourceforge.net
#
# This module is free software.
# You may distribute this module under the same terms as perl itself

package Data::Stag::ITextParser;

=head1 NAME

  ITextParser.pm     - parses stag IText format into stag events

=head1 SYNOPSIS

=cut

=head1 DESCRIPTION


=head1 AUTHOR

=cut

use Exporter;
use Carp;
use FileHandle;
use strict;
use base qw(Data::Stag::BaseGenerator Exporter);

use vars qw($VERSION);
$VERSION="0.11";

sub fmtstr {
    return 'itext';
}

sub parse_fh {
    my $self = shift;
    my $fh = shift;
    my @stack = ();
    my $txt;

    while(<$fh>) {
        chomp;
        s/\\\#/MAGIC_WORD_HASH/g;
        s/\#.*//;
        s/MAGIC_WORD_HASH/\#/g;
        next unless $_;

        # remove trailing ws
        s/\s*$//;

	my $eofre = '\<\<(\S+)' . "\$";
	if (/$eofre/) {
	    my $eof = $1;
	    s/$eofre//;
	    my $sofar = $_;
	    while (<$fh>) {
		last if /^$eof/;
		$sofar .= $_;
	    }
	    $_ = $sofar;
	}

        # get indent level
        /(\s*)(.*)/s;

        my ($indent_txt, $elt) = ($1, $2);
        my $indent = length($indent_txt);
        if ($elt =~ /^([\w\-\+\?\*]+):\s*(.*)$/s ||
	    $elt =~ /^([\@\.]):\s*(.*)$/s) {
            $elt = $1;
            my $nu_txt = $2;

            $self->pop_to_level($indent, $txt, \@stack);
            $txt = undef;
            if ($nu_txt || length($nu_txt)) {
                $txt = $nu_txt;
            }
            $self->start_event($elt);
            push(@stack, [$indent, $elt]);
        }
        else {
            # body
            $txt .= $elt if $elt;
        }
    }
    $fh->close;
    $self->pop_to_level(0, $txt, \@stack);
    return;
}

sub pop_to_level {
    my $self = shift;
    my $indent = shift;
    my $txt = shift;
    my $stack = shift;

    # if buffered pcdata, export it
    if (defined $txt) {
	# unescape :s
	$txt =~ s/^\\:/:/g;
	$txt =~ s/([^\\])\\:/$1:/g;
        $self->evbody($txt);
    }
    while (scalar(@$stack) &&
           $stack->[-1]->[0] >= $indent) {
        $self->end_event($stack->[-1]->[1]);
        pop(@$stack);
    }
    return;
}

1;
