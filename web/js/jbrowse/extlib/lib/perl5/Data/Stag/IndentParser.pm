# $Id: IndentParser.pm,v 1.4 2008/06/03 17:31:15 cmungall Exp $
#
# Copyright (C) 2002 Chris Mungall <cjm@fruitfly.org>
#
# See also - http://stag.sourceforge.net
#
# This module is free software.
# You may distribute this module under the same terms as perl itself

package Data::Stag::IndentParser;

=head1 NAME

  Data::Stag::IndentParser.pm     - parses stag indent format into stag events

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
    return 'indent';
}

our $elt_re = '^([\w\-\+\?\*\:]+)\s*(.*)';

sub parse_fh {
    my $self = shift;
    my $fh = shift;
    my @stack = ();

    while(<$fh>) {
        chomp;
        s/\\\#/MAGIC_WORD_HASH/g;
        s/\#.*//;
        s/MAGIC_WORD_HASH/\#/g;

        # remove trailing ws
        s/\s*$//;
        next unless $_;

        # get indent level
        /(\s*)(.*)/s;

        my ($indent_txt, $elt) = ($1, $2);
        my $indent = length($indent_txt);
        if ($elt =~ /^$elt_re/) {
            $elt = $1;
            my $rest = $2;
            my @atts = ();
            if ($rest =~ /^\s*\[(.*)/) {
                $rest = $1;
                while ($rest =~ /^\s*([\w:]+)\s*=\s*(.*)/) {
                    push(@atts, $1);
                    $rest = $2;
                    if ($rest =~ /(.*?)([\s?\]].*)/) {
                        push(@atts, $1);
                        $rest = $2;
                    }
                    else {
                        confess("expected SPC | ], found $rest in $_");
                    }
                }
                if ($rest =~ /^\s*\]\s*(.*)/) {
                    $rest = $1;
                }
                else {
                    confess("expected ] | TAG=VAL, found $rest");
                }
            }

            $self->pop_to_level($indent, \@stack);
            $self->start_event($elt);
            if (@atts) {
                $self->start_event('@');
                while (my ($tag,$val) = splice(@atts,0,2)) {
                    $self->event($tag,$val);
                }
                $self->end_event('@');
            }
            if ($rest =~ /^\s*\"(.*)/) {
                my $buf = $1;
                my $str = '';
                my $done = 0;
                while (!$done) {
                    while (my $c = substr($buf,0,1)) {
                        $buf = substr($buf,1);
                        if ($c eq '"') {
                            $done = 1;
                            last;
                        }
                        $str .= $c;
                    }
                    if (!$done) {
                        $buf = <$fh>;
                        if (!$buf) {
                            confess("expected \", found EOF");
                        }
                    }
                }
                $self->evbody($str);
            }
            push(@stack, [$indent, $elt]);
        }
        else {
            confess("expected INDENT-ELEMENT, found $_");
        }
    }
    $fh->close;
    $self->pop_to_level(0, \@stack);
    return;
}

sub pop_to_level {
    my $self = shift;
    my $indent = shift;
    my $stack = shift;

    while (scalar(@$stack) &&
           $stack->[-1]->[0] >= $indent) {
        $self->end_event($stack->[-1]->[1]);
        pop(@$stack);
    }
    return;
}

1;
