# $Id: SxprParser.pm,v 1.22 2008/06/03 17:31:15 cmungall Exp $
#
# Copyright (C) 2002 Chris Mungall <cjm@fruitfly.org>
#
# See also - http://stag.sourceforge.net
#
# This module is free software.
# You may distribute this module under the same terms as perl itself

package Data::Stag::SxprParser;

=head1 NAME

  SxprParser.pm     - parses Stag S-expression format

=head1 SYNOPSIS

=cut

=head1 DESCRIPTION


=head1 AUTHOR

=cut

use Exporter;
use Carp;
use FileHandle;
use strict;
use Data::Stag qw(:all);
use base qw(Data::Stag::BaseGenerator Exporter);

use vars qw($VERSION);
$VERSION="0.11";

sub fmtstr {
    return 'sxpr';
}

sub parse_fh {
    my $self = shift;
    my $fh = shift;
    my $line_no = 0;
    my $state = "init";
    my $txt = '';

    # c: comment ;
    # q: quote (double quote only)
    # o: opening tag - after round bracket, before whitespace or ()
    # 0: body, or awaiting open or close
    while (my $line = <$fh>) {
        $line_no++;

        while (length($line)) {
            my $c = substr($line,0,1,'');
            if ($state eq 'init') {
                if ($c eq '(') {
                    # at start - do nothing
                    $state = 0;
                }
                elsif ($c eq "'") {
                    # leading quote is allowed [list constructor in lisp]
                    # (good for editing in emacs)
                    next;
                    $state = 0;
                }
                elsif ($c =~ /\s/) {
                    next;
                }
                elsif ($c eq ';') {
                    $state = 'c';
                }
                else {
                    $self->throw("Line: $line_no\n$line\nExpected \"(\" at start of file");
                }
            }
            $state ne 'init' || $self->throw("assertion error: state=$state");
            
            if ($state eq 'c') { # comment
                # newline is the only char that can break out of comments
                if ($c eq "\n") {
                    $state = 0;
                }
                next;
            }
            $state ne 'c' || $self->throw("assertion error: state=$state");
            
            if ($state eq 'q') {
                if ($c eq '"') {
                    $state = 0;
                }
                else {
                    $txt .= $c;
                }
                next;
            }
            $state ne 'q' || $self->throw("assertion error: state=$state");
            
            if ($c eq '"') {
                $state = 'q';
                next;
            }
            $state ne 'q' || $self->throw("assertion error: state=$state");
            
            if ($c eq ';') {
                # can only open comments when NOT in quotes
                $state = 'c';
                next;
            }
            
            if ($c eq '(') {
                if ($state eq 'o') {
                    if (!$txt) {
                        $self->throw("Line: $line_no\b$line\ndouble open brackets (( not allowed!");
                    }
                    $self->start_event($txt);
                    $txt = '';
                }
                $state = 'o';
                next;
            }
            if ($c eq ')') {
                if ($state eq 'o') {
                    if (!$txt) {
                        $self->throw("Line: $line_no\b$line\n () not allowed!");
                    }
                    $self->start_event($txt);
                    $txt = '';
                    $self->end_event;
                }
                else {
                    if ($txt) {
                        $self->evbody($txt);
                    }
                    $txt = '';
                    $self->end_event;
                }
                $state=0;
                next;
            }
            if ($state eq 'o') {
                if ($c =~ /\s/) {
                    # reached last char of start event name
                    $self->start_event($txt);
                    $txt = '';
                    $state = 0;
                    next;
                }
                else {
                    $txt .= $c;
                    next;
                }
            }
            $state == 0 || $self->throw("assertion error: state=$state");
            if ($c =~ /\s/) {
                next;
            }
            $txt .= $c;
            next;
        }
    }
    if ($txt =~ /\S/) {
        $self->throw("text at end: $txt");
    }
}

1;
