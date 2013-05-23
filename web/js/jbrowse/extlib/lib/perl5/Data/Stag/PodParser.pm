# $Id: PodParser.pm,v 1.10 2008/06/03 17:31:15 cmungall Exp $
#
# Copyright (C) 2002 Chris Mungall <cjm@fruitfly.org>
#
# See also - http://stag.sourceforge.net
#
# This module is free software.
# You may distribute this module under the same terms as perl itself

package Data::Stag::PodParser;

=head1 NAME

  PodParser.pm     - parses perl POD documentation into stag events

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

sub trim {
    my $w = shift;
    $w = '' unless defined $w;
    $w =~ s/^\s+//g;
    $w =~ s/\s+$//g;
    return $w;
}

sub parse_fh {
    my $self = shift;
    my $fh = shift;
    my @stack = ();
    my $txt;

    my $in_pod = 0;
    my $tt = 'text';
    $self->start_event("pod");
    my $buf;
    my $section;
    while(<$fh>) {
        chomp;
        if (/^=(\S+)\s*(.*)/) {
            $tt = 'text';
            my ($ev, $data) = ($1, $2);
            $in_pod = $ev eq 'cut' ? 0 : 1;
            $data =~ s/\s*$//;
            $section = $data;
            if ($buf) {
                $self->event($tt=>trim($buf));
                $buf = undef;
            }
            $self->pop_stack_to_depth(1);
#            $self->event($ev, $data);
            if ($data) {
                $self->start_event('section');
		$self->event(type=>$ev);
                $self->event(name=>$data);
            }
        }
        else {
            if (/^\s\s\s*(.*)/) {
                if ($tt eq 'text') {
                    $tt = 'code';
                    if ($buf) {
                        $self->event(text=>trim($buf));
                        $buf = undef;
                    }
                }        
            } 
            if ($in_pod) {
                $buf .= "$_\n";
            }
       }
    }
    if ($buf) {
        $self->event($tt=>trim($buf));
    }
    $fh->close;
    $self->pop_stack_to_depth(0);
    return;
}

1;
