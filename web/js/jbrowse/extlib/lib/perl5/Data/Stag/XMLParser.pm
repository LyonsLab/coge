# $Id: XMLParser.pm,v 1.16 2008/06/03 17:31:15 cmungall Exp $
#
# Copyright (C) 2002 Chris Mungall <cjm@fruitfly.org>
#
# See also - http://stag.sourceforge.net
#
# This module is free software.
# You may distribute this module under the same terms as perl itself

package Data::Stag::XMLParser;

=head1 NAME

  XMLParser.pm     - parses XML files into stag events

=head1 SYNOPSIS

=cut

=head1 DESCRIPTION


=head1 AUTHOR

=cut

use Exporter;
use Carp;
use FileHandle;
use strict;
use XML::Parser::PerlSAX;
use Data::Stag::Util qw(rearrange);
use base qw(Data::Stag::BaseGenerator Exporter);

use vars qw($VERSION);
$VERSION="0.11";

sub fmtstr {
    return 'xml';
}

# OVERRIDE
sub parse {
    my $self = shift;
    my ($file, $str, $fh) = 
      rearrange([qw(file str fh)], @_);
    if ($str) {
	# problem with IO::String in perl5.6.1
	
	my $parser = XML::Parser::PerlSAX->new();
	my $source = {String => $str};
	my %parser_args = (Source => $source,
			   Handler => $self->handler);
	
	$parser->parse(%parser_args);
    }
    else {
	$self->SUPER::parse(-file=>$file, -str=>$str, -fh=>$fh);
    }
}

sub parse_fh {
    my $self = shift;
    my $fh = shift;
    my $parser = XML::Parser::PerlSAX->new();
    my $source = {ByteStream => $fh};
    my %parser_args = (Source => $source,
                       Handler => $self->handler);
    $parser->parse(%parser_args);
}


1;
