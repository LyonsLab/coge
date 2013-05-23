# $Id: XSLHandler.pm,v 1.2 2004/10/27 22:10:44 cmungall Exp $
#
# This GO module is maintained by Chris Mungall <cjm@fruitfly.org>
#
# see also - http://www.geneontology.org
#          - http://www.godatabase.org/dev
#
# You may distribute this module under the same terms as perl itself

=head1 NAME

  Data::Stag::XSLHandler     - 

=head1 SYNOPSIS



=cut

=head1 DESCRIPTION

=head1 PUBLIC METHODS - 

=cut

# makes objects from parser events

package Data::Stag::XSLHandler;
use base qw(Data::Stag::ChainHandler);
use FileHandle;
use XML::LibXML;
use XML::LibXSLT;

use strict;

sub xslt_file {
    my $self = shift;
    if (@_) {
        $self->{_xslt_file} = shift;
        my $fh = FileHandle->new("|xsltproc $self->{_xslt_file} -");
        $self->fh($fh);
    }
    return $self->{_xslt_file};
}

sub init {
    my $self = shift;
    $self->init_writer(@_);
    $self->SUPER::init();
    return;
}

1;
