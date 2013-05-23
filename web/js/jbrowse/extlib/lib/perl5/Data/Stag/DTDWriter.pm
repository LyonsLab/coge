package Data::Stag::DTDWriter;

=head1 NAME

  Data::Stag::DTDWriter

=head1 SYNOPSIS

  my $p = Data::Stag->parser;
  $p->handler('dtd');
  $p->parse($file);
  

=cut

=head1 DESCRIPTION


=head1 PUBLIC METHODS -

=cut

use strict;
use base qw(Data::Stag::Writer);

use vars qw($VERSION);
$VERSION="0.11";

sub end_stag {
    my $self = shift;
    my $stag = shift;
    my $dtd = $stag->autoschema->dtd;
    $self->addtext($dtd);
    return;
}

1;
