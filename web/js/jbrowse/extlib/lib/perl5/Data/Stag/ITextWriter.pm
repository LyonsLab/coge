package Data::Stag::ITextWriter;

=head1 NAME

  Data::Stag::ITextWriter

=head1 SYNOPSIS


=cut

=head1 DESCRIPTION

=head1 PUBLIC METHODS -

=cut

use strict;
use base qw(Data::Stag::Writer);

use vars qw($VERSION);
$VERSION="0.11";

sub fmtstr {
    return 'itext';
}


sub stack {
    my $self = shift;
    $self->{_stack} = shift if @_;
    return $self->{_stack};
}

sub indent_txt {
    my $self = shift;
    my $stack = $self->stack;
    return "  " x scalar(@$stack);
}

sub o {
    my $self = shift;
    my @o = @_;
    $self->addtext($self->indent_txt() . "@o");
}

sub start_event {
    my $self = shift;
    my $ev = shift || '';
    my $stack = $self->stack;
    my $tag = "$ev: ";
    if ($self->use_color) {
	$tag = color('red').$ev.color('reset').':';
    }
    $self->addtext("\n" . $self->indent_txt() . $tag);
    push(@$stack, $ev);
}
sub end_event {
    my $self = shift;
    my $ev = shift || '';
    my $stack = $self->stack;
    my $top = pop @$stack;
    use Carp;
    $top eq $ev or confess("$top ne $ev");
    return $ev;
}
sub evbody {
    my $self = shift;
    my $body = shift;
    my $str = itextesc($body);
    if ($self->use_color) {
	$str = color('white').$str;
    }    
    $self->addtext($str);
    return;
}

sub itextesc {
    my $w = shift;
    if (!defined($w)) {
	$w='';
    }
    $w =~ s/:/\\:/g;
    return $w;
}

# this should already be require'd in Writer.pm
sub color {
    Term::ANSIColor::color(@_);
}


1;
