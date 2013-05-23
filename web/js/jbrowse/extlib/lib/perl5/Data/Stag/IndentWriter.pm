package Data::Stag::IndentWriter;

=head1 NAME

  Data::Stag::IndentWriter

=head1 SYNOPSIS


=cut

=head1 DESCRIPTION

=head1 PUBLIC METHODS -

=cut

use strict;
use base qw(Data::Stag::Writer);
use Carp;

use vars qw($VERSION);
$VERSION="0.11";

sub fmtstr {
    return 'indent';
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
    $self->addtext("@_");
}

sub unclosed {
    my $self = shift;
    $self->{_unclosed} = shift if @_;
    return $self->{_unclosed};
}

sub in_attr {
    my $self = shift;
    $self->{_in_attr} = shift if @_;
    return $self->{_in_attr};
}

sub start_event {
    my $self = shift;
    my $ev = shift;
    my $stack = $self->stack;
    
    if ($ev eq '@') {
        $self->o('[');
	$self->in_attr(1);
    }
    elsif ($ev eq '.') {
	# pcdata for element with attributes
	# do nothing
    }
    elsif ($self->in_attr) {
        if ($self->in_attr > 1) {
            $self->o(' ');
        }
	$self->o("$ev=");
        $self->in_attr(2);
    }
    elsif ($ev eq '') {
        # ignore null nodes
    }
    else {
        my $tag = "$ev";
        if ($self->use_color) {
            $tag = color('red').$ev.color('reset');
        }
        my $nl = "\n";
        if (@$stack == 1) {
            $nl = "\n\n";
        }
        $self->addtext($nl . $self->indent_txt() . $tag);
        $self->unclosed(1);
    }
    push(@$stack, $ev);
}

sub end_event {
    my $self = shift;
    my $ev = shift;
    my $stack = $self->stack;
    my $top = pop @$stack;
    $top eq $ev or confess("$top ne $ev");
    if ($self->in_attr) {
	if ($ev eq '@') {
	    $self->in_attr(0);
            $self->o(']');
	}
    }
    elsif ($ev eq '.') {
	# end of pcdata for element with attributes
    }
    else {
	if (!@$stack) {
	    $self->o("\n");
	}
    }
    return $ev;
}
sub evbody {
    my $self = shift;
    my $body = shift;
    my $str = itextesc($body);
    if ($self->in_attr) {
        $self->addtext($str);
    }
    else {
        if ($self->use_color) {
            $str = color('white').$str;
        }    
        $self->addtext(' "'.$str.'"');
    }
    return;
}

sub itextesc {
    my $w = shift;
    if (!defined($w)) {
	$w='';
    }
    $w =~ s/\"/\\\"/g;
    return $w;
}

# this should already be require'd in Writer.pm
sub color {
    Term::ANSIColor::color(@_);
}


1;
