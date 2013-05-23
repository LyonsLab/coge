package Data::Stag::XMLWriter;

=head1 NAME

  Data::Stag::XMLWriter

=head1 SYNOPSIS


=cut

=head1 DESCRIPTION

=head1 PUBLIC METHODS -

=cut

use strict;
use base qw(Data::Stag::Writer);
use Data::Stag::Util qw(rearrange);

use vars qw($VERSION);
$VERSION="0.11";

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


sub ensure_closed {
    my $self = shift;
    if ($self->{_unclosed}) {
	$self->o(">") ;
	$self->{_unclosed} = 0;
    }
    return;
}

sub fmtstr {
    return 'xml';
}

sub indent_txt {
    my $self = shift;
    my $stack = $self->stack;
    return "  " x scalar(@$stack);
}

sub this_line {
    my $self = shift;
    $self->{_this_line} = shift if @_;
    return $self->{_this_line};
}

sub o {
    my $self = shift;
    my $o = "@_";
    $self->addtext( $o );
}

sub first_line {
    my $self = shift;
    $self->{_first_line} = shift if @_;
    return $self->{_first_line};
}


sub start_event {
    my $self = shift;
    my $ev = shift;
    if (!defined($ev)) {
	$ev = '';
    }
    my $stack = $self->stack;
    if (!@$stack) {
	$self->o($self->first_line || "<?xml version=\"1.0\" encoding=\"UTF-8\"?>");
    }

    if ($ev eq '@') {
	$self->in_attr(1);
    }
    elsif ($ev eq '.') {
	# pcdata for element with attributes
	# do nothing
    }
    elsif ($self->in_attr) {
	$self->o(" $ev=");
    }
    elsif ($ev eq '') {
        # ignore null nodes
    }
    else {
	$self->ensure_closed;
	$self->o("\n". $self->indent_txt . "<$ev");
	$self->unclosed(1);
    }
    push(@$stack, $ev);
}

sub end_event {
    my $self = shift;
    my $ev = shift;
    my $stack = $self->stack;
    my $popped = pop(@$stack);
    return '' unless $popped;
    if ($ev && $popped ne $ev) {
        warn("uh oh; $ev ne $popped");
    }
    if (!$ev) {
	$ev = $popped;
    }
    if ($self->in_attr) {
	if ($ev eq '@') {
	    $self->in_attr(0);
	    $self->ensure_closed;
	}
    }
    elsif ($ev eq '.') {
	# end of pcdata for element with attributes
    }
    else {
	$self->ensure_closed;
	if ($self->{_nl}) {
	    $self->o("\n" . $self->indent_txt)
	}
	$self->o("</$ev>");
	$self->{_nl} = 1;
	if (!@$stack) {
	    $self->o("\n");
	}
    }
    return $ev;
}
sub evbody {
    my $self = shift;
    my $body = shift;
    my $str= xmlesc($body);
    if ($self->in_attr) {
	$self->o("\"$str\"");
    }
    else {
	$self->ensure_closed;
	$self->{_nl} = 0;
	$self->o($str);
    }
    return;
}

our $escapes = { '&' => '&amp;',
		 '<' => '&lt;',
		 '>' => '&gt;',
		 '"' => '&quot;'
	       };

sub xmlesc {
    my $w = shift;
    if (!defined $w) {
	$w = '';
    }
    $w =~ s/([\&\<\>])/$escapes->{$1}/ge;
    $w;
}



1;
