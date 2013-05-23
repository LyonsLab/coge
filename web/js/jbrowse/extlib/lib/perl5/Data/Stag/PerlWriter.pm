package Data::Stag::PerlWriter;

=head1 NAME

  Data::Stag::PerlWriter

=head1 SYNOPSIS


=cut

=head1 DESCRIPTION


=head1 PUBLIC METHODS -

=cut

use strict;
use base qw(Data::Stag::Writer Data::Stag::Writer);

use vars qw($VERSION);
$VERSION="0.11";

sub fmtstr {
    return 'perl';
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
    my $pre = " ";

    if (($self->this_line &&
        length($self->this_line) + length($o) > 
        60) ||
        $o =~ /^\[/) {
	if ($self->indent_txt) {
	    $pre = "\n" . $self->indent_txt;
	}
	else {
	    #$pre = "'";
	}
        $self->this_line($pre.$o);
    }
    else {
        if ($o =~ /^\]/) {
            $pre = "";
        }
        $self->this_line($self->this_line . $pre.$o);
    }
    $self->addtext( $pre.$o );
}

sub start_event {
    my $self = shift;
    my $ev = shift;
    if (!defined($ev)) {
	$ev = '';
    }
    if ($self->{__opened}) {
	$self->o(',');
    }
    else {
	$self->{__opened} = 1;
	$self->o('[');
    }
    my $stack = $self->stack;
    if ($self->use_color) {
	$self->o(color('white'));
	$self->o('['.color('red').$ev);
    }
    else {
	$self->o("['$ev'=>");
    }
    push(@$stack, $ev);
}
sub end_event {
    my $self = shift;
    my $ev = shift;
    my $stack = $self->stack;
    my $popped = pop(@$stack);
    if ($ev && $popped ne $ev) {
        warn("uh oh; $ev ne $popped");
    }
    if ($self->use_color) {
	$self->o(']');
    }
    else {
	$self->o(']');
    }
    if (!@$stack) {
	$self->o("\n");
    }
    return $ev;
}
sub evbody {
    my $self = shift;
    my $body = shift;
    my $str;
    if ($self->use_color) {
	if (!defined($body)) {
	    $str = color('white').'""';
	}
	elsif ($body eq '0') {
	    $str = color('white').'"'.color('green').'0'.color('white').'"';
	}
	else {
	    $body =~ s/\[/\\\[/g;
	    $body =~ s/\]/\\\]/g;
	    $body =~ s/\"/\\\"/g;
	    $str = color('white').'"'.color('green').$body.color('white').'"';
	}
    }
    else {
	$str = perlesc($body);
    }
    $self->o($str);
    return;
}

sub perlesc {
    my $w = shift;
    return '""' unless defined $w;
    return '"0"' if $w eq '0';
    $w =~ s/\"/\\\"/g;
    return '"'.$w.'"';
}

sub color {
    Term::ANSIColor::color(@_);
}


1;
