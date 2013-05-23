package Data::Stag::SAX2Stag;

=head1 NAME

  Data::Stag::SAX2Stag

=head1 SYNOPSIS

takes in SAX events and turns them into Stag events

attributes are turned into elements

=cut

=head1 DESCRIPTION

modules for dealing with xml as nested arrays

=head1 PUBLIC METHODS -

=cut

use strict;
use base qw(Data::Stag::Base);

use vars qw($VERSION);
$VERSION="0.11";

my (%mail_args, $current_element, $message_count, $sent_count);

DEAD
DEAD
DEAD
DEAD
DEAD
DEAD
DEAD
DEAD
;
sub start_element {
    my ($self, $element) = @_;

    my $name = $element->{Name};
    my $atts = $element->{Attributes};
    foreach my $k (keys %$atts) {
        $self->event($k, $atts->{$k});
    }
    $self->start_event($name);
    $self->{Handler}->start_element($element);
    
}

sub characters {
    my ($self, $characters) = @_;
    my $char = $characters->{Data};
    my $str = $self->{__str};
    if ($char) {
        $str = "" if !defined $str;
        $str .= $char;
    }
    $self->{__str} = $str;
    $self->{Handler}->characters($characters);
}

sub end_element {
    my ($self, $element) = @_;
    my $name = $element->{Name};
    my $str = $self->{__str};
    if (defined $str) {
        $str =~ s/^\s*//;
        $str =~ s/\s*$//;
        $self->evbody($str) if $str;
    }
    $self->end_event($name);
    $self->{__str} = undef;
    $self->{Handler}->end_element($element);
}


1;
