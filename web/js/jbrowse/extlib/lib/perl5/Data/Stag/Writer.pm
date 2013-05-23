package Data::Stag::Writer;

=head1 NAME

  Data::Stag::Writer

=head1 SYNOPSIS

  # Abstract class - do not use directly
  package MyOutputter;
  use base qw(Data::Stag::Writer);

  sub e_foo {
    my ($self, $foo) = @_;
    $self->writef("data1: %s\n", $foo->get_data1);
    return;
  }


=cut

=head1 DESCRIPTION

base mixin class for all writers

=head1 INHERITANCE

    This inherits from L<Data::Stag::BaseHandler>

=head1 PUBLIC METHODS -

=head3 new

       Title: new

        Args: [fn str], [fh FILEHANDLE]
      Return: L<Data::Stag::BaseHandler>
     Example: $w = MyWriter->new(-fh=>$fh);

returns the tree that was built from all uncaught events

=head3 file

       Title: file

        Args: filename str
     Returns: filename str
     Example: $handler->file("my_output_file.txt");

Setting this will cause all output to be diverted to this file; the
file will be overwritten by default. The filehandle will not be opened
unless any events are thrown

For more fine-grained control, use $handler->fh()

=head3 fh

       Title: fh

        Args: filehandle FH
     Returns: filehandle FH
     Example: $handler->fh(\*STDOUT);

Gets/Sets the output filehandle for the writer

=head3 safe_fh

       Title: safe_fh
        Type: PROTECTED

        Args: filehandle FH
     Returns: filehandle FH
     Example: $handler->fh(\*STDOUT);

As fh(), but makes sure that the filehandle is initialized

You should use this if you are overriding this class

=head3 write

       Title: write

        Type: PROTECTED
        Args: message str
     Returns: 
     Example: $self->write($stag->get_blah);

writes output

to be used by classes that subclass this one

=head3 writef

       Title: writef

As write, analagous to printf


=cut

use strict;
use base qw(Data::Stag::BaseHandler);
use Data::Stag::Util qw(rearrange);

use vars qw($VERSION);
$VERSION="0.11";

sub init {
    my $self = shift;
    $self->init_writer(@_);
    $self->SUPER::init();
    return;
}

sub init_writer {
    my $self = shift;
    my ($file, $fh) = rearrange([qw(file fh)], @_);
    $fh = \*STDOUT unless $fh || $file;
    $self->fh($fh) if $fh;
    $self->file($file) if $file;
    return;
}

sub file {
    my $self = shift;
    if (@_) {
	$self->{_file} = shift;
	$self->{_fh} = undef;       # undo default setting of fh
    }
    return $self->{_file};
}


sub fh {
    my $self = shift;
    $self->{_fh} = shift if @_;
    return $self->{_fh};
}

sub did_i_open_fh {
    my $self = shift;
    $self->{_did_i_open_fh} = shift if @_;
    return $self->{_did_i_open_fh};
}


sub is_buffered {
    my $self = shift;
    $self->{_is_buffered} = shift if @_;
    return $self->{_is_buffered};
}

sub finish {
    my $self = shift;
    while ($self->depth) {
	$self->end_event;
    }
    if ($self->{_did_i_open_fh}) {
	$self->close_fh;
    }
    return;
}

sub close_fh {
    my $self = shift;
    my $fh = $self->fh;
    if ($fh) {
	$fh->close;
    }
}

sub safe_fh {
    my $self = shift;
    my $fh = $self->{_fh};
    my $file = $self->{_file};

    if ($file && !$fh) {
	$fh =
	  FileHandle->new(">$file") || die("cannot open file: $file");
	$self->fh($fh);
	$self->did_i_open_fh(1);
    }
    $fh;
}

sub addtext {
    my $self = shift;
    my $msg = shift;
    my $fh = $self->safe_fh;
    my $file = $self->file;

    if (!$self->is_buffered && $fh) {
	print $fh $msg;
    }
    else {
	if (!$self->{_buffer}) {
	    $self->{_buffer} = '';
	}
	$self->{_buffer} .= $msg;
    }
    return;
}
*write = \&addtext;

sub writef {
    my $self = shift;
    my $fmtstr = shift;
    $self->addtext(sprintf($fmtstr, @_));
}

sub popbuffer {
    my $self = shift;
    my $b = $self->{_buffer};
    $self->{_buffer} = '';
    return $b;
}


=head2 use_color

  Usage   -
  Returns -
  Args    -

=cut

sub use_color {
    my $self = shift;
    if (@_) {
	$self->{_use_color} = shift;
	if ($self->{_use_color}) {
	    require "Term/ANSIColor.pm";
	}
    }
    return $self->{_use_color};
}

sub DESTROY {
    my $self = shift;
    $self->finish;
}

1;
