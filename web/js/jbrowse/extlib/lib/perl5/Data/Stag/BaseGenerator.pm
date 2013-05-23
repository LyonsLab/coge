# $Id: BaseGenerator.pm,v 1.15 2004/12/21 02:26:25 cmungall Exp $
#
# Copyright (C) 2002 Chris Mungall <cjm@fruitfly.org>
#
# See also - http://stag.sourceforge.net
#
# This module is free software.
# You may distribute this module under the same terms as perl itself

package Data::Stag::BaseGenerator;

=head1 NAME

  Data::Stag::BaseGenerator     - base class for parsers and other event generators

=head1 SYNOPSIS

  # writing the parser
  package MyParser;
  use base qw(Data::Stag::BaseGenerator);
  
  sub parse_fh {
    my ($self, $fh) = shift;

    my $lnum = 0;
    $self->start_event('data');
    while (<$fh>) {
      ++$lnum;
      $self->line_no($lnum);
      # do stuff
      $self->start_event('foo');

      # ...
      $self->event(blah=>5);

      #
      if (/incorrect_line/) {
         $self->parse_err('line not in correct format');
      }

      # ...
      $self->end_event('foo');
    }
    $self->pop_stack_to_depth(0);
  }
  1;

  # using the parser
  my $p = MyParser->new;
  my $h = MyHandler->new; # see Data::Stag::BaseHandler
  my $eh = Data::Stag->makehandler;
  $p->handler($h);
  $p->errhandler($eh);
  $p->parse($file);

  # result tree
  print $h->stag->xml;

  # write parse errs on standard err
  printf \*STDERR $p->errhandler->stag->xml;

  # using the parser from the command line
  unix> stag-parse.pl -p MyParser -w xml -e err.xml > out.xml

  # using the parser from the command line via intermediate handler
  unix> stag-handle.pl -p MyParser -m MyHandler -w xml -e err.xml > out.xml

=cut

=head1 DESCRIPTION

This is the base class for all parsers and event generators

parsers/generators take some input (usually a filehandle, but a
generator could be a socket listener, for example) and fire stag
events

stag events are

=over

=item start_event NODENAME

=item evbody DATA

=item end_event NODENAME {optional}

=item event NODENAME DATA

=back

These events can be nested/hierarchical

If uncaught, these events are stacked into a stag tree, which can be
written as xml or one of the other stag formats

specialised handlers can be written to catch the events your parser
throws

For example, you may wish to write a pod parser that generates nested
events like this:

  <pod>
   <section>
     <type>head1</type>
     <name>NAME</name>
     <text>Data::Stag - Structured Tags datastructures</text>
   </section>
   ...
  </pod>

(see the source for Data::Stag::PodParser for details)

You can write handlers that take the pod-xml and generate something -
for example HTML

parsers may encounter unexpected things along the way - they may throw
an exception, and fall over - or they may choose to fire an error
event. by default, error event streams are diverted to STDERR. You can
create your own error handlers

=head1 PUBLIC METHODS

=head3 new

       Title: new

        Args: 
      Return: L<Data::Stag::BaseGenerator>
     Example: 

CONSTRUCTOR

=head3 handler

       Title: handler
    Function: GET/SET ACCESSOR METHOD
        Args: handler L<Data::Stag::BaseHandler> optional
      Return: L<Data::Stag::BaseHandler>
     Example: $p->handler(MyHandler->new);

each parser has a handler - all events generated are passed onto the
handler; the default handler simply sits there collecting events

=head3 errhandler

       Title: errhandler
    Function: GET/SET ACCESSOR METHOD
        Args: handler L<Data::Stag::BaseHandler> optional
      Return: L<Data::Stag::BaseHandler>
     Example: $p->errhandler(Data::Stag->makehandler);

each parser has an error handler - if the parser encounters things it
does not expect, it can pass errors to the errorhandler

if no errorhandler is set, an XML event handler that writes to STDERR is used

=head3 cache_errors

       Title: cache_errors
        Args: 
      Return: 
     Example: $p->cache_errors

If this is called, all errors will be cached rather than written to STDERR

The error list can be accessed like this

  $p->parse($fn);
  @errs = $p->errhandler->stag->get_error;

=head2 parse

  Example - $parser->parse($file1, $file2);
  Returns - 
  Args    - filenames str-LIST

parses a file

=head2 parse

  Example - $parser->parse_fh($fh)
  Returns - 
  Args    - fh FILEHANDLE

parses an open filehandle

=cut

=head1 PROTECTED METHODS

These methods are only of interest if you are making your own
parser/generator class

=over

=item start_event NODENAME

=item evbody DATA

=item end_event NODENAME {optional}

=item event NODENAME DATA

=back

=head1 SEE ALSO

L<Data::Stag>
L<Data::Stag::BaseHandler>

=cut

use Exporter;
@ISA = qw(Exporter);

use Carp;
use FileHandle;
use Data::Stag::Util qw(rearrange);
use Data::Stag::null;
use strict qw(subs vars refs);

# Exceptions

sub throw {
    my $self = shift;
    confess("@_");
}

sub warn {
    my $self = shift;
    warn("@_");
}

#sub last_evcall_type {
#    my $self = shift;
#    $self->{_last_evcall_type} = shift if @_;
#    return $self->{_last_evcall_type};
#}


sub stack {
    my $self = shift;
    $self->{_stack} = shift if @_;
    $self->{_stack} = [] unless $self->{_stack};
    return $self->{_stack};
}

sub stack_top {
    my $self = shift;
    $self->stack->[-1];
}

sub push_stack {
    my $self = shift;
    push(@{$self->stack}, @_);
}

sub pop_stack {
    my $self = shift;
    my $top = $self->stack_top;
    $self->end_event($top);
    $top;
    
}

sub pop_stack_to_depth {
    my $self = shift;
    my $depth = shift;
    while ($depth < $self->stack_depth) {
        $self->end_event;
    }
}

sub pop_all {
    my $self = shift;
    $self->pop_stack_to_depth(0);
}

sub stack_depth {
    my $self = shift;
    scalar(@{$self->stack});
}

*error_list = \&messages;

sub message {
    my $self = shift;
    my $msg = shift;
    unless (ref($msg)) {
        $msg =
          {msg=>$msg,
           line=>$self->line,
           line_no=>$self->line_no,
           file=>$self->file};
    }
    push(@{$self->messages},
         $msg);
}



sub new {
    my ($class, $init_h) = @_;
    my $self = {};
    $self->{handler} = Data::Stag::null->new;
    if ($init_h) {
	map {$self->{$_} = $init_h->{$_}} keys %$init_h;
    }
    bless $self, $class;
    $self->init if $self->can("init");
    $self;
}

sub load_module {

    my $self = shift;
    my $classname = shift;
    my $mod = $classname;
    $mod =~ s/::/\//g;

    if ($main::{"_<$mod.pm"}) {
    }
    else {
        require "$mod.pm";
    }
}

sub modulemap {
    my $self = shift;
    $self->{_modulemap} = shift if @_;
    return $self->{_modulemap};
}

sub handler {
    my $self = shift;
    if (@_) {
        my $h = shift;
        if ($h && !ref($h)) {
            my $base = "Data::Stag:";
	    my $mm = $self->modulemap;
	    if ($mm && $mm->{$h}) {
		$h = $mm->{$h};
	    }
            $h =~ s/^xml$/$base:XMLWriter/;
            $h =~ s/^perl$/$base:PerlWriter/;
            $h =~ s/^sxpr$/$base:SxprWriter/;
            $h =~ s/^itext$/$base:ITextWriter/;
            $h =~ s/^graph$/$base:GraphWriter/;
            $self->load_module($h);
            $h = $h->new;
	    if ($h->can("fh")) {
		$h->fh(\*STDOUT);
	    }
        }
        $self->{_handler} = $h;
        if (!$h->errhandler && $self->errhandler) {
            $h->errhandler($self->errhandler);
        }
    }
#    return $self->{_handler} || Data::Stag::null->new();
    return $self->{_handler};
}

sub cache_errors {
    my $self = shift;
    return $self->errhandler(Data::Stag->makehandler);
}

sub errhandler {
    my $self = shift;
    if (@_) {
        $self->{errhandler} = shift;
    }
    return $self->{errhandler};
}

sub err_event {
    my $self = shift;
    if (!$self->errhandler) {
	$self->errhandler(Data::Stag->getformathandler('xml'));
	$self->errhandler->fh(\*STDERR);
	
#	my $estag = Data::Stag->new(@_);
#	eval {
#	    confess;
#	};
#	$estag->set_stacktrace($@);
#	print STDERR $estag->xml;
#	exit 1;
    }
    if (!$self->errhandler->depth) {
	$self->errhandler->start_event("error_eventset");
    }
    $self->errhandler->event(@_);
    return;
}

sub err {
    my $self = shift;
    my $err = shift;
    if (ref($err)) {
	$self->throw("Bad error msg $err - must not by ref");
    }
    $self->err_event(message=>$err);
    return;
}

sub parse_err {
    my $self = shift;
    my $err = shift || '';
    if (ref($err)) {
	$self->throw("Bad error msg $err - must not by ref");
    }
    my @tags = ([message=>$err],[file=>$self->file]);
    my $line = $self->line;
    push(@tags, [line=>$line]) if defined $line;
    my $line_no = $self->line_no;
    push(@tags, [line_no=>$line_no]) if $line_no;
    my $pclass = ref($self);
    push(@tags, [parse_class=>"$pclass"]);
    $self->err_event(error=>[@tags]);
    return;
}


sub line_no {
    my $self = shift;
    $self->{_line_no} = shift if @_;
    return $self->{_line_no};
}

sub line {
    my $self = shift;
    $self->{_line} = shift if @_;
    return $self->{_line};
}
sub file {
    my $self = shift;
    $self->{_file} = shift if @_;
    return $self->{_file};
}

sub finish {
    my $self = shift;
    if ($self->errhandler && $self->errhandler->depth) {
	$self->errhandler->end_event;
    }
}


sub parse {
    my $self = shift;
    my ($file, $str, $fh) = 
      rearrange([qw(file str fh)], @_);

    $self->file($file) if $file;
    if ($str) {
        $self->load_module("IO::String");
        $fh = IO::String->new($str) || confess($str);
    }
    elsif ($file) {
        if ($file eq '-') {
            $fh = \*STDIN;
        }
        else {
            $self->load_module("FileHandle");
            $fh = FileHandle->new($file) || confess("cannot open file: $file");
        }
    }
    else {
    }
    if (!$fh) {
        confess("no filehandle");
    }
    $self->parse_fh($fh);
    # problem with IO::String closing in perl5.6.1
    unless ($str) {
        #$fh->close  || confess("cannot close file: $file");  // problem stdout
        $fh->close;
    }
    return;
}

sub handler_err {
    my $self = shift;
    $self->err_event(error=>[[message=>'handler problem'],
			     [stack=>shift]]);
}

sub errlist {
    my $self = shift;
    $self->finish;
    my $eh = $self->errhandler;
    if ($eh && $eh->stag && $eh->stag->data) {
	return ($eh->stag->get_error);
    }
    return ();
}


# MAIN EVENT HANDLING
# start/end/event/body

sub start_event {
    my $self = shift;
    push(@{$self->{_stack}}, $_[0]);
    $self->{_handler}->start_event(@_);
    return;
}
sub end_event { 
    my $self = shift; 
    my $ev = shift || $self->{_stack}->[-1];
    $self->{_handler}->end_event($ev);

    my $out = pop(@{$self->{_stack}});
    return;
}

sub event {
    my $self = shift;
    $self->{_handler}->event(@_);
    return;
}
sub evbody {
    my $self = shift;
    $self->{_handler}->evbody(@_);
#    my $lc = $self->last_evcall_type;
#    if ($lc && $lc eq 'evbody') {
#        confess("attempting to event_body illegally (body already defined)");
#    }

#    eval {
#	$self->handler->evbody(@_);
#    };
#    if ($@) {
#	$self->handler_err($@);
#    }
#    $self->last_evcall_type('evbody');
    return;
}

1;
