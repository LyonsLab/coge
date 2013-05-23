#!/usr/bin/perl -w

eval 'exec /usr/bin/perl -w -S $0 ${1+"$@"}'
    if 0; # not running under some shell

# POD docs at end

use strict;

use Data::Stag qw(:all);
use Getopt::Long;

use FileHandle;

my $exec;
my $codefile;
my $parser = '';
my $errhandler = "";
my $errf;
my $writer = '';
my %trap = ();
my $datafmt;
my $module;
my @units;
GetOptions("codefile|c=s"=>\$codefile,
	   "sub|s=s"=>\$exec,
	   "trap|t=s%"=>\%trap,
           "parser|format|p=s" => \$parser,
           "errhandler=s" => \$errhandler,
           "errf|e=s" => \$errf,
	   "writer|w=s"=>\$writer,
	   "data|d=s"=>\$datafmt,
	   "module|m=s"=>\$module,
	   "units|u=s@"=>\@units,
	   "help|h"=>sub {system("perldoc $0");exit 0},
	  );

if (!$codefile && !$exec && !$module) {
    $codefile = shift @ARGV if (@ARGV > 1);
    die "you must supply -c or -m or -s, or provide codefile" 
      unless $codefile;
}
my @files = @ARGV;

$errhandler =  Data::Stag->getformathandler($errhandler || 'xml');
if ($errf) {
    $errhandler->file($errf);
}
else {
    $errhandler->fh(\*STDERR);
}
my $catch = {};
no strict;
if ($exec) {
    $catch = eval $exec;
    if ($@) {
        die $@;
    }
    if (ref($catch) ne 'HASH') {
        print STDERR "exec \"$exec\" is not a hash";
        exit 1;
    }
}
if ($codefile) {
    $catch = do "$codefile";
    if ($@) {
        print STDERR "\n\nstag-handle error:\n";
        print STDERR "There was an error with the codefile \"$codefile\":\n\n";
        die $@;
    }
    if (ref($catch) ne 'HASH') {
        print STDERR "codefile \"$codefile\" does not return a hash";
        exit 1;
    }
}
if (%trap) {
    #    
    #    die Dumper \%trap;
    %$catch = (%$catch, %trap);
}
use strict;
my @events;
my $inner_handler;
if ($module) {
    $inner_handler = Data::Stag->makehandler($module);
} else {
    my $meth = $exec ?  $exec : $codefile;
    if (!%$catch) {
        die "method \"$meth\" did not return handler";
    }
    if (!ref($catch) || ref($catch) ne "HASH") {
        die("$meth must return hashref");
    }
    $inner_handler = Data::Stag->makehandler(%$catch);
    @events = keys %$catch;
}
if (@units) {
    @events = @units;
}
$inner_handler->errhandler($errhandler);
my $h = Data::Stag->chainhandlers([@events],
                                  $inner_handler,
                                  $writer);

while (my $fn = shift @files) {
    my $fh;
    if (!$fn || $fn eq '-') {
        $fh = \*STDIN;
        $fn = '';
    }
    else {
        $fh = FileHandle->new($fn) || die "Cannot open file: $fn";
    }
    
    my $p = Data::Stag->parser(-file=>$fn, -format=>$parser, 
                               -errhandler=>$errhandler);

    $p->handler($h);
    $p->parse_fh($fh);
    if ($datafmt) {
        print $inner_handler->stag->$datafmt();
    }
}

__END__

=head1 NAME

stag-handle.pl - streams a stag file through a handler into a writer

=head1 SYNOPSIS

  stag-handle.pl -w itext -c my-handler.pl myfile.xml > processed.itext
  stag-handle.pl -w itext -p My::Parser -m My::Handler myfile.xml > processed.itext

=head1 DESCRIPTION

will take a Stag compatible format (xml, sxpr or itext), turn the data
into an event stream passing it through my-handler.pl

=over ARGUMENTS

=item -help|h

shows this document

=item -module|m PERLMODULE

A module that is used to transform the input events
the module should inherit from L<Data::Stag::BaseHandler>

=item -unit|u NODE_NAME

(you should always use this option if you specify -m)

this is the unit that gets passed to the handler/transformer. this
will get set automatically if you use the the -c, -s or -t options

multiple units can be set

  -u foo -u bar -u boz

=item -writer|w WRITER

writer for final transformed tree; can be xml, sxpr or itext

=item -module|m MODULE

perl modules for handling events

=item -codefile|c FILE

a file containing a perlhashref containing event handlers - see below

=item -sub|s PERL

a perl hashref containing handlers 

=item -trap|t ELEMENT=SUB

=back



=head1 EXAMPLES

  unix> cat my-handler.pl
  {
    person => sub {
	my ($self, $person) = @_;
	$person->set_fullname($person->get_firstname . ' ' .
			   $person->get_lastname);
        $person;
    },
    address => sub {
	my ($self, $address) = @_;
	# remove addresses altogether from processed file
        return;
    },
  }


=cut
