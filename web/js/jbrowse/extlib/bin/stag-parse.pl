#!/usr/bin/perl -w

eval 'exec /usr/bin/perl -w -S $0 ${1+"$@"}'
    if 0; # not running under some shell

# POD docs at end

use strict;

use Carp;
use Data::Stag qw(:all);
use Getopt::Long;

my $parser = "";
my $handler = "";
my $writer = "";
my $errhandler = "";
my $errf;
my $outf;
my $mapf;
my $tosql;
my $toxml;
my $toperl;
my $debug;
my $help;
my $color;
my $root_node;

GetOptions(
           "help|h"=>\$help,
	   "root|r=s"=>\$root_node,
           "parser|format|p=s" => \$parser,
           "handler|writer|w=s" => \$writer,
           "errhandler=s" => \$errhandler,
           "errf|e=s" => \$errf,
           "out|o=s" => \$outf,
           "debug"=>\$debug,
           "colour|color"=>\$color,
          );
if ($help) {
    system("perldoc $0");
    exit 0;
}

$errhandler =  Data::Stag->getformathandler($errhandler || 'xml');
$handler = Data::Stag->getformathandler($writer);
if ($handler->can("fh")) {
    if ($outf) {
	$handler->file($outf);
    }
    else {
	$handler->fh(\*STDOUT);
    }
}

if ($color) {
    $handler->use_color(1);
}
if ($errf) {
    $errhandler->file($errf);
}
else {
    $errhandler->fh(\*STDERR);
}

my @files = @ARGV;
if (@files > 1) {
    $root_node = 'set' unless $root_node;
}
$handler->start_event($root_node) if $root_node;
foreach my $fn (@files) {

    my @pargs = (-file=>$fn, -format=>$parser, -handler=>$handler, -errhandler=>$errhandler);
    if ($fn eq '-') {
	if (!$parser) {
	    $parser = 'xml';
	}
	@pargs = (-format=>$parser, -handler=>$handler, 
		  -fh=>\*STDIN, -errhandler=>$errhandler);
    }
    my $tree = 
      Data::Stag->parse(@pargs);

    if ($toxml) {
        print $tree->xml;
    }
    if ($toperl) {
        print tree2perldump($tree);
    }
}
if ($errf) {
    $errhandler->finish;
}
$handler->end_event($root_node) if $root_node;
$handler->fh->close if $outf;
exit 0;

__END__

=head1 NAME 

stag-parse.pl - parses a file and fires events (e.g. sxpr to xml)

=head1 SYNOPSIS

  # convert XML to IText
  stag-parse.pl -p xml -w itext file1.xml file2.xml

  # use a custom parser/generator and a custom writer/generator
  stag-parse.pl -p MyMod::MyParser -w MyMod::MyWriter file.txt

=head1 DESCRIPTION

script wrapper for the Data::Stag modules

feeds in files into a parser object that generates nestarray events,
and feeds the events into a handler/writer class

=head1 ARGUMENTS

=over

=item -p|parser FORMAT

FORMAT is one of xml, sxpr or itext, or the name of a perl module

this is the class that parsers the input file(s) and generates stag
events

xml assumed as default

=item -w|writer FORMAT

FORMAT is one of xml, sxpr or itext, or the name of a perl module

this is the class that catches the events thrown by the parser; it can
be any class, but the class is typically a writer

xml assumed as default

=item -o|out FILE

the writer will use this file (defaults to STDOUT)

=item -e|errf FILE

file to store parse error handler output

=item -errhandler FORMAT/MODULE

FORMAT is one of xml, sxpr or itext, or the name of a perl module

all parse error events go to this module

=item 

=item -r|root NODE_NAME

if this is specified, NODE_NAME becomes the root of the stag tree, and
anything that was previously the root is placed below this.

this happens automatically if more than one file is parsed (because
there can only be one tree root)

=item -color

Works only if the output handler is able to provide ASCII-colors
(currently supported for itext and xml)

=back


=head1 SEE ALSO

L<Data::Stag>

This script is a wrapper for the method

  Data::Stag->parse()

=cut

