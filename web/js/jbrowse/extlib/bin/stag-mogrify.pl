#!/usr/bin/perl -w

eval 'exec /usr/bin/perl -w -S $0 ${1+"$@"}'
    if 0; # not running under some shell

=head1 NAME 

stag-mogrify.pl  - mangle stag files

=head1 SYNOPSIS

  stag-mogrify.pl  -w itext file1.xml file2.xml

=head1 DESCRIPTION

script wrapper for the Data::Stag modules

feeds in files into a parser object that generates nestarray events,
and feeds the events into a handler/writer class

=head1 ARGUMENTS

=cut



use strict;

use Carp;
use Data::Stag qw(:all);
use Getopt::Long;

my $parser = "";
my $handler = "";
my $mapf;
my $tosql;
my $toxml;
my $toperl;
my $debug;
my $help;
my @interpose = ();
my @regexps = ();
my @delete = ();
my @add = ();
my $lc;
my $merge;
GetOptions(
           "help|h"=>\$help,
           "parser|format|p=s" => \$parser,
           "handler|writer|w=s" => \$handler,
           "interpose|i=s@" => \@interpose,
           "add|a=s@" => \@add,
           "delete|d=s@" => \@delete,
	   "regexp|re|r=s@"=> \@regexps,
           "xml"=>\$toxml,
           "perl"=>\$toperl,
           "lc"=>\$lc,
           "debug"=>\$debug,
	   "merge=s"=>\$merge,
          );
if ($help) {
    system("perldoc $0");
    exit 0;
}


my @files = @ARGV;
my @mergeset = ();
foreach my $fn (@files) {

    my $tree = 
      Data::Stag->parse($fn, 
			$parser);

    if ($lc) {
	$tree->iterate(sub {
			   my $s = shift;
			   my $e = $s->element;
			   my $orig = $e;
			   $e = lc($e);
			   $s->element($e) unless $e eq $orig;
			   return
		       });
    }
    foreach my $regexp (@regexps) {
	$tree->iterate(sub {
			   my $s = shift;
			   my $e = $s->element;
			   my $orig = $e;
			   eval "\$e =~ $regexp";
			   $s->element($e) unless $e eq $orig;
			   return
		       });
    }
    foreach (@delete) {
	$tree->remove($_);
    }
    foreach (@add) {
	$tree->where($_,
		     sub {
			 my $terminal = shift;
			 my $el = $terminal->element;
			 $terminal->data([Data::Stag->new($el=>$terminal->data)]);
			 $terminal->element($el . "_data");
		     });
    }
    foreach (@interpose) {
	my ($from, $to, $inter) = split(/\,/, $_);
	if (!$to) {
	    ($from, $to) = ('', $from);
	}
	my @adds = ();
	$tree->iterate(
		       sub {
			   my $outer = shift;
			   return if $outer->isterminal;
			   if ($from && $outer->element ne $from) {
			       return;
			   }
			   my @inner = $outer->get($to);
			   return unless @inner;
			   $outer->unset($to);
			   foreach (@inner) {
			       my $ie = $inter;
			       if (!$ie) {
				   $ie = $outer->element .'_'. $to;
			       }
			       my $x = Data::Stag->new($ie=>[$_]);
#			       $outer->add($inter, $x);
			       push(@adds, [$outer, $ie, $x]);
#			       print STDERR $outer->xml;
			   }
		       });
	
	foreach (@adds) {
	    my ($outer, $inter, $x) = @$_;
	    $outer->add($inter, $x);
	}
    }

    if ($merge) {
	push(@mergeset, $tree);
    }
    else {
	if ($toxml) {
	    print $tree->xml;
	}
	else {
	    my $W = Data::Stag->getformathandler($handler);
	    $W->fh(\*STDOUT);
	    $tree->sax($W);
	}
    }
}
if ($merge) {
    my $tree = Data::Stag->new($merge=>[@mergeset]);
	if ($toxml) {
	    print $tree->xml;
	}
	else {
	    $tree->generate(-fmt=>$handler);
	}
}
