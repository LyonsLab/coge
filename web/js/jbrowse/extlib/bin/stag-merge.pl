#!/usr/bin/perl -w

eval 'exec /usr/bin/perl -w -S $0 ${1+"$@"}'
    if 0; # not running under some shell

=head1 NAME 

stag-merge.pl

=head1 SYNOPSIS

  stag-merge.pl  -xml file1.xml file2.xml

=head1 DESCRIPTION

script wrapper for the Data::Stag modules

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
my @add = ();
my $lc;
GetOptions(
           "help|h"=>\$help,
           "parser|format|p=s" => \$parser,
           "handler|writer|w=s" => \$handler,
           "interpose|i=s@" => \@interpose,
           "add|a=s@" => \@add,
	   "regexp|re|r=s@"=> \@regexps,
           "xml"=>\$toxml,
           "perl"=>\$toperl,
           "lc"=>\$lc,
           "debug"=>\$debug,
          );
if ($help) {
    system("perldoc $0");
    exit 0;
}


my @files = @ARGV;
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

    if ($toxml) {
        print $tree->xml;
    }
    else {
        $tree->generate(-fmt=>$handler);
    }
}

