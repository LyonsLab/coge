#!/usr/bin/perl -w

eval 'exec /usr/bin/perl -w -S $0 ${1+"$@"}'
    if 0; # not running under some shell

# POD docs at end


use strict;

use Data::Stag qw(:all);
use Getopt::Long;


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
my @add = ();
my $lc;
my $merge;
GetOptions(
           "help|h"=>\$help,
           "parser|format|p=s" => \$parser,
           "handler|writer|w=s" => \$handler,
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


my $funcname = shift @ARGV;
my $match = shift @ARGV;
my @files = @ARGV;

my $func;
my $agg;
if ($funcname =~ /^sub /) {
    $func = eval $funcname;
    if ($@) {
	die $@;
    }
}
elsif ($funcname eq 'sum') {
    $func = 
      sub {
	  $agg += $_ foreach @_;
      };
}
elsif ($funcname eq 'avg') {
    $func = 
      sub {
	  $agg += $_ foreach @_;
	  $agg = $agg/@_;
      };
}
elsif ($funcname eq 'cat') {
    $func = 
      sub {
	  $agg .= $_ foreach @_;
      };
}
else {
    die $funcname;
}

my @matches = ();
my $H = Data::Stag->makehandler($match => sub {
				    my $self = shift;
				    my $stag = shift;
				    my $data = $stag->data;
				    push(@matches, $data);
				});
foreach my $fn (@files) {
    Data::Stag->parse(-file=>$fn, -handler=>$H);
    $func->(@matches);
    print "$agg\n";
}
exit 0;

__END__

=head1 NAME 

stag-query.pl - aggregare queries

=head1 SYNOPSIS

  stag-query.pl avg person/age file.xml

  stag-query.pl sum person/salary file.xml

  stag-query.pl 'sub { $agg .= ", ".shift }' person/name file.xml

=head1 DESCRIPTION

Performs aggregate queries

=head1 ARGUMENTS

=cut

