#!/usr/bin/perl -w

eval 'exec /usr/bin/perl -w -S $0 ${1+"$@"}'
    if 0; # not running under some shell


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


my $e = shift @ARGV;
my @files = @ARGV;

foreach my $fn (@files) {
    my $tree = Data::Stag->parse($fn);
    my @subtrees = $tree->find($e);
    print $_->xml foreach  (@subtrees);
}

__END__

=head1 NAME 

stag-findsubtree.pl - finds nodes in a stag file

=head1 SYNOPSIS

  stag-findsubtree.pl 'person/name' file.xml

=head1 DESCRIPTION

parses in an input file and writes out subnodes

=head1 USAGE

  stag-findsubtree.pl [-p PARSER] [-w WRITER] NODE FILE

=head1 ARGUMENTS

=over

=item -p|parser FORMAT

FORMAT is one of xml, sxpr or itext, or the name of a perl module

xml assumed as default

=item -w|writer FORMAT

FORMAT is one of xml, sxpr or itext, or the name of a perl module

=item NODE

the name of the node/element 

=back

=head1 LIMITATIONS

not event based

=head1 SEE ALSO

L<Data::Stag>

=cut

