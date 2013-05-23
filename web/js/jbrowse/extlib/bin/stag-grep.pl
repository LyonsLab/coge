#!/usr/bin/perl -w

eval 'exec /usr/bin/perl -w -S $0 ${1+"$@"}'
    if 0; # not running under some shell

# POD docs at bottom of file

use strict;

use Carp;
use Data::Stag qw(:all);
use Getopt::Long;

my $fmt = "";
my $out = "";
my $mapf;
my $tosql;
my $toxml;
my $toperl;
my $debug;
my $help;
my $count;
my $ff;
my @queryl = ();
GetOptions(
           "help|h"=>\$help,
           "parser|format|p=s" => \$fmt,
           "handler|writer|w=s" => \$out,
	   "count|c" => \$count,
           "xml"=>\$toxml,
           "perl"=>\$toperl,
           "debug"=>\$debug,
	   "filterfile|f=s"=>\$ff,
	   "query|q=s@"=>\@queryl,
          );
if ($help) {
    system("perldoc $0");
    exit 0;
}

my $w = shift;
my $sub;
if ($ff) {
    $sub = do $ff;
    if ($@) {
        die $@;
    }
}
elsif (@queryl) {

    my $ev = 
      'sub { my $s=shift; '.
	join(' && ',
	     map {
		 if (
                     /([\w\/]+)\s*(==|<=|>=|<|>|=)\s*(.*)/ ||    # op
                     /([\w\/]+)\s+(\S\S)\s+(.*)/                 # lt,gt,etc
                    ) {
		     my ($var, $op, $val) = ($1, $2, $3);
                     if ($op eq '=') {
                         if ($val =~ /^\d+$/ ||
                             $val =~ /^\d*\.\d+/) {
                             $op = '==';
                         }
                         else {
                             $op = 'eq';
                             $val = "'$val'";
                         }
                     }
		     "defined \$s->get('$var') && \$s->get('$var') $op $val"
		 }
		 else {
		     die($_);
		 }
	     } @queryl).
	       '}';
    $sub = eval $ev;
    if ($@) {
	die $@;
    }
}
else {
    $sub = shift;
    $sub = eval $sub;

}
if ($@) {
    print $@;
    exit 1;
}

my $c = 0;
my @files = @ARGV;
foreach my $fn (@files) {

    my $handler = Data::Stag->makehandler($w=>sub {
					      my $self = shift;
					      my $stag = shift;
					      my $ok = $sub->($stag);
					      if ($ok) {
						  $c++;
						  return $stag;
					      }
					      else {
						  $stag->free;
						  return;
					      }
					    });
    if ($count) {
	$out = 'Data::Stag::null';
    }
    
    if (!$out) {
	$out = 'xml';
    }
    my $ch = Data::Stag->chainhandlers($w, $handler, $out);

    my @pargs = (-file=>$fn, -format=>$fmt, -handler=>$ch);
    if ($fn eq '-') {
	if (!$fmt) {
	    $fmt = 'xml';
	}
	@pargs = (-format=>$fmt, -handler=>$ch, -fh=>\*STDIN);
    }

    my $tree = 
      Data::Stag->parse(@pargs);

    if ($count) {
	print "$c\n";
    }

#    my @res =
#      $tree->where($w,
#		   $sub);
#    print $_->xml foreach @res;
}

exit 0;

__END__

=head1 NAME 

stag-grep.pl - filters a stag file (xml, itext, sxpr) for nodes of interest

=head1 SYNOPSIS

  stag-grep.pl person -q name=fred file1.xml

  stag-grep.pl person 'sub {shift->get_name =~ /^A*/}' file1.xml

  stag-grep.pl -p My::Foo -w sxpr record 'sub{..}' file2

=head1 USAGE

  stag-grep.pl [-p|parser PARSER] [-w|writer WRITER] NODE -q tag=val FILE

  stag-grep.pl [-p|parser PARSER] [-w|writer WRITER] NODE SUB FILE

  stag-grep.pl [-p|parser PARSER] [-w|writer WRITER]  NODE -f PERLFILE FILE

=head1 DESCRIPTION

parsers an input file using the specified parser (which may be a built
in stag parser, such as xml) and filters the resulting stag tree
according to a user-supplied subroutine, writing out only the
nodes/elements that pass the test.

the parser is event based, so it should be able to handle large files
(although if the node you parse is large, it will take up more memory)

=head1 ARGUMENTS

=over

=item -p|parser FORMAT

FORMAT is one of xml, sxpr or itext, or the name of a perl module

xml assumed as default

=item -w|writer FORMAT

FORMAT is one of xml, sxpr or itext, or the name of a perl module

=item -c|count

prints the number of nodes that pass the test

=item -filterfile|f

a file containing a perl subroutine (in place of the SUB argument)

=item -q|query TAG1=VAL1 -q|query TAG2=VAL2 ...  -q|query TAGN=VALN

filters based on the field TAG

other operators can be used too - eg <, <=, etc

multiple q arguments can be passed in

for more complex operations, pass in your own subroutine, see below

=item SUB

a perl subroutine. this subroutine is evaluated evry time NODE is
encountered - the stag object for NODE is passed into the subroutine.

if the subroutine passes, the node will be passed to the writer for
display

=item NODE

the name of the node/element we are filtering on

=item FILE

the file to be parser. If no parser option is supplied, this is
assumed to a be a stag compatible syntax (xml, sxpr or itext);
otherwise you should parse in a parser name or a parser module that
throws stag events

=back

=head1 SEE ALSO

L<Data::Stag>

=cut
