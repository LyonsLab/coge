#!/usr/bin/perl -w

eval 'exec /usr/bin/perl -w -S $0 ${1+"$@"}'
    if 0; # not running under some shell

# POD docs at bottom of file

use strict;

use Data::Stag qw(:all);
use Getopt::Long;

my $parser = "";
my $mapf;
my $help;
my @ignore = ();
my @report = ();
my $trace;
my $verbose;
GetOptions(
           "help|h"=>\$help,
           "parser|format|p=s"=>\$parser,
	   "ignore|s=s@"=>\@ignore,
	   "report|r=s@"=>\@report,
	   "trace|t"=>\$trace,
	   "verbose|v"=>\$verbose,
          );

our %REPORT = map {$_=>1} @report;

if ($help) {
    system("perldoc $0");
    exit 0;
}

my @files = @ARGV;
my $file1 = shift @files;
if (!@files) {
    die "you must pass in at least two files";
}
my $stag1 = Data::Stag->parse($file1, $parser);
$stag1->unset($_) foreach @ignore;
foreach my $file2 (@files) {
    my $stag2 = Data::Stag->parse($file2, $parser);
    $stag2->unset($_) foreach @ignore;
    my ($match, $reason) = match($stag1, $stag2);
    if ($match) {
	printf "SAME: $file1 $file2\n";
    }
    else {
	printf "DIFF: $file1 $file2\n";
	printf "REASON:\n";
	showreason($reason);
    }
}
exit 0;

sub showreason {
    my $reason = shift;
    my $indent = shift || 0;
    print '  ' x $indent;
    my ($msg, @children) = @$reason;
    printf $msg;
    print "\n";
    showreason($_, $indent+1) foreach @children;
}

sub match {
    my $stag1 = shift;
    my $stag2 = shift;
    if ($stag1->name ne $stag2->name) {
	return(0, mismatch("name_mismatch", [$stag1, $stag2]));
    }
    trace("comparing %s", $stag1->name, $stag2->name);
    my $t1 = $stag1->isterminal || 0;
    my $t2 = $stag2->isterminal || 0;
    if ($t1 != $t2) {
	return(0, mismatch("different_node_types", [$stag1, $stag2]));
    }
    if ($t1 && $t2) {
	if ($stag1->data eq $stag2->data) {
	    return 1;
	}
	else {
	    return(0, mismatch(sprintf("data_mismatch(%s ne %s)",
				       smallstr($stag1->data),
				       smallstr($stag2->data),
				      ), 
			       [$stag1, $stag2]));
	}
    }
    # both nodes nonterminal
    if ($t1 || $t2) {
	die "assertion error";
    }

    trace("  ..looking at kids\n");
    my @kids1 = $stag1->kids;
    my @kids2 = $stag2->kids;

    # must match exactly
    if (@kids1 != @kids2) {
	return(0, mismatch(sprintf("subelement_count_mismatch [%s <=VS=> %s]",
				   names(@kids1), names(@kids2)),
			   [$stag1, $stag2]));
    }
    # null always matches
    if (!@kids1) {
	# both must be null
	die "assertion error" unless !@kids2;
	return 1;
    }
    trace("  ..matrix:\n");
    my @filled = ();
    for (my $i=0; $i<@kids1; $i++) {
	my $kid1 = $kids1[$i];
	my $matched;
	my @reasons = ();
	for (my $j=0; $j<@kids2; $j++) {
	    next if $filled[$j];
	    my $kid2 = $kids2[$j];
	    next unless $kid1->name eq $kid2->name;
	    my ($match, $reason) = match($kid1, $kid2);
	    if ($match) {
		$filled[$j] = 1;
		$matched = 1;
		last;
	    }
	    else {
		push(@reasons, $reason);
	    }
	}
	if (!$matched) {
	    my $mismatch = 
	      mismatch("no_matching_node", [$kid1]);
	    push(@$mismatch, @reasons);
	    return(0, $mismatch);
	}
    }
    trace("  ..match!\n");
    return 1;
}

sub names {
    join(', ', map {$_->name} @_);
}

sub mismatch {
    my $msg = shift;
    my @stags = @{shift || []};
    my @names = map {$_->name} @stags;
    my $reason =
      sprintf "$msg: %s",
	join(' AND ', @names);
    if (grep {$REPORT{$_}} @names) {
	printf "$reason\n";
	if ($verbose) {
	    print $_->sxpr foreach @stags;
	}
    }
    return [$reason];
}

sub smallstr {
    my $str = shift;
    return $str if length($str) < 50;
    return substr($str, 0, 50) ."...";
}

sub trace {
    return unless $trace;
    my $fmt = shift;
    printf  $fmt, @_;
    print  "\n";
}

__END__

=head1 NAME 

stag-diff.pl - finds the difference between two stag files

=head1 SYNOPSIS

  stag-diff.pl -ignore foo-id -ignore bar-id file1.xml file2.xml

=head1 DESCRIPTION

Compares two data trees and reports whether they match. If they do not
match, the mismatch is reported.

=over ARGUMENTS

=item -help|h

shows this document

=item -ignore|i ELEMENT

these nodes are ignored for the purposes of comparison. Note that
attributes are treated as elements, prefixed by the containing element
id. For example, if you have

  <foo ID="wibble">

And you wish to ignore the ID attribute, then you would use the switch

  -ignore foo-ID

You can specify multiple elements to ignore like this

  -i foo -i bar -i baz

You can also specify paths

  -i foo/bar/bar-id

=item -parser|p FORMAT

which parser to use. The default is XML. This can also be autodetected
by the file suffix. Other alternatives are B<sxpr> and B<itext>. See
L<Data::Stag> for details.

=item -report|r ELEMENT

report mismatches as they occur on each element of type ELEMENT

multiple elements can be specified

=item -verbose|v

used in conjunction with the B<-report> switch

shows the tree of the mismatching element

=back

=head2 OUTPUT

If a mismatch is reported, a report is generated displaying the
subpart of the tree that could not be matched. This will look like
this:

REASON:
no_matching_node: annotation
  no_matching_node: feature_set
    no_matching_node: feature_span
      no_matching_node: evidence
        no_matching_node: evidence-id
          data_mismatch(:15077290 ne :15077291): evidence-id AND evidence-id

Due to the nature of tree matching, it can be difficult to specify
exactly how trees do not match. To investigate this, you may need to
use the B<-r> and B<-v> options. For the above output, I would
recommend using

  stag-diff.pl -r feature_span -v

=head2 ALGORITHM

Both trees are recursively traversed... see the actual code for how this works

The order of elements is not important; eg
  
  <foo>
    <bar>
      <baz>1</baz>
    </bar>
    <bar>
      <baz>2</baz>
    </bar>
  </foo>

matches

  <foo>
    <bar>
      <baz>2</baz>
    </bar>
    <bar>
      <baz>1</baz>
    </bar>
  </foo>

The recursive nature of this algorithm means that certain tree
comparisons will explode wrt time and memory. I think this will only
happen with very deep trees where nodes high up in the tree can only
be differentiated by nodes low down in the tree.

Both trees are loaded into memory to begin with, so it may thrash with
very large documents

=head2 AUTHOR

Chris Mungall 
cjm at fruitfly dot org

=head1 SEE ALSO

L<Data::Stag>

=cut
