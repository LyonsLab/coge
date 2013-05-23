#!/usr/bin/perl -w

eval 'exec /usr/bin/perl -w -S $0 ${1+"$@"}'
    if 0; # not running under some shell

# POD docs at end

use strict;

use Data::Stag qw(:all);
use Getopt::Long;


my @cols = ();
my $sep = "\t";
my $parser;
my $nest;
my $errhandler;
my $errf;
GetOptions(
           "parser|format|p=s" => \$parser,
           "errhandler=s" => \$errhandler,
           "errf|e=s" => \$errf,
	   "cols|c=s@"=>\@cols,
	   "help|h"=>sub { system("perldoc $0"); exit },
           "nest|n"=>\$nest,
	  );

$errhandler =  Data::Stag->getformathandler($errhandler || 'xml');
if ($errf) {
    $errhandler->file($errf);
}
else {
    $errhandler->fh(\*STDERR);
}

my $node = shift;
my $fn = shift @ARGV;
push(@cols, @ARGV);
@cols = map {split/\,/,$_} @cols;

my $np = scalar @cols;
my %idx = map {$cols[$_]=>$_} (0..$#cols);
my @vals = map {[]} @cols;
my @level_idx = ();

sub setcol {
    my ($col, $val) = @_;
    my $i = $idx{$col};
    die $col unless defined $i;
    push(@{$vals[$i]}, $val);
    return;
}

my %catch = ();
foreach my $col (@cols) {
    $catch{$col} =
      sub {
	  my ($self, $stag) = @_;
	  setcol($col, $stag->data);
	  return;
      };
}
$catch{$node} =
  sub {
      my ($self, $stag) = @_;
      push(@$_, 'NULL') foreach grep {!@$_} @vals; # "left join" - null for non-existent
      
      if ($nest) {
          print 
            join("\t", map {'{'.join(', ', @$_).'}'} @vals), 
              "\n";
      }
      else {
          my $c = 1;
          $c *= scalar(@$_) foreach @vals;
          if ($c == 1) {
              # no combinatorial explosion
              print 
                join("\t", map {$_->[0]} @vals), 
                  "\n";
              
          }
          elsif ($c == 0) {
              print STDERR $stag->xml;
              die "assertion error @vals";
          }
          else {
              # do cartesian explosion
              #
              # we COULD do this recursively but it would be slow
              my @N = map {0} @vals;
              my $done = 0;
              while (!$done) {
                  for (my $i=0; $i<@N; $i++) {
                      print "\t" if $i;
                      print $vals[$i]->[$N[$i]];
                  }
                  print "\n";
                  my $i = $#N;
                  my $carry_the_one = 1;
                  while ($i >= 0 && $carry_the_one) {
                      $N[$i] ++;
                      if ($N[$i] >= @{$vals[$i]}) {
                          $N[$i] = 0;
                          $i--;
                      }
                      else {
                          $carry_the_one = 0;
                      }
                  }
                  $done = 1 if $i<0;
              }
          }
      }
      @vals = map {[]} @cols;      
      return;
  };

my $h = Data::Stag->makehandler(%catch);
Data::Stag->parse(-file=>$fn, -format=>$parser, -handler=>$h, -errhandler=>$errhandler);
exit 0;

__END__

=head1 NAME 

stag-flatten.pl - turns stag data into a flat table

=head1 SYNOPSIS

  stag-flatten.pl -c name -c person/name dept MyFile.xml

=head1 DESCRIPTION

reads in a file in a stag format, and 'flattens' it to a tab-delimited
table format. given this data:

  (company
   (dept
    (name "special-operations")
    (person
     (name "james-bond"))
    (person
     (name "fred"))))

the above command will return a two column table

  special-operations      james-bond
  special-operations      fred

If there are multiple values for the columns within the node, then the
cartesian product will be calculated



=head1 USAGE

  stag-flatten.pl [-p PARSER] [-c COLS] [-c COLS] NODE <file>

=head1 ARGUMENTS

=over

=item -p|parser FORMAT

FORMAT is one of xml, sxpr or itext

xml assumed as default

=item -c|column COL1,COL2,COL3,..

the name of the columns/elements to write out

this can be specified either with multiple -c arguments, or with a
comma-seperated (no spaces) list of column (terminal node) names after
a single -c

=item -n|nest

if set, then the output will be a compress repeating values into the
same row; each cell in the table will be enclosed by {}, and will
contain a comma-delimited set of values

=back

=head1 SEE ALSO

L<Data::Stag>

=cut

