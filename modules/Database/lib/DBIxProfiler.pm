package DBIxProfiler;
use strict;
use warnings;

use base qw(DBIx::Class::Storage::Statistics);

use Time::HiRes qw(time);

my $start;
my $QCOUNT = 0;

sub query_start {
  my $self = shift();
  my $sql = shift();
  my @params = @_;

  print STDERR "Executing $sql: ".join(', ', @params)."\n";
  $start = time();
}

sub query_end {
  my $self = shift();
  my $sql = shift();
  my @params = @_;

  printf "Execution took %0.4f seconds.\n", time() - $start;
  $QCOUNT += 1;
  $start = undef;
}

sub DESTROY {
  print STDERR "*** QCOUNT=", $QCOUNT, "\n";
}

1;
