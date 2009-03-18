package DBIxProfiler;
use strict;
use warnings;

use base qw(DBIx::Class::Storage::Statistics);

use Time::HiRes qw(time);

my $start;
my $QCOUNT = 0;
my $TOTAL = 0;

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

  printf STDERR "Execution took %0.4f seconds.", time() - $start;
  print STDERR " QCOUNT: $QCOUNT\n";
  $QCOUNT += 1;
  $TOTAL += time() - $start;
  $start = undef;
}

sub get_query_count_and_time {
    my $self = shift();
    return $QCOUNT, $TOTAL;
}

sub DESTROY {
  #print STDERR "TOTAL QUERIES: " . $QCOUNT . "\n";
}

1;
