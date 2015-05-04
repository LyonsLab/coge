#!/usr/bin/env perl
#-------------------------------------------------------------------------------
# Purpose:  Wrapper script for Tophat for use with JEX.  TopHat's command-line
#     usage with paired-end reads will not work with JEX's workflow interface.
# Author:   mdb
# Created:  1/30/15
#------------------------------------------------------------------------------

use strict;
use warnings;

use CoGe::Accessory::Utils qw(execute);

my $read_type = shift;  # 'single' or 'paired'
my $cmd_args = shift;   # tophat executable and all command line arguments except for file lists
my @files = @ARGV;      # fastq files

if ($read_type eq 'paired') {
    my @m1 = grep { $_ =~ /\_R1/ } @files;
    my @m2 = grep { $_ =~ /\_R2/ } @files;
    die "error: invalid paired-end files" unless (@m1 and @m2);
    $cmd_args .= join(',', sort @m1) . ' ' . join(',', sort @m2);
}
else {
    $cmd_args .= join(',', @files);
}

print STDOUT "tophat.pl: $cmd_args\n";
my $rc = execute($cmd_args);

exit($rc);
