#!/usr/bin/env perl
#-------------------------------------------------------------------------------
# Purpose:  Wrapper script for Bowtie/Tophat for use with JEX.  
#     Bowtie & TopHat's command-line usage with paired-end reads will not work 
#     with JEX's workflow interface.
# Author:   mdb
# Created:  1/30/15
#------------------------------------------------------------------------------

use strict;
use warnings;

use Getopt::Long;
use CoGe::Accessory::Utils qw(execute detect_paired_end);

our ($read_type, $cmd_args, $file_list, $output_filename);

GetOptions(
    "read_type=s"  => \$read_type, # 'single' or 'paired'
    "cmd_args=s"   => \$cmd_args,  # tophat executable and all command line arguments except for file lists
    "files=s"      => \$file_list, # comma-separated list of input FASTQ files
    "output=s"     => \$output_filename, # name of output BAM file (mdb added 10/5/16 for ChIP-seq pipeline)
);

die unless ($read_type && $cmd_args && $file_list);

my @files = split(',', $file_list);

if ($read_type eq 'paired') {
    my ($m1, $m2) = detect_paired_end(\@files);
    die "error: invalid paired-end files: m1: @$m1 -- m2: @$m2" unless (@$m1 and @$m2);
    $cmd_args .= join(',', sort @$m1) . ' ' . join(',', sort @$m2);
}
else {
    $cmd_args .=  join(' ', @files);
}

print STDOUT "Running command: $cmd_args\n";
my $rc = execute($cmd_args);

# Rename output BAM file -- mdb added 10/5/16 for ChIP-seq pipeline
my $default_output_filename = 'accepted_hits.bam';
if ($rc == 0 && $output_filename && $output_filename ne $default_output_filename) {
    print STDOUT "Renaming $default_output_filename to $output_filename\n";
    rename($default_output_filename, $output_filename);
}

exit($rc);
