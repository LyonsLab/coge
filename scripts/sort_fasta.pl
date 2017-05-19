#!/usr/bin/env perl
#-------------------------------------------------------------------------------
# Purpose:  Sort fasta sections by length (descending order).  A SAMtools
#           faidx file is expected to be located in the same directory as the
#           input fasta file.  The faidx is used to prevent loading the entire
#           FASTA file into memory.
# Usage:    perl sort_fasta.pl <input_file>
# Author:   Matt Bomhoff
# Created:  10/11/16
#-------------------------------------------------------------------------------
# Parameters:
#    <none>
#-------------------------------------------------------------------------------

use strict;
use Data::Dumper;

my $input_file = shift;
die "Usage:  perl sort_fasta.pl <input_file>\n" unless (defined $input_file);

die "Cannot access input file: $input_file\n" unless (-r $input_file);

my $faidx_file = "$input_file.fai";
unless (-r $faidx_file) {
    `samtools faidx $input_file`;
}
die "Couldn't create index file" unless (-r $faidx_file);

my $index = read_faidx($faidx_file);
die unless $index;

open(my $fh, $input_file) or
    die("Error: cannot open file '$input_file'\n");

foreach my $rec (sort { $b->{length} <=> $a->{length} } @$index) {
    seek($fh, $rec->{offset}, 0);

    print '>', $rec->{name}, "\n";

    my $count = 0;
    while (my $line = <$fh>) {
        chomp $line;
        $line=~ tr/\015//d;
        $count += length $line;
        print $line, "\n";
        last if ($count >= $rec->{length});
    }
}

close($fh);
exit;

#-------------------------------------------------------------------------------
sub read_faidx {
    my $faidx_filepath = shift;
    my $index;

    # File format described here: http://www.htslib.org/doc/faidx.html
    open(my $fh, $faidx_filepath) or die("Error: cannot open file '$faidx_filepath'\n");
    while (<$fh>) {
        chomp;
        my ($name, $len, $offset, $linebases, $linewidth) = split("\t");
        next unless ($name && $len);

        push @$index, {
            name => $name,
            length => $len,
            offset => $offset,
            linebases => $linebases,
            linewidth => $linewidth
        };
    }
    close($fh);

    return $index;
}
