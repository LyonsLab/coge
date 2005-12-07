#!/usr/bin/perl -w

use strict;
use Data::Dumper;
use lib "/home/elyons/projects/pairwise/lib/";
use Pairwise;

my $pw = new Pairwise;
$pw->seqA("PELICAN");
$pw->seqB("CEELLICAN");
print join "\n",$pw->global_align();
print"\n";
$pw->print_dpm();
