#!/usr/bin/perl -w

use strict;
use CoGe::Accessory::Web;
use Data::Dumper;
use File::Spec::Functions;
use CGI;

my $q = CGI->new;
my $gid = $q->param('gid');
my $chr = $q->param('chr');
my $ws = $q->param('ws');

# Connect to the database
my $conf = CoGe::Accessory::Web->get_defaults();

my $filename = $gid . "_" . $chr . "_" . $ws . "_out.txt";

print "Content-Type: application/json\n\n";
my $path = catfile($conf->{SECTEMPDIR}, "downloads/genome", $gid);
my $file = catfile($path, $filename);
open(my $fh, $file);
my $l = <$fh>; # skip header line
my (@at, @gc, @n, @x);
while ($l = <$fh>) {
    chomp $l;
    my @tokens = split /\t/, $l;
    push @at, $tokens[3];
    push @gc, $tokens[4];
    push @n, $tokens[5];
    push @x, $tokens[6];
}
print '{"at":[';
print join(',', @at);
print '],"gc":[';
print join(',', @gc);
print '],"n":[';
print join(',', @n);
print '],"x":[';
print join(',', @x);
print ']}';
close $fh;