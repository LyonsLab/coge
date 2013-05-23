#!/usr/bin/perl -w

eval 'exec /usr/bin/perl -w -S $0 ${1+"$@"}'
    if 0; # not running under some shell

use Data::Stag qw(:all);
use Data::Stag::XMLParser;
use Data::Stag::ITextWriter;
my $p = Data::Stag::XMLParser->new;
my $h = Data::Stag::ITextWriter->new;
$p->handler($h);
foreach my $xmlfile (@ARGV) {
    $p->parse($xmlfile);
}

