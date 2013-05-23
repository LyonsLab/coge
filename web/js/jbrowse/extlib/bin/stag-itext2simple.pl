#!/usr/bin/perl -w

eval 'exec /usr/bin/perl -w -S $0 ${1+"$@"}'
    if 0; # not running under some shell

use Data::Stag qw(:all);
use Data::Stag::ITextParser;
use Data::Stag::Simple;
use Data::Stag::Base;
my $p = Data::Stag::ITextParser->new;
my $h = Data::Stag::Simple->new;
#my $h = Data::Stag::XMLWriter->new;
$p->handler($h);
foreach my $f (@ARGV) {
    $p->parse($f);
}

