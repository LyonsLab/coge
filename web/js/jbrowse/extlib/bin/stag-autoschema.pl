#!/usr/bin/perl -w

eval 'exec /usr/bin/perl -w -S $0 ${1+"$@"}'
    if 0; # not running under some shell

# POD docs at bottom of file


use strict;

use Carp;
use Data::Stag qw(:all);
use Getopt::Long;

my $parser = "";
my $handler = "sxpr";
my $mapf;
my $tosql;
my $toxml;
my $toperl;
my $debug;
my $help;
my $name;
my $defs;
my $dtd;
GetOptions(
           "help|h"=>\$help,
           "parser|format|p=s" => \$parser,
           "handler|writer|w=s" => \$handler,
           "xml"=>\$toxml,
           "perl"=>\$toperl,
           "debug"=>\$debug,
	   "name|n=s"=>\$name,
	   "defs"=>\$defs,
           "dtd"=>\$dtd,
          );
if ($help) {
    system("perldoc $0");
    exit 0;
}

$dtd = 1 if $handler eq 'dtd';

#my @hdr = ();
#if ($name) {
#    push(@hdr, (name=>$name));
#}

my @files = @ARGV;
my $tree = Data::Stag->new(null=>[]);
foreach my $fn (@files) {

    my $curr_tree = 
      Data::Stag->parse($fn, 
                        $parser);
    $tree->name($curr_tree->name);
    $tree->addkid($_) foreach $curr_tree->subnodes;
}

my $s = $tree->autoschema;
if ($defs) {
    my @sdefs = ();
    my %u = ();
    $s->iterate(sub {
                    my $stag = shift;
                    my $n = $stag->name;
                    $n =~ s/[\+\?\*]$//;
                    return if $u{$n};
                    $u{$n}=1;
                    push(@sdefs, ($n=>''));
                    return;
                });
    $s = Data::Stag->unflatten(schemadefs=>[@sdefs]);
} else {
    $s->iterate(sub {
                    my $stag = shift;
                    my $d = $stag->data;
                    if (!ref $d) {
                        $stag->data($d =~ /INT/ ? "i" : "s");
                    }
                });
}
#    my $top = 
#      Data::Stag->unflatten(schema=>[
#				     @hdr,
#				    ]);
#    $top->set_nesting($s->data);
if ($dtd) {
    print $s->dtd;
} else {
    print $s->generate(-fmt=>$handler);
}

__END__

=head1 NAME 

stag-autoschema.pl - writes the implicit stag-schema for a stag file

=head1 SYNOPSIS

  stag-autoschema.pl -w sxpr sample-data.xml

  stag-autoschema.pl -dtd sample-data.xml

=head1 DESCRIPTION

Takes a stag compatible file (xml, sxpr, itext), or a file in any
format plus a parser, and writes out the implicit underlying stag-schema

stag-schema should look relatively self-explanatory.

Here is an example stag-schema, shown in sxpr syntax:

  (db
   (person*
    (name "s"
    (address+
     (address_type "s")
     (street "s")
     (street2? "s")
     (city "s")
     (zip? "s")))))

The database db contains zero or more persons, each person has a
mandatory name and at least one address.

The cardinality mnemonics are as follows:

=over

=item +

1 or more

=item ?

0 or one

=item *

0 or more

=back

The default cardinality is 1

=head1 ARGUMENTS

=over

=item -p|parser FORMAT

FORMAT is one of xml, sxpr or itext, or the name of a perl module

xml assumed as default

=item -dtd

exports schema as DTD

=item -w|writer FORMAT

FORMAT is one of xml, sxpr or itext, or the name of a perl module, OR DTD

The default is sxpr

note that stag schemas exported as xml will be invalid xml, due to the
use of symbols *, +, ? in the node names

=back

=head1 LIMITATIONS

not event based - memory usage becomes exhorbitant on large files;
prepare a small sample beforehand

=cut


