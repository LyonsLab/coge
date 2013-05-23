#!/usr/bin/perl -w

eval 'exec /usr/bin/perl -w -S $0 ${1+"$@"}'
    if 0; # not running under some shell

# POD docs at end

use strict;

use Carp;
use Data::Stag qw(:all);
use Getopt::Long;

my $parser = "";
my $writer = "";
my $debug;
my $help;
GetOptions(
           "help|h"=>\$help,
           "parser|format|p=s" => \$parser,
           "handler|writer|w=s" => \$writer,
           "debug"=>\$debug,
          );
if ($help) {
    system("perldoc $0");
    exit 0;
}

my $joinfield = shift @ARGV;
my ($j1, $j2) = ($joinfield, $joinfield);
if ($joinfield =~ /(.*)=(.*)/) {
    $j1 = $1;
    $j2 = $2;
}
my ($jn1, $jk1) = splitnodekey($j1);
my ($jn2, $jk2) = splitnodekey($j2);


my $fn1 = shift @ARGV;
my @files = @ARGV;

print "$jn1 $jk1 $jn2 $jk2\n";
my %index = ();
my $idxhandler = 
  Data::Stag->makehandler(-NOTREE=>1,
			  $jn2=>sub {
			      my ($self,$stag) = @_;
			      my $key = $stag->sget($jk2);
			      $index{$key} = $stag->duplicate;
			      $stag;
			  });

foreach my $fn (@files) {

    Data::Stag->parse(-file=>$fn,
		      -format=>$parser,
		      -handler=>$idxhandler);
}

if ($debug) {
    printf "INDEXED KEYVALS: %d\n", scalar(keys %index);
}

my $jhandler =
  Data::Stag->makehandler($jn1=>sub {
			      my ($self,$stag) = @_;
			      my $key = $stag->sget($jk1);
			      my $val = $index{$key};
			      if ($val) {
				  $stag->add_species($val->data);
			      }
			      $stag;
			  });
my $whandler =
  Data::Stag->chainhandlers([$jn1],
			    $jhandler,
			    $writer);
Data::Stag->parse(-file=>$fn1,
		  -format=>$parser,
		  -handler=>$whandler);
exit 0;

sub splitnodekey {
    my $s = shift;
    if ($s =~ /(\w+)\/(.*)/) {
	return ($1,$2);
    }
    return ($s,'');
}

__END__

=head1 NAME 

stag-join.pl - joins two stag files together based around common key

=head1 SYNOPSIS

  stag-join.pl  -w xml country/city_id=capital/capital_id countries.xml capitals.xml

  stag-join.pl  -w itext gene/tax_id=species/tax_id genedb.itext speciesdb.itext

=head1 DESCRIPTION

Performs a relational-style INNER JOIN between two stag trees; this
effectively merges two files together, based on some kind of ID in the
file

=head1 ARGUMENTS

=item -p|parser FORMAT

FORMAT is one of xml, sxpr or itext, or the name of a perl module

xml assumed as default

=item -w|writer FORMAT

FORMAT is one of xml, sxpr or itext, or the name of a perl module

=head1 LIMITATIONS

currently not event based, so may not be memory efficicent. could be
easily rewritten to be event based

=head1 SEE ALSO

L<Data::Stag>

This script is a wrapper for the method

  Data::Stag->ijoin()
 
=cut

