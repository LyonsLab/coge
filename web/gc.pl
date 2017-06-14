#! /usr/bin/perl -w

use strict;

use CGI;

use CoGe::Accessory::Utils qw( commify );
use CoGe::Accessory::Web;
use CoGe::Core::Sequence;

my ($db, $user, $conf) = CoGe::Accessory::Web->init();
my $s = CoGe::Core::Sequence->new(29349);
my $percentages = $s->get_percentages(3, $db);
my $form = new CGI;
print $form->header, '<html><body>';
print "Total length: "
      . commify($percentages->{total})
      . " bp, GC: "
      . sprintf( "%.2f", $percentages->{gc} )
      . "%  AT: "
      . sprintf( "%.2f", $percentages->{at} )
      . "%  NX: "
      . sprintf( "%.2f", $percentages->{nx} )
      . "%";
print '</body></html>';