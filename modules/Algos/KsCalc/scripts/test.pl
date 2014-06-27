#!/usr/bin/perl -w

use strict;
use CoGe::Algos::KsCalc;
use Data::Dumper;
use File::Temp;

my $TEMPDIR = "/tmp/";
my $object = CoGe::Algos::KsCalc->new ();

$object->version(5);
$object->name1("At1g10850");
$object->name2("At1g60630");
$object->palign();
$object->print_align();
print "\n\n";
print length($object->prot1),":", length($object->gaplessP1),"\n";
print length($object->prot2),":", length($object->gaplessP2),"\n";
print length($object->dna1),":", length($object->gaplessD1),"\n";
print length($object->dna2),":", length($object->gaplessD2),"\n";
print "\n\n";
print $object->dna1,"\n";

#print $object->gaplessP1,"\n";
#print $object->gaplessP2,"\n";
#print $object->gaplessD1,"\n";
#print $object->gaplessD2,"\n";
print $object->phylip_align;
#$object->print_align(seq1=>$object->gaplessD1, seq2=>$object->gaplessD2);
#print Dumper $object;
my $aln = new File::Temp ( TEMPLATE=>'Ks__XXXXX',
			   DIR=>$TEMPDIR,
			   SUFFIX=>'.aln',
			   UNLINK=>0);
print Dumper $object->KsCalc($aln);
