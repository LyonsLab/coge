#!/usr/bin/perl -w

use strict;
use Data::Dumper;
use CoGeX;

my $connstr = 'dbi:mysql:dbname=DB;host=HOST;port=PORT';
my $coge = CoGeX->connect($connstr, 'USER', 'PASSWORD' );

#need to check: , ; = / : " ' <space> [ ] |  ( ) { }
foreach my $fn ($coge->resultset('FeatureName')->search({name=>{'like'=>'CDS%'}}))
  {
    my $name = $fn->name;
    next unless $name =~ /CDS/;
    my $feat = $fn->feature;
    print $name,"\n";
#     if (
# 	$name =~ /modifier/i ||
#        )
#       {
#	 $fn->name($name);
#	 $fn->update();
#	 $fn->delete();
#	 $feat->add_to_annotations({annotation=>$name});
#       }
    $fn->delete();
#    $feat->add_to_annotations({annotation=>$name});
#    $feat->add_to_feature_names({name=>$item})
  }
exit();
foreach my $fn ($coge->resultset('FeatureName')->search({name=>{'like'=>'%,%'}}))
  {
    my $name = $fn->name;
#    next if $name =~ /[^,] [^,]/;
    $name =~ s/(\d,) (\d)/$1$2/g;
    next if $name =~ /\S,\S/;
    my $feat = $fn->feature;
    my $name_find =0;
    if (
	$name =~ /motif/i ||
	$name =~ /synthetase/i ||
	$name =~ /domain/ ||
	$name =~ /family/ ||
	$name =~ /tree/ ||
	$name =~ /isozyme/ ||
	$name =~ /summer/ ||
	$name =~ /subunit/ ||
	$name =~ /spruce/ ||
#	$name =~ // ||
	$name =~ /similar/i
       )
      {
	$feat->add_to_annotations({annotation=>$name});
      }
    else
      {
	foreach my $item (split /,/,$name)
	  {
	    $item =~ s/ and //;
	    $item =~ s/^\s+//;
	    $item =~ s/\s+$//;
	    $item =~ s/si://;
	    next unless $item;
	    next if $item =~ /^puta/i;
	    next if $item =~ /^put$/i;
	    next if $item =~ /^biosynthetic$/i;
	    next if $item =~ /^..?$/;
	    next if $item =~ /chloropl$/;
	    next if $item =~ /^expr$/;
	    if ($item =~ /\s\S+\s/) #annotations
	      {
		$feat->add_to_annotations({annotation=>$name});
	      }
	    elsif (
		   $item =~ /type/i ||
		   $item =~ /domain/i ||
		   $item =~ /unit/i ||
		   $item =~ /family/i ||
		   $item =~ /class/i ||
		   $item =~ /kinase/i ||
		   $item =~ /catalytic/i ||
		   $item =~ /terminal/i ||
		   $item =~ /uncharacterized/i ||
		   $item =~ /tigr/i ||
		   $item =~ /depend/i ||
		   $item =~ /cyto/i ||
		   $item =~ /chloro/i ||
		   $item =~ /mitoc/i ||
		   $item =~ /NAD/i ||
		   $item =~ /-like/i ||
		   $item =~ /kDa/ ||
		   $item =~ /pseudo/i ||
		   $item =~ /similar/i ||
		   $item =~ /variant/i ||
		   $item =~ /chain/i ||
		   $item =~ /3'/i ||
		   $item =~ /5'/i ||
		   $item =~ /terminus/i ||
		   $item =~ /component/i ||
		   $item =~ /dihydr/i ||
		   $item =~ /N-t/i ||
		   $item =~ /non-/i ||
		   $item =~ /mit$/i ||
		   $item =~ /DNA-/i ||
		   $item =~ /^\S\s\S$/ ||
		   $item =~ /euka/i ||
		   $item =~ /chaper/i ||
		   $item =~ /^com$/i ||
		   $item =~ /^cat$/i ||
		   $item =~ /^chlor?$/i ||
		   $item =~ /^exp$/i ||
		   $item =~ /^chl$/i ||
		   $item =~ /eine/i ||
		   $item =~ /amine/i ||
		   $item =~ /subu/i ||
		   $item =~ /-term/i ||
		   $item =~ /protein/i ||
		   $item =~ /expr/i ||
		   $item =~ /finger/i ||
		   $item =~ /isof/i ||
		   $item =~ /hyper/i ||
		   $item =~ /transfer/i ||
		   $item =~ /transpo/i ||
		   $item =~ /homolo/i ||
		   $item =~ /contain/i ||
		   $item =~ /rotatio/i ||
		   $item =~ /partia/i ||
		   $item =~ /mito/i ||
		   $item =~ /touch/i ||
		   $item =~ /pyruvate/i ||
		   $item =~ /molyb/i ||
		   $item =~ /subfa/i ||
		   $item =~ /assemb/i ||
		   $item =~ /decarb/i ||
		   $item =~ /group/i ||
#		   $item =~ //i ||

		   $item =~/.{15}/
		  )
	      {
		$feat->add_to_annotations({annotation=>$name});
	      }
	    else
	      {
		$feat->add_to_feature_names({name=>$item});
	      }
	  }
      }
    $fn->delete();
  }
#
#error gene names with "putative" are clearly annotations
foreach my $fn ($coge->resultset('FeatureName')->search({name=>{'like'=>'%putative%'}}))
  {
    my $name = $fn->name;
    my $feat = $fn->feature;
    print $name,"\n";
    $feat->add_to_annotations({annotation=>$name});
    $fn->delete();
  }

foreach my $fn ($coge->resultset('FeatureName')->search({name=>{'like'=>"%;%"}}))
  {
    my $name = $fn->name;
    #example error: name fgenesh1_pg.C_LG_I000019; proteinId 62667; ex
    if ($name =~ /^name/)
      {
	print $name,"\n";
	$name =~ s/^name //;
	$name =~ s/;.*//;
	print $name,"\n\n";
	next unless $name;
	$fn->name($name);
	$fn->update;
      }
#
    #example error: mpt51; mpb51; fbpC1
    if ($name =~ /; \S\S\S/)
      {
	my $feat = $fn->feature;
	print $name,"\n";
	foreach my $item (split /;/,$name)
	  {
	    $item =~ s/^\s+//;
	    $item =~ s/\s+$//;
	    $feat->add_to_feature_names({name=>$item});
	    print "\t",$item,"\n";
	  }
	print "\n";
	$fn->delete();
      }
    #print $name,"\n";
  }
