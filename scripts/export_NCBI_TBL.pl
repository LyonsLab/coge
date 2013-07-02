#!/usr/bin/perl -w

use strict;
use CoGeX;
use Data::Dumper;
use Getopt::Long;
use Sort::Versions;

use vars qw($DEBUG $GO $db $user $pass $coge $dsgid);


# ./export_NCBI_TBL.pl -u coge -p 123coge321 -db coge -dsgid 16884

GetOptions ( "debug=s" => \$DEBUG,
             "go=s"    => \$GO,
             "database|db=s"=>\$db,
             "user|u=s"=>\$user,
             "password|pw|p=s"=>\$pass,
             "dsgid=i"=>\$dsgid,
           );

my $connstr = "dbi:mysql:dbname=$db;host=localhost;port=3307";
$coge = CoGeX->connect($connstr, $user, $pass );
#$coge->storage->debugobj(new DBIxProfiler());
#$coge->storage->debug(1);

my $dsg = $coge->resultset('DatasetGroup')->find($dsgid);

my %chr2ds;
foreach my $ds($dsg->datasets)
  {
    map {$chr2ds{$_}=$ds} $ds->chromosomes;
  }
foreach my $chr (sort {versioncmp($a, $b)} $dsg->chromosomes)
  {
    print ">Features $chr\n";
    foreach my $feat (sort {$a->start <=> $b->start || $a->feature_type_id <=> $b->feature_type_id} $chr2ds{$chr}->features({chromosome=>$chr}, {"order_by"=>"start ASC"}))
      {
	next if $feat->type->name eq "chromosome";
	next if $feat->type->name =~ /utr/i;
	  
	my $locs = get_locs($feat);;
	my $item = pop @$locs;
	print join ("\t", @$item, $feat->type->name),"\n";
	foreach my $loc (@$locs)
	  {
	    print join ("\t", @$loc),"\n";
	  }
	my $name_tag = get_name_tag($feat);
	my ($pri_name, @names) = $feat->names;
	print "\t\t\t",$name_tag,"\t";
	if ($name_tag =~ /_id/)
	  {
	    print "gnl|dbname|";
	  }
	print $pri_name,"\n";
	
	foreach my $name (@names)
	  {
	    print "\t\t\t", join ("\t", "alt_name", $name),"\n";
	  }
	foreach my $anno ($feat->annotations)
	  {
	    print "\t\t\t", join ("\t", $anno->type->name,$anno->annotation),"\n";
	  }
      }
  }

sub get_name_tag
  {
    my $feat = shift;
    my $tag;
    if ($feat->type->name =~ /mRNA/) {$tag = "transcript_id"}
    elsif ($feat->type->name =~ /CDS/) {$tag = "protein_id"}
    else {$tag = "locus_tag"}
    return $tag;
  }

sub get_locs
  {
    my $feat = shift;
    my @locs = map {[$_->start, $_->stop]} sort {$a->start <=> $b->start} $feat->locations;
    if ($feat->strand =~ /-/)
      {
	@locs = reverse @locs;
	foreach my $item (@locs)
	  {
	    ($item->[0], $item->[1]) = ($item->[1], $item->[0]);
	  }
      }
    return \@locs;
  }
