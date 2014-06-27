#!/usr/bin/perl -w

use CoGeX;
use Data::Dumper;
use Getopt::Long;

my ($type, $file, @anno);
GetOptions ("t|type=s"=>\$type,
	    "f|file=s"=>\$file,
	    "a|anno=s"=>\@anno,
	    );

my $connstr = 'dbi:mysql:dbname=DB;host=HOST;port=PORT';
my $coge = CoGeX->connect($connstr, 'USER', 'PASSWORD' );

my $data = process_file($file);
#print Dumper $data;
my $seq_count = 0;
org: foreach my $orgname (keys %$data)
  {
    my ($org) = $coge->resultset('Organism')->search(name=>$orgname);
    print "No org in database for $orgname!\n" unless ($org);
    my $found =0;
    ds: foreach my $ds ($coge->get_current_datasets_for_org(org=>$org->id))
      {
	feat: foreach my $feat ($ds->features({name=>$type},{join=>feature_type}))
	  {

	    my $name = $orgname.": ". join (", ", $feat->names, map {$_->annotation} $feat->annotations)."\n";
	    if (@anno)
	      {
           my $pass=0;
		foreach my $anno (@anno)
		  {
		    $pass =1 if $name =~ /$anno/i;
		  }
          next feat unless $pass;
	      }
	    $found =1;
	    print ">".$seq_count." ".$name;
	    print $feat->genomic_sequence,"\n";
	    $seq_count++;
	    next org;
	  }
      }
    print STDERR "no sequence for ",$orgname,"\n"unless $found;
  }

sub process_file
  {
    my $file = shift;
    my %data;
    open (IN, $file);
    while (<IN>)
      {
	chomp;
	my @line = split/\t/;
	my ($org, $desc) = split /  /,$line[2];
	$data{$org} = $desc;
      }
    close IN;
    return \%data;
  }
