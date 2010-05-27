#!/usr/bin/perl -w

use CoGeX;

my $connstr = 'dbi:mysql:dbname=genomes;host=biocon;port=3306';
my $coge = CoGeX->connect($connstr, 'cnssys', 'CnS' );


foreach my $feat ($coge->resultset('Feature')->search({"feature_names.name"=>"rbcL"},
						      {join=>"feature_names"},
						      ))
  {
    next unless $feat->feature_type->name eq "CDS";
    next unless $feat->dataset->organism->name =~ /chloroplast/ || $feat->dataset->organism->description =~ /chloroplast/;
    print ">",join (", ", "Names: ",$feat->names,"Organsim: ".$feat->dataset->organism->name.": ".$feat->dataset->organism->description, "Annotations: ", map {$_->annotation} $feat->annotations),"\n";
    
    print $feat->genomic_sequence,"\n";
    
  }
