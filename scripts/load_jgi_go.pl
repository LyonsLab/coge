#!/usr/bin/perl  -w

use strict;
use Data::Dumper;
use CoGeX;
use Getopt::Long;

use vars qw($GO $DEBUG $dsid $gff_file $go_file);

my $connstr = 'dbi:mysql:dbname=DB;host=HOST;port=PORT';
my $coge = CoGeX->connect($connstr, 'USER', 'PASSWORD' );
#$coge->storage->debugobj(new DBIxProfiler());
#$coge->storage->debug(1);

GetOptions ( "dsid=i" => \$dsid,
             "go=s"    => \$GO,
             "debug=s" => \$DEBUG,
             "gff_file=s"=>\$gff_file,
	     "go_file=s"=>\$go_file,
           );

$DEBUG = 1 unless defined $DEBUG;

unless ($dsid && -r $gff_file && -r $go_file)
  {
    print qq{Usage: $0 -dsid <coge database dataset_id> -gff_file <jgi gff annotation file previously loaded into CoGe> -go_file <jgi go annotation file>

The gff file needs to be previously loaded to and present in order to match jgi protein ids (arbitrary numbers from their system) to gene names in CoGe and Go annotaions

};
    exit;
  }

my $protid2name = process_gff_file($gff_file);
process_go_file($go_file);

sub process_go_file
  {
    my $file = shift;
    my %go;
    my %goid;
    my %orgs;

    open (IN, $file);
    while (<IN>)
      {
	next if /^#/;
	next if /^proteinId/;
	chomp;
	my @line = split /\t/;
	my $atg = $go{$line[3]};
	unless ($atg)
	  {
	    my $go_cat = $line[3];
	    $go_cat = "GO $go_cat";
	    $go_cat =~ s/_/ /g;
	    ($atg) = $coge->resultset('AnnotationTypeGroup')->search({name=>$go_cat});
	    $go{$line[3]} = $atg;
	  }

	print " no go cat for $line[3]\n" if $DEBUG && !$atg;
	my $at = $goid{$line[4]};
	unless ($at)
	  {
	    my $goid = $line[4];
	    $goid =~ s/GO://g;
	    foreach my $tmp ($coge->resultset('AnnotationType')->find_or_create({name=>$goid}))
	      {
		$at = $tmp unless $at;
		unless ($tmp->description)
		  {
		    warn "adding description $line[2] to $line[4]\n";
		    $tmp->description($line[2]);
		    $tmp->update;
		  }
		unless ($tmp->annotation_type_group)
		  {
		    $tmp->annotation_type_group_id($atg->id);
		    $tmp->update;
		    print "here\n";
		  }
	      }
	    $goid{$line[4]} = $at;
	  }
	foreach my $fn ($coge->resultset('FeatureName')->search({name=>$protid2name->{$line[0]}}))
	  {
	    next if $fn->feature->dataset->id ne $dsid;
	    my $ftype = $fn->feature->type->name;
	    next unless $ftype eq "gene" || $ftype eq "mRNA" || $ftype eq "CDS";
	    $orgs{$fn->feature->dataset->organism->name}=1;
	    $fn->feature->add_to_annotations({annotation_type_id=>$at->id});
	  }
      }
    close IN;
    print Dumper \%orgs;
  }

sub process_gff_file
  {
    my $file = shift;
    my %data;
    open (IN, $file);
    while (<IN>)
      {
	next if /^#/;
	chomp;
	my @line = split /\t/;
	my $name;
	foreach my $item (split /;/, $line[-1])
	  {
	    $item =~ s/"//g;
	    $item =~ s/^\s+//;
	    $item =~ s/\s+$//;
	    my ($type, $thing) = split/\s/,$item,2;
	    $name = $thing if $type eq "name";
	    $data{$thing} = $name if $type eq "proteinId";
	  }
      }
    close IN;
    return \%data;
  }
