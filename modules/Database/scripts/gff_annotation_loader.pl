#!/usr/bin/perl -w

use strict;
use Data::Dumper;
use CoGeX;
use Getopt::Long;

my $GO = 0;
my $DEBUG = 1;
my $dsid;
my $add_gene =0;
my @names;
my $connstr = 'dbi:mysql:dbname=genomes;host=biocon;port=3306';
my$coge = CoGeX->connect($connstr, 'cnssys', 'CnS' );
#$coge->storage->debugobj(new DBIxProfiler());
#$coge->storage->debug(1);


GetOptions ( "dsid=i" => \$dsid,
	     "go=s"    => \$GO,
	     "debug=s" => \$DEBUG,
	     "name=s" => \@names,
	     "add_gene_feature" => \$add_gene,
	   );

my $ds = $coge->resultset('Dataset')->find($dsid);

unless ($dsid)
  {
    warn "No dataset id specified\nUsage: $0 -dsid dataset_id -go 1 < gff_annotation_file.gff\n";
    exit;
  }

unless ($ds)
  {
    warn "unable to find a valid dataset entry for $dsid\n";
    exit;
  }
push @names, "mRNA" unless @names;

warn "-go flag is not true, nothing will be added to the database.\n" unless $GO;
my %data;
my %annos;
my %feat_types; #store feature type objects
my ($anno_type) = $coge->resultset('AnnotationType')->search({name=>"note"}); #generic annotation type
my $prev_type;
while (<>)
  {
    next if /^#/;
    chomp;
    my @line = split /\t/;
    next if $line[2] eq "clone";
    next if $line[2] eq "mRNA";
    next if $line[2] eq "intron";
    my $chr = $line[0];
    $chr =~ s/chromosome//i;
    $chr =~ s/chr//i;
    $chr =~ s/^_//i;
    my %names;
    my $name;
    foreach my $item (split /;/, $line[-1])
      {
	my $tmp;
	$item =~ s/"//g;
	my ($type, $info) = $item =~ /=/ ? split (/=/,$item,2) : (split / /,$item,2);
	foreach my $namecheck (@names)
	  {
	    $type = "name" if $type eq $namecheck;
	  }
	if ($type eq "name")
	  {
	    $tmp = $info;
	  }
	next if $type eq "exonNumber";
	next if $type eq "transcriptId";
	next if $type eq "proteinId";
	next if $tmp && $tmp =~ /intron/;
	next if $tmp && $tmp =~ /exon/;
	next if $tmp && $tmp =~ /cds/;
	$name = $tmp unless $name;
	$names{$tmp}=1 if $tmp;
	next unless $name; #No name, don't know what to do!
	$annos{$name}{$info}=1 if $type eq "Description";
	print join ("\t",@line),"\n" if ($type eq "biotype" && !$name);
	$annos{$name}{$info}=1 if $type eq "biotype";
      }
    next unless $name; #No name, don't know what to do!
    my $strand = $line[6] =~ /-/ ? -1 :1;
    $line[2] = "mRNA" if $line[2] eq "exon";
    push @{$data{$name}{$line[2]}{loc}}, {
					  start=>$line[3],
					  stop=>$line[4],
					  strand=>$strand,
					  chr=>$chr,
					 };
    map {$data{$name}{$line[2]}{names}{$_}=1} keys %names;
#    print Dumper \%data;
  }
if ($add_gene)
  {
    name: foreach my $name (keys %data)
      {
	my $start;
	my $stop;
	my $strand;
	my $chr;
	my %names;
	foreach my $type (keys %{$data{$name}})
	  {
	    map {$names{$_}=1} keys %{$data{$name}{$type}{names}};
	    foreach my $loc (@{$data{$name}{$type}{loc}})
	      {
		next name if $type eq "gene";
		$start = $loc->{start} unless $start;
		$start = $loc->{start} if $loc->{start} < $start;
		$stop = $loc->{stop} unless $stop;
		$stop = $loc->{stop} if $loc->{stop} > $stop;
		$strand = $loc->{strand};
		$chr = $loc->{chr};
	      }
	  }
	$data{$name}{gene}{loc}=[{
				  start=>$start,
				  stop=>$stop,
				  strand=>$strand,
				  chr=>$chr,
				 }];
	$data{$name}{gene}{names}= \%names;
      }
  }
#print Dumper \%data;
#print Dumper \%annos;
#exit;


foreach my $name (keys %data)
  {
    foreach my $feat_type (keys %{$data{$name}})
      {
	my ($start) = sort {$a<=>$b} map {$_->{start}} @{$data{$name}{$feat_type}{loc}};
	my ($stop) = sort {$b<=>$a} map {$_->{stop}} @{$data{$name}{$feat_type}{loc}};
	my ($strand) = map {$_->{strand}} @{$data{$name}{$feat_type}{loc}};
	my ($chr) = map {$_->{chr}} @{$data{$name}{$feat_type}{loc}};
	$feat_types{$feat_type} = $coge->resultset('FeatureType')->find_or_create( { name => $feat_type } ) if $GO && !$feat_types{$feat_type};
	my $feat_type_obj = $feat_types{$feat_type};
	
	print "Creating feature of type $feat_type\n" if $DEBUG;
	
	my $feat = $ds->add_to_features({
					 feature_type_id => $feat_type_obj->id,
					 start=>$start,
					 stop=>$stop,
					 chromosome=>$chr,
					 strand=>$strand,
					}) if $GO;
	my $featid = $feat ? $feat->id : "no_go";
	foreach my $loc (@{$data{$name}{$feat_type}{loc}})
	  {
	    print "Adding location $chr:(".$loc->{start}."-".$loc->{stop}.", $strand)\n" if $DEBUG;
	    my $loc_tmp = $feat->add_to_locations(
						  {
						   start      => $loc->{start},
						   stop       => $loc->{stop},
						   strand     => $loc->{strand},
						   chromosome => $loc->{chr}
						  }
						 ) if $GO;
	  }
	foreach my $tmp (keys %{$data{$name}{$feat_type}{names}})
	  {
	    print "Adding name $tmp to feature ", $featid ,"\n" if $DEBUG;
	    my $feat_name = $feat->add_to_feature_names({
							 name=>$tmp,
							 #				   feature_id=>$featid,
							}) if $GO ;
	  }
	if ($DEBUG && $annos{$name})
	  {
	    foreach my $anno (keys %{$annos{$name}})
	      {
	    print "Adding annotation $anno\n" if $DEBUG;
		my $annoo = $feat->add_to_annotations({annotation=>$anno, annotation_type_id => $anno_type->id}) if $GO && $anno;
	      }
	  }
      }
  }
