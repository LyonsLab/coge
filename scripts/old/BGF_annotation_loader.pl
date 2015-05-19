#!/usr/bin/perl -w

use strict;
use Data::Dumper;
use CoGeX;
use Getopt::Long;

my $GO = 0;
my $DEBUG = 1;
my $dsid;
my $add_gene =0;
my $connstr = 'dbi:mysql:dbname=DB;host=HOST;port=PORT';
my $coge = CoGeX->connect($connstr, 'USER', 'PASSWORD' );
#$coge->storage->debugobj(new DBIxProfiler());
#$coge->storage->debug(1);

GetOptions ( "dsid=i" => \$dsid,
	     "go=s"    => \$GO,
	     "debug=s" => \$DEBUG,
	     "add_gene_feature" => \$add_gene,
	   );

my $ds = $coge->resultset('Dataset')->find($dsid);

unless ($ds)
  {
    warn "unable to find a valid dataset entry for $dsid\n";
    exit;
  }

warn "-go flag is not true, nothing will be added to the database.\n" unless $GO;
my %data;
my %annos;
while (<>)
  {
    next if /^GENE_ID/;
    s/^\s+//;
    chomp;
    my @line = split /\t/;
    my $chr = $line[0];
    $chr =~ s/chromosome//i;
    $chr =~ s/chr//i;
    $chr =~ s/^_//i;
    $chr =~ s/^0//;
    my $name = $line[2];
    my $strand = $line[10] eq "+" ? 1 : -1;
    #add gene
    ($line[11],$line[12]) = ($line[12],$line[11]) if $line[11]>$line[12];
    push @{$data{$name}{gene}}, {
				     start=>$line[11],
				     stop=>$line[12],
				     strand=>$strand,
				     chr=>$chr,
				    };
    my @start = split/,/,$line[-2];
    my @stop = split/,/,$line[-1];
    for (my $i=0; $i< @start; $i++)
      {
	my ($start, $stop) = ($start[$i], $stop[$i]);
	($start, $stop) = ($stop, $start) if $start > $stop;
	push @{$data{$name}{mRNA}}, {
				     start=>$start,
				     stop=>$stop,
				     strand=>$strand,
				     chr=>$chr,
				    };

	push @{$data{$name}{CDS}}, {
				     start=>$start,
				     stop=>$stop,
				     strand=>$strand,
				     chr=>$chr,
				    };

      }
#    push @{$data{$name}{$line[2]}}, {
#				     start=>$line[3],
#				     stop=>$line[4],
#				     strand=>$strand,
#				     chr=>$chr,
#				    };
  }
if ($add_gene)
  {
    name: foreach my $name (keys %data)
      {
	my $start;
	my $stop;
	my $strand;
	my $chr;
	foreach my $type (keys %{$data{$name}})
	  {
	    foreach my $loc (@{$data{$name}{$type}})
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
	$data{$name}{gene}=[{
			     start=>$start,
			     stop=>$stop,
			     strand=>$strand,
			     chr=>$chr,
			    }];
      }
  }
#print Dumper \%data;
#print Dumper \%annos;
#exit;

my ($anno_type) = $coge->resultset('AnnotationType')->search({name=>"note"});
foreach my $name (keys %data)
  {
    foreach my $feat_type (keys %{$data{$name}})
      {
	my ($start) = sort {$a<=>$b} map {$_->{start}} @{$data{$name}{$feat_type}};
	my ($stop) = sort {$b<=>$a} map {$_->{stop}} @{$data{$name}{$feat_type}};
	my ($strand) = map {$_->{strand}} @{$data{$name}{$feat_type}};
	my ($chr) = map {$_->{chr}} @{$data{$name}{$feat_type}};
	my $feat_type_obj = $coge->resultset('FeatureType')->find_or_create( { name => $feat_type } ) if $GO;
	print "Creating feature of type $feat_type\n" if $DEBUG;

	my $feat = $ds->add_to_features({
					 feature_type_id => $feat_type_obj->id,
					 start=>$start,
					 stop=>$stop,
					 chromosome=>$chr,
					 strand=>$strand,
					}) if $GO;
	my $featid = $feat->id if $feat;
	foreach my $loc (@{$data{$name}{$feat_type}})
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
	print "Adding name $name to feature ", $featid ,"\n" if $DEBUG;
	my $feat_name = $feat->add_to_feature_names({
						     name=>$name,
						     #				   feature_id=>$featid,
						    }) if $GO ;
	if ($DEBUG && $annos{$name})
	  {
	    print "Adding annotation $annos{$name}\n" ;
	    foreach my $anno (@{$annos{$name}})
	      {
		my $annoo = $feat->add_to_annotations({annotation=>$anno, annotation_type_id => $anno_type->id}) if $GO && $anno;
	      }
	  }
      }
  }
