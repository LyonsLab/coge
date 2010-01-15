#!/usr/bin/perl -w

use strict;
use Data::Dumper;
use CoGeX;
use DBIxProfiler;

use Getopt::Long;

my ($coge, $fasta_name, $name_re, @datasets, $debug, $dsg);

GetOptions(
	   "fasta_name=s" => \$fasta_name,
	   "name_search|name_re|name|n=s"=>\$name_re,
	   "dataset|ds=s"=>\@datasets,
	   "dataset_group|dsg=s"=>\$dsg,
	   "debug" => \$debug,
	   );
  unless (@datasets || $dsg)
  {
    print qq#
Welcome t0 $0

Usage:  $0 -dataset 24 -dataset NC_000001  -name_search regex_search -fasta_name output.faa
    or for multiple datasets
        $0 -ds 40504 -ds 40500 -ds 40495 -ds 40499 -ds 40503 -fasta_name arabidopsis_v9.fasta -name_re 'AT\\dG\\d{5}\$' > arabidopsis_v9.gff

Options:

 -dataset | -ds                        Datasets to retrieve, either by name or by database id

 -dataset | -ds                        Dataset_group to retrieve, either by name or by database id

 -name_search | -name_re | -name | -n  OPTIONAL: Regular expression to search for specific feature names

 -fasta_name                           OPTIONAL:  create a fasta file of the sequences

 -debug                                OPTIONAL:  print debugging messages

#;
    exit();
  }

$coge = CoGeX->dbconnect();
#$coge->storage->debugobj(new DBIxProfiler());
#$coge->storage->debug(1);

# gt sketch -seqid 1 -addintrons yes -start 1000 -style default.style -force -end 70000  out.png grape.gff3

push @datasets, get_datasets($dsg) if $dsg;
my $chrs = get_locs(\@datasets);
my $schrs = get_sequence(ds=>\@datasets, file_name=>$fasta_name) if $fasta_name;

foreach my $k (keys %$chrs){
    if (! $schrs->{$k}){
        print STDERR "NOTFOUND in seqs:" . $k . "\n" if $debug;
    }
}
foreach my $k (keys %$schrs){
    if (! $chrs->{$k}){
        print STDERR "NOTFOUND in locs:" . $k . "\n" if $debug;
    }
}

sub get_datasets
  {
    my $dsgid = shift;
    my $dsg = $coge->resultset('DatasetGroup')->resolve($dsgid);
    my @ds = $dsg->datasets;
    return map {$_->id} @ds;
  }

sub get_locs {

    my $datasets = shift;
    my $SEP = "\t";
    my %chrs;
    foreach my $ds (@$datasets){
        my ($dso) = $coge->resultset('Dataset')->resolve($ds);
        foreach my $chr ($dso->get_chromosomes){
            $chrs{$chr} = $dso->last_chromosome_position($chr);
        }
    }
    my @chrs = sort { $a cmp $b } keys %chrs;

    print "##gff-version\t3\n";
    foreach my $chr (@chrs){
        print "##sequence-region $chr 1 " . $chrs{$chr} . "\n";
    }
    my %fids = ();
    my $count=0;
    my %types;
    foreach my $chr (@chrs){
        my %seen = ();
        #if ($i++ > 2){ print STDERR "*" x 1000 . "\nENDING EARLY";  last; }
        
        my $feat_rs = $coge->resultset('Feature')->search( {
                  'me.dataset_id' => { 'IN' => $datasets },
                  'me.chromosome' => $chr,
#                  'me.feature_type_id'  =>  1
                } , { 
                   'prefetch'           => [ 'feature_type', 'feature_names'] 
                  ,'order_by'           => [ 'me.start', 'me.feature_type_id'] #go by order in genome, then make sure that genes (feature_type_id == 1) is first
                });

        #gff: chr  organization feature_type  start stop strand . name
        print STDERR "dataset_ids: " . join(",", @$datasets) . ";  chr: $chr\n" if $debug;
        while(my $feat = $feat_rs->next()){
            if ($fids{$feat->feature_id}){ next; }
            my @feat_names;
            if($name_re){
                @feat_names = grep { $_ =~ /$name_re/i } $feat->names(); 
		next unless @feat_names;
            }
            else {
                @feat_names = $feat->names(); 
            }
            my $strand = $feat->strand == 1 ? '+' : '-';
	    my ($names) = join (",", @feat_names);
            my $attrs = "ID=$count";
	    $attrs .= ";Name=$names" if $names;
	    my $annos =join (",", map {$_->annotation_type->name.": ".$_->annotation} $feat->annotations);
            my $gstr = join("\t", ($chr, 'coge', $feat->feature_type->name, $feat->start, $feat->stop, ".", $strand, ".", $attrs));
	    $gstr .= ";Note=$annos" if $annos;
            if($seen{$gstr}){ next; }
            $seen{$gstr} = 1;
            print $gstr . "\n";
	    $types{$feat->feature_type->name}++;
            $fids{$feat->feature_id} = 1; #feat_id has been used
            my $parent = $count;
	    $count++;
	    next unless $feat->feature_type_id == 1 && @feat_names; #if not a gene, don't do the next set of searches.
            my $mrna_rs = $coge->resultset('Feature')->search( {
                  'me.dataset_id' => { 'IN' => $datasets },
                  'me.chromosome' => $chr,
                  'feature_names.name' =>  {'IN'=>[@feat_names]},
                  'me.feature_type_id'  =>  2
                } , { 
                    'join' => 'feature_names',
                   'prefetch'           => [ 'feature_type', 'locations'] 
                  ,'order_by'           => [ 'me.start', 'locations.start']
                });

            while(my $f = $mrna_rs->next()){
	      if($fids{$f->feature_id}){ next; }
		my $mrna_names = join (",", $f->names);
                my $mrna_attrs = "Parent=$parent";
		$mrna_attrs .= ";Name=$mrna_names" if $mrna_names;
		my $mrna_annos =join (",", map {$_->annotation_type->name.": ".$_->annotation} $f->annotations);
                foreach my $loc ($f->locations()){
		  next if $loc->start > $feat->stop || $loc->stop < $feat->start; #outside of genes boundaries;  Have to count it as something else
		  my $gstr = join("\t", ($f->chr, 'coge', $f->feature_type->name, $loc->start, $loc->stop, ".", $strand, ".", $mrna_attrs));
		  $gstr.=";Note=$mrna_annos" if $mrna_annos;
		  if($seen{$gstr}){ next; }
		  $seen{$gstr} = 1;
		  print $gstr . "\n";
		  $fids{$f->feature_id} = 1; #feat_id has been used;
		  $types{$f->feature_type->name}++;
                }
            }
            my $sub_rs = $coge->resultset('Feature')->search( {
                  'me.dataset_id' => { 'IN' => $datasets },
                  'me.chromosome' => $chr,
                  'feature_names.name' =>  {'IN'=>[@feat_names]},
                  , 'me.feature_type_id'  =>  { 'NOT IN' => [1,2] }
                } , { 
                   'join'               => [ 'feature_names'],
                   'prefetch'           => [ 'feature_type', 'locations'] 
                  ,'order_by'           => [ 'me.chromosome', 'me.start']
                });
            while(my $f = $sub_rs->next()){
	      if($fids{$f->feature_id}){ next; }
	      #my $locs = $f->locations({}, {'order_by' => 'start'});
	      my $other_names = join (",", $f->names);
	      my $other_attrs = "Parent=$parent";
	      $other_attrs .= ";Name=$other_names" if $other_names;
	      my $other_annos =join (",", map {$_->annotation_type->name.": ".$_->annotation} $f->annotations);
	      foreach my $loc ($f->locations()){
		next if $loc->start > $feat->stop || $loc->stop < $feat->start; #outside of genes boundaries;  Have to count it as something else
		my $gstr = join("\t", ($f->chr, 'coge', $f->feature_type->name, $loc->start, $loc->stop, ".", $strand, ".", $other_attrs));
		$gstr.=";Note=$other_annos" if $other_annos;

		if($seen{$gstr}){ next; }
		$seen{$gstr} = 1;
		print $gstr . "\n";
		$fids{$f->feature_id} = 1; #feat_id has been used;
		$types{$f->feature_type->name}++;
	      }
            }
        }
    }
    if ($debug)
      {
	print STDERR "Feature_types:\n\t";
	print STDERR join ("\n\t", map {$_."\t".$types{$_}} sort keys %types),"\n";
	  
      }
    return \%chrs;
}


sub get_sequence {
  my %opts = @_;
  my $datasets = $opts{ds};
  my $file_name = $opts{file_name};
  open(FA, ">", $file_name);
    my %chrs;
    foreach my $ds (@$datasets){
        my $ds = $coge->resultset('Dataset')->resolve($ds);
        my %seen;
        foreach my $chr ($ds->get_chromosomes){
            # TODO: this will break with contigs.
            next if $chr =~ /random/;
            $chrs{$chr} = 1;
            print FA ">$chr\n";
            print FA $ds->get_genomic_sequence(chromosome => $chr) . "\n";
        }
    }
    close FA;
    return \%chrs;
}
