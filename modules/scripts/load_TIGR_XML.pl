#!/usr/bin/perl 

## Written by Eric Lyons, UC Berkeley 2004
## Contact:  elyons@nature.berkeley.edu


use CNS::TIGRXMLLIB::TIGR::TIGR_XML_parser;
use CNS::TIGRXMLLIB::TIGR::Gene_obj;
use CNS::TIGRXMLLIB::BSBC::Feature;
use CNS::TIGRXMLLIB::BSBC::FeatureDB;
use strict;
use Getopt::Long;
use Data::Dumper;
use DBI;
use Carp;


use CoGe::Genome;

my ($file, $verbose, $help, $chr, $empty_db, $version);

GetOptions('file=s'=>\$file,
	   'verbose|v'=>\$verbose,
	   'help|h' => \$help,
	   'chromosome|chr=s' => \$chr,
	   'empty_db' =>\$empty_db,
	   'version=s' => \$version,
	  );

my $org_name = "Rice";
my $org_desc = "Oryza sativa (japonica)";
my $data_source_name = "TIGR";
my $data_source_desc = "The Institute for Genomic Research";
my $data_source_link = "www.tigr.org";
my ($data_info_name) = $file=~ /([^\/]*$)/;
my $data_info_desc = "Rice Version $version Chromosome $chr, TIGR XML format";
my $data_info_link = "ftp://ftp.tigr.org/pub/data/Eukaryotic_Projects/o_sativa/annotation_dbs/pseudomolecules/version_3.0";
my $genomic_seq_len = 10000;

unless (-r $file)
  {
    $help = 1;
    print "WARNING: You need a valid XML file.\n";
  }

unless ($chr)
  {
    $help = 1;
    print "WARNING: You need to specify a chromosome!\n\n";
  }
unless ($version)
  {
    $help = 1;
    print "WARNING: You need to specify a data source version!\n\n";
  }


show_help() if $help; 

my $pwd = `pwd`;
chomp $pwd;
unless ($file =~ /\//) {
    #full path not specified.
    $file = "$pwd/$file";
}

my ($genome_db, $organism, $data_source, $data_info) = initialize_genome_db($org_name, $org_desc, $data_source_name, $data_source_desc, $data_source_link, $data_info_name, $data_info_desc, $data_info_link, $version);

empty_db($data_info) if $empty_db;


my $TIGRparser = new CNS::TIGRXMLLIB::TIGR::TIGR_XML_parser();
$TIGRparser->capture_genes_from_assembly_xml("$file");
load_genomic_sequence($TIGRparser->get_assembly_sequence(), $genomic_seq_len);

#print $TIGRparser->toString();

my $feature_count = 0;
my $gene_count = 0;
my $public_locus_count = 0;
my %count;
foreach my $gene ($TIGRparser->get_genes()) {
  $gene->refine_gene_object();
  $gene->create_all_sequence_types(\$TIGRparser->{assembly_seq}) if $TIGRparser->{assembly_seq};
  my $entry = $gene->toString;
  my $locus = new CNS::TIGRXMLLIB::BSBC::Feature(verbose=>$verbose);
#  print $entry;
  $gene_count++;
  $public_locus_count++ if $gene->{'pub_locus'};
  foreach my $item (split /ISOFORM:/, $entry)
    {
      $locus->process_entry($gene, $item);
      $count{$locus->gene_type}++;
#      print Dumper $locus;
      $feature_count++;
      my @feats;
      #gene feature
      my $gene_feat = create_feature("gene");
      load_locations($gene_feat, $locus->gene_location, $locus->strand, $chr);
      load_xref($gene_feat, $locus->supporting_cdna());
      push @feats, $gene_feat;
      if  ($locus->gene_type() =~ /rna/i)
	{
	  my $type = $locus->gene_type();
	  $type =~ s/rna/RNA/;
	  my $rna_feat = create_feature($type);
	  load_locations($rna_feat, $locus->exon_locations, $locus->strand, $chr);
	  push @feats, $rna_feat;
	}
      if ($locus->gene_type() eq "protein-coding")
   	{
	  #mRNA feature
	  if ($locus->exon_locations) #untested!
	    {
	      my $mrna_feat = create_feature("mRNA");
	      load_locations($mrna_feat, $locus->exon_locations, $locus->strand, $chr);
	      push @feats, $mrna_feat;
	    }
   	  #CDS feature
	  if ($locus->CDS_locations) #untested
	    {
	      my $cds_feat = create_feature("CDS");
	      load_locations($cds_feat, $locus->CDS_locations, $locus->strand, $chr);
	      load_protein($cds_feat, $locus->protein_seq());
	      push @feats, $cds_feat;
	    }
	}
      #names
      foreach my $feat (@feats)
	{
	  #names
	  foreach my $name (
			    $locus->public_locus(),
			    $locus->public_locus_model(),
			    $locus->TU_feat_name(),
			    $locus->model_feat_name(),
			    $locus->alt_locus(),
			    $locus->locus(),
			    $locus->gene_synonyms(),
			    split/;\s*/, $locus->common_name(),
			   )
	    {
	      next unless $name;
	      load_name($feat, $name);
	    }
	  #GO
	  load_go($feat, $locus->GO);
	  #additional annotations
	  load_annotation($feat, "secondary products", $locus->secondary_products()) if $locus->secondary_products;
	  load_annotation($feat, "public comment", $locus->public_comment()) if $locus->public_comment;
	  
	}
    }
}

print "Total processed genes: $gene_count\n";
print "  genes with a public locus:  $public_locus_count\n";
print "Total processed entries (a gene may have multiple isoforms): $feature_count\n";
print "Types processed:";
print Dumper (\%count);

sub load_xref
  {
    my ($feat, $xrefs) = @_;
    return 0 unless ref $xrefs =~ /ARRAY/i;
    foreach my $xref (@$xrefs)
      {
	load_annotation($feat, $xref->{dbsource}, $xref->{accn}, "xref");
	load_name($feat, $xref->{accn});
	load_annotation($feat, "GI", $xref->{gi}, "xref");
      }
  }

sub load_go
  {
    my ($feat, $gos) = @_;
    return 0 unless ref ($gos) =~ /ARRAY/i;
    foreach my $go (@$gos)
      {
	load_annotation($feat, $go->{goid}, $go->{goterm}, "GO ".$go->{gocat});
      }
  }

sub load_annotation
  {
    my ($feat, $type, $annotation, $group) = @_;
    return 0 unless $annotation && $type;
    my %type_opts = (name=>$type);
    if ($group)
      {
	my $anno_type_group = $genome_db->get_annotation_type_group_obj->find_or_create({name=>$group}); 
	$type_opts{annotation_type_group_id}=$anno_type_group->id();
	  
      }
    my $anno_type = $genome_db->get_annotation_type_obj->find_or_create({
									 %type_opts
									});
    $genome_db->get_annotation_obj->create({
					    annotation=>$annotation,
					    feature_id=>$feat->id(),
					    annotation_type_id=>$anno_type->id(),
					   });
  }

sub load_protein
  {
    my ($feat, $seq) = @_;
    return 0 unless $seq;
    my $seq_type = $genome_db->get_sequence_type_obj->find_or_create({
								      name=>"protein",
								      description=>"translation",
								     });
    my $sequence = $genome_db->get_sequence_obj->create({
							 sequence_type_id=>$seq_type->id(),
							 sequence_data=>$seq,
							 feature_id=>$feat->id(),
							});
  }

sub load_name
  {
    my ($feat, $name, $desc) = @_;
    return 0 unless $name;
    my $feat_name = $genome_db->get_feature_name_obj->create({
							      name=>$name,
							      description=>$desc,
							      feature_id=>$feat->id
							     });
  }

sub load_locations
  {
    my ($feat, $locs,$strand, $chr) = @_;
    return 0 unless ref ($locs) =~ /ARRAY/i;
    return 0 unless $strand;
    return 0 unless $chr;
    my @locs;
    foreach my $loc (@$locs)
      {
	my $location = $genome_db->get_location_obj->create({
							     feature_id=>$feat->id,
							     start=>$loc->{sbegin},
							     stop=> $loc->{send},
							     strand=>$strand, #expect this to be "1" or "-1"
							     chromosome=>$chr, #make sure that chromosomes are in real numbers (not roman numerals [IX])
							    });
	push @locs, $location;
      }
    return \@locs;
  }


sub create_feature
  {
    my ($type) = @_;
    return 0 unless $type;
    my $feat_type = $genome_db->get_feature_type_obj->find_or_create({name=>$type});

    my $db_feature = $genome_db->get_feature_obj->create({feature_type_id=>$feat_type->id,
							  data_information_id=>$data_info->id,
							  organism_id=>$organism->id,
							 });
    return $db_feature;
  }

sub load_genomic_sequence
  {
    my ($seq, $len) = @_;
    my $seqlen = length $seq;
    print STDERR "Loading genomic sequence ($seqlen nt)\n";
    $seq =~ s/\s//g;
    $seq =~ s/\n//g;
    my $i = 0;
    my $gso = $genome_db->get_genomic_sequence_obj;
    while ($i < $seqlen)
      {
	my $str = substr($seq,$i,$len);
	my $start = $i+1;
	my $stop = $i + length $str;
	$gso->create({start=>$start,
                      stop=>$stop,
                      chromosome=>$chr,
                      sequence_data=>$str,
                      organism_id=>$organism->id,
                      data_information_id=>$data_info->id,
                  });
	$i += $len;
      }
  }

sub initialize_genome_db
  {
    my ($org_name, $org_desc, $data_source_name, $data_source_desc, $data_source_link, $data_info_name, $data_info_desc, $data_info_link, $version) = @_;
    my $genome = Comp_Genomics::Genome->new();
    my $organism = $genome->get_organism_obj->find_or_create({name=>$org_name, #this adds "Arabidopsis thaliana"
							      description=>$org_desc, #this adds the complete kingdom to species list
							     }); 
    my $data_source = $genome->get_data_source_obj->find_or_create({name=>$data_source_name,
								    description=>$data_source_desc,
								    link=>$data_source_link,
								   });
    my ($file_name) = $file=~ /([^\/]*$)/;
    my $data_information = $genome->get_data_information_obj->find_or_create({name=>$data_info_name,
									      description=>$data_info_desc,
									      link=>$data_info_link,
									      data_source_id=>$data_source->id,
									      version=>$version,
									     });
    return ($genome, $organism, $data_source, $data_information);

  }

sub empty_db
  {
    my $data_information = shift;
    $data_information->delete;
    exit;
  }


sub show_help
  {
    print qq{
Welcome to $0!

This program loads TIGR XML data into Berkeley's department of BSBC
arabidopsis CNS database.


Usage:
$0 --file <xml_file> --chr <chromosome_number>


Options:

--file           Name of XML file to parse and load

--verbose|v      Verbose mode, prints out more messages of data loading process

--help|h         Prints out this help page

--chromosome|chr Chromosome to which the data belongs

--empty_db       if set, database is emptied of source and program exits.

--version        version of data
};
exit;
  }

