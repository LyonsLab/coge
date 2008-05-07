#!/usr/bin/perl 

## Written by Eric Lyons, UC Berkeley 2004
## Contact:  elyons@nature.berkeley.edu

use lib "/home/elyons/projects/TIGR_XML/parser/TIGRXMLLIB/";
#use TIGR_XML_parser;
#use Gene_obj;
#use CNS::TIGRXMLLIB::BSBC::Feature;
#use CNS::TIGRXMLLIB::BSBC::FeatureDB;
use XML::Simple;
use strict;
use Getopt::Long;
use Data::Dumper;
use DBI;
use Carp;
use CoGeX;


#use CoGe::Genome;

my ($file, $verbose, $help, $chr, $empty_db, $org_name, $org_desc, $org_id, $data_source_name, $data_source_desc, $data_source_link, $data_source_id, $dataset_id, $dataset_name, $dataset_desc, $dataset_link, $dataset_version, $coge, $add_genomic_seq, $GO);

GetOptions('file|f=s'=>\$file,
	   'verbose|debug'=>\$verbose,
	   'help|h' => \$help,
	   'chromosome|chr|c=s' => \$chr,
	   'empty_db|e' =>\$empty_db,
	   'version|ds_version|v=s' => \$dataset_version,
	   'org_name|on=s' => \$org_name,
	   'org_desc|od=s' => \$org_desc,
	   'org_id|oid=s'=>\$org_id,
	   'data_source_name|dsn=s' => \$data_source_name,
	   'data_source_desc|dsd=s' => \$data_source_desc,
	   'data_source_link|dsl=s' => \$data_source_link,
	   'data_source_id|dsid=s' => \$data_source_id,
	   'dataset_name|din=s' => \$dataset_name,
	   'dataset_desc|did=s' => \$dataset_desc,
	   'dataset_link|dil=s' => \$dataset_link,
	   'dataset_id|did=s'=>\$dataset_id,
	   'add_genomic_seq' => \$add_genomic_seq,
	   'go'=>\$GO,
	  );

$org_name = "Arabidopsis thaliana" unless $org_name;
$org_desc = "" unless $org_desc;
$data_source_name = "TAIR" unless $data_source_name;
$data_source_desc = "The Arabidopsis Information Resource" unless $data_source_desc;
$data_source_link = "www.arabidopsis.org" unless $data_source_link;
($dataset_name) = $file=~ /([^\/]*$)/ unless $dataset_name;
$dataset_desc = "TAIR Version $dataset_version Chromosome $chr, TIGR XML format" unless $dataset_desc;
$dataset_link = "ftp://ftp.arabidopsis.org/home/tair/Genes/TAIR8_genome_release/Tair8_XML/" unless $dataset_link;


unless (-r $file)
  {
    $help = 1;
    print "WARNING: You need a valid XML file.\n";
  }

unless ($dataset_version)
  {
    $help = 1;
    print "WARNING: You need to specify a dataset version!\n\n";
  }


my $xml = XMLin($file);
$chr = $xml->{PSEUDOCHROMOSOME}{ASSEMBLY}{CHROMOSOME} unless $chr;
$chr =~ s/^0+//; #remove those preceeding 0

my $connstr = 'dbi:mysql:dbname=genomes;host=biocon;port=3306';
my $coge = CoGeX->connect($connstr, 'cnssys', 'CnS' );

my $organism = process_header ($xml->{PSEUDOCHROMOSOME}{ASSEMBLY}{HEADER}) if $xml->{PSEUDOCHROMOSOME}{ASSEMBLY}{HEADER};

#my $organism = get_organism(name=>$org_name, desc=>$org_desc, id=>$org_id);
my $data_source = get_data_source(name=>$data_source_name, desc=>$data_source_desc, link=>$data_source_link, id=>$data_source_id);


my $genomic_seq_len = 10000;


show_help() if $help; 

my $pwd = `pwd`;
chomp $pwd;
unless ($file =~ /\//) {
    #full path not specified.
    $file = "$pwd/$file";
}

my ($dataset);
$dataset = $coge->resultset('Dataset')->find($dataset_id) if $dataset_id;
$dataset = $coge->resultset('Dataset')->find_or_create(
							    {
							     name                => $dataset_name,
							     description         => $dataset_desc,
							     link                => $dataset_link,
							     organism_id         => $organism->id,
							     data_source_id      => $data_source->id(),
							     version             => $dataset_version,
							    })  if $GO && !$dataset;
empty_db($dataset) if $empty_db;


#Dump of $xml->{PSEUDOCHROMOSOME}{ASSEMBLY}
#$VAR1 = [
#          'DATABASE',
#          'ASMBL_ID',
#          'COORDSET',
#          'CHROMOSOME', ->chr number
#          'CURRENT_DATE',
#          'CLONE_ID',
#          'GENE_LIST', ->important
#          'MISC_INFO',
#          'ASSEMBLY_SEQUENCE',  -> main sequence
#          'HEADER' -> organims name and lineage desc
#        ];

load_genomic_sequence($xml->{PSEUDOCHROMOSOME}{ASSEMBLY}{ASSEMBLY_SEQUENCE}, $genomic_seq_len) if $add_genomic_seq && $xml->{PSEUDOCHROMOSOME}{ASSEMBLY}{ASSEMBLY_SEQUENCE};

#print Dumper $xml;
my $gene_list = $xml->{PSEUDOCHROMOSOME}{ASSEMBLY}{GENE_LIST} ? $xml->{PSEUDOCHROMOSOME}{ASSEMBLY}{GENE_LIST} : $xml->{ASSEMBLY}{GENE_LIST};
#print Dumper $gene_list;
process_genes($gene_list);

exit;

sub process_genes
  {
    my $genes = shift;
    my $feature_count = 0;
    my $gene_count = 0;
    my $public_locus_count = 0;
    my %count;
    foreach my $type (keys %$genes)
      {
	print "Process genes of type $type\n";
#	next unless $type =~ /RNA_GENES/; #tmp for time
	foreach my $sub_type (keys %{$genes->{$type}})
	  {
	    print "\tProcess subtype $sub_type\n";
	    foreach my $item (@{$genes->{$type}{$sub_type}})
	      {
#		if ($sub_type eq "TU" || $sub_type eq "TRANSPOSABLE_ELEMENT_GENE")
		  {

		    $feature_count++;
		    $gene_count++;
		    $count{$type}{$sub_type}++;
		    my $gene = "gene";
#		    $gene = "Transposable Element" if $sub_type eq "TRANSPOSABLE_ELEMENT_GENE";
		    $gene = "pseudo".$gene if $item->{GENE_INFO}{IS_PSEUDOGENE};
		    
		    my $models;
		    $models = $item->{MODEL} if ref ($item->{MODEL}) =~ /array/i;
		    push @$models, $item->{MODEL} if ref ($item->{MODEL}) =~ /hash/i;
		    my %go; #gene ontologies
		    if ($item->{GENE_INFO}{GENE_ONTOLOGY})
		      {
			my $go_stuff;
			if ($item->{GENE_INFO}{GENE_ONTOLOGY} && ref ($item->{GENE_INFO}{GENE_ONTOLOGY}{GO_ID}) =~ /array/i)
			  {
			    $go_stuff = $item->{GENE_INFO}{GENE_ONTOLOGY}{GO_ID};
			  }
			elsif ($item->{GENE_INFO}{GENE_ONTOLOGY})
			  {
			    push @$go_stuff, $item->{GENE_INFO}{GENE_ONTOLOGY}{GO_ID} if $item->{GENE_INFO}{GENE_ONTOLOGY}{GO_ID};
			  }
			foreach my $go_item (@$go_stuff)
			  {
			    $go{$go_item->{ASSIGNMENT}} = {
							   goid=>$go_item->{ASSIGNMENT}=~/(\d+)/,
							   goterm=>$go_item->{GO_TERM},
							   gocat=>$go_item->{GO_TYPE},
							  };
			  }
		      }

		    foreach my $model (@$models)
		      {
			my %names;
			my %anno;
			my %locs;
			
			$names{$item->{FEAT_NAME}}=1 if $item->{FEAT_NAME};
			$names{$item->{GENE_INFO}{PUB_LOCUS}}=2 if $item->{GENE_INFO}{PUB_LOCUS};
			$names{$item->{GENE_INFO}{LOCUS}}=2 if $item->{GENE_INFO}{LOCUS};
			$names{$model->{FEAT_NAME}}=1 if $model->{FEAT_NAME};
			$names{$model->{PUB_LOCUS}}=1 if $model->{PUB_LOCUS};

			$anno{$item->{COMMENT}}=1 if $item->{COMMENT};
			$anno{$item->{GENE_INFO}{PUB_COMMENT}} =1 if $item->{GENE_INFO}{PUB_COMMENT};
			if (ref ($item->{GENE_INFO}{COM_NAME}) =~ /array/i)
			  {
			    foreach my $cname (@{$item->{GENE_INFO}{COM_NAME}})
			      {
				$anno{$cname->{content}}=1 if $cname->{content};
			      }
			  }
			else
			  {
			    $anno{$item->{GENE_INFO}{COM_NAME}{content}}=1;
			  }
			my $prot_seq = $model->{PROTEIN_SEQUENCE};
			$prot_seq =~ s/\*$//;
			$gene = "pseudogene" if $prot_seq =~ /\*/; #found that genes with nonsense mutations were not labeled as pseudogenes in some datasets
			my $strand = $item->{COORDSET}{END5} < $item->{COORDSET}{END3} ? 1 : -1;
			$locs{$gene} = [[$item->{COORDSET}{END5}, $item->{COORDSET}{END3}, $strand, $chr]];
			
			if (ref ($model->{EXON}) =~ /array/i)
			  {
			    foreach my $exon (@{$model->{EXON}})
			      {
				push @{$locs{mRNA}}, [$exon->{COORDSET}{END5},$exon->{COORDSET}{END3}, $strand, $chr] if $exon->{COORDSET};
				push @{$locs{CDS}}, [$exon->{CDS}{COORDSET}{END5},$exon->{CDS}{COORDSET}{END3}, $strand, $chr] if $exon->{CDS};
			      }
			  }
			elsif ($model->{EXON})
			  {
			    my $exon = $model->{EXON};
			    push @{$locs{mRNA}}, [$exon->{COORDSET}{END5},$exon->{COORDSET}{END3}, $strand, $chr] if $exon->{COORDSET};
			    push @{$locs{CDS}}, [$exon->{CDS}{COORDSET}{END5},$exon->{CDS}{COORDSET}{END3}, $strand, $chr] if $exon->{CDS};
			  }
			if ($sub_type eq "PRE-TRNA")
			  {
			    $names{$item->{TRNA}{FEAT_NAME}}=1 if $item->{TRNA} && $item->{TRNA}{FEAT_NAME};
			    $names{$item->{TRNA}{'RNA-EXON'}{FEAT_NAME}}=1 if $item->{TRNA} && $item->{TRNA}{'RNA-EXON'} && $item->{TRNA}{'RNA-EXON'}{FEAT_NAME};
			    if ($item->{GENE_INFO}{PUB_COMMENT} =~ /anticodon: (\w+)/)
			      {
				my $anticodon = uc($1);
				$anno{anticodon} = $anticodon;
				my $codon = reverse $anticodon;
				$codon =~ tr/ATCG/TAGC/;
				$anno{codon}=$codon;
			      }
			    $locs{"pre-tRNA"} = $locs{mRNA};
			    delete $locs{mRNA};
			  }
			elsif($sub_type eq "TRANSPOSABLE_ELEMENT_GENE")
			  {
			    $locs{"Transposable Element"} = $locs{mRNA};
			    delete $locs{mRNA};
			  }
			elsif($sub_type eq "OTHER_RNA")
			  {
			    $locs{"Other RNA"} = $locs{mRNA};
			    delete $locs{mRNA};
			  }
			elsif($sub_type eq "MIRNA")
			  {
			    $locs{"miRNA"} = $locs{mRNA};
			    delete $locs{mRNA};
			  }
			elsif($sub_type eq "SNORNA")
			  {
			    $locs{"snoRNA"} = $locs{mRNA};
			    delete $locs{mRNA};
			  }
			elsif($sub_type eq "SNRNA")
			  {
			    $locs{"snRNA"} = $locs{mRNA};
			    delete $locs{mRNA};
			  }
			elsif($sub_type eq "RRNA")
			  {
			    $locs{"rRNA"} = $locs{mRNA};
			    delete $locs{mRNA};
			  }
			else
			  {
			    unless ($sub_type eq "TU")
			      {
				print "WARNING:  unknown sub_type $sub_type:\n";
				print Dumper $item, \%names, \%anno, \%go, \%locs;
				next;
			      }
			    
			  }
			load_info(
				  locs=>\%locs,
				  annos=>\%anno,
				  go=>\%go,
				  names=>\%names,
				  prot_seq=>$prot_seq,
				 );

		      }
		  }
		}
	  }
      }

    print "Total processed genes: $gene_count\n";
    print "  genes with a public locus:  $public_locus_count\n";
    print "Total processed entries (a gene may have multiple isoforms): $feature_count\n";
    print "Types processed:";
    print Dumper (\%count);
  }

# foreach my $gene ($TIGRparser->get_genes()) {
#   $gene->refine_gene_object();
#   $gene->create_all_sequence_types(\$TIGRparser->{assembly_seq}) if $TIGRparser->{assembly_seq};
#   print Dumper $gene;
#   next;
#   my $entry = $gene->toString;
# #  my $locus = new CNS::TIGRXMLLIB::BSBC::Feature(verbose=>$verbose);
#   $gene_count++;
#   $public_locus_count++ if $gene->{'pub_locus'};
#   foreach my $item (split /ISOFORM:/, $entry)
#     {
#       print $item,"\n";
#       my $locus = gen_locus($item);

#       next;
# #      $locus->process_entry($gene, $item);
#       $count{$locus->gene_type}++;
#       $feature_count++;
#       my @feats;
#       #gene feature
#       my $gene_feat = $locus->pseudogene? create_feature("pseudogene"): create_feature("gene");
#       next unless load_locations($gene_feat, $locus->gene_location, $locus->strand, $chr);
#       load_xref($gene_feat, $locus->supporting_cdna());
#       push @feats, $gene_feat;
#       if  ($locus->gene_type() =~ /rna/i)
# 	{
# 	  my $type = $locus->gene_type();
# 	  $type =~ s/rna/RNA/;
# 	  my $rna_feat = create_feature($type);
# 	  my $locs = load_locations($rna_feat, $locus->exon_locations, $locus->strand, $chr);
# 	  push @feats, $rna_feat if $locs;
# 	}
#       if ($locus->gene_type() eq "protein-coding")
#    	{
# 	  #mRNA feature
# 	  if ($locus->exon_locations) #untested!
# 	    {
# 	      my $mrna_feat = create_feature("mRNA");
# 	      my $locs = load_locations($mrna_feat, $locus->exon_locations, $locus->strand, $chr);
# 	      push @feats, $mrna_feat if $locs;
# 	    }
#    	  #CDS feature
# 	  if ($locus->CDS_locations) #untested
# 	    {
# 	      my $cds_feat = create_feature("CDS");
# 	      my $locs = load_locations($cds_feat, $locus->CDS_locations, $locus->strand, $chr);
# 	      load_protein($cds_feat, $locus->protein_seq()) if $locus->protein_seq();
# 	      push @feats, $cds_feat if $locs;
# 	    }
# 	}
#       #names
#       foreach my $feat (@feats)
# 	{
# 	  next if ref ($feat) =~ /Deleted/;
# 	  #names
# 	  foreach my $name (
# 			    $locus->public_locus(),
# 			    $locus->public_locus_model(),
# 			    $locus->TU_feat_name(),
# 			    $locus->model_feat_name(),
# 			    $locus->alt_locus(),
# 			    $locus->locus(),
# 			    $locus->gene_synonyms(),
# 			    split/;\s*/, $locus->common_name(),
# 			   )
# 	    {
# 	      next unless $name;
# 	      load_name($feat, $name);
# 	    }
# 	  #GO
# 	  load_go($feat, $locus->GO);
# 	  #additional annotations
# 	  load_annotation($feat, "secondary products", $locus->secondary_products()) if $locus->secondary_products;
# 	  load_annotation($feat, "public comment", $locus->public_comment()) if $locus->public_comment;
	  
# 	}
#     }
# }


sub load_info
    {
      my %opts = @_;
      my $locs = $opts{locs};
      my $annos = $opts{annos};
      my $go = $opts{go};
      my $names = $opts{names};
      my $prot_seq = $opts{prot_seq};
      foreach my $type (keys %$locs)
	{
	  my $locs = $locs->{$type};
	  my $feat = create_feature($type);
	  
	  load_locations($feat, $locs);
	  foreach my $name (keys %$names)
	    {
	      my $primary = $names->{$name} == 2 ? 1 : 0; #"2" is used to designate primary feature names;
	      load_name(feat=>$feat, 
			name=>$name,
			primary=>$primary);
	    }
	  foreach my $anno (keys %$annos)
	    {
	      if (ref ($annos->{$anno}) =~ /hash/)
		{
		  my $type = $anno;
		  foreach my $item (map {split/;/,$_} keys %{$annos->{$type}})
		    {
		      $item =~ s/^\s+//;
		      $item =~ s/\s+$//;
		      next unless $item;
		      load_annotation($feat, $item, $type);
		    }
		}
	      else
		{
		  load_annotation($feat, $anno);
		}
	    }
	  if ($type eq "CDS" && $prot_seq)
	    {
	      load_protein($feat, $prot_seq);
	    }
	  load_go ($feat, [values %$go]) if ($go && keys %$go);
	}
    }

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
	load_annotation($feat, $go->{goterm}, $go->{goid}, "GO ".$go->{gocat});
      }
  }

sub load_annotation
  {
    my ($feat, $annotation, $type, $group) = @_;
#    print Dumper \@_;
    $type = "annotation" unless $type;
    return 0 unless $annotation && $type;
    my %type_opts = (name=>$type);
    if ($group)
      {
	my $anno_type_group = $coge->resultset('AnnotationTypeGroup')->find_or_create({name=>$group}); 
	$type_opts{annotation_type_group_id}=$anno_type_group->id() if $GO;
	  
      }
    my $anno_type = $coge->resultset('AnnotationType')->find_or_create({
									%type_opts
								       }) if $GO;
    $feat->add_to_annotations({
			       annotation=>$annotation,
			       annotation_type_id=>$anno_type->id(),
			      }) if $GO;
  }

sub load_protein
  {
    my ($feat, $seq) = @_;
    return 0 unless $seq;
    my $seq_type = $coge->resultset('SequenceType')->find_or_create({
								     name=>"protein",
								     description=>"translation",
								    }) if $GO;
    my $sequence = $feat->add_to_sequences({
					    sequence_type_id=>$seq_type->id(),
					    sequence_data=>$seq,
					   }) if $GO;
  }

sub load_name
  {
#    my ($feat, $name, $desc) = @_;
    my %opts = @_;
    my $feat = $opts{feat};
    my $name = $opts{name};
    my $desc = $opts{desc};
    my $primary = $opts{primary};
    return 0 unless $name;
    my $feat_name = $feat->add_to_feature_names({
						 name=>$name,
						 description=>$desc,
						 primary_name=>$primary,
						}) if $GO;
  }

sub load_locations
  {
    my ($feat, $locs) = @_;
    return 0 unless ref ($locs) =~ /ARRAY/i;
    my @locs;
    my $start;
    my $stop;
    my $strand;
    my $chr;
    foreach my $loc (@$locs)
      {
	my ($st, $sp);
	($st, $sp, $strand, $chr) = @$loc;
	$strand = 1 if $strand =~ /\+/;
	$strand = -1 if $strand =~ /-/;
	return unless $chr && $strand;
	($st, $sp) = ($sp, $st) if $sp < $st;
	$start = $st unless $start;
	$start = $st if $st < $start;
	$stop = $sp unless $stop;
	$stop = $sp if $sp > $stop;
	my $location = $feat->add_to_locations({
						start=>$st,
						stop=> $sp,
						strand=>$strand, #expect this to be "1" or "-1"
						chromosome=>$chr, #make sure that chromosomes are in real numbers (not roman numerals [IX])
					       }) if $GO;
	push @locs, $location;
      }
    if (@locs && $GO)
      {
	$feat->start($start);
	$feat->stop($stop);
	$feat->strand($strand);
	$feat->chromosome($chr);
	$feat->update;
	return \@locs;
      }
    elsif ($GO)
      { #no locations, no feature
	$feat->delete;
	return 0;
      }
  }


sub create_feature
  {
    my ($type) = @_;
    return 0 unless $type;
    $type = "tRNA" if $type eq "TRNA";
    my $feat_type = $coge->resultset('FeatureType')->find_or_create({name=>$type}) if $GO;

    my $db_feature = $coge->resultset('Feature')->create(
							 {feature_type_id=>$feat_type->id,
							  dataset_id=>$dataset->id,
							 }) if $GO;
    return $db_feature;
  }

sub load_genomic_sequence
  {
    my ($seq, $len) = @_;
    $seq =~ s/\n//g;
    $seq =~ s/\r//g;
    $seq =~ s/\s+//g;
    my $seqlen = length $seq;
    
    print STDERR "Loading genomic sequence ($seqlen nt)\n";
    $seq =~ s/\s//g;
    $seq =~ s/\n//g;
    my $i = 0;
    while ($i < $seqlen)
      {
	my $str = substr($seq,$i,$len);
	my $start = $i+1;
	my $stop = $i + length $str;
	$dataset->add_to_genomic_sequences({start => $start,
					    stop => $stop,
					    chromosome => $chr,
					    sequence_data => $str,
					   }) if $str && $GO;
	$i += $len;
      }
  }

sub empty_db
  {
    my $data_information = shift;
    print "Clearning database of ", $data_information->name,"\n";
    $data_information->delete;
    exit;
  }

sub process_header
    {
      my $header = shift;
      my $name = $header->{ORGANISM};
      my $desc = $header->{LINEAGE};
      return get_organism(name=>$name, desc=>$desc);
    }

sub get_organism
{
  my %opts = @_;
  my $name = $opts{name};
  my $desc = $opts{desc};
  my $id = $opts{id};
  return $coge->resultset('Organism')->find($id) if $id;
  return $coge->resultset('Organism')->find_or_create(
                                                   {
                                                    name=>$name,
                                                    description=>$desc,
                                                   }) if $GO;
}

sub get_data_source
{
  my %opts = @_;
  my $name = $opts{name};
  my $desc = $opts{desc};
  my $link = $opts{link};
  my $id = $opts{id};
  return $coge->resultset('DataSource')->find($id) if $id;
  return $coge->resultset('DataSource')->find_or_create(
                                                     {
                                                      name=>$name,
                                                      description=>$desc,
                                                      link=>$link,
                                                     }) if $GO;
}


sub gen_locus
  {
    my $item = shift;
    my %locus;
    my $loc =0;
    
    foreach my $line(split/\n/,$item)
      {
#	$loc = 0 if $
	if ($loc)
	  {
	  }
	else
	  {
	    my ($type, $stuff) = split /:\s*/, $line,2;
	    $type =~ s/^\s+//;
	    $type =~ s/\s+$//;
	    $stuff =~ s/^\s+//;
	    $stuff =~ s/\s+$//;
	    if ($type =~ /_coords/)
	      {
		
	      }
	    $locus{$type}{$stuff}=1;
	    $loc = 0;
	  }
      }
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

