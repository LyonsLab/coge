#!/usr/bin/perl -w
use strict;
use Comp_Genomics::Genome;

my $chr = shift;
my $file = shift;

my $org_name = "Arabidopsis thaliana";
my $org_desc = "";
my $data_source_name = "TAIR";
my $data_source_desc = "The Arabidopsis Information Resource";
my $data_source_link = "www.arabidopsis.org";
my ($data_info_name) = $file=~ /([^\/]*$)/;
my $data_info_desc = "TAIR Version 6 Chromosome $chr, TIGR XML format";
my $data_info_link = "ftp.arabidopsis.org/home/tair/Genes/TAIR6_genome_release/tair6_xml/";
my $version = 6;

my ($gdb, $org, $data_source, $data_info) = initialize_genome_db($org_name, $org_desc, $data_source_name, $data_source_desc, $data_source_link, $data_info_name, $data_info_desc, $data_info_link, $version);

my $gso = $gdb->get_genomic_sequence_obj();

my $len = 10000;
while (<>)
  {
    chomp;
    s/<.*?>//g;
    s/\s+//g;
    next unless length $_ > $len;
    my $i = 0;
    while ($i < length $_)
      {
	my $str = substr($_,$i,$len);
#	print $str,"\n";
#	next;
	my $start = $i+1;
	my $stop = $i + length $str;
	$gso->create({start=>$start,
                      stop=>$stop,
                      chromosome=>$chr,
                      sequence_data=>$str,
                      organism_id=>$org->id,
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
