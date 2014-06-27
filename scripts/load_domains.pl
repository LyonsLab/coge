#!/usr/bin/perl -w

use strict;
use CoGe::Genome;
use Data::Dumper;
use Getopt::Long;

my ($file, $verbose, $help, $version);

GetOptions('file|f=s'=>\$file,
	   'verboes|debug'=>\$verbose,
	   'help|h'=>\$help,
	   'version|ver|v=i'=>\$version,
	  );

$help = 1 unless $file && $version;
show_help() if $help;

my $db = CoGe::Genome->new;
my ($data_source) = $db->get_data_source_obj->search({data_source_id=>1});
my ($org) = $db->get_org_obj->search({organism_id=>1});
my $data_info = $db->get_data_info_obj->find_or_create({data_source_id=>$data_source->id, organism_id=>$org->id, name=>"all.domains.txt", description=>"TAIR version 6 list of protein domains determined by a variety of algorithms.", version=>6, link=>"ftp://ftp.arabidopsis.org/home/tair/Proteins/Domains/all.domains.txt"});
my $data = process_file($file);
my $feat_type = $db->get_feature_type_obj->find_or_create({name=>"Functional Domains", description=>"Protein domains detected by a variety of algorithms"});
process_data($data);
#$feat_type->delete;
sub process_domains
  {
    my %opts = @_;
    my $feat = $opts{feat};
    my $doms = $opts{doms};
    foreach my $dom_id (keys %$doms)
      {
	my @locs;
	my ($app, $d_id, $d_name, $ip_id, $ip_name, $eval);
	foreach my $dom (@{$doms->{$dom_id}})
	  {
	    $app = $dom->{application};
	    $d_id = $dom->{domain_id};
	    $d_name = $dom->{domain_desc};
	    $eval = $dom->{eval};
	    $ip_id = $dom->{interpro_id};
	    $ip_name = $dom->{interpro_name};
	    foreach my $seq ($feat->sequences)
	      {
		my $dstart = $seq->get_genomic_position($dom->{start});
		my $dstop = $seq->get_genomic_position($dom->{stop});
		push @locs, $seq->get_genomic_locations(start=>$dom->{start}, stop=>$dom->{stop});
	      }
	  }
	my $new_feat = $db->get_feat_obj->create({feature_type_id=>$feat_type->id(), data_information_id=>$data_info->id});
	my $atg = $db->get_anno_type_group_obj->find_or_create({name=>$app, description=>"Domain finding algorithm"});
	my $at = $db->get_anno_type_obj->find_or_create({annotation_type_group_id=>$atg->id, name=>$ip_id, description=>$ip_name});
	my $a = $db->get_anno_obj->find_or_create({annotation_type_id=>$at->id, feature_id=>$new_feat->id, annotation=>$ip_name." ($d_id, eval: $eval)"});
	$db->get_feat_name_obj->find_or_create({feature_id=>$new_feat->id, name=>$ip_name, description=>$ip_id});
	$db->get_feat_name_obj->find_or_create({feature_id=>$new_feat->id, name=>$ip_id, description=>$ip_name});
	$db->get_feat_name_obj->find_or_create({feature_id=>$new_feat->id, name=>$d_name, description=>$d_id});
	$db->get_feat_name_obj->find_or_create({feature_id=>$new_feat->id, name=>$d_id, description=>$d_name});
	foreach my $loc (@locs)
	  {
	    $db->get_loc_obj->create({feature_id=>$new_feat->id,
				      start => $loc->start,
				      stop => $loc->stop,
				      strand => $loc->strand,
				      chromosome=> $loc->chr,
				     });
	    delete $loc->{__Changed};
	  }

      }
  }

sub process_data
  {
    my $data = shift;
    my $total = scalar keys %$data;
    my $i = 0;
    foreach my $accn (keys %$data)
      {
	$i++;
	print $i,"/",$total,":\t", $accn,"\n";
	foreach my $feat ($db->get_features_by_name_and_version(name=>$accn, version=>$version))
	  {
	    next unless $feat->type->name =~ /CDS/i;
	    process_domains(feat=>$feat, doms=>$data->{$accn});
	  }
      }
  }

sub process_file
  {
    my $file = shift;
    my %data;
    open (IN, $file) || die "Can't open $file for reading: $!";
    while (<IN>)
      {
	chomp;
	my @line = split /\t/;
	next unless $line[10];
	next if $line[10] =~ /NULL/;
#	next unless $line[2] =~ /HMM/;
	push @{$data{$line[0]}{$line[3]}}, {
				  application=>$line[2],
				  domain_id=>$line[3],
				  domain_desc=>$line[4],
				  start=>$line[5],
				  stop=>$line[6],
				  eval=>$line[7],
				  interpro_id=>$line[9],
				  interpro_name=>$line[10],
				 }
      }
    close IN;
    return \%data;
  }

sub show_help
  {
    print qq{
Welcome to $0!

This program loads tab-delimited domain information.  It was originally designed for the domain information file from arabidopsis.org for version 6 of the arabidopsis genome.

--file|f          Name of tab delimited file with domain data

--verbose|debug   Verbose and debugging mode

--version|ver|v   Version of the data which to add the domain features

--help            prints this page
};
    exit;
  }
