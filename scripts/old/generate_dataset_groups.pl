#!/usr/bin/perl -w

use strict;
use Data::Dumper;
use CoGeX;

my $coge = CoGeX->dbconnect();
my $DEBUG= 0;

print join ("\t", qw{#VERSION ORGANISM_ID GENOMIC_SEQUENCE_TYPE_ID DATASET_IDS}),"\n";
foreach my $org ($coge->resultset('Organism')->all())
#foreach my $org ($coge->resultset('Organism')->search({name=>{like=>'%Gluconacetobacter diazotrophicus PAl 5%'}}))
#foreach my $org ($coge->resultset('Organism')->search({name=>{like=>'%coli%'}}))
  {
    print "#",$org->name,"\n" if $DEBUG;
    my %chr_data;
    my %contig_data;
    my $contig_flag = 0;
    dataset: foreach my $ds ($org->datasets)
      {
#	next unless $ds->sequence_type->name =~ /masked/i;
	foreach my $chr ($ds->chromosomes)
	  {
	    my $ver = $ds->version;
	    print "#WARNING: ",$ds->name," does not have a version\n" unless $ver;
	    if ($chr =~ /contig/i || $chr =~ /scaffold/i || $chr =~ /super/i || $org->id == 333) #if this is set, we actually don't need to go through the rest of this stuff
	      #funky orgs:
	      #  333 => zea mayze
	      {
		push @{$contig_data{$ds->datasource->name}{$ds->sequence_type->name}{$ds->version}},$ds;
		$contig_flag = 1;
		next dataset;
	      }
	    if ($chr_data{$ds->datasource->name}{$ds->sequence_type->name}{$chr}{$ver})
	      {
		print "#WARNING CHROMOSOME EXISTS FOR THIS DATA SOURCE AT THIS VERSION, ELEVATING VERSION:  current: ", $ds->name,"  previous:  ", $chr_data{$ds->datasource->name}{$ds->sequence_type->name}{$chr}{$ver}->name,"\n";
		my $tmp = $chr_data{$ds->datasource->name}{$ds->sequence_type->name}{$chr};
		my ($tmp_ver) = sort {$b<=>$a} keys %$tmp;
		$tmp_ver++;
		$ver = $tmp_ver;
	      }
	    $chr_data{$ds->datasource->name}{$ds->sequence_type->name}{$chr}{$ver}=$ds;
	  }
      }
    dump_contig(\%contig_data) if keys %contig_data;
    dump_data(\%chr_data) if keys %chr_data;
    #    print Dumper \%data;
#    <STDIN>;
  }

sub dump_data
  {
    my $data = shift;
    foreach my $source (keys %$data)
      {
	foreach my $seq_type (keys %{$data->{$source}})
	  {
	    my %chrs = map {$_=>[]} keys %{$data->{$source}{$seq_type}};
	    my $ver_count;
	    foreach my $chr (keys %{$data->{$source}{$seq_type}})
	      {
		my $tmp_ver_count;
		foreach my $ver (sort {$a <=> $b} keys %{$data->{$source}{$seq_type}{$chr}})
		  {
		    push @{$chrs{$chr}},$ver;
		    $tmp_ver_count++;
		  }
		$ver_count = $tmp_ver_count unless $ver_count;
		$ver_count = $tmp_ver_count if $tmp_ver_count > $ver_count;
	      }
	    for (my $i=0; $i<$ver_count; $i++)
	      {
		print "#",$source,"\t",$seq_type,"\n" if $DEBUG;
		print "#",join ("\t",qw{CHR VER}),"\n" if $DEBUG;
		my @ds;
		foreach my $chr (sort keys %chrs)
		  {
		    my $ver;
		    $ver = $chrs{$chr}[$i] if defined $chrs{$chr}[$i];
		    $ver = $chrs{$chr}[-1] unless defined $ver;
		    print "#",$chr,"\t",$ver,"\n" if $DEBUG;
		    unless ($data->{$source}{$seq_type}{$chr}{$ver})
		      {
			print "#ERROR $source $seq_type $chr $ver\n";
		      }
		    push @ds, $data->{$source}{$seq_type}{$chr}{$ver}; #we need to get our dataset group
		  }
		print "#","-"x20,"\n" if $DEBUG;
		if (@ds)
		  {
		    @ds = sort { $b->version <=> $a->version } @ds;
		    my %dsids = map {$_->id,1}@ds;
#		    print join ("\n", map{$_->name} @ds),"\n";
		    print join ("\t", $ds[0]->version, $ds[0]->organism->id, $ds[0]->sequence_type->id, join ("::", keys %dsids)),"\n";
		  }
	      }
	  }
      }
  }

sub dump_contig
  {
#    $data->{$ds->datasource->name}{$ds->sequence_type->name}{$ds->version}{$ds->id}=$ds;
    my $data = shift;
    foreach my $source (keys %$data)
      {
	foreach my $type (keys %{$data->{$source}})
	  {
	    foreach my $ver (keys %{$data->{$source}{$type}})
	      {
		my $ds = $data->{$source}{$type}{$ver};
		print "#",join ("\t", (map {$_->name} @$ds), $ds->[0]->organism->name),"\n" if $DEBUG;
		print join ("\t", $ds->[0]->version, $ds->[0]->organism->id, $ds->[0]->sequence_type->id, join ("::",map{$_->id} @$ds )),"\n";
	      }
	  }
      }
  }
