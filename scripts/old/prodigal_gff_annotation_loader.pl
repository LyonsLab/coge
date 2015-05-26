#!/usr/bin/perl -w

use strict;
use Data::Dumper;
use CoGeX;
use Getopt::Long;
use URI::Escape;

my $GO = 0;
my ($dsid, $ds_name, $ds_desc, $ds_link, $ds_version, $source_name, $source_desc, $source_link, $source_id, $gff_file, $anno_file, $DEBUG);
my $add_gene =0;
my $add_cds =0;
my $add_type_to_name =0;
my @names;
my @skip_types;
my @anno_names;
my $connstr = 'dbi:mysql:dbname=coge;host=localhost;port=PORT';
my$coge = CoGeX->connect($connstr, 'USER', 'PASSWORD' );
#$coge->storage->debugobj(new DBIxProfiler());
#$coge->storage->debug(1);

#Phytozome example:
#scaffold_13     Ptrichocarpav2_0        mRNA    6802038 6807521 .       -       .       mRNA POPTR_0013s07870.1;Name POPTR_0013s07870.1;PACid 17328357
#scaffold_13     Ptrichocarpav2_0        5'-UTR  6807412 6807521 .       -       .       mRNA POPTR_0013s07870.1;PACid 17328357
#scaffold_13     Ptrichocarpav2_0        5'-UTR  6807257 6807320 .       -       .       mRNA POPTR_0013s07870.1;PACid 17328357
#scaffold_13     Ptrichocarpav2_0        5'-UTR  6806896 6806922 .       -       .       mRNA POPTR_0013s07870.1;PACid 17328357
#scaffold_13     Ptrichocarpav2_0        CDS     6806864 6806895 .       -       0       mRNA POPTR_0013s07870.1;PACid 17328357
#scaffold_13     Ptrichocarpav2_0        CDS     6805536 6805595 .       -       2       mRNA POPTR_0013s07870.1;PACid 17328357
#scaffold_13     Ptrichocarpav2_0        CDS     6805335 6805404 .       -       2       mRNA POPTR_0013s07870.1;PACid 17328357
#scaffold_13     Ptrichocarpav2_0        CDS     6804457 6804526 .       -       0       mRNA POPTR_0013s07870.1;PACid 17328357
#scaffold_13     Ptrichocarpav2_0        CDS     6803629 6803655 .       -       1       mRNA POPTR_0013s07870.1;PACid 17328357
#scaffold_13     Ptrichocarpav2_0        3'-UTR  6802038 6802490 .       -       .       mRNA POPTR_0013s07870.1;PACid 17328357
#scaffold_13     Ptrichocarpav2_0        CDS     6802491 6802540 .       -       1       mRNA POPTR_0013s07870.1;PACid 17328357

#names for features:  mRNA Name

#link to annotation page: http://www.phytozome.net/genePage.php?search=1&detail=1&crown&method=0&searchText=transcriptid%3A17335774
#last digits (17335774) are in GFF under PACid

#data transformations needed for CoGe:
# mRNA->gene
#5'-UTR->mRNA
#3'-UTR->mRNA
#CDS replicated to mRNA

#pain in my ass!
my $sgd_link = "http://www.yeastgenome.org/cgi-bin/locus.fpl?dbid=";

GetOptions (
	     "source_name=s" => \$source_name, # datasource
	     "source_desc=s" => \$source_desc,
	     "source_link=s" => \$source_link,
	     "source_id=s"   => \$source_id,
	    "ds_name=s" => \$ds_name,# datasetid
	    "ds_desc=s" => \$ds_desc,
	    "ds_link=s" => \$ds_link,
	    "ds_version=s" => \$ds_version,
	    "dsid=i" => \$dsid,
	     "go=s"    => \$GO,
	     "debug=s" => \$DEBUG,
	     "name=s" => \@names,
	     "anno_name=s" => \@anno_names,
	     "add_gene_feature" => \$add_gene,
	     "add_cds_feature" => \$add_cds,
	     "add_type_to_name"=>\$add_type_to_name, #adds type (column 2) to name
	     "skip_type=s"=>\@skip_types,
	    "gff_file=s"=>\$gff_file,
	    "anno_file=s"=>\$anno_file,
	   );
$DEBUG = 1 unless defined $DEBUG; #turn on by default

if ($source_name)
  {
    my $source = $coge->resultset("DataSource")->find_or_create({name=>$source_name,description=>$source_desc, link=>$source_link});
    $source_id = $source->id;
  }

my $ds = generate_ds(ds_name => $ds_name,
		     ds_desc => $ds_desc,
		     ds_link => $ds_link,
		     ds_version => $ds_version,
		     ds_id =>$dsid,
		     source_id=>$source_id,
		    );
#  $coge->resultset('Dataset')->find($dsid);

unless ($ds)
  {
    warn "unable to find or create a valid dataset entry";
    exit;
  }
print "Working on dataset: ", $ds->name. " (".$ds->id.")\n";
#some defaults to check for in names and annotations
push @names, "ID";
push @names, "Name";
push @names, "Alias";
push @names, "gene";
push @names, "Parent";
push @names, "Locus_id";
push @names, "ID_converter";
push @names, "Gene_symbols";
push @skip_types, "Link_to";
push @skip_types, "References";
push @skip_types, "Sequence_download";

push @anno_names, "Note";
push @anno_names, "NIAS_FLcDNA";
push @anno_names, "Comment";
push @anno_names, "GO";
push @anno_names, "ORF_evidence"; #can we link to SGD?
push @anno_names, "Transcript_evidence";
push @anno_names, "Status";
push @anno_names, "InterPro";
push @anno_names, "start_type";
push @anno_names, "rbs_motif";
push @anno_names, "rbs_spacer";
my %anno_names = map {$_,1} @anno_names if @anno_names;
my %check_names = map {$_,1} @names;
my %skip_types = map {$_,1} @skip_types;

warn "-go flag is not true, nothing will be added to the database.\n" unless $GO;
my %data;
my %annos;
my %anno_name_lookup;
my %feat_types; #store feature type objects
#my ($anno_type) = $coge->resultset('AnnotationType')->find_or_create({name=>"Phytozome link"}); #generic annotation type
my $prev_type;
my %master_names;

#add information from annotation_file
#process_annotation_file(file=>$anno_file, annos=>\%annos, names=>\%anno_name_lookup) if $anno_file;
open (IN, $gff_file);

my %seen_types;
my %added_types;
my %seen_chr;
while (<IN>)
  {
    next if /^#/;
    chomp;
    next unless $_;
    my @line = split /\t/;
    next if $line[2] eq "clone";
    next if $line[2] eq "intron";
    next if $line[2] eq "chromosome";
    my $chr;
    $chr = $line[0];
    $chr =~ s/Prodigal_Seq_//;
    $chr =~ s/%.*//;
    $chr =~ s/chromosome//i;
    $chr =~ s/chr//i;
    $chr =~ s/^_//i;
    $chr =~ s/^0//g;

    ($chr) = split /\s+/,$chr;

    my %names;
    my $name;
    foreach my $item (split /;/, $line[-1])
      {
	my $tmp;
	$item =~ s/"//g;
	$item =~ s/^\s+//;
	$item =~ s/\s+$//;
	my ($type, $info) = $item =~ /=/ ? split (/=/,$item,2) : (split / /,$item,2);
	$seen_types{$type}++;
	next if $skip_types{$type};
	if ($check_names{$type})
	  {
	    foreach my $item (split /,/, $info)
	      {
		$names{$item} =1;
		$name = "Prodigal_prediction_".$item unless $name;
		if ($item =~ /\.\d+$/)
		  {
		    my $tmp = $item;
		    $tmp =~ s/\.\d+$//;
		    $names{$tmp}=1;
		  }
	      }
	  }
	next unless $name; #No name, don't know what to do!
	$info = uri_unescape($info); #remove URL formatting
	$annos{$name}{$info}={type=>$type} if $anno_names{$type};
	#add phytozome link:
	if ($type eq "dbxref")
	  {
	    $info =~ s/SGD://;
	    $annos{$name}{"View annotation at SGD"}={link=>$sgd_link."$info", type=>"SGD link"}
	  }
      }
    next unless $name; #No name, don't know what to do!

    #add in additional names from the annotation file
    foreach my $tmp (keys %names)
      {
	next unless $anno_name_lookup{$tmp};
	foreach my $tmp2 (keys %{$anno_name_lookup{$tmp}})
	  {
	    $names{$tmp2}=1;
	  }
      }
    my @names = keys %names;
    foreach my $n (@names)
      {
	foreach my $n2 (@names)
	  {
	    $master_names{$n}{$n2}=1;
	  }
      }

    my $strand = $line[6] =~ /-/ ? -1 :1;
    my $type = ($line[2]);
    my @type = ($type);
#    push @type, "CDS" if $add_cds && $type eq "mRNA";

    #phytozome replications of CDS to mRNA
    push @type, "mRNA" if $type eq "CDS";
    push @type, "gene" if $type eq "CDS";

    foreach my $type (@type)
      {
	$added_types{$type}++;
	$seen_chr{$chr}++;
	push @{$data{$line[1]}{$chr}{$name}{$type}{loc}}, {
						     start=>$line[3],
						     stop=>$line[4],
						     strand=>$strand,
						     chr=>$chr,
						    };
	foreach my $n (keys %names)
	  {
	    $data{$line[1]}{$chr}{$name}{$type}{names}{$n}=1;
	    foreach my $n1 ( keys %{$master_names{$n}})
	      {
		$data{$line[1]}{$chr}{$name}{$type}{names}{$n1}=1
	      }
	  }
      }
#    print Dumper \%data;
  }
if ($add_gene)
  {
    print "HERE!\n";
    foreach my $source (keys %data)
      {
	foreach my $chr_loc (keys %{$data{$source}})
	  {
	  name: foreach my $name (keys %{$data{$source}{$chr_loc}})
	      {
		my $start;
		my $stop;
		my $strand;
		my $chr;
		my %names;
		foreach my $type (keys %{$data{$source}{$chr_loc}{$name}})
		  {
		    map {$names{$_}=1} keys %{$data{$source}{$chr_loc}{$name}{$type}{names}};
		    foreach my $loc (@{$data{$source}{$chr_loc}{$name}{$type}{loc}})
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
		$data{$source}{$chr_loc}{$name}{gene}{loc}=[{
							 start=>$start,
							 stop=>$stop,
							 strand=>$strand,
							 chr=>$chr,
							}];
		$data{$source}{$chr_loc}{$name}{gene}{names}= \%names;
	      }
	  }
      }
  }
#print Dumper \%data;
#print Dumper \%annos;
#exit;

#time to load information into database

my %anno_types; #hash to store annotation type objects

foreach my $source (keys %data)
  {
    foreach my $chr_loc (keys %{$data{$source}})
      {

	foreach my $name (keys %{$data{$source}{$chr_loc}})
	  {
	    foreach my $feat_type (keys %{$data{$source}{$chr_loc}{$name}})
	      {
		print "\n" if $DEBUG;
		my ($start) = sort {$a<=>$b} map {$_->{start}} @{$data{$source}{$chr_loc}{$name}{$feat_type}{loc}};
		my ($stop) = sort {$b<=>$a} map {$_->{stop}} @{$data{$source}{$chr_loc}{$name}{$feat_type}{loc}};
		my ($strand) = map {$_->{strand}} @{$data{$source}{$chr_loc}{$name}{$feat_type}{loc}};
		my ($chr) = map {$_->{chr}} @{$data{$source}{$chr_loc}{$name}{$feat_type}{loc}};
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
		my %seen_locs;
		foreach my $loc (@{$data{$source}{$chr_loc}{$name}{$feat_type}{loc}})
		  {
		    next if $seen_locs{$loc->{start}}{$loc->{stop}};
		    $seen_locs{$loc->{start}}{$loc->{stop}}=1;
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
		my $names = $data{$source}{$chr_loc}{$name}{$feat_type}{names};
		my %names;
		foreach my $name (keys %$names)
		  {
		    $names{$name}=1;
		    foreach my $name2 (keys %{$master_names{$name}})
		      {
			$names{$name2}=1;
		      }
		  }
		foreach my $tmp (keys %names)
		  {

		    print "Adding name $tmp to feature ", $featid ,"\n" if $DEBUG;
		    my $feat_name = $feat->add_to_feature_names({
								 name=>$tmp,
								 #				   feature_id=>$featid,
								}) if $GO ;
		    if ($annos{$tmp})
		      {
			foreach my $anno (keys %{$annos{$tmp}})
			  {
			    next unless $anno;
			    my $type_name = $annos{$tmp}{$anno}{type} || "Note";
			    my ($anno_type) = $anno_types{$type_name};
			    unless ($anno_type)
			      {
				($anno_type) = $coge->resultset('AnnotationType')->find_or_create({name=>$type_name});
				$anno_types{$type_name} = $anno_type;
			      }
			    my $link = $annos{$tmp}{$anno}{link};
			    print "Adding annotation ($type_name): $anno\n" if $DEBUG;
			    print "\tlink: $link\n" if $DEBUG && $link;
			    my $anno = $feat->add_to_annotations({annotation=>$anno, link=>$link, annotation_type_id => $anno_type->id}) if $GO && $anno;
			  }
		      }
		  }
	      }
	  }
      }
  }

print "Completed working on dataset: ", $ds->name. " (".$ds->id.")\n";
close IN;
print "\n";
print "Seen data types:\n";
print join ("\n", map {$_."\t".$seen_types{$_}} sort keys %seen_types),"\n";

print"\n";
print "Added data types:\n";
print join ("\n", map {$_."\t".$added_types{$_}} sort keys %added_types),"\n";

print "\n";
print "Seen chromosomes:\n";
print join ("\n", map {$_."\t".$seen_chr{$_}} sort keys %seen_chr),"\n";

sub generate_ds
  {
    my %opts = @_;
    my $ds_name = $opts{ds_name};
    my $ds_desc = $opts{ds_desc};
    my $ds_link = $opts{ds_link};
    my $ds_version = $opts{ds_version};
    my $ds_id = $opts{ds_id};
    my $source_id = $opts{source_id};
    unless ($ds_name || $ds_id)
      {
	warn "no dataset name or database id specified\n";
	return;
      }
    my $ds = $ds_id ? $coge->resultset('Dataset')->find($ds_id) :
      $coge->resultset('Dataset')->find_or_create({
						   name                => $ds_name,
						   description         => $ds_desc,
						   link                => $ds_link,
						   data_source_id      => $source_id,
						   version=>$ds_version,
						  });;
    return $ds;

  }

sub process_annotation_file
  {
    my %opts = @_;
    my $file = $opts{file};
    my $annos = $opts{annos};
    my $anno_names = $opts{names};
    open (IN, $file);
    while (<IN>)
      {
	chomp;
	next unless $_;
	my @line = split /\t/;
	my $name = $line[0];
	my $name2 = $line[0];
	$name2 =~ s/\.\d+$//; #get rid of trailing version number if present

	#extra sorghum names
	unless ($line[1] =~ /no.*symbols/)
	  {
	    foreach my $tmp (split/,/,$line[1])
	      {

		$anno_names->{$name}{$tmp}=1;
		$anno_names->{$name2}{$tmp}=1;
	      }
	  }

	#sorghum annotations
	unless ($line[2] =~ /no.*defline/)
	  {
#	    print $line[2],"\n";
	    foreach my $tmp (split/;/,$line[2])
	      {
		$tmp =~ s/^\s+//;
		$tmp =~ s/\s+$//;
		$tmp =~ s/^\[\s*(.*?)\s*\]$/$1/;
#		$tmp =~ s/\s*\]$//;
		$tmp =~ s/\s*,$//;
#		print $name,"\t", $tmp,"\n";
		$annos->{$name}{$tmp}={};
#		$anno_names->{$name2}{$tmp}=1;
	      }
	  }
	unless ($line[3] =~ /no.*ids/)
	  {
	    foreach my $tmp (split/,/,$line[3])
	      {
#		print $name,"\t",q$tmp,"\n";
		my $link = "http://pfam.sanger.ac.uk/family/".$tmp;
		$annos->{$name}{"View at Pfam ($tmp)"}={link=>$link, type=>"Pfam link"};
	      }
	  }
	unless ($line[4] =~ /no.*ids/)
	  {
	    foreach my $tmp (split/,/,$line[4])
	      {
#		print $name,"\t",$tmp,"\n";
		my $link = "http://www.pantherdb.org/panther/familyList.do?searchType=basic&fieldName=all&searchAll=true&listType=6&fieldValue=".$tmp;
		$annos->{$name}{"View at PantherDB ($tmp)"}={link=>$link, type=>"PantherDB link"};
	      }
	  }
	unless ($line[5] =~ /no.*ids/)
	  {
	    foreach my $tmp (split/,/,$line[4])
	      {
#		print $name,"\t",$tmp,"\n";
#		my $link = "http://www.pantherdb.org/panther/familyList.do?searchType=basic&fieldName=all&searchAll=true&listType=6&fieldValue=".$tmp;
		$annos->{$name}{$tmp}={type=>"KOG"};
	      }
	  }

	unless ($line[6] =~ /no.*id/)
	  {
	    foreach my $tmp (split/,/,$line[6])
	      {
		my $link = "http://www.expasy.org/enzyme/".$tmp;
#		print $name,"\t",$tmp,"\t", $link, "\n";
		$annos->{$name}{"View at ExPASy ($tmp)"}={link=>$link, type=>"EC link"};
	      }
	  }
	unless ($line[7] =~ /no.*orthology/)
	  {
	    foreach my $tmp (split/,/,$line[7])
	      {
		my $link = "http://www.genome.jp/dbget-bin/www_bget?ko:".$tmp;
#		print $name,"\t",$tmp,"\t", $link, "\n";
		$annos->{$name}{"View at KEGG ($tmp)"}={link=>$link, type=>"KEGG link"};
	      }
	  }
	unless ($line[8] =~ /no.*hit/)
	  {
	    foreach my $tmp (split/,/,$line[8])
	      {
		my $link = "/CoGe/FeatView.pl?accn=".$tmp;
#		print $name,"\t",$tmp,"\t", $link, "\n";
		$annos->{$name}{"$tmp"}={link=>$link, type=>"Best Arabidopsis Match"};
	      }
	  }
#	unless (!$line[9] || $line[9] =~ /no.*symbols/) ## skipping
#	  {
#	    foreach my $tmp (split/,/, $line[9])
#	      {
#		my $link = "/CoGe/FeatView.pl?accn=".$tmp;
#		print $name,"\t",$tmp,"\t", $link, "\n";
#		$annos->{$name}{"$tmp"}={link=>$link, type=>"Best Arabidopsis Match"};
#	      }
#	  }
	unless (!$line[10] || $line[10] =~ /no.*hit/)
	  {
	    foreach my $tmp (split/;/,$line[10])
	      {
		$tmp =~ s/^\s+//;
		$tmp =~ s/\s+$//;
#		my $link = "/CoGe/FeatView.pl?accn=".$tmp;
#		print $name,"\t",$tmp, "\n";
		$annos->{$name}{"$tmp"}={type=>"Arabidopsis annotation"};
	      }
	  }

      }
    close IN;
  }
