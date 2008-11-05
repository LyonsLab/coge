#!/usr/bin/perl

use strict;
use CoGe::Accessory::GenBank;
use LWP::Simple;
use Roman;
use Data::Dumper;
use Getopt::Long;
use CoGeX;

# variables
my ($DEBUG, $GO, $ERASE, $DELETED, @files, $dir, @accns, $tmpdir, $help, $chromosome, $ds_link, $test);
my $genomic_seq_len = 10000; 		# length to break up genomic sequence
my $connstr = 'dbi:mysql:dbname=genomes;host=biocon;port=3306';
my $coge = CoGeX->connect($connstr, 'cnssys', 'CnS' );
#$coge->storage->debugobj(new DBIxProfiler());
#$coge->storage->debug(1);


$| =1;

# parse command line options
GetOptions ( 
	    "debug" 			=> \$DEBUG,
	    "go"    			=> \$GO,
	    "erase|e"  		=> \$ERASE, 
	    "file|f=s" 		=> \@files,
	    "dir|d=s"  		=> \$dir,
	    "accn|a=s"			=> \@accns,
	    "temp_dir|td=s"	=> \$tmpdir,
	    "help|h"			=> \$help,
	    "chromosome|chr=s"=>\$chromosome,
	    "dataset_link=s"=>\$ds_link,
	    "test"=>\$test, #to add the name test to dataset name for testing purposes
	   );
$chromosome = 1 unless $chromosome;
$tmpdir = "/tmp/gb" unless $tmpdir;	# set default directory to /tmp
$DEBUG = 0 unless defined $DEBUG;   # set to 1 to enable debug printing
$GO = 0 unless defined $GO; 		# set to 1 to actually make db calls (instead of testing)
$ERASE = 0 unless defined $ERASE;
$DELETED = 0 if $ERASE;
#my $ERASE = 0; 					# set to 1 to clear the database of entries created for $dataset
#my @files = @ARGV;
help () if $help;

# read in each file into a GBFile object
# if the directory is given, process from dir
# if accns is given, process by accns and take precedence over directory
#push @files, process_dir($dir) if $dir;
#push @files, process_accns(@accns) if scalar @accns;

print "Go = $GO \n" if $DEBUG;

# loop through all the files in input
my $increment_chr_version =0;
my $file_count =0;
my $orig_chr;

foreach my $accn (@accns)
  {
    print "Processing accn $accn... \n";
    $file_count++;
    # open our input data file...
    my $genbank = new CoGe::Accessory::GenBank();
    $genbank->get_genbank_from_ncbi(file=> "$tmpdir/$accn.gbk", accn=>$accn);
    # remove path from file name
#    my ($file) = $longfile=~ /([^\/]*$)/;
#    chomp $file;

    my ($organism, $data_source, $dataset);
    my $ds_count = 0;
    $chromosome = $genbank->chromosome if $genbank->chromosome;
    my $EXIT =0;
    $chromosome = $orig_chr if $orig_chr;
    $chromosome = $orig_chr.".".$ds_count if ($orig_chr && $orig_chr eq $chromosome);
    $orig_chr = $chromosome unless defined $orig_chr;
    $chromosome = $orig_chr.".".$file_count if $increment_chr_version;
    $chromosome =~ s/^\s+//;
    $chromosome =~ s/\s+$//;
    $chromosome =~ s/\s+/_/g;
    $ds_count++;
    print "Working on genbank ".$genbank->accession,", chromosome $chromosome","\n";
    unless ($organism && $data_source)
      {
	($organism) = get_organism($genbank);
	if ($organism)
	  {
	    print "Organism info:";
	    print "\t",$organism->id,": ";
	    print "\t",$organism->name,": ";
	    print "\t",$organism->description,"\n";
	  }
	$data_source = get_data_source();
      }
    my $dataset_desc = "LOCUS: "     . $genbank->locus();
    $dataset_desc   .= ", ACCESSION: " . $genbank->accession();
    $dataset_desc   .= ", VERSION: "   . $genbank->version();
    # testing to see if already in database
    my $version = $genbank->version();
    my $name = $genbank->accession;
    $name .= ".gbk" unless $name =~ /gbk$/;
    $name .= ".test" if $test;
    my @dataset_test = $coge->resultset('Dataset')->search(
							   {
							    name => $name,
							    version=>$genbank->version,
							    description         => $dataset_desc,
							    organism_id         => $organism->id,
							    data_source_id      => $data_source->id(),
							   }) if $GO;
    # loop through results and see if we have any matches
    foreach my $test (@dataset_test)
      {
	    # test to see if match
	    if (($test->name == $accn || $test->name == $genbank->accession) && $test->version == $genbank->version)
	      {
		# stop executing if a match is found
		$EXIT = 1;
		print "Genome is already inserted in the database...\n";
		
		# if match and erase flag is set, erase the entry
		if ($ERASE)
		  {
		    print "Clearing database of entries associated with ".$test->name." (version ".$test->version.")...";
		    $test->delete();
		    print "Deleted!\n";
		    $DELETED = 1;
		    next;
		  }
	      }
	  }
	
	# if delete is set but nothing has been deleted
	if ($ERASE && !$DELETED)
	  {
	    $EXIT = 1;
	    print "No match found. Nothing deleted.\n";
	  }
	next if $EXIT;
	
	# actually create a dataset object if 
	$ds_link = "http://www.ncbi.nlm.nih.gov/entrez/viewer.fcgi?db=nucleotide&dopt=gbwithparts&list_uids=".$accn unless $ds_link;
	$dataset = $coge->resultset('Dataset')->find_or_create(
							       {
								name                => $name,
								description         => $dataset_desc,
								link                => $ds_link,
								organism_id         => $organism->id,
								data_source_id      => $data_source->id(),
								version             => $genbank->version,
							       })  if $GO;
	
    
	### Main Feature Processing Loop ###
    
    my @gbs;
    if (@{$genbank->wgs_data})
      {
	@gbs = @{$genbank->wgs_data};
	
      }
    else
      {
	push @gbs, $genbank;
      }
    foreach my $entry (@gbs)
      {
	$chromosome = "contig_".$entry->accession if @{$genbank->wgs_data};
	print "Processing features for ".$entry->accession."...\n" unless $EXIT;
	foreach my $feature (@{$entry->features()})
	  {
	    # search the db for that feature
	    unless ($feature->type())
	      {
		print "Feature has no feature type name \$feature->type():\n";
		print Dumper $feature;
		next;
	      }
	    if ($feature->type() =~ /source/i)
	      {
		next;
	      }
	    my $feat_type = $coge->resultset('FeatureType')->find_or_create({ name => $feature->type() })  if $GO;
	    # create a db_feature for to link this feature with the dataset table
	    my ($start, $stop, $strand) = get_feature_location($feature);
	    my $db_feature = $coge->resultset('Feature')->create({
								  feature_type_id     => $feat_type->id,
								  dataset_id => $dataset->id,
								  chromosome=> $chromosome,
								  strand=>$strand,
								  start=>$start,
								  stop=>$stop,
								 }) if $GO;
	    
	    # expect first feature to be the source feature!
	    if ($feature->type() =~ /source/i) #source is now skipped.
	      {
		# generate name based on organism name and chromosome
		my $feat_name = $coge->resultset('FeatureName')->create(
									{
									 name => $organism->name,
									 description => "Chromosome " . $chromosome,
									 feature_id  => $db_feature->id
									}) if $GO;
		# generate name for accession
		$feat_name = $coge->resultset('FeatureName')->create(
								     {
								      name => $entry->accession,
								      feature_id => $db_feature->id()
								     }) if $GO;
		# generate name for version
		$feat_name = $coge->resultset('FeatureName')->create(
								     {
								      name => $entry->accession.".".$entry->version,
								      feature_id => $db_feature->id
								     }) if $GO;
		
		# generate name for GI
		$feat_name = $coge->resultset('FeatureName')->create(
								     {
								      name => $entry->gi,
								      description=>"GI number",
								      feature_id => $db_feature->id
								     }) if $GO;
	      }
	    
	    # add a location entry
	    my $loc_string = $feature->location;
	    $loc_string =~ s/complement//g;
	    $loc_string =~ s/join//;
	    $loc_string =~ s/order//;
	    $loc_string =~ s/\(|\)//g;
	    
	    # loop through the locations
	    foreach my $loc (split /,/,$loc_string)
	      {
		$loc =~ s/<|>//g;
		my ($start, $stop) = split /\.\./, $loc;
		$start =~ s/\^.*//;
		$stop =~ s/\^.*//;
		$stop = $start unless $stop;
		die "problem with start $start or stop $stop\n" unless $start =~ /^\d+$/ && $stop =~ /^\d+$/;
		my $location = $db_feature->add_to_locations(
							     {
							      start      => $start,
							      stop       => $stop,
							      strand     => $feature->strand,
							      chromosome => $chromosome
							     }) if $GO;
	      }
	    # now work through the qualifiers for this feature
	    # start by getting the hashref of qualifiers
	    my $annot = $feature->qualifiers();
	    
	    my %names;
	    #	      print Dumper $annot;
	    foreach  my $anno (keys %{$annot})
	      {
		my $stuff = $annot->{$anno};
		# deal with db_xref: (taxon:3702) (GeneID:821318) (GI:18379324)
		if ($anno =~ /xref/i)
		  {
		    my $anno_type_group = $coge->resultset('AnnotationTypeGroup')->find_or_create( { name => $anno } )  if $GO;
		    # go through each of the entries in the db_xref qualifier values and split on ':', then add entries individually
		    foreach my $xref (@{$stuff})
		      {
			my @inner = split(/:/, $xref );
			# first add the annot_type_obj
			my $anno_type = $coge->resultset('AnnotationType')->find_or_create(
											   {
											    name => $inner[0],
											    annotation_type_group_id => $anno_type_group->id(),
											   }) if $GO;
			
			# now create the row for the data value of the xref
			my $sub_anno = $db_feature->add_to_annotations(
								       {
									annotation         => $inner[1],
									annotation_type_id => $anno_type->id()
								       })  if $GO;
		      }
		  }
		elsif ($anno =~ /locus_tag/i
		       || $anno =~ /transcript_id/i
		       || $anno =~ /protein_id/i
		       || $anno =~ /gene/i
		       || $anno =~ /standard_name/i
		       || $anno =~ /synonym/i # synonyms are embedded in the /note= tag! these are names
		       || $anno eq "names")    
		  {
		    foreach my $item (@{$stuff})
		      {
			foreach my $thing (split/;/,$item)
			  {
			    $thing =~ s/^\s+//;
			    $thing =~ s/\s+$//;
			    $names{$thing}=1;
			  }
		      }
		  }
		elsif ($anno =~ /translation/i) # this needs to be entered into the sequence table
		  {
		    my $seq_type = $coge->resultset('SequenceType')->find_or_create(
										    {
										     name        => "protein",
										     description => "translation"
										    }) if $GO;
		    foreach my $item (@{$stuff})
		      {
			$item =~ s/\s+//g;
			my $sequence = $db_feature->add_to_sequences({
								      sequence_type_id => $seq_type->id(),
								      sequence_data    => $item,
								     }) if $GO;
		      }
		  }
		elsif ($anno eq "note")
		  {
		    # if go annot are present, they'll be in the note qualifier,
		    # so process is specifically
		    foreach my $item (@{$stuff})
		      {
			my $leftover = "";
			my @temp = split( /;/, $item );
			foreach my $go_raw (@temp)
			  {
			    if ($go_raw =~ /go_/)
			      {
				while ($go_raw =~ /(go_.*?):\s+(.*?)\[goid G?O?:?(.*?)\]/g)
				  {
				    # example:
				    # go_function: nucleic acid binding [goid 0003676]
				    my $anno_type_group = $coge->resultset('AnnoationTypeGroup')->find_or_create( { name => $1 } ) if $GO;
				    # $1 should be "go_function"
				    my $anno_type = $coge->resultset('AnnotationType')->find_or_create(
												       {
													name => $3,    #this should be "0003676"
													annotation_type_group_id => $anno_type_group->id(),
												       }
												      ) if $GO;
				    my $sub_anno = $db_feature->add_to_annotations(
										   {
										    annotation => $2,    #this should be "nucleic acid binding"
										    annotation_type_id => $anno_type->id
										   }
										  )  if $GO;
				  }
			      } else {
				$leftover .= " " . $go_raw if $go_raw;
			      }
			    # now just add the note remainder
			    $leftover =~ s/^\s+//;
			    $leftover =~ s/\s+$//;
			    if ($leftover)
			      {
				my $anno_type = $coge->resultset('AnnotationType')->find_or_create( { name => $anno } ) if $GO;
				my $sub_anno = $db_feature->add_to_annotations(
									       {
										annotation         => $leftover,
										annotation_type_id => $anno_type->id(),
									       }
									      ) if $GO;
			      }
			  }
		      }
		  }
		else           ##everything else
		  {
		    foreach my $item (@{$stuff})
		      {
			my $anno_type = $coge->resultset('AnnotationType')->find_or_create( { name => $anno } ) if $GO;
			my $sub_anno = $db_feature->add_to_annotations(
								       {
									annotation         => $item,
									annotation_type_id => $anno_type->id(),
								       }
								      ) if $GO;
		      }
		  }
	      }
	    foreach my $name (keys %names)
	      {
		$name =~ s/\s+$//g;
		$name =~ s/^\s+//g;
		my $feat_name = $db_feature->add_to_feature_names({
								   name       => $name,
								  }) if $GO;
	      }
	  }
	print "Processing Genomic Sequence. . .\n" unless $EXIT;
	load_genomic_sequence(len=> $genomic_seq_len, ds=>$dataset, seq=>$entry->sequence, chr=>$chromosome);
      }
    print "completed parsing $accn!\n" if $DEBUG;
    
    # cleanup
    #	$organism->delete() if $GO;
  }

sub load_genomic_sequence
  {
    my %opts = @_;
    my $seq = $opts{seq};
    my $len = $opts{len};
    my $chr = $opts{chr};
    my $ds = $opts{ds};
    my $seqlen = length $seq;
    print "Loading genomic sequence ($seqlen nt)\n" if $DEBUG;
    $seq =~ s/\s//g;
    $seq =~ s/\n//g;
    my $i = 0;
    while ($i < $seqlen)
      {
	my $str = substr($seq,$i,$len);
	my $start = $i+1;
	my $stop = $i + length $str;
	$ds->add_to_genomic_sequences({start => $start,
				      stop => $stop,
				      chromosome => $chr,
				      sequence_data => $str,
				     }) if $GO;
	$i += $len;
      }
  }

sub process_dir
{
	my $dir = shift;
	my @files;
	opendir (DIR, $dir) || die "Can't open $dir for reading: $!";
	while (my $f = readdir(DIR))
	{
		next if -d "$dir/$f";
		push @files, "$dir/$f" if -r "$dir/$f";
	}
	return @files;
}

sub get_genbank_from_ncbi
{
	my $accn = shift;
	my $url = "http://www.ncbi.nlm.nih.gov/entrez/viewer.fcgi?db=nucleotide&qty=1&c_start=1&dopt=gbwithparts&send=Send&sendto=t&from=begin&to=end&extrafeatpresent=1&ef_CDD=8&ef_MGC=16&ef_HPRD=32&ef_STS=64&ef_tRNA=128&ef_microRNA=256&list_uids=$accn";
	print "Fetching sequence for $accn from $url\n . . .\n" if $DEBUG;
	my $content = get ($url);
	my $message = $content ? "Have content for $accn.\n" : "Did not get a sequence for $accn\n";
	print "$message" if $DEBUG;
	
	return $content;
}

sub process_accns
{
	my @accns = @_;
	my @files;
	id: foreach my $id (@accns)
	{
#	  $id =~ s/NZ_//; #remove this from Whole Shotgun Sequence as the entry without this may have more assembly (for whatever reason)
	  my $found = 0;
	  print "Searching for $id. . .";
	  foreach my $ds ($coge->resultset('Dataset')->search(name=>"$id.gbk"))
	    {
	      print "Found: ", $ds->name,": ",$ds->description,"\n";
	      $found =1;
	      next;
	    }
	  next id if $found &! $ERASE;
	  mkdir ($tmpdir) unless -d $tmpdir;
	  my $file = "$tmpdir/$id.gbk";
	  unless (-r $file)
	    {
	      print "getting $id from ncbi\n";
	      my $gb = get_genbank_from_ncbi($id);
	      unless ($gb)
		{
		  warn "no content was generated for accession $id.  Skipping.\n";
		  next;
		}
	      open (OUT, ">$file") || die "can't open $file: $!";
	      print OUT $gb;
	      close OUT;
	    }
	  push @files, $file;
	}
	return wantarray ? @files : \@files;
}

sub help
{
	print qq
	{
		Welcome to $0!  This program loads genbank entries into the CoGe genomes database

		Options:
	};
}

sub get_feature_location
  {
    my $feat = shift;
    my $loc_string = $feat->location;
    my $strand = $feat->location =~ /complement/ ? "-1" : 1;
    $loc_string =~ s/complement//g;
    $loc_string =~ s/join//;
    $loc_string =~ s/order//;
    $loc_string =~ s/\(|\)//g;
    my ($rstart, $rstop);
    foreach my $loc (split /,/,$loc_string)
      {
	my ($start, $stop) = split /\.\./, $loc;
	($start) = $start =~ /(\d+)/;
	($stop) = $stop =~ /(\d+)/;
	$start =~ s/\^.*//;
	$stop =~ s/\^.*//;
	$stop = $start unless $stop;
	$rstart = $start unless $rstart;
	$rstop = $stop unless $stop;
	$rstart = $start if $start < $rstart;
	$rstop = $stop if $stop > $rstop;
      }
    return ($rstart, $rstop, $strand);
  }

sub get_organism
{
  my ($entry) = shift;
  my $name = $entry->data_source();
  $name =~ s/'//g;
#  print $name,"\n";
#  print $entry->organism,"\n";
  my $org = $coge->resultset('Organism')->find({name=>$name}) if $GO;
  unless ($org)
    {
      $org = $coge->resultset('Organism')->find_or_create(
							  {
							   name=>$name,
							   description=>$entry->organism()
							  }) if $GO;
    }
  return $org;
}
sub get_data_source
{
  return $coge->resultset('DataSource')->find_or_create(
						     {
						      name=>'NCBI',
						      description=>"National Center for Biotechnology Information",
						      link=>'www.ncbi.nih.gov'
						     }) if $GO;
}

