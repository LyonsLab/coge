#!/usr/bin/perl

use strict;
use CoGe::Accessory::GBlite;
use Roman;
use Data::Dumper;
use Getopt::Long;
use LWP::Simple;
use CoGeX;

# variables
my ($DEBUG, $GO, $ERASE, $DELETED, @files, $dir, @accns, $tmpdir, $help, $chromosome,$di_id);
my $genomic_seq_len = 10000; 		# length to break up genomic sequence
my $connstr = 'dbi:mysql:dbname=DB;host=HOST;port=PORT';
my $coge = CoGeX->connect($connstr, 'USER', 'PASSWORD' );
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
	    "di_id=s"=>\$di_id,
	   );
$chromosome = 1 unless $chromosome;
$tmpdir = "/tmp/gb" unless $tmpdir;	# set default directory to /tmp
$DEBUG = 0 unless defined $DEBUG;   # set to 1 to enable debug printing
$GO = 0 unless defined $GO; 		# set to 1 to actually make db calls (instead of testing)
$ERASE = 0 unless defined $ERASE;
$DELETED = 0 if $ERASE;

help () if $help;
return unless $di_id;
my $dataset = $coge->resultset('Dataset')->find($di_id);
exit unless $dataset;

# read in each file into a GBFile object
# if the directory is given, process from dir
# if accns is given, process by accns and take precedence over directory
push @files, process_dir($dir) if $dir;
push @files, process_accns(@accns) if scalar @accns;

print "Go = $GO \n" if $DEBUG;

# loop through all the files in input
foreach my $longfile (@files)
  {
    # open our input data file...
    my $genbank = new CoGe::Accessory::GBlite( $longfile );
    # remove path from file name
    my ($file) = $longfile=~ /([^\/]*$)/;
    chomp $file;
    print "Processing $file... \n";
    my $ds_count = 0;
    my $orig_chr;
    while(my $entry = $genbank->nextEntry)
	{
	  my $EXIT =0;
	  $ds_count++;
	  print "Working on entry ".$entry->accession,"\n";
	  $chromosome = $entry->accession;
	  my $name = $entry->accession;
	  $name .= ".gbk" unless $name =~ /gbk$/;
	  ### Main Feature Processing Loop ###
	  print "Processing features for ".$entry->accession."...\n" unless $EXIT;
	  foreach my $feature (@{$entry->features()})
	    {
	      # search the db for that feature
	      next if $feature->key() eq "misc_feature";
	      unless ($feature->key())
		{
		  print "Feature has no feature type name \$feature->key():\n";
		  print Dumper $feature;
		  next;
		}

	      my $feat_type = $coge->resultset('FeatureType')->find_or_create({ name => $feature->key() })  if $GO;
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
	      if ($feature->key() =~ /source/i)
		{
		  # generate name based on organism name and chromosome
		  my $feat_name = $coge->resultset('FeatureName')->create(
									  {
									   name => $dataset->organism->name,
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
								       }) if $entry->version && $GO;

		  # generate name for GI
		  $feat_name = $coge->resultset('FeatureName')->create(
								       {
									name => $entry->gi,
									feature_id => $db_feature->id
								       }) if $entry->gi && $GO;
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
	      foreach  my $anno (keys %{$annot})
		{
		  # deal with db_xref: (taxon:3702) (GeneID:821318) (GI:18379324)
		  if ($anno =~ /xref/i)
		    {
		      my $anno_type_group = $coge->resultset('AnnotationTypeGroup')->find_or_create( { name => $anno } )  if $GO;
		      # go through each of the entries in the db_xref qualifier values and split on ':', then add entries individually
		      my @xrefs = split(/ /, $annot->{$anno}); # split values on \s
		      foreach my $xref (@xrefs)
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
			 || $anno =~ /synonym/i )    # synonyms are embedded in the /note= tag! these are names
		    {
		      my $feat_name = $db_feature->add_to_feature_names({
									 name       => $annot->{$anno},
									}) if $GO;
		    }
		  elsif ($anno =~ /translation/i) # this needs to be entered into the sequence table
		    {
		      my $seq_type = $coge->resultset('SequenceType')->find_or_create(
										      {
										       name        => "protein",
										       description => "translation"
										      }) if $GO;
		      my $sequence = $db_feature->add_to_sequences({
								    sequence_type_id => $seq_type->id(),
								    sequence_data    => $annot->{$anno},
								   }) if $GO;
		    }
		  elsif ($anno eq "note")
		    {
		      # if go annot are present, they'll be in the note qualifier,
		      # so process is specifically
		      my $leftover = "";
		      my @temp = split( /;/, $annot->{$anno} );
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
		  else           ##everything else
		    {
		      my $anno_type = $coge->resultset('AnnotationType')->find_or_create( { name => $anno } ) if $GO;
		      my $sub_anno = $db_feature->add_to_annotations(
								     {
								      annotation         => $annot->{$anno},
								      annotation_type_id => $anno_type->id(),
								     }
								    ) if $GO;
		    }
		}
	    }
	  print "Processing Genomic Sequence. . .\n" unless $EXIT;
	  load_genomic_sequence(len=> $genomic_seq_len, ds=>$dataset, seq=>$entry->sequence, chr=>$chromosome);
	}
	print "completed parsing $file!\n" if $DEBUG;

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
    my ($feat_type) = $coge->resultset('FeatureType')->search({name=>'contig'});
    my ($feat) = $coge->resultset('Feature')->create({
						      feature_type_id     => $feat_type->id,
						      dataset_id => $dataset->id,
						      chromosome=> $chr,
						      strand=>1,
						      start=>1,
						      stop=>$seqlen,
						     }) if $GO;
    $feat->add_to_locations(
			    {
			     start      => 1,
			     stop       => $seqlen,
			     strand     => 1,
			     chromosome => $chr
			    }) if $GO;
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
	  my $found = 0;
	  print "Searching for $id. . .\n";
	  foreach my $ds ($coge->resultset('Dataset')->search(name=>"$id.gbk"))
	    {
	      print "\t", "Found: ", $ds->name,": ",$ds->description,"\n";
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
