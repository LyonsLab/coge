#!/usr/bin/perl -w

use strict;
use CoGe::Accessory::GenBank;
use LWP::Simple;
use Data::Dumper;
use Getopt::Long;
use CoGeX;
use File::Path;

# variables
my ($DEBUG, $GO, $ERASE,$autoupdate,@accns, $tmpdir, $help, $user_chr, $ds_link, $test);
my $connstr = 'dbi:mysql:dbname=coge;host=biocon.berkeley.edu;port=3306';
my $coge = CoGeX->connect($connstr, 'cnssys', 'CnS' );
#$coge->storage->debugobj(new DBIxProfiler());
#$coge->storage->debug(1);

$| =1;

# parse command line options
GetOptions (
	    "debug" 			=> \$DEBUG,
	    "go"    			=> \$GO,
	    "erase|e"  		=> \$ERASE,
	    "accn|a=s"			=> \@accns,
	    "temp_dir|td=s"	=> \$tmpdir,
	    "help|h"			=> \$help,
	    "user_chr|chr=s"=>\$user_chr,
	    "dataset_link=s"=>\$ds_link,
	    "test"=>\$test, #to add the name test to dataset name for testing purposes
	    "autoupdate"=>\$autoupdate,
	   );
$user_chr = 1 unless $user_chr;
$tmpdir = "/tmp/gb" unless $tmpdir;	# set default directory to /tmp

mkpath ($tmpdir) unless -d $tmpdir;
$DEBUG = 0 unless defined $DEBUG;   # set to 1 to enable debug printing
$GO = 0 unless defined $GO; 		# set to 1 to actually make db calls (instead of testing)
$ERASE = 0 unless defined $ERASE;
$autoupdate = 0 unless $autoupdate;

my $formatdb =  "/usr/bin/formatdb -p F -o T"; #path to blast's formatdb program
my %dsg;
#my $ERASE = 0; 					# set to 1 to clear the database of entries created for $dataset
help () if $help;

print "Go = $GO \n" if $DEBUG;

my $data_source = get_data_source(); #for NCBI

# loop through all the accessions
my %previous_datasets; #storage for previously loaded accessions
my $dsg; #storage for coge dataset_group_obj
my %dsg_update;#storage for previous dataset groups that had a dataset deleted and needs to be reloaded
accn: foreach my $accn (@accns)
  {
    print "Working on $accn...";
    my $previous = check_accn($accn);
    foreach my $item (@$previous)
      {
	if (!$item->{version_diff} && !$item->{length_diff})
	  {
#	    push @previous_datasets, $item->{ds};
	    $previous_datasets{$item->{ds}->id} = $item->{ds};
	    print "previously loaded\n";
	    next accn;
	  }
	elsif (!$item->{version_diff} && $item->{length_diff})
	  {
	    print "Detected a difference in total genomic length between CoGe (".$item->{coge_length}.") and NCBI(".$item->{ncbi_length}.").  Would you like to delete and reload? (y/n)";
	    my $ans = <STDIN> unless $autoupdate;
	    if ($autoupdate || $ans =~ /y/i)
	      {
	       print "Autoupdate flag set to true.  Automatically reloading dataset.\n" if $autoupdate;
		my $ds = $item->{ds};
		foreach my $item ($ds->dataset_groups)
		  {
		    $dsg_update{$item->id}{dsg}=$item;
		    push @{$dsg_update{$item->id}{accn}}, $accn;
#		    delete_dataset_group($item);
		  }
		$ds->delete;
	      }
	    else
	      {
		$previous_datasets{$item->{ds}->id} = $item->{ds};
		next accn;
	      }
	  }
      }
    if ($ERASE)
      {
	print "skipping loading due to ERASE flag being set\n";
	next;
      }

    my $genbank = new CoGe::Accessory::GenBank();
    $genbank->debug(1);
    $genbank->get_genbank_from_ncbi(file=> "$tmpdir/$accn.gbk", accn=>$accn,);
    my $EXIT =0;
    my $chromosome;
    $chromosome = $user_chr if $user_chr;
    if ($genbank->chromosome)
      {
	$chromosome = $genbank->chromosome ;
	print "#"x20,"\n";
	print $genbank->data_source();
	print "GBChr:  $chromosome\n";
	print "#"x20,"\n";
      }
    $chromosome =~ s/chromosome//i;
    $chromosome =~ s/chr//i;
    $chromosome =~ s/^\s+//;
    $chromosome =~ s/\s+$//;
    $chromosome =~ s/\s+/_/g;
    $chromosome =~ s/\.0$//;
    print "\tchromosome: $chromosome","\n";
    my ($organism, $dataset);
    unless ($organism)
      {
	($organism) = get_organism($genbank);
	if ($organism)
	  {
	    print "Organism info:";
	    print "\t",$organism->id,": ";
	    print "\t",$organism->name,"\n";
	    print "\t",$organism->description,"\n";
	  }
	else
	  {
	    print "WARNING:  Unable to retrieve an organism object for $accn.  Probably a problem with the genbank object file parsing\n" unless !$GO;
	    next if $GO;
	  }
      }
    my $dataset_desc = "LOCUS: "     . $genbank->locus();
    $dataset_desc   .= ", ACCESSION: " . $genbank->accession();
    $dataset_desc   .= ", VERSION: "   . $genbank->version();
    # testing to see if already in database
    my $version = $genbank->version();
    my $name = $genbank->accession;
    $name .= ".gbk" unless $name =~ /gbk$/;
    $name .= ".test" if $test;
    # actually create a dataset object if
    $ds_link = "http://www.ncbi.nlm.nih.gov/entrez/viewer.fcgi?db=nucleotide&dopt=gbwithparts&list_uids=".$accn unless $ds_link;
    $dataset = $coge->resultset('Dataset')->create(
						   {
						    name                => $name,
						    description         => $dataset_desc,
						    link                => $ds_link,
						    data_source_id      => $data_source->id(),
						    version             => $genbank->version,
						   })  if $GO;
    ### Main Feature Processing Loop ###
    my @gbs;
    if (@{$genbank->wgs_data})
      {
	entry: foreach my $entry (@{$genbank->wgs_data})
	  {
	    my $accn = $entry->accession;
	    while (!$accn)
	      {

		$entry->get_genbank_from_ncbi(reload=>1);
		print "#"x20,"\n";
		print "Warning.  Didn't not retrieve a valid accession for entry.\n";
		print "Requested id: ".$entry->requested_id."\n";
		print "Trying to retrieve valid entry.\n";
		print "#"x20,"\n";
		$accn=$entry->accession;
	      }
	    print "Checking WGS $accn...";
	    my $previous = check_accn($accn);
	    foreach my $item (@$previous)
	      {
		if (!$item->{version_diff} && !$item->{length_diff})
		  {
		    $previous_datasets{$item->{ds}->id} = $item->{ds};
		    print "previously loaded\n";
		    next entry;
		  }
		elsif (!$item->{version_diff} && $item->{length_diff})
		  {
		    print "Detected a difference in total genomic length between CoGe (".$item->{coge_length}.") and NCBI(".$item->{ncbi_length}.").  Would you like to delete and reload? (y/n)";
		    my $ans = <STDIN> unless $autoupdate;
		    if ($autoupdate || $ans =~ /y/i)
		      {
			print "Autoupdate flag set to true.  Automatically reloading dataset.\n" if $autoupdate;
			my $ds = $item->{ds};
			my $dsg_flag;
			foreach my $item ($ds->dataset_groups)
			  {
			    delete_dataset_group($item);
			  }
			$ds->delete unless $dsg_flag;
		      }
		    else {
		      $previous_datasets{$item->{ds}->id} = $item->{ds};
		      next entry;
		    }
		  }
		else
		  {
		    print "not present.  Will be loaded.\n";
		  }
	      }
	    push @gbs, $entry;
	  }
#	@gbs = @{$genbank->wgs_data};

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
		#change source to chromosome
		$feature->type('chromosome');
#		next;
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
	    if ($feature->type() =~ /chromosome/i)
	      {
		# generate name based on chromosome
		my $feat_name = $coge->resultset('FeatureName')->create(
									{
									 name => $chromosome,
#									 description => "Chromosome " . $chromosome,
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

	    # loop through the locatons
	    foreach my $loc (split /,/,$loc_string)
	      {
		$loc =~ s/<|>//g;
		my ($start, $stop) = split /\.\./, $loc;
		$stop = $start unless $stop;
		$start =~ s/\^.*//;
		$stop =~ s/\^.*//;

		die "problem with $accn start $start or stop $stop\n" unless $start =~ /^\d+$/ && $stop =~ /^\d+$/;
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
	$dsg = generate_dsg(version=>$version, org_id=>$organism->id, gst_id=>1) if $organism && !$dsg; #gst_id 1 is for unmasked sequence data
	load_genomic_sequence(dsg=>$dsg, seq=>$entry->sequence, chr=>$chromosome);
	if ($GO)
	  {
	    my $load = 1;
	    foreach my $dsc ($dsg->dataset_connectors) #check to see if there is a prior link to the dataset -- this will happen when loading whole genome shotgun sequence
	      {
		$load = 0 if $dsc->dataset_id == $dataset->id;
	      }
	    $dsg->add_to_dataset_connectors({dataset_id=>$dataset->id}) if $load;
	    if ($dataset->version > $dsg->version)
	      {
		$dsg->version($dataset->version);
		$dsg->update;
	      }
	  }

	#format blastable db

      }
    print "completed parsing $accn!\n";# if $DEBUG;
  }
#need to add previous datasets if new dataset was added with a new dataset group
if ($GO)
  {

    if ($dsg && keys %previous_datasets)
      {
	my $ver; #need a higher version number than previous
	foreach my $ds (values %previous_datasets)
	  {
	    if ($dsg->dataset_connectors({dataset_id=>$ds->id}))
		{
		  my $name = $ds->name;
		  print "$name has been previously added to this dataset group.  Skipping\n";
		  next;
		}
	    foreach my $item ($ds->dataset_groups)
	      {
		$ver = $item->version unless $ver;
		$ver = $item->version if $item->version > $ver;
	      }
	    foreach my $chr($ds->chromosomes)
	      {
		load_genomic_sequence(dsg=>$dsg, seq=>$ds->genomic_sequence(chr=>$chr), chr=>$chr);
	      }
	    $dsg->add_to_dataset_connectors({dataset_id=>$ds->id});
	  }
	#incement and update dsg version if new version number is higher
	$ver++;
	if ($ver > $dsg->version)
	  {
	    $dsg->version($ver);
	    $dsg->update();
	  }
      }
  }

if (keys %dsg_update)
  {
    foreach my $item (values %dsg_update)
      {
	my $dsg = $item->{dsg};
	print "#"x20,"\n";
	print "Dataset Group ".$dsg->name." (".$dsg->id.") had a dataset deleted.  You will need to remove or update this dataset group.\n";
	print "accessions that were deleted and reloaded:\n";
	foreach my $accn (@{$item->{accn}})
	  {
	    print "\t", $accn,"\n";
#	    my $previous = check_accn($accn);
#	    my ($ds) = sort {$b->version<=>$a->version} @$previous;
#	    $dsg->add_to_dataset_connectors({dataset_id=>$ds->id});
	  }
	if ($autoupdate)
	  {
	    print "Autoupdate flag has been set to true.  This Dataset_group is being deleted from the database.\n";
	    delete_dataset_group($dsg);
	  }
	print "#"x20,"\n";
      }
  }

if ($GO && $dsg)
  {
    print "Creating blastable database\n";
    my $cmd = $formatdb." -i ".$dsg->file_path;
    print "\tFormatdb running $cmd\n";
    `$cmd`;
  }

if ($ERASE)
  {
    print "You are going to erase data!  Press any key to continue.  Control-c to abort!\n";
    <STDIN>;

    my %dsgs;
    $dsgs{$dsg->id}=$dsg if $dsg;
    foreach my $ds (values %previous_datasets)
      {
	my $dsg_flag;
	foreach my $item ($ds->dataset_groups)
	  {
	    $dsg_flag=1;
	    $dsgs{$item->id}=$item;
	  }
	$ds->delete unless $dsg_flag;
      }
    foreach my $item (values %dsgs)
      {
	delete_dataset_group($item);
      }
  }

sub load_genomic_sequence
  {
    my %opts = @_;
    my $seq = $opts{seq};
    my $chr = $opts{chr};
#    my $ds = $opts{ds};
    my $dsg = $opts{dsg};
    return unless $dsg;
    my $seqlen = length $seq;
    if (my ($item) = $dsg->genomic_sequences({chromosome=>$chr}))
	{
	  my $prev_length = $item->sequence_length;
	  print "$chr has previously beed added to this dataset_group.  Previous length: $prev_length.  Currently length: $seqlen.  Skipping.\n";
	  return;
	}
    print "Loading genomic sequence ($seqlen nt)\n";# if $DEBUG;
    $seq =~ s/\s//g;
    $seq =~ s/\n//g;
    $dsg->add_to_genomic_sequences({sequence_length=>$seqlen,
				    chromosome=>$chr,
				   }) if $GO;
    my $path = $dsg->file_path;
    $path =~ s/\/[^\/]*$/\//;
    mkpath($path);
    mkpath($path."/chr");
    #append sequence ot master file for dataset group
    open (OUT, ">>".$path."/".$dsg->id.".faa");
    my $head = $chr =~ /^\d+$/ ? ">gi" : ">lcl";
    $head .= "|".$chr;
    print OUT "$head\n$seq\n";
    close OUT;
    #create individual file for chromosome
    open (OUT, ">".$path."/chr/$chr");
    print OUT $seq;
    close OUT;
    #must add a feature of type chromosome to the dataset so the dataset "knows" its chromosomes
#    my $feat = get_feature(type=>"chromosome", name=>"chromosome $chr", ds=>$ds, chr=>$chr, start=>1, stop=>length($seq));
#    add_location(chr=>$chr, start=>1, stop=>length($seq), feat=>$feat, strand=>1);
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
	$stop = $start unless $stop;
	($start) = $start =~ /(\d+)/;
	($stop) = $stop =~ /(\d+)/;
	$start =~ s/\^.*//;
	$stop =~ s/\^.*//;
	$rstart = $start unless $rstart;
	$rstop = $stop unless $rstop;
	$rstart = $start if $start < $rstart;
	$rstop = $stop if $stop > $rstop;
      }
    return ($rstart, $rstop, $strand);
  }

sub get_organism
{
  my ($entry) = shift;
  my $name = $entry->data_source();
  if ($entry->strain)
    {
      my $strain = $entry->strain;
      $strain =~ s/\(/;/g;
      $strain =~ s/\)//g;
      $strain =~ s/strain://g;
      $strain =~ s/\\//g;
      $strain =~ s/=/;/g;
      $name =~ s/strain//;
      $name =~ s/\sstr\.?\s/ /;
      my @strains = split /;/,$strain;
      my @parse_strains;
      foreach (sort @strains)
	{
	  s/str\.\s//;
	  s/^\s+//;
	  s/\s+$//;
	  next unless $_;
	  $name =~ s/$_//;
	  push @parse_strains, $_;
	}
      $name =~ s/\s\s+/ /g;
      $name =~ s/\s+$//;
      my $add = join ("; ", sort @parse_strains);
      $name .= " strain" unless $add =~ /strain/;
      $name .= " $add";
    }
  if ($entry->substrain)
    {
      my $sstrain = $entry->substrain;
      $sstrain =~ s/\(/\\\(/g;
      $sstrain =~ s/\)/\\\)/g;
      $sstrain =~ s/substrain://g;
      $name =~ s/\ssubstr\.?\s/ /;
      $sstrain =~ s/\\//g;
      $sstrain =~ s/=/;/g;
      my @sstrains = split /;/,$sstrain;
      my @parse_sstrains;
      foreach (sort @sstrains)
	{
	  s/substr\.\s//;
	  s/^\s+//;
	  s/\s+$//;
	  next unless $_;
	  $name =~ s/$_//;
	  push @parse_sstrains, $_;
	}
      $name =~ s/\s\s+/ /g;
      $name =~ s/\s+$//;
      my $add = join ("; ", sort @parse_sstrains);
      $name .= " substrain" unless $add =~ /substrain/;
      $name .= " $add";
    }
  $name =~ s/'//g;
  $name =~ s/\(\s*\)//g;
  $name =~ s/\s\s+/ /g;
  $name =~ s/^\s+//;
  $name =~ s/\s+$//;
#  print $name,"\n";
#  print $entry->organism,"\n";
  my $desc = $entry->organism();
  $desc =~ s/^.*?::\s*//;
  $desc =~ s/\.$//;
  print qq{
Organism Information from Genbank Entry:
  $name
  $desc
} if $DEBUG;
  my $org = $coge->resultset('Organism')->find({name=>$name, description=>$desc});
  unless ($org)
    {
      unless ($name)
	{
	  print "WARNING: ", $entry->accession, " has no organism name\n";
	  return;
	}
      $org = $coge->resultset('Organism')->find_or_create(
							  {
							   name=>$name,
							   description=>$desc,
							  }) if $GO && $name;
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
						     });
}

sub generate_dsg
  {
    my %opts = @_;
    my $name = $opts{name};
    my $desc = $opts{desc};
    my $version = $opts{version};
    my $org_id = $opts{org_id};
    my $gst_id = $opts{gst_id};
    my $dsg_id = $opts{dsg_id};
    my $dsg = $dsg_id ? $coge->resultset('DatasetGroup')->find($dsg_id) :
      $coge->resultset('DatasetGroup')->create({name=>$name,
					       description=>$desc,
					       version=>$version,
					       organism_id=>$org_id,
					       genomic_sequence_type_id=>$gst_id,
					      }) if $GO;
    return unless $dsg;
    unless ($dsg->file_path)
      {
	my $path = "/opt/apache/CoGe/data/genomic_sequence/".$dsg->get_path."/".$dsg->id.".faa";
	print "Path to sequence: ", $path,"\n";
	$dsg->file_path($path);
	$dsg->update;
      }
    return $dsg;
  }

sub check_accn
  {
    my $accn = shift;
    my $gi = get_gi($accn);
    print "gi|".$gi."...";
    my $summary = get_gi_summary($gi);
#    print $summary,"\n";
    my ($version) = $summary =~ /gi\|\d+\|ref\|.*?\.(\d+)\|/i;
    my ($length) = $summary =~ /<Item Name="Length" Type="Integer">(\d+)<\/Item>/i;
    my ($taxaid) = $summary =~ /<Item Name="TaxId" Type="Integer">(\d+)<\/Item>/i;

    my @results;
    foreach my $ds ($coge->resultset('Dataset')->search({name=>$accn.".gbk"}))
      {
	my $version_diff = $ds->version eq $version ? 0 : 1;
	my $length_diff;
	my $cogelength;

	foreach my $feat ($ds->features({"feature_type.name"=>"chromosome"},
					{"join"=>"feature_type"},
					))
	  {
	    $cogelength += $feat->length;
	  }
	$length_diff = $cogelength eq $length ? 0 : 1;
	push @results, {ds=>$ds,
			version_diff=>$version_diff,
			length_diff=>$length_diff,
			coge_length =>$cogelength,
			ncbi_length => $length,

		       };
      }
    return \@results;
  }

sub delete_dataset_group
  {
    my $dsg = shift;
    my $path = $dsg->file_path;
    $path =~ s/[^\/]*$//;
    my $cmd = "rm -rf $path";
    print "Removing genomic sequence:  running: $cmd";
    `$cmd`;
    print "Deleting dataset group: ".$dsg->name,"\n";
    $dsg->delete();
  }

###NCBI eutils stuff

sub get_gi
  {
    my $accn = shift;
    my $esearch = "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=nucleotide&term=";
    my $result = get($esearch."$accn");
    my ($id) = $result =~ /<id>(.*?)<\/id>/i;
    return $id;
  }

sub get_gi_summary
  {
    my $gi = shift;
    my $esummary = "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=nucleotide&rettype=gbwithparts&retmode=text&complexity=0&id=";
    my $result = get($esummary.$gi);
    return $result;
  }

sub help
{
	print qq
	{
		Welcome to $0!  This program loads genbank entries into the CoGe genomes database

		Options:
	};
}
