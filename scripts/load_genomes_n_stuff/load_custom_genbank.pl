#!/usr/bin/perl -w

use strict;
use CoGe::Accessory::GenBank2;
use DBI;
use CoGeX;
use Roman;
use Data::Dumper;
use Getopt::Long;
use File::Path;

use vars qw($DEBUG $coge $GENOMIC_SEQ_LEN $GO $ERASE);

my ($gb_file, $gb_dir, $org_name, $org_desc, $org_id, $org_restricted, $source_name, $source_desc, $source_link, $source_id, $ds_name, $ds_desc, $ds_link, $ds_version, $ds_id, $chr, $seq_type_name, $seq_type_desc, $seq_type_id, $chr_basename, $add_chr_name, $use_fasta_header, $dsg_name, $dsg_desc, $dsg_version, $dsg_id, $restricted, $db, $user, $pass, $seq_dir, $tmp_dir, $chr_prefix);

##Example usage:
#./fasta_genome_loader.pl -org_name "Allenigales" -source_id 24 -ds_name oldsuper2.fasta -ds_version 2 -use_fasta_header -nt ~/projects/genome/data/Selaginella_moellendorffii/pre-v2/oldsuper2.fasta

GetOptions ( "debug=s" => \$DEBUG,
             "go=s"    => \$GO,
             "erase|e" => \$ERASE,
             "gb_file|file|gb=s" => \$gb_file,
             "gb_dir|dir=s"=>\$gb_dir,
             "org_name=s" => \$org_name,
             "org_desc=s" => \$org_desc,
             "org_id=s"   => \$org_id,
             "org_restricted" => \$org_restricted,#set flag to make organism restriced
             "source_name=s" => \$source_name, # datasource
             "source_desc=s" => \$source_desc,
             "source_link=s" => \$source_link,
             "source_id=s"   => \$source_id,
             "ds_name=s" => \$ds_name,# datasetid
             "ds_desc=s" => \$ds_desc,
             "ds_link=s" => \$ds_link,
             "ds_version=s" => \$ds_version,
             "ds_id=s"=>\$ds_id,
             "dsg_name=s" => \$dsg_name,
             "dsg_desc=s" => \$dsg_desc,
             "dsg_version=s" => \$dsg_version,
             "dsg_id=s"=>\$dsg_id,
             "restricted=i"=>\$restricted,

             "chr=s"=>\$chr,
             "seq_type_name=s" => \$seq_type_name,
             "seq_type_desc=s" => \$seq_type_desc,
             "seq_type_id=i"=>\$seq_type_id, # masked50 == id 2
             "chr_basename=s"=>\$chr_basename,
             "add_chr_name=s"=>\$add_chr_name,
             "use_fasta_header"=>\$use_fasta_header,
             "database|db=s"=>\$db,
            "user|u=s"=>\$user,
            "password|pw=s"=>\$pass,
             "seq_dir|sd=s"=>\$seq_dir, #the base level for where CoGe's genome sequences are stored
	     "tmp_dir|td=s"=>\$tmp_dir, #place to write temp files
	     "chr_prefix|cp=s"=>\$chr_prefix, #text to add before the chromosome name.  e.g. "contig_"
           );

$tmp_dir = "/tmp" unless $tmp_dir;
$DEBUG = 1 unless defined $DEBUG; # set to '1' to get updates on what's going on
$GO = 0 unless defined $GO; #set to 1 to actually make db calls.
($ds_name) = $gb_file=~ /([^\/]*)$/ unless ($ds_name);

$restricted = 0 unless defined $restricted;
print STDERR "Running $0\n";

my $formatdb =  "/usr/bin/formatdb -p F -o T"; #path to blast's formatdb program

my $connstr = "dbi:mysql:dbname=$db;host=localhost;port=3307" if $db;
$coge = CoGeX->connect($connstr, $user, $pass ) if $db;

if ($org_name && $GO)
  {
    my $org = $coge->resultset("Organism")->find_or_create({name=>$org_name,description=>$org_desc}) if $GO;
    if ($org_restricted && $GO)
      {
        $org->restricted(1);
        $org->update;
      }
    $org_id = $org->id;
  }

if ($source_name && $GO)
  {
    my $source = $coge->resultset("DataSource")->find_or_create({name=>$source_name,description=>$source_desc, link=>$source_link});
    $source_id = $source->id;
  }

unless (($org_id && $source_id))
  {
    print "Need a valid organism id and data source id in order to create dataset and dataset_group objects.  Loading aborted.\n\n";
    exit unless !$GO; #exit if this is a real run, otherwise continue and assume someone is testing the run
    print "Running in test mode!\n";
  }

if ($seq_type_name && $GO)
  {
    my $gst = $coge->resultset('GenomicSequenceType')->find_or_create({name=>$seq_type_name,
                                                                       description=>$seq_type_desc,
                                                                      });
    $seq_type_id = $gst->id;
  }
$seq_type_id = 1 unless $seq_type_id;  #default to unmasked sequence data

my $ds = generate_ds(ds_name => $ds_name,
                     ds_desc => $ds_desc,
                     ds_link => $ds_link,
                     ds_version => $ds_version,
                     ds_id =>$ds_id,
                     source_id=>$source_id,
                     restricted=>$restricted,
                    ) if $GO;

unless ($ds || !$GO)
  {
    warn "dataset object not initialized.  Exiting.";
    exit;
  }

$dsg_version = $ds->version unless $dsg_version || !$GO;
my $dsg = generate_dsg(name=>$dsg_name,
                       desc=>$dsg_desc,
                       version=>$dsg_version,
                       dsg_id=>$dsg_id,
                       org_id=>$org_id,
                       gst_id=>$seq_type_id,
                       restricted=>$restricted,
                       );
#link dataset and dataset_group

$coge->resultset('DatasetConnector')->find_or_create({dataset_id=>$ds->id, dataset_group_id=>$dsg->id}) if $GO;

if ($ERASE)
  {
    print "Clearing database of entries associated with ".$ds->name.". . .";
    $ds->delete();
    $dsg->delete();
    print "finished!\n";
    exit;
  }
print "dataset_id: ".$ds->id,"\n" if $ds;
print "dataset_group_id: ".$dsg->id,"\n\n\n" if $dsg;

process_gb (file=>$gb_file, ds=>$ds, dsg=>$dsg, dir=>$gb_dir, chr=>$chr, chr_prefix=>$chr_prefix) if $gb_file || $gb_dir;
print "dataset_id: ".$ds->id,"\n" if $ds;
print "dataset_group_id: ".$dsg->id,"\n" if $dsg;

sub process_gb
  {
    my %opts = @_;
    my $file =$opts{file};
    print $file,"\n";
    my $ds = $opts{ds};
    my $dsg = $opts{dsg};
    my $dir = $opts{dir};
    my $chr = $opts{chr};
    my $chr_prefix = $opts{chr_prefix};
    my @files;
    push @files, $file if $file && -r $file;
    if ($dir && -d $dir)
      {
        opendir (DIR, $dir);
        while (my $item = readdir(DIR))
          {
            next if $item =~ /^\.\.?$/;
            push @files, "$dir/$item" if -r "$dir/$item";
          }
        closedir(DIR);
      }
    foreach my $file (@files)
      {
        check_gb_file (file=>$file, ds=>$ds, chr=>$chr, dsg=>$dsg, chr_prefix=>$chr_prefix);
      }
    if ($GO)
      {
	my $cmd = $formatdb." -i ".$dsg->file_path;#."/".$dsg->id.".faa";
	print "\tFormatdb running $cmd\n";
	`$cmd`;
      }
  }

sub check_gb_file
  {
    #need to determine if this is a multi-entry gb file
    my %opts = @_;
    my $file =$opts{file};
    my $ds = $opts{ds};
    my $dsg = $opts{dsg};
    my $chr = $opts{chr};
    my $chr_prefix = $opts{chr_prefix};
    my $cmd = qq{grep "//" $file | wc -l};
    print "running $cmd\n";
    my $res = `$cmd`;
    chomp $res;
    my $files=[];
    if ($res > 1)
      {
	print "detected multi gb entries in file (count: $res).  Dividing into tmp files.\n";
	$files = divide_gb_file(file=>$file);
	print join ("\n", map {"\t".$_} @$files),"\n";
      }
    else
      {
	push @$files, $file;
      }
    foreach $file (@$files)
      {
	process_gb_file (file=>$file, ds=>$ds, chr=>$chr, dsg=>$dsg, chr_prefix=>$chr_prefix);
      }

  }

sub process_gb_file
  {
    my %opts = @_;
    my $file =$opts{file};
    my $ds = $opts{ds};
    my $dsg = $opts{dsg};
    my $chr = $opts{chr};
    my $chr_prefix = $opts{chr_prefix};
    my $gb = new CoGe::Accessory::GenBank2;
    $gb->parse_genbank_file(file=>$file);
    $chr = $gb->chromosome unless $chr;
    $chr = $gb->locus unless $chr;
    $chr = $chr_prefix.$chr if $chr_prefix;
    print "\n";
    print "Working on chromosome $chr (file: $file)\n";
    load_features(features=>$gb->features, chr=>$chr, ds=>$ds, dsg=>$dsg);
    load_genomic_sequence(dsg=>$dsg, seq=>$gb->sequence, chr=>$chr);

    if ($GO)
      {
	my $load = 1;
	foreach my $dsc ($dsg->dataset_connectors) #check to see if there is a prior link to the dataset -- this will happen when loading whole genome shotgun sequence
	  {
	    $load = 0 if $dsc->dataset_id == $ds->id;
	  }
	$dsg->add_to_dataset_connectors({dataset_id=>$ds->id}) if $load;
	if ($ds->version > $dsg->version)
	  {
	    $dsg->version($ds->version);
	    $dsg->update;
	  }
      }

  }

sub load_features
  {
    my %opts = @_;
    my $chr = $opts{chr};
    my $features = $opts{features};
    my $ds = $opts{ds};
    my $dsg = $opts{dsg};

    print "Processing features for $chr\n";
    my $count =1;
    foreach my $feature (@$features)
      {
	unless ($feature->type())
	  {
	    print "Feature has no feature type name \$feature->type():\n";
	    print Dumper $feature;
	    next;
	  }
	$count++;
	if ($feature->type() =~ /source/i)
	  {
	    #change source to chromosome
	    $feature->type('chromosome');
	    #		next;
	  }
	print "\tFeature Type: ".$feature->type,"\n" if $DEBUG;
	my $feat_type = $coge->resultset('FeatureType')->find_or_create({ name => $feature->type() })  if $GO;
	# create a db_feature for to link this feature with the dataset table
	my ($start, $stop, $strand) = get_feature_location($feature);
	my $db_feature = $coge->resultset('Feature')->create({
							      feature_type_id     => $feat_type->id,
							      dataset_id => $ds->id,
							      chromosome=> $chr,
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
								     name => $chr,
								     feature_id  => $db_feature->id
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

	    die "problem with $chr start $start or stop $stop\n" unless $start =~ /^\d+$/ && $stop =~ /^\d+$/;
	    print "\tAdding Location: $start - $stop\n" if $DEBUG;
	    my $location = $db_feature->add_to_locations(
							 {
							  start      => $start,
							  stop       => $stop,
							  strand     => $feature->strand,
							  chromosome => $chr
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
		my $master =1; #make first one master
		foreach my $item (@{$stuff})
		  {
		    foreach my $thing (split/;/,$item)
		      {
			$thing =~ s/^\s+//;
			$thing =~ s/\s+$//;
			$names{$thing}=0 unless defined $names{$thing};
			$names{$thing}=1 if $anno =~ /locus_tag/i; #primary_name;
		      }
		  }
	      }
	    elsif ($anno =~ /translation/i) # this needs to be entered into the sequence table
	      {
		next; #skip this.  Protein sequences are translated on the fly from DNA sequence
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
		    next unless $item;
		    print "\tAdding annotation $anno: $item\n" if $DEBUG;
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
	#add the names
	unless (keys %names)
	  { #no names!
	    my $name;
	    if ($dsg)
	      {
		foreach my $item (split /\s+/, $dsg->organism->name)
		  {
		    $item =~ /(^.)/g;
		    $name .= $1;
		  }
		$name .= "_";
	      }
	    $name .=$chr."_".$count;
	    $names{$name}=1;
	  }
	foreach my $name (keys %names)
	  {
	    my $master = $names{$name} ? 1 : 0;
	    $name =~ s/\s+$//g;
	    $name =~ s/^\s+//g;
	    print "\tAdding Name: $name\n" if $DEBUG;
	    my $feat_name = $db_feature->add_to_feature_names({
							       name       => $name,
							       primary_name => $master,
							      }) if $GO;
	  }
      }
    print "Processed $count features\n";
  }

sub divide_gb_file
  {
    my %opts = @_;
    my $file = $opts{file};
    my @files;
    $/ = "//\n";
    open (IN, $file);
    my $count = 1;
    my $base = int(rand(10000000));
    mkdir("$tmp_dir/$base");
    while (<IN>)
      {
	my $file = "$tmp_dir/$base/$count.gbk";
	push @files, $file;
	open (OUT, ">$file");
	print OUT $_;
	close OUT;
	$count++;
      }
    close IN;
    return \@files;
  }

sub generate_ds
  {
    my %opts = @_;
    my $ds_name = $opts{ds_name};
    my $ds_desc = $opts{ds_desc};
    my $ds_link = $opts{ds_link};
    my $ds_version = $opts{ds_version};
    my $ds_id = $opts{ds_id};
    my $source_id = $opts{source_id};
    my $restriced = $opts{restricted};
    unless ($ds_name || $ds_id)
      {
	warn "no dataset name or database id specified\n";
	return;
      }
    my $ds = $ds_id ? $coge->resultset('Dataset')->find($ds_id) :
      $coge->resultset('Dataset')->create({
					   name                => $ds_name,
					   description         => $ds_desc,
					   link                => $ds_link,
					   data_source_id      => $source_id,
					   restricted          => $restricted,
					   version=>$ds_version,
					  })  if $GO;
    return $ds;

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
    my $restricted = $opts{restricted} || 0;
    my $dsg = $dsg_id ? $coge->resultset('DatasetGroup')->find($dsg_id) :
      $coge->resultset('DatasetGroup')->create({name=>$name,
                                                description=>$desc,
                                                version=>$version,
                                                organism_id=>$org_id,
                                                genomic_sequence_type_id=>$gst_id,
                                                restricted=>$restricted,
                                               }) if $GO;
    return unless $dsg;
    unless ($dsg->file_path)
      {
        my $path = "$seq_dir/".$dsg->get_path."/".$dsg->id.".faa";
        print $path,"\n";
        $dsg->file_path($path);
        $dsg->update;
      }
    return $dsg;
  }

sub load_genomic_sequence
  {
    my %opts = @_;
    my $seq = $opts{seq};
    my $chr = $opts{chr};
    my $dsg = $opts{dsg};
    return unless $dsg;
    my $seqlen = length $seq;
    if (my ($item) = $dsg->genomic_sequences({chromosome=>$chr}))
	{
	  my $prev_length = $item->sequence_length;
	  print "$chr has previously been added to this dataset_group.  Previous length: $prev_length.  Currently length: $seqlen.  Skipping.\n";
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
