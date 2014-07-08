#!/usr/bin/perl

# script name: genomes_dbfill.pl
# function: fills the genomes database with info from a genbank
# genbank file, uses GFDB::GBlite to parse, and Eric's db api to fill
# the db

# bct (01/15/06)
# updated and modified by ehl (Mar 06)

use DBI;
use strict;
use CoGe::Accessory::GBlite;
use CoGe::Genome;
use Roman;
use Data::Dumper;
use Getopt::Long;

my ($HELP, $DEBUG, $GO, $ERASE, @files);
GetOptions(
	   'help|h'=>\$HELP,
	   'debug|d|verbose|v'=>\$DEBUG,
	   'file|f=s'=>\@files,
	   'go|g'=>\$GO,
	   'delete|erase|e'=>\$ERASE,
	  );
$HELP = 1 unless @files;
$GO = 1 if $ERASE; #set to 1 to actually make db calls.
help() if $HELP;

# vars
my($statement) = "";
my $genomic_seq_len = 10000; #length to break up genomic sequence

my $genome = new CoGe::Genome();

# read in each file into a GBFile object
foreach my $longfile ( @files ) {

  # open our input data file...
  my $genbank = new CoGe::Accessory::GBlite( $longfile );
  # now remove the path from the file name
  my $file = `basename $longfile`;
  chomp $file;
  print "Processing $file... \n" if $DEBUG;
  my ($organism, $data_source, $data_information);

  # end of file-specific procession, not start on gb entries.  In most
  # cases, there will be just one entry per file.  This supports
  # multiple entries in one file.
  while(my $entry = $genbank->nextEntry) {
#    print Dumper $entry;
    unless ($organism && $data_source && $data_information)
      {
	# first do the file-specific stuff
	$organism = $genome->get_organism_obj->find_or_create(
								 {
								  name=>$entry->source(),
#								  name=>$entry->organism(),
								  description=>$entry->organism_long()
								 })  if $GO;
	if ( $DEBUG ) {
	  print "organism_obj FoC: name - ", $entry->organism(), "\n";
	  print "organism_obj FoC: description - ", $entry->organism_long(), "\n";
	  print "-="x30, "\n";
	  print "organism_obj returned: ", $organism->id(), "\n" if $organism;
	}

	$data_source = $genome->get_data_source_obj()->find_or_create(
									 {
									  name=>'NCBI',
									  description=>"National Center for Biotechnology Information",
									  link=>'www.ncbi.nih.gov'
									 }) if $GO;
	if ( $DEBUG ) {
	  print "data_source_obj FoC: name - NCBI\n";
	  print "data_source_obj FoC: description - National Center for Biotechnology Information\n";
	  print "data_source_obj FoC: link - www.ncbi.nih.gov\n";
	  print "-="x30, "\n";
	  print "data_source_obj returned: ", $data_source->id(), "\n" if $data_source;
	}

	my $data_information_desc = "LOCUS: "     . $entry->locus();
	$data_information_desc   .= ", ACCESSION: " . $entry->accession();
	$data_information_desc   .= ", VERSION: "   . $entry->version();

	$data_information = $genome->get_data_information_obj()->find_or_create(
										{
										 name                => $file,
										 description                => $data_information_desc,
										 link                => $longfile, # run this from the same level
										 # as the file you're parsing,
										 # the link will be correctly set
										 organism_id         => $organism->id,
										 data_source_id      => $data_source->id(),
										 version=>$entry->version,
										})  if $GO;
	if ( $DEBUG ) {
	  print "data_information_obj FoC: name - ", $file, "\n";
	  print "data_information_obj FoC: description - ", $data_information_desc, "\n";
	  print "data_information_obj FoC: link - ", $longfile, "\n";
	  print "data_information_obj FoC: organism_id - ", $organism->id(), "\n" if $organism;
	  print "data_information_obj FoC: data_source_id - ", $data_source->id(), "\n" if $data_source;
	  print "-="x30, "\n";
	  print "data_information_obj returned: ", $data_information->id(), "\n" if $data_information;
	}
	if ($ERASE)
	  {
	    print "Clearing database of entries associated with ".$data_information->name.". . .";
	    $data_information->delete();
	    print "finished!\n";
	    exit;
	  }
      }
    my $chromosome = "1";

    ### Main Feature Processing Loop ###
    foreach my $feature ( @{ $entry->features() } ) {
      # start with getting a hook in the db for this feature type obj
      # $feature->key() should contain something like "mRNA" or "CDS"
      my $feat_type = $genome->get_feature_type_obj->find_or_create( { name => $feature->key() } )  if $GO;
      if ( $DEBUG ) {
        print "feature_type_obj FoC: name - ", $feature->key(), "\n";
        print "-="x30, "\n";
        print "feature_type_obj returned: ", $feat_type->id(), "\n" if $feat_type;
      }

      # create a db_feature for to link this feature with the
      # data_information table
      my $db_feature = $genome->get_feature_obj->create(
        {
          feature_type_id     => $feat_type->id,
          data_information_id => $data_information->id,
        }
      ) if $GO;
      if ( $DEBUG ) {
        print "feature_obj C: feature_type_id - ", $feat_type->id(), "\n" if $feat_type;
        print "feature_obj C: data_information_id - ", $data_information->id(), "\n" if $data_information;
        print "-="x30, "\n";
        print "feature_obj returned: ", $db_feature->id(), "\n" if $db_feature;
      }
      #expect first feature to be the source feature!
      if ( $feature->key() =~ /source/i ) {
        my $quals = $feature->qualifiers(); # get the source qualifiers

        if ( exists $quals->{chromosome} ) {
          $chromosome = $quals->{chromosome};
	  $chromosome =~ s/\s//g;
          $chromosome = arabic( $chromosome ) if isroman($chromosome); #if roman, convert to arabic
        }
        print "FEATURE: source - chromosome = $chromosome\n" if $DEBUG;

        #generate name based on organism name and chromosome
        my $feat_name = $genome->get_feature_name_obj->create(
							      {
							       name        => $organism->name,
							       description => "Chromosome " . $chromosome,
							       feature_id  => $db_feature->id
							      }
							     ) if $GO;
        if ( $DEBUG ) {
          print "feature_name_obj C: name - ", $organism->name(), "\n" if $organism;
          print "feature_name_obj C: description - ", "Chromosome $chromosome", "\n";
          print "feature_name_obj C: feature_id - ", $db_feature->id(), "\n" if $db_feature;
          print "-="x30, "\n";
          print "feature_name_obj returned: ", $feat_name->id(), "\n" if $feat_name;
        }

        #generate name for accession
        $feat_name = $genome->get_feature_name_obj->create(
							   {
							    name       => $entry->accession,
							    feature_id => $db_feature->id()
							   }
							  ) if $GO;
        if ( $DEBUG ) {
          print "feature_name_obj C: name - ", $entry->accession(), "\n";
          print "feature_name_obj C: feature_id - ", $db_feature->id(), "\n" if $db_feature;
          print "-="x30, "\n";
          print "feature_name_obj returned: ", $feat_name->id(), "\n" if $feat_name;
        }

        #generate name for version
        $feat_name = $genome->get_feature_name_obj->create(
							   {
							    name       => $entry->accession.".".$entry->version,
							    feature_id => $db_feature->id
							   }
							  ) if $GO;
        if ( $DEBUG ) {
          print "feature_name_obj C: name - ", $entry->accession.".".$entry->version(), "\n";
          print "feature_name_obj C: feature_id - ", $db_feature->id(), "\n" if $db_feature;
          print "-="x30, "\n";
          print "feature_name_obj returned: ", $feat_name->id(), "\n" if $feat_name;
        }
        #generate name for GI
        $feat_name = $genome->get_feature_name_obj->create(
							   {
							    name       => "GI:".$entry->gi,
							    feature_id => $db_feature->id
							   }
							  ) if $GO;
        if ( $DEBUG ) {
          print "feature_name_obj C: name - ", "GI:".$entry->gi(), "\n";
          print "feature_name_obj C: feature_id - ", $db_feature->id(), "\n" if $db_feature;
          print "-="x30, "\n";
          print "feature_name_obj returned: ", $feat_name->id(), "\n" if $feat_name;
        }
      }

      # add a location entry
      my $loc_string = $feature->location;
      $loc_string =~ s/complement//;
      $loc_string =~ s/join//;
      $loc_string =~ s/\(|\)//g;
      foreach my $loc (split /,/,$loc_string)
	{
	  $loc =~ s/<|>//g;
	  my ($start, $stop) = split /\.\./, $loc;
	  $stop = $start unless $stop;
	  die "problem with start $start or stop $stop\n" unless $start =~ /^\d+$/ && $stop =~ /^\d+$/;
	  my $location = $genome->get_location_obj->create(
							   {
							    feature_id => $db_feature->id,
							    start      => $start,
							    stop       => $stop,
							    strand     => $feature->strand,
							    chromosome => $chromosome
							   }
							  ) if $GO;
	  if ( $DEBUG ) {
	    print "location_obj C: feature_id - ", $db_feature->id(), "\n" if $db_feature;
	    print "location_obj C: start - ", $start, "\n";
	    print "location_obj C: stop - ", $stop, "\n";
	    print "location_obj C: strand - ", $feature->strand(), "\n";
	    print "location_obj C: chromosome - ", $chromosome, "\n";
	    print "-="x30, "\n";
	    print "location_obj returned: ", $location->id(), "\n" if $location;
	  }
	}
      # now work through the qualifiers for this feature
      # start by getting the hashref of qualifiers
      my $annot = $feature->qualifiers();
      foreach  my $anno ( keys %{ $annot } ) {
        #deal with db_xref: (taxon:3702) (GeneID:821318) (GI:18379324)
        if ( $anno =~ /xref/i ) {
          print "Processing db_xref qualifier...\n" if $DEBUG;
          my $anno_type_group =
            $genome->get_annotation_type_group_obj->find_or_create( { name => $anno } )  if $GO;
          if ( $DEBUG ) {
            print "annotation_type_group_obj FoC: name - ", $anno, "\n";
            print "-="x30, "\n";
            print "annotation_type_group_obj returned: ", $anno_type_group->id(), "\n" if $anno_type_group;
          }

          # go through each of the entries in the db_xref qualifier
          # values and split on ':', then add entries individually
          my @xrefs = split(/ /, $annot->{$anno}); # split values on \s
          foreach my $xref ( @xrefs ) {
            my @inner = split(/:/, $xref );
            # first add the annot_type_obj
            my $anno_type = $genome->get_annotation_type_obj->find_or_create(
              {
                name                     => $inner[0],
                annotation_type_group_id => $anno_type_group->id(),
              }
            )  if $GO;
            if ( $DEBUG ) {
              print "annotation_type_obj FoC: name - ", $inner[0], "\n";
              print "annotation_type_obj FoC: annotation_type_group_id - ",
                      $anno_type_group->id(), "\n" if $anno_type_group;
              print "-="x30, "\n";
              print "annotation_type_obj returned: ", $anno_type->id(), "\n" if $anno_type;
            }

            # now create the row for the data value of the xref
            my $sub_anno = $genome->get_annotation_obj->create(
              {
                annotation         => $inner[1],
                feature_id         => $db_feature->id(),
                annotation_type_id => $anno_type->id()
              })  if $GO;
            if ( $DEBUG ) {
              print "annotation_obj C: annotation - ", $inner[1], "\n";
              print "annotation_obj C: feature_id - ", $db_feature->id(), "\n" if $db_feature;
              print "-="x30, "\n";
              print "annotation_obj returned: ", $sub_anno->id(), "\n" if $sub_anno;
            }
          } # end of processing db_xref
        }
        elsif (
             $anno =~ /locus_tag/i
          || $anno =~ /transcript_id/i
          || $anno =~ /protein_id/i
          || $anno =~ /gene/i
          || $anno =~ /synonym/i    ##synonyms are embedded in the /note= tag!
          )                         #these are names
        {
          my $feat_name = $genome->get_feature_name_obj->create(
            {
              name       => $annot->{$anno},
              feature_id => $db_feature->id(),
            }
          ) if $GO;
          if ( $DEBUG ) {
            print "feature_name_obj C: name - ", $annot->{$anno}, "\n";
            print "feature_name_obj C: feature_id - ", $db_feature->id(), "\n" if $db_feature;
            print "-="x30, "\n";
            print "feature_name_obj returned: ", $feat_name->id(), "\n" if $feat_name;
          }
        }
        elsif ( $anno =~ /translation/i ) #this needs to be entered into the sequence table
        {
          my $seq_type = $genome->get_sequence_type_obj->find_or_create(
            {
              name        => "protein",
              description => "translation",
            }
          ) if $GO;
	  if ( $DEBUG )
	    {
	      print "seq_type_obj C: name - protein\n";
	      print "seq_type_obj C: description - translation\n";
	      print "-="x30, "\n";
	      print "seq_type_obj returned: ", $seq_type->id(), "\n" if $seq_type;
	    }
          my $sequence = $genome->get_sequence_obj->create(
            {
              sequence_type_id => $seq_type->id(),
              sequence_data    => $annot->{$anno},
              feature_id       => $db_feature->id(),
            }
          ) if $GO;
	  if ( $DEBUG )
	    {
	      print "seq_obj C: data = ".$annot->{$anno}."\n";
	      print "seq_obj C: sequence_type_id ", $seq_type->id(), "\n" if $seq_type;
	      print "seq_obj C: feature_id - ", $db_feature->id(), "\n" if $db_feature;
	      print "-="x30, "\n";
	      print "seq_obj returned: ", $sequence->id(), "\n" if $sequence;
	    }
        }
        elsif ( $anno eq "note")
        {
          # if go annot are present, they'll be in the note qualifier,
          # so process is specifically
          my $leftover = "";
          my @temp = split( /;/, $annot->{$anno} );
          foreach my $go_raw ( @temp )
          {
            #                  $1        $2            $3
            if ( $go_raw =~ /go_/ ) {
	      while ($go_raw =~ /(go_.*?):\s+(.*?)\[goid G?O?:?(.*?)\]/g) {
		#example:
		#go_function: nucleic acid binding [goid 0003676]
		my $anno_type_group =
		  $genome->get_annotation_type_group_obj->find_or_create( { name => $1 } ) if $GO;
		if ( $DEBUG )
		  {
		    print "annotation_type_group_obj FoC: name - ", $1, "\n";
		    print "-="x30, "\n";
		    print "annotation_type_group_obj returned: ", $anno_type_group->id(), "\n" if $anno_type_group;
		  }

		# $1 should be "go_function"
		my $anno_type = $genome->get_annotation_type_obj->find_or_create(
										 {
										  name => $3,    #this should be "0003676"
										  annotation_type_group_id => $anno_type_group->id(),
										 }
										) if $GO;
		if ( $DEBUG ) {
		  print "annotation_type_obj FoC: name - ", $3, "\n";
		  print "annotation_type_obj FoC: annotation_type_group_id - ",
		    $anno_type_group->id(), "\n" if $anno_type_group;
		  print "-="x30, "\n";
		  print "annotation_type_obj returned: ", $anno_type->id(), "\n" if $anno_type;
		}

		my $sub_anno = $genome->get_annotation_obj->create(
								   {
								    annotation => $2,    #this should be "nucleic acid binding"
								    feature_id         => $db_feature->id,
								    annotation_type_id => $anno_type->id
								   }
								  )  if $GO;
		if ( $DEBUG ) {
		  print "annotation_obj C: annotation - ", $2, "\n";
		  print "annotation_obj C: feature_id - ", $db_feature->id(), "\n" if $db_feature;
		  print "-="x30, "\n";
		  print "annotation_obj returned: ", $sub_anno->id(), "\n" if $sub_anno;
		}
	      }
            } else {
              $leftover .= " " . $go_raw if $go_raw;
            }
            # now just add the note remainder
	    $leftover =~ s/^\s+//;
	    $leftover =~ s/\s+$//;
	    if ($leftover)
	      {
		my $anno_type =
                  $genome->get_annotation_type_obj->find_or_create( { name => $anno } ) if $GO;
		if ( $DEBUG ) {
		  print "annotation_type_obj FoC: name - ", $anno, "\n";
		  print "-="x30, "\n";
		  print "annotation_type_obj returned: ", $anno_type->id(), "\n" if $anno_type;
		}

		my $sub_anno = $genome->get_annotation_obj->find_or_create(
									   {
									    annotation         => $leftover,
									    feature_id         => $db_feature->id(),
									    annotation_type_id => $anno_type->id(),
									   }
									  ) if $GO;
		if ( $DEBUG ) {
		  print "annotation_obj C: annotation - ", $leftover, "\n";
		  print "annotation_obj C: feature_id - ", $db_feature->id(), "\n" if $db_feature;
		  print "-="x30, "\n";
		  print "annotation_obj returned: ", $sub_anno->id(), "\n" if $sub_anno;
		}
	      }
          }
        }
        else           ##everything else
        {
          my $anno_type = $genome->get_annotation_type_obj->find_or_create( { name => $anno, } )  if $GO;
	  if ( $DEBUG ) {
	    print "annotation_type_obj FoC: name - ", $anno, "\n";
	    print "-="x30, "\n";
	    print "annotation_type_obj returned: ", $anno_type->id(), "\n" if $anno_type;
	  }
	  my $sub_anno = $genome->get_annotation_obj->create(
							     {
							      annotation         => $annot->{$anno},
							      feature_id         => $db_feature->id(),
							      annotation_type_id => $anno_type->id(),
							     }
							    ) if $GO;
	  if ( $DEBUG ) {
	    print "annotation_obj C: annotation - ", $annot->{$anno}, "\n";
	    print "annotation_obj C: feature_id - ", $db_feature->id(), "\n" if $db_feature;
	    print "-="x30, "\n";
	    print "annotation_obj returned: ", $sub_anno->id(), "\n" if $sub_anno;
	  }
	}
      }
      print "\n" if $DEBUG;
    } #end "feature" while
    print "Processing Genomic Sequence. . .\n" if $DEBUG;
    load_genomic_sequence(len=> $genomic_seq_len, di=>$data_information, seq=>$entry->sequence, chr=>$chromosome);
    print "\n";
  } #end "entry" while
  print "completed parsing $file!\n" if $DEBUG;
} # end "file" while

sub load_genomic_sequence
  {
    my %opts = @_;
    my $seq = $opts{seq};
    my $len = $opts{len};
    my $chr = $opts{chr};
    my $di = $opts{di};
    my $seqlen = length $seq;
    print "Loading genomic sequence ($seqlen nt)\n" if $DEBUG;
    $seq =~ s/\s//g;
    $seq =~ s/\n//g;
    my $i = 0;
    my $gso = $genome->get_genomic_sequence_obj;
    while ($i < $seqlen)
      {
	my $str = substr($seq,$i,$len);
	my $start = $i+1;
	my $stop = $i + length $str;
	$gso->create({start=>$start,
                      stop=>$stop,
                      chromosome=>$chr,
                      sequence_data=>$str,
                      data_information_id=>$di->id,
                  }) if $GO;
	$i += $len;
      }
  }

sub help
  {
    print qq{
Welcome to $0;

This program load genbank files into the CoGe Genomes database.  WARNING: You may
need to customize this program for your particular genbank file in order to properly
handle  specific annotation types such as Geneontology or functional domains.

Options                (Valid values)

-help | -h             (0|1) Print this message

-file | -f             (String) Path to genebank genome file for loading into the
                       CoGe system

-go   | -g             (0|1) You must set this to 1 in order for the data to be loaded
                       (NOTE:  this is a failsafe switch to make sure you are ready!)

-delete | -erase | -e  (0|1) Erase the data instead of loading the data.  You will want
                       to set this to 1 if you had an abortive data-loading run.

-debug   | -d          (0|1) Verbose debugging output

-verbose | -v          (0|1) Verbose debugging output
};
    exit;
  }

__END__
