#!/usr/bin/perl -w

use strict;
use Data::Dumper;
use CoGeX;
use Getopt::Long;
use URI::Escape;

my $GO = 0;
my ($dbname, $dbuser, $dbpassword, $dsid, $ds_name, $ds_desc, $ds_link, $ds_version,
	$source_name, $source_desc, $source_link, $source_id, $gff_file,
	$anno_file, $DEBUG, $auto);
my $add_gene         = 0;
my $add_cds          = 0;
my $add_type_to_name = 0;
my @names;
my @skip_types;
my @anno_names;

#data transformations needed for CoGe:
#gene:  full extent of a gene model.  Beginning to end
#mRNA: UTRs and Codons
#CDS: coding sequence (codons)
#whatever else is fine

GetOptions(
	"db=s"				=> \$dbname,
	"user=s"			=> \$dbuser,
	"password=s"		=> \$dbpassword,
	"source_name=s"    	=> \$source_name,         # datasource
	"source_desc=s"    	=> \$source_desc,
	"source_link=s"    	=> \$source_link,
	"source_id=s"      	=> \$source_id,
	"ds_name=s"        	=> \$ds_name,             # datasetid
	"ds_desc=s"        	=> \$ds_desc,
	"ds_link=s"        	=> \$ds_link,
	"ds_version=s"     	=> \$ds_version,
	"dsid=i"           	=> \$dsid,
	"go=s"             	=> \$GO,
	"debug=s"          	=> \$DEBUG,
	"name=s"           	=> \@names,
	"anno_name=s"      	=> \@anno_names,
	"add_gene_feature" 	=> \$add_gene,
	"add_cds_feature"  	=> \$add_cds,
	"add_type_to_name" 	=> \$add_type_to_name,    #adds type (column 2) to name
	"skip_type=s"      	=> \@skip_types,
	"gff_file=s"       	=> \$gff_file,
	"anno_file=s"      	=> \$anno_file,
	"auto"             	=> \$auto,
);
$DEBUG = 1 unless defined $DEBUG;                #turn on by default
warn "-go flag is not true, nothing will be added to the database.\n" unless $GO;

my $connstr = "dbi:mysql:dbname=$dbname;host=localhost;port=PORT";
my $coge    = CoGeX->connect( $connstr, $dbuser, $dbpassword );
#$coge->storage->debugobj(new DBIxProfiler());
#$coge->storage->debug(1);

if ($source_name) {
	my $source = $coge->resultset("DataSource")->find_or_create(
		{
			name        => $source_name,
			description => $source_desc,
			link        => $source_link
		}
	);
	$source_id = $source->id;
}

my $ds = generate_ds(
	ds_name    => $ds_name,
	ds_desc    => $ds_desc,
	ds_link    => $ds_link,
	ds_version => $ds_version,
	ds_id      => $dsid,
	source_id  => $source_id,
);

unless ($ds) {
	warn "unable to find or create a valid dataset entry";
	exit;
}
print "Working on dataset: ", $ds->name . " (" . $ds->id . ")\n";

#let's get a list of chromosomes
print STDERR "Building list of valid chromosomes. . .";
my %valid_chrs = map { $_ => 1 } $ds->chromosomes;
print STDERR "finished.\n";

#some defaults to check for in names and annotations
push @names,      "ID";
push @names,      "name";
push @names,      "Name";
push @names,      "Alias";
push @names,      "gene";
push @names,      "Parent";
push @names,      "Locus_id";
push @names,      "ID_converter";
push @names,      "Gene_symbols";
push @names,      "gene_id";
push @skip_types, "Link_to";
push @skip_types, "References";
push @skip_types, "Sequence_download";
push @skip_types, "transcript_id";

push @anno_names, "Source";
push @anno_names, "Note";
push @anno_names, "NIAS_FLcDNA";
push @anno_names, "Comment";
push @anno_names, "GO";
push @anno_names, "ORF_evidence";          #can we link to SGD?
push @anno_names, "Transcript_evidence";
push @anno_names, "Status";
push @anno_names, "InterPro";
push @anno_names, "Description";
push @anno_names, "Function";
push @anno_names, "Derives_from";

my @skip_names_re = qw(
  :five_prime
  :three_prime
  :exon
  \.exon
  :utr
  \.utr
  :cds
  \.cds
  cds\.
  :hsp
  \.hsp

  intron
  _E\d
);

my %anno_names  = map { $_, 1 } @anno_names if @anno_names;
my %check_names = map { $_, 1 } @names;
my %skip_types  = map { $_, 1 } @skip_types;

my %data;
my %annos;
my %anno_name_lookup;
my %feat_types;    #store feature type objects

#my ($anno_type) = $coge->resultset('AnnotationType')->find_or_create({name=>"Phytozome link"}); #generic annotation type
my $prev_type;

#add information from annotation_file
process_annotation_file(
	file  => $anno_file,
	annos => \%annos,
	names => \%anno_name_lookup
  )
  if $anno_file;
open( IN, $gff_file ) || die "$!";

my %seen_types;

my $line_count = 0;
my $last_RNA   = "mRNA"; #storage for the last RNA type seen.  For converting exons to appropriate RNA type.
while (<IN>) {
	next if /^#/;
	next if /^Error/;
	chomp;
	next unless $_;
	$line_count++;
	print STDERR "processed $line_count lines\n" unless $line_count % 10000;

	#last if $line_count > 50000;
	my @line = split /\t/;
	next if $line[2] eq "clone";
	next if $line[2] eq "intron";
	next if $line[2] eq "chromosome";
	next if $line[2] eq "start_codon";
	next if $line[2] eq "stop_codon";

	#process and check chromosomes
	my $chr;
	$chr = $line[0];

	#    $chr =~ s/ig_//;
	$chr =~ s/%.*//;
	$chr =~ s/chromosome//i;
	$chr =~ s/chr//i;
	$chr =~ s/^_//i;
	$chr =~ s/^0//g;
	($chr) = split(/\s+/, $chr);
	unless ( $valid_chrs{$chr} ) {
		print "WARNING:  Chromosome '$chr' does not exist in the dataset! Pausing. . . (hit return to continue)\n";
		<STDIN> unless $auto;
		next;
	}

	#    next if $line[2] eq "gene";
	#    next if $line[2] eq "";
	next if $line[2] eq "transcript";
	next if $line[2] eq "protein";

#    $line[2] = "mRNA" if $line[2] eq "transcript";
#in many GFF files, the mRNA is what CoGe calls a Gene (the full extent of the transcribed sequence including introns and exons.  Instead, what the GFF calls an exon is really the transcribed mRNA.  In this cases, we want to hold the mRNA information to link to Parents and whatever annotation it contains, but don't want to actually add the location.  We will change the feature type to something weird that can be handled downstream correctly -- specifically the locations
	if ( $line[2] =~ /(.*RNA.*)/ ) {
		$last_RNA = $line[2];
		$line[2] = "$1_no_locs";
	}
	$line[2] = $last_RNA if $line[2] eq "exon";
	$line[2] = $last_RNA if $line[2] eq "five_prime_UTR";
	$line[2] = $last_RNA if $line[2] eq "three_prime_UTR";

	my %names;
	my $name;
	my ( $parent, $id );    #two important things in a GFF line to track;
	foreach my $item ( split(/;/, $line[-1]) ) {
		my $tmp;
		$item =~ s/"//g;
		$item =~ s/^\s+//;
		$item =~ s/\s+$//;
		next unless $item;
		my ( $type, $info ) = ( split(/[\s=]/, $item, 2) );
		$seen_types{$type}++;
		$parent = $info if $type eq "Parent";
		$id     = $info if $type eq "ID";
		next if $skip_types{$type};

		if ( $check_names{$type} ) {
		  outer: foreach my $item ( split(/,/, $info) ) {
				$names{$item} = 1;

 #these nexts will skip from using the primary name as the ID to the Parent name
				foreach my $re (@skip_names_re) {
					next outer if $item =~ /$re/i;
				}
				$name = $item unless $name;
				if ( $item =~ /^LOC_/ ) {
					my $tmp = $item;
					$tmp =~ s/^LOC_//;
					$names{$tmp} = 1;
				}
			}
		}
		next unless $name;    #No name, don't know what to do!
		$info = uri_unescape($info);    #remove URL formatting
		$annos{$name}{$info} = { type => $type } if $anno_names{$type};

		#add phytozome link:
		if ( $type eq "dbxref" ) {

#	    $info =~ s/SGD://;
#	    $annos{$name}{"View annotation at SGD"}={link=>$sgd_link."$info", type=>"SGD link"}
		}
	}

	next unless $name;                  #No name, don't know what to do!

	#add in additional names from the annotation file
	foreach my $tmp ( keys %names ) {
		next unless $anno_name_lookup{$tmp};
		foreach my $tmp2 ( keys %{ $anno_name_lookup{$tmp} } ) {
			$names{$tmp2} = 1;
		}
	}
	my $strand;
	$strand = -1 if $line[6] =~ /-/;
	$strand = 1  if $line[6] =~ /\+/;
	$strand = 0  if $line[6] =~ /\./;
	my @type = ( $line[2] );

	#    push @type, "CDS" if $add_cds && $type eq "mRNA";

	#phytozome replications of CDS to mRNA
	#    push @type, "mRNA" if $type eq "CDS";
	#    push @type, "mRNA" if $type =~ /UTR/;
	#replicate mRNA to gene
	#    push @type, "gene" if $type eq "mRNA";

	foreach my $tmp (@type) {
		my $tmp_name = $name;
		my $type     = $tmp;
		$type =~ s/_no_locs//;
		$tmp_name = $parent if $type eq "gene" && $parent;    #ugly hack

		#	print join ("\t", $name, $tmp_name),"\n";
		#initialize data structure

		$data{ $line[1] }{$chr}{$tmp_name}{$type} = {} unless $data{ $line[1] }{$chr}{$tmp_name}{$type};
		foreach my $n ( keys %names ) {
			$data{ $line[1] }{$chr}{$tmp_name}{$type}{names}{$n} = 1;
		}
		next if $tmp =~ /_no_locs/;  #skip adding locations for things like mRNA
		push @{ $data{ $line[1] }{$chr}{$tmp_name}{$type}{loc} },
		  {
			start  => $line[3],
			stop   => $line[4],
			strand => $strand,
			chr    => $chr,
		  };
	}

	#    print Dumper \%data;
}
if ($add_gene) {
	foreach my $source ( keys %data ) {
		foreach my $chr_loc ( keys %{ $data{$source} } ) {
		  name: foreach my $name ( keys %{ $data{$source}{$chr_loc} } ) {
				my $start;
				my $stop;
				my $strand;
				my $chr;
				my %names;
				foreach my $type ( keys %{ $data{$source}{$chr_loc}{$name} } ) {
					map { $names{$_} = 1 } keys %{ $data{$source}{$chr_loc}{$name}{$type}{names} };
					foreach my $loc ( @{ $data{$source}{$chr_loc}{$name}{$type}{loc} } )
					{
						next name if $type eq "gene";
						$start = $loc->{start} unless $start;
						$start = $loc->{start} if $loc->{start} < $start;
						$stop = $loc->{stop} unless $stop;
						$stop   = $loc->{stop}   if $loc->{stop} > $stop;
						$strand = $loc->{strand} if $loc->{strand};
						$chr    = $loc->{chr};
					}
					foreach my $loc ( @{ $data{$source}{$chr_loc}{$name}{$type}{loc} } )
					{
						$loc->{strand} = $strand;
					}
				}
				$data{$source}{$chr_loc}{$name}{gene}{loc} = [
					{
						start  => $start,
						stop   => $stop,
						strand => $strand,
						chr    => $chr,
					}
				];
				$data{$source}{$chr_loc}{$name}{gene}{names} = \%names;
			}
		}
	}
}

#print Dumper \%data;
#print Dumper \%annos;
#exit;

#time to load information into database

my %anno_types;    #hash to store annotation type objects

foreach my $source ( keys %data ) {
	foreach my $chr_loc ( sort { $a cmp $b } keys %{ $data{$source} } ) {
		foreach my $name ( sort { $a cmp $b } keys %{ $data{$source}{$chr_loc} } )
		{
			foreach my $feat_type ( sort { $a cmp $b } keys %{ $data{$source}{$chr_loc}{$name} } )
			{
				sleep 0.1;
				print "\n" if $DEBUG;
				my ($start) = sort { $a <=> $b } map  { $_->{start} } @{ $data{$source}{$chr_loc}{$name}{$feat_type}{loc} };
				my ($stop) = sort { $b <=> $a } map  { $_->{stop} } @{ $data{$source}{$chr_loc}{$name}{$feat_type}{loc} };
				my ($strand) = map { $_->{strand} } @{ $data{$source}{$chr_loc}{$name}{$feat_type}{loc} };
				my ($chr) = map { $_->{chr} } @{ $data{$source}{$chr_loc}{$name}{$feat_type}{loc} };
				$feat_types{$feat_type} = $coge->resultset('FeatureType')->find_or_create( { name => $feat_type } )
				  if $GO && !$feat_types{$feat_type};
				my $feat_type_obj = $feat_types{$feat_type};

				print "Creating feature of type $feat_type\n" if $DEBUG;

				my $feat = $ds->add_to_features(
					{
						feature_type_id => $feat_type_obj->id,
						start           => $start,
						stop            => $stop,
						chromosome      => $chr,
						strand          => $strand,
					}
				  )
				  if $GO;
				my $featid = $feat ? $feat->id : "no_go";
				my %seen_locs;
				my $loc_count = 0;

				foreach my $loc ( sort { $a->{start} <=> $b->{start} } @{ $data{$source}{$chr_loc}{$name}{$feat_type}{loc} } )
				{
					$loc_count++;
					next if $feat_type eq "gene" && $loc_count > 1; #only use the first one as this will be the full length of the gene.  Stupid hack

					next if $seen_locs{ $loc->{start} }{ $loc->{stop} };
					$seen_locs{ $loc->{start} }{ $loc->{stop} } = 1;
					print "Adding location $chr:(" . $loc->{start} . "-" . $loc->{stop} . ", $strand)\n" if $DEBUG;
					my $loc_tmp = $feat->add_to_locations(
						{
							start      => $loc->{start},
							stop       => $loc->{stop},
							strand     => $loc->{strand},
							chromosome => $loc->{chr}
						}
					  )
					  if $GO;
				}
				my $names = $data{$source}{$chr_loc}{$name}{$feat_type}{names};
				my %names;
				foreach my $name ( keys %$names ) {
					$names{$name} = 1;
				}
				my %seen_annos; #hash to store annotations so duplicates aren't added
			  master_names: foreach my $tmp ( keys %names ) {
					foreach my $re (@skip_names_re) {
						next master_names if $tmp =~ /$re/i;
					}
					my $master = 0;
					$master = 1 if $tmp eq $name;
					print "Adding name $tmp to feature ", $featid if $DEBUG;
					print " (MASTER)" if $master;
					print "\n"        if $DEBUG;

					my $feat_name = $feat->add_to_feature_names(
						{	name         => $tmp,
							primary_name => $master,
						}
					  )
					  if $GO;
					if ( $annos{$tmp} ) {
						foreach my $anno ( keys %{ $annos{$tmp} } ) {
							next unless $anno;
							next if $seen_annos{$anno};
							$seen_annos{$anno} = 1;
							my $type_name = $annos{$tmp}{$anno}{type} || "Note";
							my ($anno_type) = $anno_types{$type_name};
							unless ($anno_type) {
								($anno_type) = $coge->resultset('AnnotationType')->find_or_create( { name => $type_name } );
								$anno_types{$type_name} = $anno_type;
							}
							my $link = $annos{$tmp}{$anno}{link};
							print "Adding annotation ($type_name): $anno\n" if $DEBUG;
							print "\tlink: $link\n" if $DEBUG && $link;
							my $anno = $feat->add_to_annotations(
								{
									annotation         => $anno,
									link               => $link,
									annotation_type_id => $anno_type->id
								}
							  )
							  if $GO && $anno;
						}
					}
				}
			}
		}
	}
}

print "Completed working on dataset: ", $ds->name . " (" . $ds->id . ")\n";
close IN;
print "Seen data types:\n";
print join( "\n", map { $_ . "\t" . $seen_types{$_} } sort keys %seen_types ), "\n";

sub generate_ds {
	my %opts       = @_;
	my $ds_name    = $opts{ds_name};
	my $ds_desc    = $opts{ds_desc};
	my $ds_link    = $opts{ds_link};
	my $ds_version = $opts{ds_version};
	my $ds_id      = $opts{ds_id};
	my $source_id  = $opts{source_id};
	unless ( $ds_name || $ds_id ) {
		warn "no dataset name or database id specified\n";
		return;
	}
	my $ds = $ds_id
	  ? $coge->resultset('Dataset')->find($ds_id)
	  : $coge->resultset('Dataset')->find_or_create(
		{
			name           => $ds_name,
			description    => $ds_desc,
			link           => $ds_link,
			data_source_id => $source_id,
			version        => $ds_version,
		}
	  );
	return $ds;

}

sub process_annotation_file {
	print STDERR "Processing annotation file\n";
	my %opts       = @_;
	my $file       = $opts{file};
	my $annos      = $opts{annos};
	my $anno_names = $opts{names};
	open( IN, $file );
	while (<IN>) {
		chomp;
		next unless $_;
		my @line  = split /\t/;
		my $name  = $line[0];
		my $name2 = $line[0];
		$name2 =~ s/\.\d+$//;    #get rid of trailing version number if present

		#pfam
		unless ( !$line[1] || $line[1] =~ /no.*ids/ ) {
			foreach my $tmp ( split(/,/, $line[1]) ) {
				#print $name,"\t",q$tmp,"\n";
				my $link = "http://pfam.sanger.ac.uk/family/" . $tmp;
				$annos->{$name}{"View at Pfam ($tmp)"} =
				  { link => $link, type => "Pfam link" };
			}
		}

		#sorghum annotations
		#	unless ($line[2] =~ /no.*defline/)
		#	  {
		#	    foreach my $tmp (split/;/,$line[2])
		#	      {
		#		$tmp =~ s/^\s+//;
		#		$tmp =~ s/\s+$//;
		#		$tmp =~ s/^\[\s*(.*?)\s*\]$/$1/;
		#		$tmp =~ s/\s*,$//;
		#		$annos->{$name}{$tmp}={};
		#	      }
		#	  }
		unless ( !$line[2] || $line[2] =~ /no.*ids/ ) {
			foreach my $tmp ( split(/,/, $line[2]) ) {
				#print $name,"\t",$tmp,"\n";
				my $link = "http://www.pantherdb.org/panther/familyList.do?searchType=basic&fieldName=all&searchAll=true&listType=6&fieldValue=" . $tmp;
				$annos->{$name}{"View at PantherDB ($tmp)"} =
				  { link => $link, type => "PantherDB link" };
			}
		}
		unless ( !$line[3] || $line[3] =~ /no.*ids/ ) {
			foreach my $tmp ( split(/,/, $line[3]) ) {

#		print $name,"\t",$tmp,"\n";
#		my $link = "http://www.pantherdb.org/panther/familyList.do?searchType=basic&fieldName=all&searchAll=true&listType=6&fieldValue=".$tmp;
				$annos->{$name}{$tmp} = { type => "KOG" };
			}
		}
		unless ( !$line[4] || $line[4] =~ /no.*orthology/ ) {
			foreach my $tmp ( split(/,/, $line[4]) ) {
				my $link = "http://www.genome.jp/dbget-bin/www_bget?ko:" . $tmp;
				$annos->{$name}{"View at KEGG ($tmp)"} =
				  { link => $link, type => "KEGG link" };
			}
		}

		unless ( !$line[5] || $line[5] =~ /no.*orthology/ ) {
			foreach my $tmp ( split(/,/, $line[5]) ) {
				my $link = "http://www.genome.jp/dbget-bin/www_bget?ko:" . $tmp;
				$annos->{$name}{"View at KEGG ($tmp)"} =
				  { link => $link, type => "KEGG link" };
			}
		}

		#	unless ($line[6] =~ /no.*id/)
		#	  {
		#	    foreach my $tmp (split/,/,$line[6])
		#	      {
		#		my $link = "http://www.expasy.org/enzyme/".$tmp;
##		print $name,"\t",$tmp,"\t", $link, "\n";
	 #		$annos->{$name}{"View at ExPASy ($tmp)"}={link=>$link, type=>"EC link"};
	 #	      }
	 #	  }
		unless ( !$line[6] || $line[6] =~ /no.*hit/ ) {
			foreach my $tmp ( split(/,/, $line[6]) ) {
				my $link = "FeatView.pl?accn=" . $tmp;

				#print $name,"\t",$tmp,"\t", $link, "\n";
				$annos->{$name}{"$tmp"} =
				  { link => $link, type => "Best Arabidopsis Match" };
			}
		}
		unless ( !$line[7] || $line[7] =~ /no.*symbols/ )    ## skipping
		{
			foreach my $tmp ( split(/,/, $line[7]) ) {
				my $link = "FeatView.pl?accn=" . $tmp;

				#print $name,"\t",$tmp,"\t", $link, "\n";
				$annos->{$name}{"$tmp"} =
				  { link => $link, type => "Arabidopsis Symbol" };
			}
		}
		unless ( !$line[8] || $line[8] =~ /no.*hit/ ) {
			foreach my $tmp ( split(/;/, $line[8]) ) {
				$tmp =~ s/^\s+//;
				$tmp =~ s/\s+$//;

				#my $link = "/CoGe/FeatView.pl?accn=".$tmp;
				#print $name,"\t",$tmp, "\n";
				$annos->{$name}{"$tmp"} = { type => "Arabidopsis annotation" };
			}
		}
		unless ( !$line[9] || $line[9] =~ /no.*hit/ ) {
			foreach my $tmp ( split(/,/, $line[9]) ) {
				$tmp =~ s/LOC_//;
				my $link = "FeatView.pl?accn=" . $tmp;

				#print $name,"\t",$tmp,"\t", $link, "\n";
				$annos->{$name}{"$tmp"} =
				  { link => $link, type => "Best Rice Match" };
			}
		}
		unless ( !$line[10] || $line[10] =~ /no.*symbols/ )    ## skipping
		{
			foreach my $tmp ( split(/,/, $line[10]) ) {
				my $link = "FeatView.pl?accn=" . $tmp;

				#print $name,"\t",$tmp,"\t", $link, "\n";
				$annos->{$name}{"$tmp"} =
				  { link => $link, type => "Rice Symbol" };
			}
		}
		unless ( !$line[11] || $line[11] =~ /no.*hit/ ) {
			foreach my $tmp ( split(/;/, $line[8]) ) {
				$tmp =~ s/^\s+//;
				$tmp =~ s/\s+$//;

				#my $link = "/CoGe/FeatView.pl?accn=".$tmp;
				#print $name,"\t",$tmp, "\n";
				$annos->{$name}{"$tmp"} = { type => "Rice annotation" };
			}
		}

	}
	close IN;
}
