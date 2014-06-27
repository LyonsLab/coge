#!/usr/bin/perl -w

use strict;
use Data::Dumper;
use CoGeX;
use Getopt::Long;

my $GO = 0;
my $DEBUG = 1;
my ($dsid, $ds_name, $ds_desc, $ds_link, $ds_version, $source_name, $source_desc, $source_link, $source_id);
my $add_gene =0;
my $add_cds =0;
my $add_type_to_name =0;
my @names;
my @skip_anno_types;
my @skip_feat_types;
my @anno_names;
my $connstr = 'dbi:mysql:dbname=coge;host=HOST;port=PORT';
my$coge = CoGeX->connect($connstr, 'USER', 'PASSWORD' );
#$coge->storage->debugobj(new DBIxProfiler());
#$coge->storage->debug(1);

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
	     "skip_anno_type=s"=>\@skip_anno_types,
	     "skip_feat_type=s"=>\@skip_feat_types,
	   );
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
push @names, "mRNA" unless @names;
push @anno_names, "Description";
push @anno_names, "biotype";

my %anno_names = map {$_,1} @anno_names if @anno_names;
my %check_names = map {$_,1} @names;
my %skip_anno_types = map {$_,1} @skip_anno_types;
my %skip_feat_types = map {$_,1} @skip_feat_types;

warn "-go flag is not true, nothing will be added to the database.\n" unless $GO;
my %data;
my %annos;
my %master_names;
my %feat_types; #store feature type objects
my ($anno_type) = $coge->resultset('AnnotationType')->search({name=>"note"}); #generic annotation type
my $prev_type;
while (<>)
  {
    next if /^#/;
    chomp;
    next unless $_;
    my @line = split /\t/;
    next if $line[2] eq "clone";
#    next if $line[2] eq "mRNA";
    next if $line[2] eq "intron";
    next if $line[2] eq "chromosome";
    my $chr;
    $chr = $line[0];
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
	$info =~ s/\%([A-Fa-f0-9]{2})/pack('C', hex($1))/seg;
	next if $skip_anno_types{$type};
	if ($check_names{$type})
	  {
	    $names{$info} =1;
	    $name = $info unless $name;
	    if ($info =~ /\.\d+$/)
	      {
		my $tmp = $info;
		$tmp =~ s/\.\d$//;
		$names{$tmp}=1;
	      }
	    if ($info =~ /^LOC_/)
	      {
		my $tmp = $info;
		$tmp =~ s/^LOC_//;
		$names{$tmp}=1;
		$tmp =~ s/\.\d$//;
		$names{$tmp}=1;
	      }
	  }
	$line[2] = "Transposable Element" if $info =~ /transpos/i;

	next unless $name; #No name, don't know what to do!
	$annos{$name}{"$type: $info"}=1 if $anno_names{$type};
      }
    foreach my $i (keys %names)
      {
	foreach my $j (keys %names)
	  {
	    $master_names{$i}{$j}=1;
	    $master_names{$j}{$i}=1;
	  }
      }
    next unless $name; #No name, don't know what to do!
    next if $skip_feat_types{$line[2]};
    my $strand = $line[6] =~ /-/ ? -1 :1;
    my $type = ($line[2]);
    $type = "mRNA" if $type =~ /^exon$/i;
    $type = "mRNA" if $type =~ /^five_prime_UTR$/i;
    $type = "mRNA" if $type =~ /^three_prime_UTR$/i;
    my @type = ($type);
    push @type, "CDS" if $add_cds && $type eq "mRNA";
####add
    push @type, "mRNA" if $type eq "CDS";

    foreach my $type (@type)
      {
	push @{$data{$line[1]}{$name}{$type}{loc}}, {
						     start=>$line[3],
						     stop=>$line[4],
						     strand=>$strand,
						     chr=>$chr,
						    };
	map {$data{$line[1]}{$name}{$type}{names}{$_}=1} keys %names;
      }
#    print Dumper \%data;
  }
if ($add_gene)
  {
    foreach my $source (keys %data)
      {
      name: foreach my $name (keys %{$data{$source}})
	  {
	    my $start;
	    my $stop;
	    my $strand;
	    my $chr;
	    my %names;
	    foreach my $type (keys %{$data{$source}{$name}})
	      {
		map {$names{$_}=1} keys %{$data{$source}{$name}{$type}{names}};
		foreach my $loc (@{$data{$source}{$name}{$type}{loc}})
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
	    $data{$source}{$name}{gene}{loc}=[{
					       start=>$start,
					       stop=>$stop,
					       strand=>$strand,
					       chr=>$chr,
					      }];
	    $data{$source}{$name}{gene}{names}= \%names;
	  }
      }
  }
#print Dumper \%data;
#print Dumper \%annos;
#exit;

foreach my $source (keys %data)
  {
    foreach my $name (keys %{$data{$source}})
      {
	foreach my $feat_type (keys %{$data{$source}{$name}})
	  {
	    my ($start) = sort {$a<=>$b} map {$_->{start}} @{$data{$source}{$name}{$feat_type}{loc}};
	    my ($stop) = sort {$b<=>$a} map {$_->{stop}} @{$data{$source}{$name}{$feat_type}{loc}};
	    my ($strand) = map {$_->{strand}} @{$data{$source}{$name}{$feat_type}{loc}};
	    my ($chr) = map {$_->{chr}} @{$data{$source}{$name}{$feat_type}{loc}};
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
	    foreach my $loc (@{$data{$source}{$name}{$feat_type}{loc}})
	      {
		print "\tAdding location $chr:(".$loc->{start}."-".$loc->{stop}.", $strand)\n" if $DEBUG;
		my $loc_tmp = $feat->add_to_locations(
						      {
						       start      => $loc->{start},
						       stop       => $loc->{stop},
						       strand     => $loc->{strand},
						       chromosome => $loc->{chr}
						      }
						     ) if $GO;
	      }
	    my $names = $data{$source}{$name}{$feat_type}{names};
 	    my @names = keys %$names;
 	    foreach my $tmp (@names)
 	      {
 		foreach my $item (keys %{$master_names{$tmp}})
 		  {
 		    $names->{$item}=1;
 		  }
 	      }
	    foreach my $tmp (keys %{$names})
	      {
		print "\tAdding name $tmp to feature ", $featid ,"\n" if $DEBUG;
		my $feat_name = $feat->add_to_feature_names({
							     name=>$tmp,
							     #				   feature_id=>$featid,
							    }) if $GO ;
		if ($annos{$tmp})
		  {
		    foreach my $anno (keys %{$annos{$tmp}})
		      {
			print "\tAdding annotation $anno\n" if $DEBUG;
			my $annoo = $feat->add_to_annotations({annotation=>$anno, annotation_type_id => $anno_type->id}) if $GO && $anno;
		      }
		  }
	      }
	  }
      }
  }

print "Completed working on dataset: ", $ds->name. " (".$ds->id.")\n";

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
