#!/usr/bin/perl

use DBI;
use strict;
use CoGeX;
use CoGeX::DataSource;
use CoGeX::Dataset;
use CoGeX::Organism;
use CoGeX::GenomicSequence;
use Roman;
use Data::Dumper;
use Getopt::Long;

use vars qw($DEBUG $coge $GENOMIC_SEQ_LEN $GO $ERASE);

my ($nt_file, $nt_dir, $org_name, $org_desc, $org_id, $ds_name, $ds_desc, $ds_link, $ds_id, $di_name, $di_desc, $di_link, $di_version, $di_id, $use_contigs_as_features, $chr, $seq_type_name, $seq_type_desc, $seq_type_id, $chr_basename, $add_chr_name, $use_fasta_header);

GetOptions ( "debug=s" => \$DEBUG,
	     "go=s"    => \$GO,
	     "erase|e" => \$ERASE,
	     "fasta_file|fasta|faa|nt_file|nt=s" => \$nt_file,
	     "fasta_dir|nt_dir|dir=s"=>\$nt_dir,
	     "org_name=s" => \$org_name,
	     "org_desc=s" => \$org_desc,
	     "org_id=s"   => \$org_id,
	     "ds_name=s" => \$ds_name, # datasource
	     "ds_desc=s" => \$ds_desc,
	     "ds_link=s" => \$ds_link,
	     "ds_id=s"   => \$ds_id,
	     "di_name=s" => \$di_name,# datasetid
	     "di_desc=s" => \$di_desc,
	     "di_link=s" => \$di_link,
	     "di_version=s" => \$di_version,
	     "di_id=s"=>\$di_id,
	     "chr=s"=>\$chr,
	     "use_contigs_as_features=i" => \$use_contigs_as_features,
	     "seq_type_name=s" => \$seq_type_name,
	     "seq_type_desc=s" => \$seq_type_desc,
	     "seq_type_id=i"=>\$seq_type_id, # masked50 == id 2
	     "chr_basename=s"=>\$chr_basename,
	     "add_chr_name=s"=>\$add_chr_name,
	     "use_fasta_header"=>\$use_fasta_header
	   );

$DEBUG = 1 unless defined $DEBUG; # set to '1' to get updates on what's going on
$GO = 1 unless defined $GO; #set to 1 to actually make db calls.

$use_contigs_as_features =0 unless defined $use_contigs_as_features;

my $GENOMIC_SEQ_LEN = 10000; #length to break up genomic sequence
my $connstr = 'dbi:mysql:dbname=genomes;host=HOST;port=PORT';
$coge = CoGeX->connect($connstr, 'USER', 'PASSWORD' );
#$coge->storage->debugobj(new DBIxProfiler());
#$coge->storage->debug(1);

my $di = $di_id ? $coge->resultset('Dataset')->find($di_id) : generate_di(
									  org_name => $org_name,
									  org_desc => $org_desc,
									  org_id=>$org_id,
									  ds_name => $ds_name,
									  ds_desc => $ds_desc,
									  ds_link => $ds_link,
									  ds_id=>$ds_id,
									  di_name => $di_name,
									  di_desc => $di_desc,
									  di_link => $di_link,
									  di_version => $di_version,
									 );
unless ($di)
  {
    warn "dataset object not initialized.  Exiting.";
    exit;
  }

my $gstype = generate_gstype(id=>$seq_type_id, name=>$seq_type_name, desc=>$seq_type_desc);
unless ($gstype)
  {
    warn "genomic sequence type object not initialized.  Exiting.";
    exit;
  }

if ($ERASE)
  {
#    my $di = $DB->get_dataset_obj->resolve_dataset($di;
    print "Clearing database of entries associated with ".$di->name.". . .";
    $di->delete();
    print "finished!\n";
    exit;
  }

process_nt (file=>$nt_file, di=>$di, dir=>$nt_dir, chr=>$chr) if $nt_file || $nt_dir;
sub add_location
  {
    my %opts = @_;
    my $chr=$opts{chr};
    my $start=$opts{start};
    my $stop=$opts{stop};
    my $strand=$opts{strand};
    my $feat=$opts{feat};
    print "Adding location $chr:($start-$stop, $strand)\n" if $DEBUG;
    my $loc = $coge->resultset('Location')->create(
						   {
						    feature_id => $feat->id,
						    start      => $start,
						    stop       => $stop,
						    strand     => $strand,
						    chromosome => $chr
						   }
						  ) if $GO;
    return $loc;
  }

sub get_feature
  {
    my %opts = @_;
    my $type = $opts{type};
    my $name = $opts{name};
    my $chr = $opts{chr};
    my $start = $opts{start};
    my $stop = $opts{stop};
    print "Working on $name of type $type\n" if $DEBUG;
    my $di = $opts{di};
    my $feat_type = $coge->resultset('FeatureType')->find_or_create( { name => $type } ) if $GO;
    print "Creating feature of type $type\n" if $DEBUG;
    my $feat = $di->add_to_features ({
				   feature_type_id => $feat_type->id,
				   #dataset_id => $di->id,
				   start=>$start,
				   stop=>$stop,
				   chromosome=>$chr,
				   strand=>1,
				  }) if ($GO);
    if ($DEBUG)
      {
	print "Creating feature_name $name for feat ";
	print $feat->id if $feat;
	print "\n";
      }
    my $feat_name = $coge->resultset('FeatureName')->create({
							     name=>$name,
							     feature_id=>$feat->id,
							    }) if $GO;
    return $feat;
  }

sub process_nt
  {
    my %opts = @_;
    my $file =$opts{file};
    print $file,"\n";
    my $di = $opts{di};
    my $dir = $opts{dir};
    my $chr = $opts{chr};
    my @files;
    push @files, $file if $file && -r $file;
    if (-d $dir)
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
	process_nt_file (file=>$file, di=>$di, chr=>$chr);
      }
  }

sub process_nt_file
  {
    my %opts = @_;
    my $file =$opts{file};
    my $di = $opts{di};
    my $chr = $opts{chr};
    return unless $file;
    print "Processing $file\n";
    $/ = "\n>";
    open (IN, $file) || die "can't open $file for reading: $!";
    while (<IN>)
      {
	s/>//g;
	my ($name, $seq) = split /\n/, $_,2;
	$seq =~ s/\n//g;

	my $chrtmp;
	$chrtmp = $chr if defined $chr;
	$chrtmp = $name if $use_fasta_header;
	($chrtmp) = $name=~/^(\S+)/ unless $chrtmp || $add_chr_name;
	$chrtmp = $name unless defined $chrtmp;
	$chrtmp =~ s/chromosome//;
	$chrtmp =~ s/chr//;
	$chrtmp =~ s/^0+//;
	$chrtmp =~ s/^_+//;
	$chrtmp =~ s/^\s//;
	$chrtmp =~ s/\s$//;
	$chrtmp =0 unless $chrtmp;
	$chrtmp = $chr_basename.$chrtmp if $chr_basename;
        #($chrtmp) = $name =~ /(^\S+)/;
	load_genomic_sequence(chr=>$chrtmp,
			      seq=>$seq,
			      di=>$di,
			      );
      }
    close IN;
  }

sub load_genomic_sequence
  {
    my %opts = @_;
    my $seq = $opts{seq};
    my $len = $opts{len};
    my $chr = $opts{chr};
    my $di = $opts{di};
    $len = $GENOMIC_SEQ_LEN unless $len;
    $seq =~ s/\s//g;
    $seq =~ s/\n//g;
    my $seqlen = length $seq;
    print "Loading genomic sequence for chr $chr ($seqlen nt)\n" if $DEBUG;
    my $i = 0;
#    my $gso = $coge->get_genomic_sequence_obj;
    while ($i < $seqlen)
      {
	my $str = substr($seq,$i,$len);
	my $start = $i+1;
	my $stop = $i + length $str;
	$i += $len;
	unless ($str)
	  {
	    print join ("\t", $seqlen, $start, $stop, $i),"\n";
	    exit;
	  }
	$di->add_to_genomic_sequences({start=>$start,
				       stop=>$stop,
				       chromosome=>$chr,
				       sequence_data=>$str,
				       genomic_sequence_type_id=>$gstype->id,
				      }) if $GO;
	#	print "\t$start-$stop\n";
      }
    if ($use_contigs_as_features)
      {
	my $feat = get_feature(type=>"contig", name=>$chr, di=>$di, chr=>$chr, start=>1, stop=>length($seq));
	add_location(chr=>$chr, start=>1, stop=>length($seq), feat=>$feat, strand=>1);

      }
  }

sub generate_di
  {
    my %opts = @_;
    my $org_id = $opts{org_id};
    my $org_name = $opts{org_name};
    my $org_desc = $opts{org_desc};
    my $ds_id = $opts{ds_id};
    my $ds_name = $opts{ds_name};
    my $ds_desc = $opts{ds_desc};
    my $ds_link = $opts{ds_link};
    my $di_name = $opts{di_name};
    my $di_desc = $opts{di_desc};
    my $di_link = $opts{di_link};
    my $di_version = $opts{di_version};
    unless ($org_id || $org_name)
      {
	warn "no organism id or name specified\n";
	return;
      }
    my $org = $coge->resultset('Organism')->find($org_id) if $org_id;
    $org = $coge->resultset('Organism')->find_or_create({
							 name=>$org_name,
							 description=>$org_desc,
							})  if $GO && !$org;
    unless ($ds_id || $ds_name)
      {
	warn "no data source id or name specified\n";
	return;
      }

    my $ds = $coge->resultset('DataSource')->find($ds_id) if $ds_id;
    $ds = $coge->resultset('DataSource')->find_or_create({
							  name=>$ds_name,
							  description=>$ds_desc,
							  link=>$ds_link,
							 }) if $GO && !$ds;
    unless ($di_name)
      {
	warn "no dataset name specified\n";
	return;
      }
    my $di = $coge->resultset('Dataset')->find_or_create({
									   name                => $di_name,
									   description         => $di_desc,
									   link                => $di_link,
									   organism_id         => $org->id,
									   data_source_id      => $ds->id(),
									   version=>$di_version,
									  })  if $GO;
    return $di;

  }

sub generate_gstype
  {
    my %opts = @_;
    my $name = $opts{name};
    my $desc = $opts{desc};
    my $id = $opts{id};
    my $gstype;
    $id = 1 unless ($id || $name);
    if ($id)
      {
	$gstype = $coge->resultset('GenomicSequenceType')->find($id);
      }
    else
      {
	$gstype = $coge->resultset('GenomicSequenceType')->find_or_create({name=>$name,
									     description=>$desc,
									    });
      }
    return $gstype;
  }

__END__
