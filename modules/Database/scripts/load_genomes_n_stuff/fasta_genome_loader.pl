#!/usr/bin/perl

use DBI;
use strict;
use CoGeX;
use Roman;
use Data::Dumper;
use Getopt::Long;
use File::Path;

use vars qw($DEBUG $coge $GENOMIC_SEQ_LEN $GO $ERASE);

my (
	$nt_file,       $nt_dir,           $org_name,    $org_desc,
	$org_id,        $org_restricted,   $source_name, $source_desc,
	$source_link,   $source_id,        $ds_name,     $ds_desc,
	$ds_link,       $ds_version,       $ds_id,       $chr,
	$seq_type_name, $seq_type_desc,    $seq_type_id, $chr_basename,
	$add_chr_name,  $use_fasta_header, $dsg_name,    $dsg_desc,
	$dsg_version,   $dsg_id,           $restricted,  $db,
	$user,          $pass,             $seq_dir
);

##Example usage:
#./fasta_genome_loader.pl -org_name "Allenigales" -source_id 24 -ds_name oldsuper2.fasta -ds_version 2 -use_fasta_header -nt ~/projects/genome/data/Selaginella_moellendorffii/pre-v2/oldsuper2.fasta

GetOptions(
	"debug=s"                           => \$DEBUG,
	"go=s"                              => \$GO,
	"erase|e"                           => \$ERASE,
	"fasta_file|fasta|faa|nt_file|nt=s" => \$nt_file,
	"fasta_dir|nt_dir|dir=s"            => \$nt_dir,
	"org_name=s"                        => \$org_name,
	"org_desc=s"                        => \$org_desc,
	"org_id=s"                          => \$org_id,
	"org_restricted" => \$org_restricted,   #set flag to make organism restriced
	"source_name=s"  => \$source_name,      # datasource
	"source_desc=s"  => \$source_desc,
	"source_link=s"  => \$source_link,
	"source_id=s"    => \$source_id,
	"ds_name=s"      => \$ds_name,          # datasetid
	"ds_desc=s"      => \$ds_desc,
	"ds_link=s"      => \$ds_link,
	"ds_version=s"   => \$ds_version,
	"ds_id=s"        => \$ds_id,
	"dsg_name=s"     => \$dsg_name, # FIXME mdb 9/7/12, need to add alias for "genome" rename
	"dsg_desc=s"     => \$dsg_desc,
	"dsg_version=s"  => \$dsg_version,
	"dsg_id=s"       => \$dsg_id,
	"restricted=i"   => \$restricted,

	"chr=s"            => \$chr,
	"seq_type_name=s"  => \$seq_type_name,
	"seq_type_desc=s"  => \$seq_type_desc,
	"seq_type_id=i"    => \$seq_type_id,        # masked50 == id 2
	"chr_basename=s"   => \$chr_basename,
	"add_chr_name=s"   => \$add_chr_name,
	"use_fasta_header" => \$use_fasta_header,
	"database|db=s"    => \$db,
	"user|u=s"         => \$user,
	"password|pw=s"    => \$pass,
	"seq_dir|sd=s"     => \$seq_dir,
);

$DEBUG = 1 unless defined $DEBUG; # set to '1' to get updates on what's going on
$GO    = 1 unless defined $GO;    #set to 1 to actually make db calls.
($ds_name) = $nt_file =~ /([^\/]*)$/ unless ($ds_name);

$restricted = 0 unless defined $restricted;
print STDERR "Running $0\n";

my $formatdb = "/usr/bin/formatdb -p F -o T";  #path to blast's formatdb program

my $connstr = "dbi:mysql:dbname=$db;host=localhost;port=3307";
$coge = CoGeX->connect( $connstr, $user, $pass );

#$coge->storage->debugobj(new DBIxProfiler());
#$coge->storage->debug(1);

if ($org_name) {
	my $org =
	  $coge->resultset("Organism")
	  ->find_or_create( { name => $org_name, description => $org_desc } );
	if ($org_restricted) {
		$org->restricted(1);
		$org->update;
	}
	$org_id = $org->id;
}

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

unless ( $org_id && $source_id ) {
	print "Need a valid organism id and data source id in order to create dataset and genome objects.  Loading aborted.\n";
	exit;
}

if ($seq_type_name) {
	my $gst = $coge->resultset('GenomicSequenceType')->find_or_create(
		{
			name        => $seq_type_name,
			description => $seq_type_desc,
		}
	);
	$seq_type_id = $gst->id;
}
$seq_type_id = 1 unless $seq_type_id;    #default to unmasked sequence data

my $ds = generate_ds( ds_name => $ds_name, ds_desc => $ds_desc, ds_link => $ds_link, ds_version => $ds_version, ds_id => $ds_id, source_id => $source_id, restricted => $restricted );

unless ($ds) {
	warn "dataset object not initialized.  Exiting.";
	exit;
}

$dsg_version = $ds->version unless $dsg_version;
my $dsg = generate_dsg( name => $dsg_name, desc => $dsg_desc, version => $dsg_version, dsg_id => $dsg_id, org_id => $org_id, gst_id => $seq_type_id, restricted => $restricted );

#link dataset and genome

$coge->resultset('DatasetConnector')
  ->find_or_create( { dataset_id => $ds->id, genome_id => $dsg->id } );

if ($ERASE) {
	print "Clearing database of entries associated with " . $ds->name . ". . .";
	$ds->delete();
	$dsg->delete();
	print "finished!\n";
	exit;
}
print "dataset_id: " . $ds->id,        "\n";
print "genome_id: " . $dsg->id, "\n\n\n";

process_nt( file => $nt_file, ds => $ds, dsg => $dsg, dir => $nt_dir, chr => $chr ) if $nt_file || $nt_dir;
print "dataset_id: " . $ds->id,        "\n";
print "genome_id: " . $dsg->id, "\n";

sub add_location {
	my %opts   = @_;
	my $chr    = $opts{chr};
	my $start  = $opts{start};
	my $stop   = $opts{stop};
	my $strand = $opts{strand};
	my $feat   = $opts{feat};
	print "Adding location $chr:($start-$stop, $strand)\n" if $DEBUG;
	my $loc = $coge->resultset('Location')->find_or_create(
		{
			feature_id => $feat->id,
			start      => $start,
			stop       => $stop,
			strand     => $strand,
			chromosome => $chr
		}
	  )
	  if $GO;
	return $loc;
}

sub get_feature {
	my %opts  = @_;
	my $type  = $opts{type};
	my $name  = $opts{name};
	my $chr   = $opts{chr};
	my $start = $opts{start};
	my $stop  = $opts{stop};
	print "Working on $name of type $type\n" if $DEBUG;
	my $ds        = $opts{ds};
	my $feat_type =
	  $coge->resultset('FeatureType')->find_or_create( { name => $type } )
	  if $GO;
	print "Creating feature of type $type\n" if $DEBUG;
	my $feat = $coge->resultset('Feature')->find_or_create(
		{
			dataset_id      => $ds->id,
			feature_type_id => $feat_type->id,
			start           => $start,
			stop            => $stop,
			chromosome      => $chr,
			strand          => 1,
		}
	  )
	  if ($GO);

	if ($DEBUG) {
		print "Creating feature_name $name for feat ";
		print $feat->id if $feat;
		print "\n";
	}
	my $feat_name = $coge->resultset('FeatureName')->find_or_create(
		{
			name       => $name,
			feature_id => $feat->id,
		}
	  )
	  if $GO;
	return $feat;
}

sub process_nt {
	my %opts = @_;
	my $file = $opts{file};
	print $file, "\n";
	my $ds  = $opts{ds};
	my $dsg = $opts{dsg};
	my $dir = $opts{dir};
	my $chr = $opts{chr};
	my @files;
	push @files, $file if $file && -r $file;

	if ( -d $dir ) {
		opendir( DIR, $dir );
		while ( my $item = readdir(DIR) ) {
			next if $item =~ /^\.\.?$/;
			push @files, "$dir/$item" if -r "$dir/$item";
		}
		closedir(DIR);
	}
	foreach my $file (@files) {
		process_nt_file( file => $file, ds => $ds, chr => $chr, dsg => $dsg );
	}
	my $cmd = $formatdb . " -i " . $dsg->file_path;    #."/".$dsg->id.".faa";
	print "\tFormatdb running $cmd\n";
	`$cmd`;
}

sub process_nt_file {
	my %opts = @_;
	my $file = $opts{file};
	my $ds   = $opts{ds};
	my $dsg  = $opts{dsg};
	my $chr  = $opts{chr};
	return unless $file;
	print "Processing $file\n";
	$/ = "\n>";
	open( IN, $file ) || die "can't open $file for reading: $!";

	while (<IN>) {
		s/>//g;
		my ( $name, $seq ) = split /\n/, $_, 2;
		$seq =~ s/\n//g;

		my $chrtmp;
		$chrtmp = $chr  if defined $chr;
		$chrtmp = $name if $use_fasta_header;
		($chrtmp) = $name =~ /^(\S+)/ unless $chrtmp || $add_chr_name;
		$chrtmp = $name unless defined $chrtmp;
		$chrtmp =~ s/^lcl\|//;
		$chrtmp =~ s/chromosome//i;
		$chrtmp =~ s/^chr//i;
		$chrtmp =~ s/^0+//;
		$chrtmp =~ s/^_+//;
		$chrtmp =~ s/\s+/ /;
		$chrtmp =~ s/^\s//;
		$chrtmp =~ s/\s$//;
		$chrtmp = 0 unless $chrtmp;
		$chrtmp = $chr_basename . $chrtmp if $chr_basename;

		#($chrtmp) = $name =~ /(^\S+)/;
		load_genomic_sequence(
			chr => $chrtmp,
			seq => $seq,
			ds  => $ds,
			dsg => $dsg,
		);
	}
	close IN;
}

sub load_genomic_sequence {
	my %opts = @_;
	my $seq  = $opts{seq};
	my $chr  = $opts{chr};
	my $ds   = $opts{ds};
	my $dsg  = $opts{dsg};
	$seq =~ s/\s//g;
	$seq =~ s/\n//g;
	my $seqlen = length $seq;

	unless ($seqlen) {
		print "Chromosome $chr has no sequencing.  Skipping. . .\n";
		return;
	}
	print "Loading genomic sequence for chromosome: $chr ($seqlen nt)\n"
	  if $DEBUG;

	$dsg->add_to_genomic_sequences(
		{
			sequence_length => $seqlen,
			chromosome      => $chr,
		}
	  )
	  if $GO;
	my $path = $dsg->file_path;
	$path =~ s/\/[^\/]*$/\//;
	mkpath($path);
	mkpath( $path . "/chr" );

	#append sequence ot master file for dataset group
	open( OUT, ">>" . $path . "/" . $dsg->id . ".faa" );
	my $head = $chr =~ /^\d+$/ ? ">gi" : ">lcl";
	$head .= "|" . $chr;
	print OUT "$head\n$seq\n";
	close OUT;

	#create individual file for chromosome
	open( OUT, ">" . $path . "/chr/$chr" );
	print OUT $seq;
	close OUT;

#must add a feature of type chromosome to the dataset so the dataset "knows" its chromosomes
	my $feat = get_feature(
		type  => "chromosome",
		name  => "chromosome $chr",
		ds    => $ds,
		chr   => $chr,
		start => 1,
		stop  => length($seq)
	);
	add_location(
		chr    => $chr,
		start  => 1,
		stop   => length($seq),
		feat   => $feat,
		strand => 1
	);
}

sub generate_ds {
	my %opts       = @_;
	my $ds_name    = $opts{ds_name};
	my $ds_desc    = $opts{ds_desc};
	my $ds_link    = $opts{ds_link};
	my $ds_version = $opts{ds_version};
	my $ds_id      = $opts{ds_id};
	my $source_id  = $opts{source_id};
	my $restriced  = $opts{restricted};
	unless ( $ds_name || $ds_id ) {
		warn "no dataset name or database id specified\n";
		return;
	}
	my $ds = $ds_id
	  ? $coge->resultset('Dataset')->find($ds_id)
	  : $coge->resultset('Dataset')->create(
		{
			name           => $ds_name,
			description    => $ds_desc,
			link           => $ds_link,
			data_source_id => $source_id,
			restricted     => $restricted,
			version        => $ds_version,
		}
	  )
	  if $GO;
	return $ds;

}

sub generate_dsg {
	my %opts       = @_;
	my $name       = $opts{name};
	my $desc       = $opts{desc};
	my $version    = $opts{version};
	my $org_id     = $opts{org_id};
	my $gst_id     = $opts{gst_id};
	my $dsg_id     = $opts{dsg_id};
	my $restricted = $opts{restricted} || 0;
	my $dsg        = $dsg_id
	  ? $coge->resultset('Genome')->find($dsg_id)
	  : $coge->resultset('Genome')->create(
		{
			name                     => $name,
			description              => $desc,
			version                  => $version,
			organism_id              => $org_id,
			genomic_sequence_type_id => $gst_id,
			restricted               => $restricted,
		}
	  )
	  if $GO;
	return unless $dsg;

	unless ( $dsg->file_path ) {
		my $path = "$seq_dir/" . $dsg->get_path . "/" . $dsg->id . ".faa";
		print $path, "\n";
		$dsg->file_path($path);
		$dsg->update;
	}
	return $dsg;
}

__END__
