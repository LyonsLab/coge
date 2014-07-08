#!/usr/bin/perl -w

use strict;
use CoGeX;
use CoGe::Accessory::Web;
use Data::Dumper;
use Getopt::Long;
use File::Path;

use vars
  qw($DEBUG $GO $conf_file $coge @gids $restricted $P $org_name $version);

GetOptions(
    "debug=s"        => \$DEBUG,
    "go=s"           => \$GO,
    "conf_file|cf=s" => \$conf_file,
    "gid=i"          => \@gids,
    "restricted|r=i" => \$restricted,
    "org_name|on=s"  => \$org_name,
    "version|v=s"    => \$version,
);
$org_name = "Merged" unless $org_name;
$version  = 1        unless defined $version;
$P = CoGe::Accessory::Web::get_defaults($conf_file);

unless ( $P && $P->{DBNAME} ) { usage(); }

my $TEMPDIR = $P->{TEMPDIR} . "copy_genome";
mkpath( $TEMPDIR, 1, 0777 );

my $DBNAME = $P->{DBNAME};
my $DBHOST = $P->{DBHOST};
my $DBPORT = $P->{DBPORT};
my $DBUSER = $P->{DBUSER};
my $DBPASS = $P->{DBPASS};

my $connstr =
  "dbi:mysql:dbname=" . $DBNAME . ";host=" . $DBHOST . ";port=" . $DBPORT;
$coge = CoGeX->dbconnect(
    db_connection_string => $connstr,
    db_name              => $DBUSER,
    db_passwd            => $DBPASS
);

#$coge->storage->debugobj(new DBIxProfiler());
#$coge->storage->debug(1);

my $fasta_genome_loader   = $P->{SCRIPTDIR} . '/load_genome.pl';
my $replicate_annotations = $P->{SCRIPTDIR} . '/copy_genome/replicate_annotations.pl';

unless (@gids) {
    usage();
}

unless ( -r $fasta_genome_loader ) {
    print STDERR qq{
Error:  unable to read fasta_genome_loader: $fasta_genome_loader
};
    exit;
}

my @genomes;
foreach my $gid (@gids) {
    my ($dsg) = $coge->resultset('Genome')->find($gid);
    $restricted = $dsg->restricted
      unless defined $restricted
    ; #use whatever is specified by the original genome unless user specifies otherwise
    unless ($dsg) {
        print STDERR qq{
Error: Unable to create a datasetgroup object for $gid
};
        next;
    }
    push @genomes, $dsg;
}

unless (@genomes) {
    print STDERR qq{No valid genome objects were created.
};
    exit;
}

my $chr_headers = get_chr_headers_for_genomes( \@genomes );
print Dumper $chr_headers;
my $new_faa =
  get_and_clean_sequence( genomes => \@genomes, chr_headers => $chr_headers );

my $stid = get_sequence_type_id_for_merged();

print STDERR "Loading genome\n";
my $new_gid = load_genome(
    faa      => $new_faa,
    genomes  => \@genomes,
    stid     => $stid,
    org_name => $org_name,
    version  => $version
);
print STDERR "NEW DSGID: $new_gid\n" if $DEBUG;

unless ($new_gid) {
    print STDERR
"Error: unable to capture a new genome id for the masked genome.  Exiting.\n";
    exit if $GO;
}
print STDERR "adding annotations\n";
add_annotations(
    new_gid     => $new_gid,
    genomes     => \@genomes,
    chr_headers => $chr_headers
);
print STDERR "finished\n";

sub add_annotations {
    my %opts        = @_;
    my $new_gid     = $opts{new_gid};
    my $genomes     = $opts{genomes};
    my $chr_headers = $chr_headers;
    my $new_genome  = $coge->resultset('Genome')->find($new_gid);
    my $new_ds      = map_chr_to_ds($new_genome);

    foreach my $genome (@genomes) {
        my $ds = map_chr_to_ds($genome);
        foreach my $chr ( sort keys %$ds ) {
            my $new_chr = $chr_headers->{ $genome->id } . "_$chr";
            unless ( $new_ds->{$new_chr} ) {
                print STDERR
"Chromosome $new_chr does not exist in new genome.  Please check the chromosome names if you think this is in error.  Skipping to next chromosome\n";
                next;
            }
            foreach
              my $f1 ( $ds->{$chr}->features->search( { chromosome => $chr } ) )
            {
                my ($f2) = $new_ds->{$new_chr}->features->search(
                    {
                        chromosome      => $new_chr,
                        feature_type_id => $f1->feature_type_id,
                        start           => $f1->start,
                        stop            => $f1->stop,
                        strand          => $f1->strand,
                    }
                );
                if ($f2) {
                    print "Feature already exists.\n" if $DEBUG;
                }
                else {
                    ($f2) = replicate_feature(
                        f1      => $f1,
                        ds2     => $new_ds->{$new_chr},
                        new_chr => $new_chr
                    );
                }
                next unless $f2;    #can't do much without the annotation
                replicate_names( f1 => $f1, f2 => $f2 );
                replicate_annotations( f1 => $f1, f2 => $f2 );
            }
        }
    }
}

sub replicate_annotations {
    my %opts = @_;
    my $f1   = $opts{f1};
    my $f2   = $opts{f2};
    print "Processing annotations\n" if $DEBUG;
    foreach my $anno ( $f1->annotations ) {
        my (@annos) = $f2->annotations->search(
            {
                annotation         => $anno->annotation,
                annotation_type_id => $anno->annotation_type_id
            }
        );
        if ( scalar @annos ) {
            print "\tAnnotation exists\n" if $DEBUG;
            next;
        }
        print "\tAdding annotation\n" if $DEBUG;
        $f2->add_to_feature_annotations(
            {
                annotation         => $anno->annotation,
                annotation_type_id => $anno->annotation_type_id
            }
        ) if $GO;
    }
}

sub replicate_names {
    my %opts = @_;
    my $f1   = $opts{f1};
    my $f2   = $opts{f2};
    print "Processing names\n" if $DEBUG;
    foreach my $name ( $f1->feature_names ) {
        my @names = $f2->feature_names->search(
            { name => $name->name, description => $name->description } );
        if ( scalar @names ) {
            print "\tName exists\n" if $DEBUG;
            next;
        }
        print "\tAdding name\n" if $DEBUG;
        $f2->add_to_feature_names(
            {
                name         => $name->name,
                description  => $name->description,
                primary_name => $name->primary_name
            }
        ) if $GO;
    }
}

sub replicate_feature {
    my %opts    = @_;
    my $f1      = $opts{f1};
    my $ds2     = $opts{ds2};
    my $new_chr = $opts{new_chr};
    print "Creating feature\n" if $DEBUG;
    my ($f2) = $ds2->add_to_features(
        {
            chromosome      => $new_chr,
            feature_type_id => $f1->feature_type_id,
            start           => $f1->start,
            stop            => $f1->stop,
            strand          => $f1->strand,
        }
    ) if $GO;
    foreach my $loc ( $f1->locations ) {
        print "\tadding location\n" if $DEBUG;
        $f2->add_to_locations(
            {
                start      => $loc->start,
                stop       => $loc->stop,
                strand     => $loc->strand,
                chromosome => $loc->chromosome,
            }
        ) if $GO;
    }
    return $f2;
}

sub map_chr_to_ds {
    my $dsg = shift;
    my %ds;
    foreach my $ds ( $dsg->datasets ) {
        foreach my $chr ( $ds->chromosomes ) {
            $ds{$chr} = $ds;
        }
    }
    return \%ds;
}

sub load_genome {
    my %opts     = @_;
    my $faa      = $opts{faa};
    my $genomes  = $opts{genomes};
    my $stid     = $opts{stid};
    my $org_name = $opts{org_name};
    my $version  = $opts{version};
    my $cmd      = $fasta_genome_loader;
    $cmd .= " -org_name '$org_name'";

    #    $cmd .= " -genome_name '". $dsg->name."'" if $dsg->name;
    my $genome_name =
      join( "_" . map { $_->id } sort { $a->id <=> $b->id } @$genomes );
    my $message = "Merged genome of: ";
    $message .= join(
        ", ",
        map {
                "<a href=OrganismView.pl?gid="
              . $_->id
              . " target=_new>"
              . $_->organism->name . "</a>"
          } sort { $a->id <=> $b->id } @$genomes
    );
    $cmd .= " -genome_message '$message'";
    $cmd .= " -source_name 'CoGe Merged Genome'";
    $cmd .= " -ds_version '$version'";
    $cmd .= " -ds_name 'CoGe Merged Genome'";
    $cmd .= " -seq_type_id " . $stid;
    $cmd .= " -u " . $DBUSER;
    $cmd .= " -pw " . $DBPASS;
    $cmd .= " -db " . $DBNAME;
    $cmd .= " -restricted " . $restricted if $restricted;
    my ($dsg) = $genomes->[0];
    my $seq_dir = $dsg->file_path;
    $seq_dir =~ s/[^\/]*$//;
    $seq_dir =~ s/\d+\///g;
    $cmd .= " -sd " . $seq_dir;
    $cmd .= " -nt " . $faa;
    print STDERR "Running: ", $cmd, "\n";
    my $dsgid;

    if ($GO) {
        open( CMD, "$cmd |" );
        while (<CMD>) {
            if (/genome_id:\s+(\d+)/) {
                $dsgid = $1;
                print STDERR "Captured dsgid: $dsgid\n" if $DEBUG;
            }
            print STDERR $_ if $DEBUG;
        }
        close CMD;
    }
    return $dsgid;
}

sub get_and_clean_sequence {
    my %opts        = @_;
    my $genomes     = $opts{genomes};
    my $chr_headers = $opts{chr_headers};
    my $rnd         = int( rand(1000000000000) );
    my $new_faa     = $TEMPDIR . "/" . $rnd . ".merged.faa";
    open( OUT, ">$new_faa" );
    foreach my $genome (@$genomes) {
        my $header = $chr_headers->{ $genome->id };
        unless ($header) {
            print STDERR "Warning: No header for "
              . $genome->id
              . ".  Skipping\n";
            next;
        }
        open( IN, $genome->file_path );
        while (<IN>) {
            if (/^>/) {
                s/>//g;
                s/lcl\|//g;
                s/gi\|//g;
                print OUT ">" . $header . "_" . $_;
            }
            else {
                print OUT $_;
            }
        }
        close IN;
    }
    close OUT;
    return $new_faa;
}

sub get_sequence_type_id_for_merged {
    my $st =
      $coge->resultset('GenomicSequenceType')
      ->find_or_create( { name => "Megered genomes" } );
    return $st->id;
}

sub get_chr_headers_for_genomes {
    my $genomes = shift;
    my %data1;
    my %out;
    foreach my $genome (@$genomes) {
        my $name = $genome->organism->name;
        my $new_name;
        my $count = 0;
        foreach my $item ( split /\s+/, $name ) {
            last if $count >= 2;
            my ($fl) = $item =~ /^(\S)/;
            $new_name .= uc($fl);
            $count++;
        }
        if ( $data1{$new_name} ) {
            my $id = $data1{$new_name};
            $new_name = $genome->id;
            $out{$id} = $id;
        }
        $data1{$new_name} = $genome->id;
        $out{ $genome->id } = $new_name;
    }
    return \%out;
}

sub usage {
    print STDERR qq{
Welcome to $0

Purpose:  take a list genome IDs from coge and generates a merged genome from them

Usage:  $0  -conf_file <coge configuration file> -gid <coge database id> -gid <coge database id>  -go 1

Options:
   -go 1                  |      Make the database calls.  Default 0

   -conf_file | cf        |      CoGe conf file

   -gid                   |      CoGe genome id

   -org_name              |      Name for the hybrid genome.  Default: Merged

   -version               |      Version of genome  (Default: 1)

   -restricted            |      mark genome as restricted (Default: 0)

Copyright Eric Lyons 2013

};
    exit;
}
