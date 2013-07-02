#! /usr/bin/perl -w
use v5.10;
use strict;
no warnings 'redefine';
umask(0);

use Benchmark;
use DBI;
use Getopt::Long;
use Parallel::ForkManager;

use CoGeX;
use CoGe::Accessory::LogUser;
use CoGe::Accessory::Web;

our (
    $cogeweb, $basename, $gid,     $feature, $fasta,
    $coge,    $P,        $TEMPDIR, $NWALIGN, $DBNAME,
    $DBHOST,  $DBPORT,   $DBUSER,  $DBPASS,  $CONFIG);

GetOptions(
    "genome_id|gid=s"   => \$gid,
    "feature_type|ft=s" => \$feature,
    "fasta|f=s"         => \$fasta,
    "config|cfg=s"      => \$CONFIG,);

$P = CoGe::Accessory::Web::get_defaults($CONFIG);
$ENV{PATH} = join ":",
  (
    $P->{COGEDIR}, $P->{BINDIR}, $P->{BINDIR} . "SynMap",
    "/usr/bin", "/usr/local/bin");
$TEMPDIR = $P->{TEMPDIR} . "SynMap";
$NWALIGN = $P->{NWALIGN};

$DBNAME = $P->{DBNAME};
$DBHOST = $P->{DBHOST};
$DBPORT = $P->{DBPORT};
$DBUSER = $P->{DBUSER};
$DBPASS = $P->{DBPASS};

my $connstr =
  "dbi:mysql:dbname=" . $DBNAME . ";host=" . $DBHOST . ";port=" . $DBPORT;
$coge = CoGeX->connect($connstr, $DBUSER, $DBPASS);

$cogeweb = CoGe::Accessory::Web::initialize_basefile(
    basename => $basename,
    tempdir  => $TEMPDIR);
$| = 1;    # Disable buffering

gen_fasta(gid => $gid, feature_type => $feature, fasta => $fasta);

sub gen_fasta
{
    my %opts         = @_;
    my $gid          = $opts{gid};
    my $feature_type = $opts{feature_type};
    my $file         = $opts{fasta};
    my ($genome) = $coge->resultset('Genome')->search({"me.genome_id" => $gid},
        {join => 'genomic_sequences', prefetch => 'genomic_sequences'});

    open(OUT, ">$file") || die "Can't open $file for writing: $!";
    if ($feature_type eq "CDS" || $feature_type eq "protein") {
        my $count = 1;
        my @res   = $coge->resultset('Feature')->search({
                feature_type_id => [3, 5, 8],
                genome_id       => $gid},
            {
                join     => [{dataset => 'dataset_connectors'}],
                prefetch => ['feature_names']});

        my @feats =
          sort {$a->chromosome cmp $b->chromosome || $a->start <=> $b->start}
          @res;

        CoGe::Accessory::Web::write_log(
            "Getting sequence for "
              . scalar(@feats)
              . " features of types CDS, tRNA, and rRNA.",
            $cogeweb->logfile);
        foreach my $feat (@feats) {
            my ($chr) = $feat->chromosome;    #=~/(\d+)/;
            my $name;
            foreach my $n ($feat->names) {
                $name = $n;
                last unless $name =~ /\s/;
            }
            unless ($name) {
                #    print STDERR "Error:  missing valid name for feature_id ".$feat->id."\n";
                $name = $feat->id;
            }

            $name =~ s/\s+/_/g;
            my $title = join("||",
                $chr, $feat->start, $feat->stop, $name, $feat->strand,
                $feat->type->name, $feat->id, $count);
            if ($feature_type eq "CDS") {
                my $seq = $feat->genomic_sequence(dsgid => $genome->id);
                next unless $seq;

                #skip sequences that are only 'x' | 'n';
                next unless $seq =~ /[^x|n]/i;
                print OUT ">" . $title . "\n";
                print OUT $seq, "\n";
                $count++;
            } elsif ($feature_type eq "protein") {
                next unless $feat->feature_type_id == 3;
                my (@seqs) = $feat->protein_sequence(dsgid => $genome->id);
                next unless scalar @seqs;
                next
                  if scalar @seqs > 1;   #didn't find the correct reading frame;
                next unless $seqs[0] =~ /[^x]/i;
                $title = ">" . $title . "\n";

                #       print OUT $title, join ($title, @seqs),"\n";
                print OUT $title, $seqs[0], "\n";
                $count++;
            }
        }
    } else {
        my @chr = sort $genome->get_chromosomes;
        CoGe::Accessory::Web::write_log(
            "Getting sequence for "
              . scalar(@chr)
              . " chromosomes (genome sequence)",
            $cogeweb->logfile);
        $file = $genome->file_path;

        #   foreach my $chr (@chr)
        #     {
        #       my $seq = $genome->get_genomic_sequence(chr=>$chr);
        #       next unless $seq;
        #       print OUT ">".$chr."\n";
        #       print OUT $seq,"\n";
        #     }
    }
    close OUT;
    return 1 if -r $file;
    CoGe::Accessory::Web::write_log("Error with fasta file creation",
        $cogeweb->logfile);

    return 0;
}
