#! /usr/bin/perl -w
use v5.10;
use strict;
no warnings 'redefine';
umask(0);

use Benchmark;
use Carp;
use DBI;
use Getopt::Long;
use Parallel::ForkManager;

use CoGeX;
use CoGe::Accessory::Web;
use DBIxProfiler;

$| = 1;
our (
    $cogeweb, $basename, $gid,     $feature, $fasta,  $coge,
    $P,       $TEMPDIR,  $NWALIGN, $DBNAME,  $DBHOST, $DBPORT,
    $DBUSER,  $DBPASS,   $CONFIG,  $debug
);

GetOptions(
    "genome_id|gid=s"   => \$gid,
    "feature_type|ft=s" => \$feature,
    "fasta|f=s"         => \$fasta,
    "config|cfg=s"      => \$CONFIG,
    "debug"             => \$debug
);

$P = CoGe::Accessory::Web::get_defaults($CONFIG);
$ENV{PATH} = join(":",
  ( $P->{COGEDIR}, $P->{BINDIR}, $P->{BINDIR} . "SynMap",
    "/usr/bin", "/usr/local/bin" )
);
$TEMPDIR = $P->{TEMPDIR} . "SynMap";
$NWALIGN = get_command_path('NWALIGN');

$DBNAME = $P->{DBNAME};
$DBHOST = $P->{DBHOST};
$DBPORT = $P->{DBPORT};
$DBUSER = $P->{DBUSER};
$DBPASS = $P->{DBPASS};

my $connstr = "dbi:mysql:dbname=" . $DBNAME . ";host=" . $DBHOST . ";port=" . $DBPORT;
$coge = CoGeX->connect( $connstr, $DBUSER, $DBPASS );

#$coge->storage->debugobj(new DBIxProfiler());
#$coge->storage->debug(1);

$cogeweb = CoGe::Accessory::Web::initialize_basefile(
    basename => $basename,
    tempdir  => $TEMPDIR
);
$| = 1;    # Disable buffering

gen_fasta( gid => $gid, feature_type => $feature, fasta => $fasta );

sub gen_fasta {
    my %opts         = @_;
    my $gid          = $opts{gid};
    my $feature_type = $opts{feature_type};
    my $file         = $opts{fasta};
    my ($genome)     =
#      $coge->resultset('Genome')->search( { "me.genome_id" => $gid },
#        { join => 'genomic_sequences', prefetch => 'genomic_sequences' } );
      $coge->resultset('Genome')->search( { "me.genome_id" => $gid } );

    my $output;

    my %datasets;      #storage for dataset objects based on dataset id
    my %feat_types;    #storage for feature type objects

    if ( $feature_type eq "CDS" || $feature_type eq "protein" ) {
        my $count = 1;
        my $t1 = new Benchmark if $debug;
        my %chr_sequence; # = map{$_=>$genome->get_genomic_sequence(chr=>$_)} $genome->get_chromosomes;
        my $t2 = new Benchmark if $debug;
        my $d0 = timestr( timediff( $t2, $t1 ) ) if $debug;
        say STDERR "Populate genomic sequence time: $d0" if $debug;
        my @res = $coge->resultset('Feature')->search(
            {
                feature_type_id => [ 3, 5, 8 ],
                genome_id       => $gid
            },
            {
                join => [ { dataset => 'dataset_connectors' } ],
                prefetch => [ 'feature_names', 'locations' ]
            }
        );
        my $t3 = new Benchmark if $debug;
        my $d1 = timestr( timediff( $t3, $t2 ) ) if $debug;
        say STDERR "Query time: $d1"             if $debug;

        my @feats = sort { $a->chromosome cmp $b->chromosome || $a->start <=> $b->start } @res;

        CoGe::Accessory::Web::write_log(
            "Getting sequence for "
              . scalar(@feats)
              . " features of types CDS, tRNA, and rRNA.",
            $cogeweb->logfile
        );
        foreach my $feat (@feats) {
            my $t3 = new Benchmark if $debug;

            # mdb added 3/25/14, issue 338 - handle start/stop reversed
            my $start = $feat->start;
            my $stop = $feat->stop;
            ($start, $stop) = ($stop, $start) if ($start > $stop);

            my ($chr) = $feat->chromosome;#=~/(\d+)/;
            my $name;
            foreach my $n ( $feat->names ) {
                $name = $n;
                last unless $name =~ /\s/;
            }
            unless ($name) {
                #print STDERR "Error:  missing valid name for feature_id ".$feat->id."\n";
                $name = $feat->id;
            }

            $name =~ s/\s+/_/g;
            my $feat_type_name;
            if ( $feat_types{ $feat->feature_type_id } ) {
                $feat_type_name = $feat_types{ $feat->feature_type_id }->name;
            }
            else {
                my $feat_type = $feat->type;
                $feat_types{ $feat_type->id } = $feat_type;
                $feat_type_name = $feat_type->name;
            }
            my $title = join( "||",
                $chr, $start, $stop, $name, $feat->strand,
                $feat_type_name, $feat->id, $count );

            #print STDERR "getting sequence\n";
            if ( $feature_type eq "CDS" ) {
                my $dataset;
                if ( $datasets{ $feat->dataset_id } ) {
                    $dataset = $datasets{ $feat->dataset_id };
                }
                else {
                    $dataset = $feat->dataset;
                    $datasets{ $dataset->id } = $dataset;
                }
                $chr_sequence{ $feat->chromosome } = $genome->get_genomic_sequence(chr => $feat->chromosome) unless $chr_sequence{ $feat->chromosome };
                next unless $chr_sequence{ $feat->chromosome };

                my $start_index;
                # Check that the start location is greater than 1
                if ($start < 1) {
                    carp "feature start location is less than 1, start=$start";
                    $start_index = 0;
                } else {
                    $start_index = $start - 1;
                }

                my $subseq = substr(
                    $chr_sequence{ $feat->chromosome },
                    $start_index,
                    $stop - $start + 1
                );
                #unless ($subseq) {
                #print STDERR join ("\t", $feat->chromosome, $feat->start, $feat->stop),"\n";
                #}
                #print STDERR $subseq,"\n",$genome->get_genomic_sequence(start=>$feat->start, stop=>$feat->stop, chr=>$feat->chromosome),"\n","#"x20,"\n";

                my $seq = $feat->genomic_sequence(
                    dsgid   => $genome,
                    dataset => $dataset,
                    seq     => $subseq
                );
                my $t4 = new Benchmark if $debug;
                my $d2 = timestr( timediff( $t4, $t3 ) ) if $debug;
                say STDERR "Seq $count time: $d2" if $debug;

                next unless $seq;

                #skip sequences that are only 'x' | 'n';
                next unless $seq =~ /[^x|n]/i;

                #print OUT ">" . $title . "\n";
                if ($seq) {
                    $output .= ">" . $title . "\n";

                    #print OUT $seq, "\n";
                    $output .= $seq . "\n";
                    $count++;
                }
            }
            elsif ( $feature_type eq "protein" ) {
                next unless $feat->feature_type_id == 3;
                my (@seqs) = $feat->protein_sequence( dsgid => $genome->id );
                next unless scalar @seqs;
                next
                  if scalar @seqs > 1;   #didn't find the correct reading frame;
                next unless $seqs[0] =~ /[^x]/i;
                $title = ">" . $title . "\n";

                #print OUT $title, $seqs[0], "\n";
                $output .= $title . $seqs[0] . "\n";
                $count++;
            }
        }
    }
    else {
        my @chr = sort $genome->get_chromosomes;
        CoGe::Accessory::Web::write_log(
            "Getting sequence for "
              . scalar(@chr)
              . " chromosomes (genome sequence)",
            $cogeweb->logfile
        );

        $file = $genome->file_path;

        #foreach my $chr (@chr)
        #{
        #    my $seq = $genome->get_genomic_sequence(chr=>$chr);
        #    next unless $seq;
        #    print OUT ">".$chr."\n";
        #    print OUT $seq,"\n";
        #}
    }

    if ($output) {
        open( OUT, ">$file" ) || die "Can't open $file for writing: $!";
        print OUT $output;
        close OUT;
    }
    return 1 if -r $file;

    CoGe::Accessory::Web::write_log( "Error with fasta file creation", $cogeweb->logfile );
    return 0;
}
