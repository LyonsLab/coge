#! /usr/bin/perl -w
use v5.10;
use strict;
no warnings 'redefine';
umask(0);

use Benchmark;
use DBI;
use Getopt::Long;
use Parallel::ForkManager;
use Data::Dumper;
use File::Spec;
use CoGeX;
use CoGe::Accessory::LogUser;
use CoGe::Accessory::Web;
use CoGe::Algos::KsCalc;


our (
    $cogeweb, $basename, $infile,  $dbfile,   $blockfile, $coge,
    $P,       $TEMPDIR,  $NWALIGN, $MAX_PROC, $DBNAME,    $DBHOST,
    $DBPORT,  $DBUSER,   $DBPASS,  $CONFIG);

GetOptions(
    "basename=s"   => \$basename,
    "infile=s"     => \$infile,
    "dbfile=s"     => \$dbfile,
    "blockfile=s"  => \$blockfile,
    "config|cfg=s" => \$CONFIG,);

$P = CoGe::Accessory::Web::get_defaults($CONFIG);
$ENV{PATH} = join ":",
  (
    $P->{COGEDIR}, $P->{BINDIR}, $P->{BINDIR} . "SynMap",
    "/usr/bin", "/usr/local/bin");
$ENV{HOME} = $P->{COGEDIR};
my $config = File::Spec->catdir($ENV{HOME},"coge.conf");
#print STDERR Dumper \%ENV;


$TEMPDIR  = $P->{TEMPDIR} . "SynMap";
$MAX_PROC = $P->{MAX_PROC};
$NWALIGN  = $P->{NWALIGN};

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

run();

sub run
{
    my ($outfile) = $infile =~ /^(.*?CDS-CDS)/;
    ($outfile) = $infile =~ /^(.*?protein-protein)/ unless $outfile;
    exit 1 unless $outfile;

    CoGe::Accessory::Web::write_log("#" x (20), $cogeweb->logfile);
    CoGe::Accessory::Web::write_log(
        "Running CodeML for synonymous/nonsynonmous rate calculations",
        $cogeweb->logfile);

    my $ret = gen_ks_db(infile => $infile, outfile => $dbfile);
    die if $ret;

    CoGe::Accessory::Web::write_log("#" x (20), $cogeweb->logfile);
    CoGe::Accessory::Web::write_log("", $cogeweb->logfile);

    CoGe::Accessory::Web::write_log("#" x (20), $cogeweb->logfile);
    CoGe::Accessory::Web::write_log("Generate Ks Blocks File",
        $cogeweb->logfile);

    CoGe::Accessory::Web::write_log("#" x (20), $cogeweb->logfile);
    CoGe::Accessory::Web::write_log("", $cogeweb->logfile);
    $ret = gen_ks_blocks_file(
        infile  => $infile,
        dbfile  => $dbfile,
        outfile => $blockfile);
    die if $ret;
}

sub gen_ks_db
{
    my %opts    = @_;
    my $infile  = $opts{infile};
    my $outfile = $opts{outfile};
    CoGe::Accessory::Web::write_log("Generating ks data.", $cogeweb->logfile);
    CoGe::Accessory::Web::write_log("initializing ks database $outfile",
        $cogeweb->logfile);
    my $create = qq{
CREATE TABLE ks_data
(
id INTEGER PRIMARY KEY,
fid1 integer,
fid2 integer,
dS varchar,
dN varchar,
dN_dS varchar,
protein_align_1,
protein_align_2,
DNA_align_1,
DNA_align_2
)
};

    my $dbh = DBI->connect("dbi:SQLite:dbname=$outfile", "", "");
    $dbh->do($create) if $create;
    $dbh->do('create INDEX fid1 ON ks_data (fid1)');
    $dbh->do('create INDEX fid2 ON ks_data (fid2)');
    $dbh->disconnect;

    my $ksdata = get_ks_data(db_file => $outfile);
    $/ = "\n";
    open(IN, $infile);
    my @data;
    while (<IN>) {

        next if /^#/;
        chomp;
        my @line  = split /\t/;
        my @item1 = split /\|\|/, $line[1];
        my @item2 = split /\|\|/, $line[5];
        unless ($item1[6] && $item2[6]) {
            warn "Line does not appear to contain coge feature ids:  $_\n";
            next;
        }
        next unless $item1[5] eq "CDS" && $item2[5] eq "CDS";
        next if $ksdata->{$item1[6]}{$item2[6]};
        push @data, [$line[1], $line[5], $item1[6], $item2[6]];
    }
    close IN;
    my $ports =
      initialize_nwalign_servers(start_port => 3000, procs => $MAX_PROC);
    my $pm = new Parallel::ForkManager($MAX_PROC * 4);
    my $i  = 0;

    my $start_time = new Benchmark;
    $dbh = DBI->connect("dbi:SQLite:dbname=$outfile", "", "");

    $dbh->do("PRAGMA synchronous=OFF");
    $dbh->do("PRAGMA cache_size = 25000");
    $dbh->do("PRAGMA count_changes=OFF");
    $dbh->do("PRAGMA journal_mode=MEMORY");
    $dbh->do("PRAGMA temp_store=MEMORY");

    foreach my $item (@data) {
        $i++;
        $i = 0 if $i == $MAX_PROC * 4;
        $pm->start and next;
        my ($fid1) = $item->[2] =~ /(^\d+$)/;
        my ($fid2) = $item->[3] =~ /(^\d+$)/;
        my ($feat1) = $coge->resultset('Feature')->find($fid1);
        my ($feat2) = $coge->resultset('Feature')->find($fid2);
        my $max_res;
        my $ks = new CoGe::Algos::KsCalc(config=>$config);
        $ks->nwalign_server_port($ports->[$i]);
        $ks->feat1($feat1);
        $ks->feat2($feat2);

        #   for (1..5)
        #     {
        my $res = $ks->KsCalc(config=>$config);    #send in port number?
        $max_res = $res unless $max_res;
        $max_res = $res
          if $res->{dS} && $max_res->{dS} && $res->{dS} < $max_res->{dS};

        #     }
        unless ($max_res) {
            $max_res = {};
        }
        my ($dS, $dN, $dNS) = ("NA", "NA", "NA");
        if (keys %$max_res) {
            $dS  = $max_res->{dS}      if defined $max_res->{dS};
            $dN  = $max_res->{dN}      if defined $max_res->{dN};
            $dNS = $max_res->{'dN/dS'} if defined $max_res->{'dN/dS'};
        }

        #alignments
        my $dalign1 = $ks->dalign1;    #DNA sequence 1
        my $dalign2 = $ks->dalign2;    #DNA sequence 2
        my $palign1 = $ks->palign1;    #PROTEIN sequence 1
        my $palign2 = $ks->palign2;    #PROTEIN sequence 2

        my $insert = qq{
INSERT INTO ks_data (fid1, fid2, dS, dN, dN_dS, protein_align_1, protein_align_2, DNA_align_1, DNA_align_2) values ($fid1, $fid2, "$dS", "$dN", "$dNS", "$palign1", "$palign2", "$dalign1", "$dalign2")
};

        my $insert_success = 0;
        while (!$insert_success) {
            $insert_success = $dbh->do($insert);
            unless ($insert_success) {

                #               print STDERR $insert;
                sleep .1;
            }
        }

        $pm->finish;
    }
    $pm->wait_all_children();
    $dbh->disconnect();
    my $finished_time = new Benchmark;
    my $completion_time = timestr(timediff($finished_time, $start_time));
    say STDERR "Completed in: $completion_time";

    CoGe::Accessory::Web::write_log("Completed generating ks data.",
        $cogeweb->logfile);

    return 0;
}

sub get_ks_data
{
    my %opts    = @_;
    my $db_file = $opts{db_file};
    my %ksdata;
    unless (-r $db_file) {
        return \%ksdata;
    }
    CoGe::Accessory::Web::write_log("\tconnecting to ks database $db_file",
        $cogeweb->logfile);
    my $select = "select * from ks_data";
    my $dbh    = DBI->connect("dbi:SQLite:dbname=$db_file", "", "");
    my $sth    = $dbh->prepare($select);
    $sth->execute();
    CoGe::Accessory::Web::write_log(
        "\texecuting select all from ks database $db_file",
        $cogeweb->logfile);
    my $total   = 0;
    my $no_data = 0;

    while (my $data = $sth->fetchrow_arrayref) {
        $total++;
        if ($data->[3] eq "") {
            $no_data++;

            #       next; #uncomment to force recalculation of missing data
        }

        $ksdata{$data->[1]}{$data->[2]} =
          $data->[3]
          ? {
            dS      => $data->[3],
            dN      => $data->[4],
            'dN/dS' => $data->[5]}
          : {};    # unless $data->[3] eq "";
    }
    CoGe::Accessory::Web::write_log("\tgathered data from ks database $db_file",
        $cogeweb->logfile);
    $sth->finish;
    undef $sth;
    $dbh->disconnect();
    CoGe::Accessory::Web::write_log("\tdisconnecting from ks database $db_file",
        $cogeweb->logfile);
    return \%ksdata;
}

sub gen_ks_blocks_file
{
    my %opts    = @_;
    my $infile  = $opts{infile};
    my $dbfile  = $opts{dbfile};
    my $outfile = $opts{outfile};

    my $ksdata = get_ks_data(db_file => $dbfile);
    $/ = "\n";
    open(IN,  $infile);
    open(OUT, ">" . $outfile);
    print OUT
      "#This file contains synonymous rate values in the first two columns:\n";
    my @block;
    my $block_title;

    while (<IN>) {
        if (/^#/) {
            my $output = process_block(
                ksdata => $ksdata,
                block  => \@block,
                header => $block_title) if $block_title;
            print OUT $output if $output;
            @block       = ();
            $block_title = $_;

            #beginning of a block;
        } else {
            push @block, $_;
        }
    }
    close IN;
    my $output = process_block(
        ksdata => $ksdata,
        block  => \@block,
        header => $block_title) if $block_title;
    print OUT $output;
    close OUT;

    return 0;
}

sub initialize_nwalign_servers
{
    my %opts       = @_;
    my $start_port = $opts{start_port};
    my $procs      = $opts{procs};
    my @ports;
    use IO::Socket;

    for (1 .. $procs) {
        unless (
            IO::Socket::INET->new(
                PeerAddr => 'localhost',
                PeerPort => $start_port)
          ) {
            system("$NWALIGN --server $start_port &");
        }

        push @ports, $start_port;
        $start_port++;
    }
    sleep 1;    #let them get started;
    return \@ports;
}

sub process_block
{
    my %opts   = @_;
    my $ksdata = $opts{ksdata};
    my $block  = $opts{block};
    my $header = $opts{header};
    my $output;
    my @ks;
    my @kn;
    foreach my $item (@$block) {
        my @line = split /\t/,   $item;
        my @seq1 = split /\|\|/, $line[1];
        my @seq2 = split /\|\|/, $line[5];
        my $ks   = $ksdata->{$seq1[6]}{$seq2[6]};
        if (defined $ks->{dS}) {
            unshift @line, $ks->{dN};
            unshift @line, $ks->{dS};
            push @ks, $ks->{dS} unless $ks->{dS} eq "NA";
            push @kn, $ks->{dN} unless $ks->{dN} eq "NA";
        } else {
            unshift @line, "undef";
            unshift @line, "undef";
        }
        $output .= join "\t", @line;
    }
    my $mean_ks = 0;
    if (scalar @ks) {
        map {$mean_ks += $_} @ks;
        $mean_ks = sprintf("%.4f", $mean_ks / scalar @ks);
    } else {
        $mean_ks = "NA";
    }
    my $mean_kn = 0;
    if (scalar @kn) {
        map {$mean_kn += $_} @kn;
        $mean_kn = sprintf("%.4f", $mean_kn / scalar @kn);
    } else {
        $mean_kn = "NA";
    }
    chomp $header;
    $header .= "  Mean Ks:  $mean_ks\tMean Kn: $mean_kn\n";
    $header .= join("\t",
        "#Ks",
        qw(Kn a<db_genome_id>_<chr> chr1||start1||stop1||name1||strand1||type1||db_feature_id1||percent_id1 start1 stop1 b<db_genome_id>_<chr> chr2||start2||stop2||name2||strand2||type2||db_feature_id2||percent_id2 start2 stop2 eval block_score GEVO_link)
    ) . "\n";
    return $header . $output;
}
