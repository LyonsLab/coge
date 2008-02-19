#! /usr/bin/perl -w

use strict;
use Data::Dumper;
use CoGeX;
use Getopt::Long;

my $connstr = 'dbi:mysql:dbname=genomes;host=biocon;port=3306';
my$coge = CoGeX->connect($connstr, 'cnssys', 'CnS' );

my $DEBUG = 1;
my $GO=0;
my $ds1; my $ds2;
my $add_gene = 0;
GetOptions ( 
             "ds1=i"            => \$ds1,
             "ds2=i"            => \$ds2,
	         "go=s"             => \$GO,
	         "add_gene_feature" => \$add_gene,
	   );

exit() unless $ds1 && $ds2;

my $dsq = $coge->resultset('Dataset')->find($ds1);
my $dss = $coge->resultset('Dataset')->find($ds2);
my ($anno_type) = $coge->resultset('AnnotationType')->search({name=>"note"});
my $cds_id = 202;

print STDERR $dsq->dataset_id . "\n";
print STDERR $dss->dataset_id . "\n";

warn "-go flag is not true, nothing will be added to the database.\n" unless $GO;

my %data;
my %annos;
my %qseen;
my %sseen;
# line like:
# OS02G04490,SB10G029285, 1986128,1986156,59117428,59117400|1986133,1986163,59117429,59117399|1992129,1992147,59096365,59096347|1984218,1984236,59096368,59096350|1992128,1992143,59096365,59096350|1982831,1982849,59113429,59113411|1995931,1995957,59090974,59090948|1993267,1993287,59103979,59103959|2001316,2001330,59109300,59109286|1982780,1982794,59109486,59109472|1978893,1978907,59110202,59110188|2000840,2000857,59114956,59114939|1979357,1979377,59116046,59116026|1986130,1986153,59117429,59117406
while (my $line = <>) {
    next if $line =~ /^#/;
    chomp $line;
    $line =~ s/,\s/,/g;
    my ($qname, $sname, @line)= split(/[,|\|]/, $line);
    my ($qchr) = $qname =~ /(\d+)G/g;
    my ($schr) = $sname =~ /(\d+)G/g;
    $qchr =~ s/^0//;
    $schr =~ s/^0//;

    my $tot = 0;
    for(my $i=0; $i<scalar(@line); $i+=4){
        my $strand = 1;
        $tot++;
        my ($qstart, $qstop, $sstart, $sstop) = @line[$i .. $i + 4];
        print STDERR "$qname, $sname, $qchr, $schr, $qstart, qstop: $qstop, $sstart, $sstop\n";
        if ($sstart > $sstop){ 
            $strand = -1;
            ($sstart, $sstop) = ($sstop, $sstart);
        }

        # dont do repeats.
        if ($qseen{ $qstart . "," . $qstop }){ next; }
        if ($sseen{ $sstart . "," . $sstop }){ next; }
        $qseen{ $qstart . "," . $qstop } = 1;
        $sseen{ $sstart . "," . $sstop } = 1;
        
  	    my $qfeat = $dsq->add_to_features({
					 feature_type_id => $cds_id,
					 start=>$qstart, stop=>$qstop,
					 chromosome=>$qchr, strand=>$strand,
        }) if $GO;
  	    my $sfeat = $dss->add_to_features({
					 feature_type_id => $cds_id,
					 start=>$sstart, stop=>$sstop,
					 chromosome=>$schr, strand=>$strand,
        }) if $GO;

	    my $qloc = $qfeat->add_to_locations( {
						   start      => $qstart,
						   stop       => $qstop,
						   strand     => $strand,
						   chromosome => $qchr
						  }) if $GO;
	    my $sloc = $sfeat->add_to_locations( {
						   start      => $sstart,
						   stop       => $sstop,
						   strand     => $strand,
						   chromosome => $schr
						  }) if $GO;

        my $qfeat_name = $qfeat->add_to_feature_names({
						     name=>$qname . "_CNS_" . $i/4,
						    }) if $GO ;
        my $sfeat_name = $sfeat->add_to_feature_names({
						     name=>$sname . "_CNS_" . $i/4,
						    }) if $GO ;
        
        print $qfeat->feature_id . "\n" if $qfeat;
        print $sfeat->feature_id . "\n" if $qfeat;
      }
      #if ($tot > 40){ exit(); }
}
