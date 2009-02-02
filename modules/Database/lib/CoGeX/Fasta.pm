package CoGeX::Fasta;

use strict;
use Data::Dumper;

sub get_seq {
    my $module = shift;
    print Dumper $module;
    my $opts = shift;
    print Dumper $opts;
    print $opts->{'blastdb'} . "\n";
    my $fastacmd = $opts->{fastacmd} || 'fastacmd';
    my $blastdb = $opts->{blastdb} || $opts->{db};
    my $seqid   = $opts->{seqid} || $opts->{chr} || $opts->{chromosome}; # chr 
    print $blastdb . "\n";

    my $cmd = "$fastacmd -d $blastdb -s $seqid ";
    if($opts->{start} && ($opts->{stop} || $opts->{end})){
        $cmd .= "-L " . $opts->{start} . "," . ($opts->{stop} || $opts->{end} ). " ";
    }
    if($opts->{reverse_complement} || $opts->{rc}){
        $cmd .= "-S 2";
    }
    print $cmd . "\n";
    return `$cmd`;

}
1;
