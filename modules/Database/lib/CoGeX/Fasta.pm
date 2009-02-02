package CoGeX::Fasta;

use strict;

sub get_seq {
    my $module = shift;
    my $opts = shift;
    my $fastacmd = $opts->{fastacmd} || 'fastacmd';
    my $blastdb = $opts->{blastdb} || $opts->{db};
    my $seqid   = $opts->{seqid} || $opts->{chr} || $opts->{chromosome}; # chr 

    my $cmd = "$fastacmd -d $blastdb -s $seqid -l 999999999 ";
    if($opts->{start} && ($opts->{stop} || $opts->{end})){
        $cmd .= "-L " . $opts->{start} . "," . ($opts->{stop} || $opts->{end} ). " ";
    }
    if($opts->{reverse_complement} || $opts->{rc}){
        $cmd .= "-S 2";
    }
    #print $cmd . "\n";
    open(FASTA, $cmd . "|") || die "cant run $cmd";
    # get rid of the header line...
    <FASTA>;
    my $seq = <FASTA>;
    return $seq;

}
1;
