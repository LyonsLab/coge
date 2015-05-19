# -*- perl -*-
use strict;

use CoGeX;
my $connstr = 'dbi:mysql:dbname=DB;host=HOST;port=PORT';
my $s = CoGeX->connect($connstr, 'USER', 'PASSWORD' );

# example for getting 10,000-mer chunks.
my $dataset = [565];
get_10kmers( 'sorghum', $dataset );

# example for getting accns from a file:
# only gets OS01G12345, not OS01G12345.1
# gets the CDS if available, else gene else *RNA, psuedogene
# prints to an output file $org.fasta.
# needs a file named $org_accns.txt which is a list of accn names.
# can be made with an sql query like:

#    SELECT DISTINCT(fn.name) from feature_name fn, feature f where
#    fn.feature_id = f.feature_id AND f.dataset_id BETWEEN 582 AND 593 AND
#    fn.name LIKE 'Os%%g%' AND fn.name NOT LIKE '%.%' ORDER BY fn.name;

my $name_re  = '^OS(\d\d)G\d{5}$';
my $name_len = 10;
my $org      = 'rice';
my @datasets = 582 .. 593;
get_accn_locs( $name_re, $name_len, $org, \@datasets );

sub get_accn_locs {
    my ( $name_re, $name_len, $org, $datasets ) = @_;
    my %seen;
    my %order;
    open( ACCN, "<", $org . "_accns.txt" )
      or die "must have a list of accns in $org" . "_accns.txt";
    open( LOC, ">", $org . ".fasta" );
    while ( my $accn = <ACCN> ) {
        chomp $accn;
        next unless $accn =~ /$name_re/i;
        foreach my $feat (
            $s->resultset('Feature')->search(
                {
                    'feature_names.name' => $accn,
                    'me.dataset_id'      => { 'IN' => $datasets }
                },
                {
                    prefetch => [ 'feature_names', 'feature_type' ],
                    order_by => ['feature_type.name']
                }

            )
          )
        {
            my @feat_names = map { $_->name } $feat->feature_names();
            if (
                scalar(
                    grep { $seen{ substr( uc($_), 0, $name_len ) } } @feat_names
                )
              )
            {
                next;
            }
            my @names = sort grep { $_ =~ /$name_re/i } @feat_names;
            my $name = uc( @names[0] );
            if ( !$name ) {
                @names = sort grep { $_ =~ /$name_re/i } @feat_names;
                $name = @names[0];    # . "|" . $feat->feature_id;
                if ( !$name ) {
                    print STDERR join( "\t", @feat_names ) . "\n";
                    exit();
                }
            }
            my ($chr) = $name =~ /$name_re/;
            $chr = sprintf( "%02i", $chr );
            map { $seen{ substr( uc($_), 0, $name_len ) } = 1 } @feat_names;
            my $start = sprintf( "%09i", $feat->start );
            my $stop  = sprintf( "%09i", $feat->stop );
            my $strand = ( $feat->strand =~ /-/ ) ? -1 : 1;

            my $header =
                $chr . "||"
              . $name . "||"
              . $start . "||"
              . $stop . "||"
              . $strand . "||"
              . uc( $feat->feature_type->name ) . "||"
              . $feat->feature_id;

            $order{$header} = 1;
            print LOC ">" . $header . "\n";

            print LOC $feat->genomic_sequence() . "\n";
        }
    }
    open( ORDER, ">", $org . ".order" );
    map { print ORDER $_ . "\n" } sort keys %order;
    close(ORDER);
}

sub get_10kmers {
    my ( $org, $datasets ) = @_;
    my %files;
    my %order;
    foreach my $gs (
        $s->resultset('GenomicSequence')->search(
            {
                'dataset_id' => { 'IN'       => $datasets },
                'chromosome' => { 'NOT LIKE' => '%super%' }
            },
            { order_by => ['start'] }
        )
      )
    {
        my ($chr) = $gs->chromosome =~ /(\d+)/;
        my $file = $org . "_" . $chr . ".fasta";
        my $FH;
        if ( !$files{$file} ) {

            open( $FH, ">", $file );
            $files{$file} = $FH;
        }
        else {
            $FH = $files{$file};
        }
        my $header = $chr . "||NONE||" . sprintf( "%09i", $gs->start() );
        print $FH ">" . $header . "\n";
        print $FH $gs->sequence_data() . "\n";

        $order{$header} = 1;
    }
    open( ORDER, ">", $org . ".order" );
    map { print ORDER $_ . "\n" } sort keys %order;
    close(ORDER);

}
