# -*- perl -*-
use strict;

use CoGeX;
my $connstr = 'dbi:mysql:dbname=DB;host=HOST;port=PORT';
my $s = CoGeX->connect($connstr, 'USER', 'PASSWORD' );
my %seen;

#$s->storage->debug(1);
open( ACCN, "<", "TAIR_sequenced_genes" );
while ( my $line = <ACCN> ) {
    chomp $line;
    my @accn = split( "\t", $line );
    my $accn = $accn[0];
    next unless $accn =~ /^AT\dG\d{5}$/i;
    foreach my $feat (
        $s->resultset('Feature')->search(
            {
                'feature_names.name' => { 'like' => "$accn%" },
                'feature_type.name'  => { 'like' => [ 'CDS', '%' ] },
                'dataset.version' => 7
            },
            {
                prefetch => [ 'feature_names', 'feature_type' ],
                join     => 'dataset',
                order_by => ['feature_type.name']
            }

        )
      )
    {
        if (
            scalar(
                grep { $seen{ substr( uc( $_->name ), 0, 9 ) } }
                  $feat->feature_names
            )
          )
        {
            next;
        }
        my @names =
          sort grep { $_ =~ /^AT\dG\d{5}\.\d$/i }
          map { $_->name } $feat->feature_names;
        my $name = @names[0];
        if ( !$name ) {
            @names =
              sort grep { $_ =~ /^AT\dG\d{5}$/i }
              map { $_->name } $feat->feature_names;
            $name = @names[0];    # . "|" . $feat->feature_id;
            if ( !$name ) {
                print STDERR join( "\t", map { $_->name } $feat->feature_names )
                  . "\n";
                exit();
            }
        }
        map { $seen{ substr( uc( $_->name ), 0, 9 ) } = 1 }
          $feat->feature_names;

        #map { print $_  . "\t" } @names;
        print ">" . $name . "\t" . $feat->feature_type->name . "\n";
        print $feat->genomic_sequence() . "\n";

    }
}
