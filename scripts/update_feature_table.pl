#!/usr/bin/perl -w

use CoGeX;
use strict;
use Data::Dumper;
$| = 1;

my $connstr = 'dbi:mysql:dbname=DB;host=HOST;port=PORT';
my $s = CoGeX->connect($connstr, 'USER', 'PASSWORD' );

$s->storage->debug(0);

my $FCOUNT = $s->resultset('Feature')->count();

my $count = 1;
while ( $count > 0 ) {
    my $rs = $s->resultset('Feature')->search(
        undef,
        {
            page => $count,
            rows => 1
        }
    );
    my $f = $rs->next();

    #  print "id = ", $f->id(), "\n";
    #  print "start = ", $f->start(), "\n";
    #  print "stop = ", $f->stop(), "\n";
    #  print "strand = ", $f->strand(), "\n";
    #  print "chromosome = ", $f->chromosome(), "\n";

    $f->update(
        {
            fstart      => $f->start(),
            fstop       => $f->stop(),
            fchromosome => $f->chromosome(),
            fstrand     => $f->strand()
        }
    );
    $count += 1;
}

#my $feature = $rc->next();
#while ( my $feature = $rc->next() ) {
#  print STDERR $feature, "\n";
#  $feature->update( start      => $feature->start(),
#                    stop       => $feature->stop(),
#                    chromosome => $feature->chromosome(),
#                    strand     => $feature->strand()
#                  );
#}

#sub esearch {
#
#    my @results;
#    my $rs = $s->resultset('Feature')->esearch( {
#                        'feature_type.name' =>   'gene' ,
#                        'feature_names.name' => ['-and' ,
#                                { '-like' => 'At2g2%'} ,
#                                { '-not_like' => ['%.%'] }
#                            ]
#                        }
#                        );
#
#    #print "got resultset\n";
#    while (my $feat = $rs->next()){
#        my $fn = $feat->feature_names;
#        foreach my $name ($fn->next()){
#            push(@results,  $name->name );
#        }
#        push(@results, $feat->dataset->organism->name );
#        map { push(@results, $_->annotation)  } $feat->annotations;
#    }
#}
