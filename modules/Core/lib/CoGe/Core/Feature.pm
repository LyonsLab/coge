package CoGe::Core::Feature;

use strict;
use warnings;

use Data::Dumper;

BEGIN {
    our ( @EXPORT, @EXPORT_OK, @ISA, $VERSION );
    require Exporter;

    $VERSION = 0.1;
    @ISA = qw( Exporter );
    @EXPORT = qw( search_features );
}

sub search_features {
    my $db = shift;
    my $search_term = shift;
    return unless ($db and $search_term);
    
    my $prefetch = {
        prefetch => [
            {
                'feature' => {
                    'dataset' =>
                      { 'dataset_connectors' => { 'genome' => 'organism' } }
                }
            }
        ]
    };
    
    my @results = $db->resultset('FeatureName')->search( { 'me.name' => $search_term }, $prefetch )->search_literal( 'MATCH(me.name) AGAINST (?)', $search_term );
    #print STDERR Dumper \@results, "\n";
    
    my %features;
    foreach my $fn (@results) {
        my $fid = $fn->feature_id;
        my $feat = $fn->feature;
        my $genome = $feat->dataset->first_genome;
        next unless $genome;
        
        $features{$fid} = {
            id => $fid,
            type => $feat->type->name,
            start => $feat->start,
            stop => $feat->stop,
            strand => $feat->strand,
            chromosome => $feat->chromosome,
            genome => {
                id => $genome->id,
                version => $genome->version,
                
            }
        };
        #$features{$fid}{primary_name} = $feat->primary_name->name if $feat->primary_name;
        
        my $feat_names = $feat->names;
        my $num_names = scalar @$feat_names;
        if ($num_names > 10) {
            push @{$features{$fid}{warnings}}, "Names truncated, too many names to list ($num_names)";
            #$feat_names = splice($feat_names, 0, 10);
        }
        else {
            $features{$fid}{names} = $feat_names;
        }
    }
    print STDERR Dumper \%features, "\n";

    return wantarray ? values %features : [ values %features ];
}

1;
