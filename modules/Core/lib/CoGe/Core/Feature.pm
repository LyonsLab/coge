package CoGe::Core::Feature;

use strict;
use warnings;

use Data::Dumper;

BEGIN {
    our ( @EXPORT, @EXPORT_OK, @ISA, $VERSION );
    require Exporter;

    $VERSION = 0.1;
    @ISA = qw( Exporter );
    @EXPORT = qw( search_features get_feature );
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
        $features{$fid} = get_feature( feature => $fn->feature, fast => 1 );
    }

    return wantarray ? values %features : [ values %features ];
}

sub get_feature {
    my %opts = @_;
    my $db = $opts{db};
    my $fid = $opts{fid};
    my $feature = $opts{feature};
    my $fast = $opts{fast};
    
    unless ($feature) {
        $feature = $db->resultset("Feature")->find($fid);
        return unless (defined $feature);
    }
    
    my $genome = $feature->dataset->first_genome;
    unless ($genome) {
        return { error => { Error => "No genome for this feature" } };
    }
    
    my %data = (
        id => int($feature->id),
        type => $feature->type->name,
        start => $feature->start,
        stop => $feature->stop,
        strand => $feature->strand,
        chromosome => $feature->chromosome,
        genome => {
            id => $genome->id,
            version => $genome->version,
            organism => $genome->organism->name
        }
    );
    
    my @names = $feature->names;
    my $num_names = scalar @names;
    if ($num_names > 10) {
        push @{$data{warning}}, "Names truncated, too many names to list ($num_names)";
        #@names = splice(@names, 0, 10);
    }
    else {
        $data{names} = \@names;
    }
    $data{primary_name} = $feature->primary_name->name if $feature->primary_name;
    
    if (!$fast) {
        $data{sequence} = $feature->genomic_sequence;
        my @annos;
        foreach ($feature->annos) {
            my %anno;
            $anno{text} = $_->annotation;
            $anno{type} = $_->type->name;
            $anno{category} = $_->type->group->name if ($_->type->group);
            push @annos, \%anno;
        }
        $data{annotations} = \@annos;
    }
    
    return \%data;
}

1;
