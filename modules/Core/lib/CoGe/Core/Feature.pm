package CoGe::Core::Feature;

use strict;
use warnings;

use Data::Dumper;

BEGIN {
    our ( @EXPORT, @EXPORT_OK, @ISA, $VERSION );
    require Exporter;

    $VERSION = 0.1;
    @ISA = qw( Exporter );
    @EXPORT = qw( search_features get_feature get_genome_for_feature );
}

sub search_features {
    my %opts = @_;
    my $db = $opts{db};
    my $user = $opts{user};
    my $search_term = $opts{search_term};
    return unless ($db and $search_term);
    my $exact = $opts{exact}; # find exact matches, default is fuzzy matching
    
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
    
    my @results;
    if ($exact) {
        @results = $db->resultset('FeatureName')->search({'me.name' => $search_term}, $prefetch);
    }
    else { # fuzzy match (full-text search)
        @results = $db->resultset('FeatureName')->search(undef, $prefetch)->search_literal( "MATCH(me.name) AGAINST ('\"$search_term\"')" );
    }
    
    my %features;
    foreach my $fn (@results) {
        my $feature = $fn->feature;
        next unless $feature;

        my $fid = $fn->feature_id;
        my $f = get_feature( user => $user, feature => $feature, fast => 1 );
        $features{$fid} = $f unless $f->{error};
        #TODO add proper error handling and response
    }

    return wantarray ? values %features : [ values %features ];
}

sub get_feature {
    my %opts = @_;
    my $db = $opts{db};
    my $user = $opts{user};
    my $fid = $opts{fid};
    my $feature = $opts{feature};
    my $fast = $opts{fast};
    
    unless ($feature) {
        $feature = $db->resultset("Feature")->find($fid);
        return unless (defined $feature);
    }
    
    my $genome = get_genome_for_feature(feature => $feature, user => $user);
    unless ($genome) {
        return { error => { message => "Access Denied" } };
    }
    
    my %data = (
        id => int($feature->id),
        type => $feature->type->name,
        start => int($feature->start),
        stop => int($feature->stop),
        strand => int($feature->strand),
        chromosome => $feature->chromosome,
        genome => {
            id => int($genome->id),
            version => $genome->version,
            organism => $genome->organism->name
        }
    );
    
    if (!$fast) {
        my @names = $feature->names;
        my $num_names = scalar @names;
        if ($num_names > 50) {
            push @{$data{warning}}, "Names truncated, too many names to list ($num_names)";
            #@names = splice(@names, 0, 10);
        }
        else {
            $data{names} = \@names;
            #$data{primary_name} = $feature->primary_name->name if $feature->primary_name;
        }
        
        my $seq = $feature->genomic_sequence;
        $data{sequence} = $seq;
        my $l = length($seq);
        if ($l > 0) {
            my $gc = $seq =~ tr/gcGC//;
            my $at = $seq =~ tr/atAT//;
            my $nx = $seq =~ tr/nxNX//;
            $data{percentages} = {
                gc => $gc / $l * 100,
                at => $at / $l * 100,
                nx => $nx / $l * 100
            };
        }
        
        my @annos;
        foreach ($feature->annos) {
            my %anno;
            $anno{value} = $_->annotation;
            $anno{type} = $_->type->name;
            $anno{category} = $_->type->group->name if ($_->type->group);
            push @annos, \%anno;
        }
        $data{annotations} = \@annos if (@annos);
        
        my @locations;
        foreach ($feature->locations) {
            my %loc;
            $loc{start}  = int($_->start);
            $loc{stop}   = int($_->stop);
            $loc{strand} = int($_->strand);
            push @locations, \%loc;
        }
        $data{locations} = \@locations if (@locations);
    }
    
    return \%data;
}

sub get_genome_for_feature { # can be used to test whether user has access to the feature
    my %opts = @_;
    my $feature = $opts{feature};
    my $user = $opts{user};
    return unless $feature;
    
    my $genome;
    foreach ($feature->genomes) {
        next if ($_->deleted);
        next unless ( !$_->restricted || (defined $user && $user->has_access_to_genome($_)) );
        $genome = $_;
        last;
    }
    
    return $genome;
}

1;
