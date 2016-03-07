package CoGe::Services::API::JBrowse::GCcontent;

use Mojo::Base 'Mojolicious::Controller';
use Mojo::JSON;
use CoGe::Core::Storage qw( get_genome_seq );
use CoGe::Services::Auth qw(init);
use URI::Escape qw(uri_unescape);
use List::Util qw( min max );

sub stats_global {
    my $self = shift;
    print STDERR "JBrowse::GCcontent::stats_global\n";
    $self->render(json => { 
        "featureDensity" => 0.02,
        "scoreMin" => 0,
        "scoreMax" => 1,
        "featureDensityByRegion" => 50000,
    });
}

sub features {
    my $self = shift;
    my $gid  = $self->stash('id');
    my $chr  = $self->stash('chr');
    $chr = uri_unescape($chr) if (defined $chr);
    my $start = $self->param('start');
    my $end   = $self->param('end');
    my $scale = $self->param('scale');
    my $basesPerSpan = $self->param('basesPerSpan');
    my $len   = $end - $start;
    print STDERR "JBrowse::GC_content::features gid=$gid chr=$chr start=$start end=$end\n";

    # Check params
    my $null_response = $self->render(json => { "features" => [] });
    if ( $end <= 0 ) {
        return $null_response;
    }

    # Authenticate user and connect to the database
    my ($db, $user) = CoGe::Services::Auth::init($self);

    # Retrieve genome
    my $genome = $db->resultset('Genome')->find($gid);
    return unless $genome;

    # Adjust location - note: incoming coordinates are interbase!
    $start = 0 if ( $start < 0 );
    my $chrLen = $genome->get_chromosome_length($chr);
    $start = min( $start, $chrLen );
    $end   = min( $end,   $chrLen );
    if ( $start == $end ) {
        return $null_response;
    }

    # Check permissions
    if ( $genome->restricted
        and ( not defined $user or not $user->has_access_to_genome($genome) ) )
    {
        return $null_response;
    }

    # Extract requested piece of sequence file
    my $seq = get_genome_seq(
        gid   => $gid,
        chr   => $chr,
        start => $start + 1, # convert from interbase to base
        stop  => $end
    );
    
    # Create bins and measure %gc content per bin
    my $score = 0;
    $basesPerSpan = 1 if $basesPerSpan < 0;
    my $numBins = min(100, int($len / $basesPerSpan));
    my $binSz = max(1, int($len / $numBins));
    print STDERR "matt: numBins=$numBins binSz=$binSz\n";
    my @bins;
    my $offset = 0;
    while ($offset < $len) {
        my $subSeq = substr($seq, $offset, $binSz);
        
        # Count occurrences of bases
        my %unique;
        map { $unique{lc($_)}++ } $subSeq =~ /(.)/ig;
        
        # Determine most frequent base
        my $mostFreqBase;
        map { 
            $mostFreqBase = $_ 
                if ( !$mostFreqBase || $unique{$_} > $unique{$mostFreqBase} )
        } keys %unique;
        
        # Count G's and C's
        my $gc_count =  $unique{'g'} // 0;
           $gc_count += $unique{'c'} // 0;
           
        push @bins, {
            start => $start + $offset,
            end => $start + $offset + $binSz,
            score => $gc_count / $binSz,
            nucleotide => $mostFreqBase
        };
        
        $offset += $binSz;
    }

    $self->render(json => { "features" => \@bins });
}

1;
