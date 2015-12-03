# mdb 12/2/15: what is this for?  It's commented out in service.pl.
package CoGe::Services::Data::Feature;
use v5.14;
use base 'CGI::Application';

use CoGeX;
use CoGe::Accessory::Web;
use CoGe::Core::Storage qw( get_experiment_data );
use File::Spec;
use File::Slurp;
use Data::Dumper;
use List::Util qw(max);
use JSON qw(encode_json);
use POSIX qw(floor ceil);

our ($FEATURE, $P);

sub setup {
    my $self = shift;
    $self->run_modes( 'get' => 'get' );
    $self->mode_param('rm');
}

sub valid_experiment {
    my $exp = shift;

    return if $exp->restricted;

    #XXX: Make sure this is correct
    #FIXME: Don't use hardcoded value for data type
    # Exclude polymorphic and alignment data
    return unless $exp->data_type == 1;

    return 1;
}

sub parse_data {
    my $exp = shift;
    my $output = shift;
    my $eid = $exp->id;
    my (@s1, @s2);

    foreach (@$output) {
        chomp;

        if (/^\"/) {    #if (/^\"$chr\"/) { # potential result line
            s/"//g;
            my @items = split(/,\s*/);
            next if ( @items != 6 );
            # || $items[0] !~ /^\"?$chr/); # make sure it's a row output line
            #
            for ( my $i = 0 ; $i < @items ; $i++ ) {
                $items[$i] = 1 if $items[$i] !~ /\w/; # what's this for?
            }

            my ( $chr, $start, $end, $strand, $value1, $value2 ) = @items;

            $end = $start + 1 if ( $end == $start ); #FIXME revisit this
            $strand = -1 if ( $strand == 0 );
            $value1 = $strand * $value1;

            push @s1, $value1 + 0.0;
            push @s2, $value2 + 0.0;
        }
    }

    my @score1 = sort { $a<=>$b } @s1;
    my @score2 = sort { $a<=>$b } @s2;

    return [$eid, {}] unless @score1 and @score2;

    my $results = {};

    # Correctly set upper and lower confidence bands
    #
    $results->{score1} = {
        value => median(@score1),
        upper => $score1[-1],
        lower => $score1[0],
    };

    $results->{score2} = {
        value => median(@score2),
        upper => $score2[-1],
        lower => $score2[0],
    };

    return [$eid, $results];
}

sub fetch_data {
    return $_, get_experiment_data(
        chr => $FEATURE->chr,
        start => $FEATURE->start,
        stop => $FEATURE->stop,
        eid => $_->id,
        data_type => $_->data_type);
}

sub experiment_info {
    my $info = {
        id   => $_->id,
        name => $_->name,
        description => $_->desc,
        source => $_->source->name,
    };

    # Add source link or experiment link
    if ($_->source->link) {
        $info->{source_link} = $_->source->link;
    } else {
        my $base = $P->{SERVER};
        $base =~ s/\/$//;

        $info->{source_link} = $base . "/ExperimentView.pl?eid=" . $_->id;
    }

    return $info;
}

sub median {
    # assumes @a is sorted
    my @a = @_;

    # Compute median
    my $middle = (@a+1)/2 - 1;
    my $a1 = $a[floor($middle)];
    my $a2 = $a[ceil($middle)];
    my $med = ($a1+$a2) / 2.0;

    return $med;
}

sub get {
    my $self = shift;
    my ($experiments, $data) = {} x 2;
    my $fid = $self->param('fid');

    # Connect to the database
    my ( $db, $user, $config ) = CoGe::Accessory::Web->init();

    # Set global variables
    $P = CoGe::Accessory::Web::get_defaults($config);
    $FEATURE = $db->resultset("Feature")->find($fid);

    return encode_json({}) unless $FEATURE;

    my @genomes = $FEATURE->genomes;
    my @experiments = grep { valid_experiment $_ } $genomes[0]->experiments;
    my @data = map { parse_data(fetch_data($_)) } @experiments;

    # Find all feature experimental data
    my @filtered =  grep { defined @$_[1]->{score1}; } @data;

    # Convert to hash ref (#experiment id, scores)
    my $filtered_data = { map {
                @$_[0], @$_[1];
        } @filtered
    };

    # Filter experiment info by features found
    my @filtered_info = grep {
        exists $filtered_data->{$_->{id}}
    } map(&experiment_info, @experiments);

    return encode_json({
        feature_id  => int($fid),
        gene        => $FEATURE->name,
        experiments => \@filtered_info,
        data        => $filtered_data
    });
}

1;
