package CoGe::Core::Sequence;

use strict;
use warnings;

use File::Slurp;
use JSON qw(decode_json encode_json);

use CoGe::Core::Chromosomes;
use CoGe::Core::Storage qw(get_genome_file get_genome_path);

##################################################

sub new {
	my ($class, $gid) = @_;
	my $self = {};
    $self->{gid} = $gid;
	return bless $self, $class;
}

##################################################

sub cache_features {
    my ($self, $features, $type_id, $db) = @_;
    my $chromosomes = CoGe::Core::Chromosomes->new($self->{gid})->hash;
    open FILE, '>' . $features;
    my $datasets = $db->storage->dbh->selectall_arrayref('SELECT dataset_id FROM dataset_connector WHERE genome_id=' . $self->{gid});
    for my $dataset (@$datasets) {
        my $sth = $db->storage->dbh->prepare('SELECT location.chromosome,location.start,location.stop,location.strand FROM feature JOIN location USING (feature_id) WHERE dataset_id=' . $dataset->[0] . ' AND feature_type_id=' . $type_id . ' AND feature.chromosome=location.chromosome AND feature.strand=location.strand AND location.start>=feature.start AND location.start<=feature.stop AND location.stop>=feature.start AND location.stop<=feature.stop');
        $sth->execute;
        my $loc = $sth->fetchrow_arrayref;
        while ($loc) {
            my $chr = $chromosomes->{$loc->[0]};
            my $offset = $chr->{offset};
            print FILE ($offset + $loc->[1]) . "\t" . ($loc->[2] - $loc->[1] + 1) . "\t" . $loc->[2] . "\n";
            $loc = $sth->fetchrow_arrayref;
        }
    }
    close FILE;
}

##################################################

sub calc_percentages {
    my ($self, $type_id, $db) = @_;
    my $seq = read_file(get_genome_file($self->{gid}));
    my $features = get_genome_path($self->{gid}) . 'feat' . $type_id . '.fi';
    $self->cache_features($features, $type_id, $db) unless -r $features;
    my $gc = 0;
    my $at = 0;
    my $nx = 0;
    my $num = 0;
    open FILE, $features;
    while (<FILE>) {
        my @a = split /\t/, $_;
        my $s = substr($seq, $a[0], $a[1]);
        $gc += $s =~ tr/gcGC//;
        $at += $s =~ tr/atAT//;
        $nx += $s =~ tr/nxNX//;
        $num++;
        warn $num if $num % 1000 == 0;
    }
    close FILE;
    my $total = $gc + $at + $nx;
    if ($total) {
        $gc = $gc / $total * 100;
        $at = $at / $total * 100;
        $nx = $nx / $total * 100;
    }
    return { gc => $gc, at => $at, nx => $nx, total => $total };
}

##################################################

sub calc_wobble {
    my ($self, $type_id, $db) = @_;
    my $seq = read_file(get_genome_file($self->{gid}));
    my $features = get_genome_path($self->{gid}) . 'feat' . $type_id . '.fi';
    $self->cache_features($features, $type_id, $db) unless -r $features;
    my $gc = 0;
    my $at = 0;
    my $nx = 0;
    my $num = 0;
    open FILE, $features;
    while (<FILE>) {
        my @a = split /\t/, $_;
        my $s = substr($seq, $a[0], $a[1]);
        $gc += $s =~ tr/gcGC//;
        $at += $s =~ tr/atAT//;
        $nx += $s =~ tr/nxNX//;
        $num++;
        warn $num if $num % 1000 == 0;
    }
    close FILE;
    my $total = $gc + $at + $nx;
    if ($total) {
        $gc = $gc / $total * 100;
        $at = $at / $total * 100;
        $nx = $nx / $total * 100;
    }
    return { gc => $gc, at => $at, nx => $nx, total => $total };
}

##################################################

sub get_percentages {
    my ($self, $type_id, $db) = @_;
    my $percentages_file = get_genome_path($self->{gid}) . $type_id . '_percentages' . '.json';
    return decode_json(read_file($percentages_file)) if -r $percentages_file;

    my $percentages = $self->calc_percentages($type_id, $db);
    open FILE, '>' . $percentages_file;
    print FILE encode_json($percentages);
    close FILE;
    return $percentages;
}

1;
