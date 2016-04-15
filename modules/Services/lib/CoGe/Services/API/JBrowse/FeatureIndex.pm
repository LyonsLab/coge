package CoGe::Services::API::JBrowse::FeatureIndex;

use strict;
use warnings;

use CoGeDBI qw( get_dataset_ids );
use Path::Class;

sub new {
	my ($class, $gid, $type, $chromosome, $conf) = @_;
	my $self = {
		chromosome => $chromosome,
		conf => $conf,
		gid => $gid,
		type => $type
	};
	return bless $self, $class;
}

sub _add_features {
    my ($chr, $type_id, $dsid, $hits, $dbh) = @_;
    my $query = 'SELECT chromosome,start,stop FROM feature WHERE dataset_id=' . $dsid;
    if ($chr) {
    	$query .= " AND chromosome='" . $chr . "'";
    }
    $query .= ' AND feature_type_id=' . $type_id;
    $query .= ' ORDER BY chromosome,start';
    my $sth = $dbh->prepare($query);
    $sth->execute();
    while (my @row = $sth->fetchrow_array) {
    	push @$hits, \@row;
    	# my $chromosome = $hits->{$row[0]};
    	# if (!$chromosome) {
    	# 	$hits->{$row[0]} = [\@row];
    	# }
    	# else {
	    # 	push @$chromosome, \@row;
    	# }
    }
}

sub get_features {
	my $self = shift;
	my $dbh = shift;
#	my $dir = dir($self->{conf}{CACHEDIR}, $self->{gid}, 'features');
#	my $file = file($dir, $self->{chromosome} . '_' . $self->{type} . '.loc');
#	if (-e $file) {
#		return _read_index($file);
#	}
	my $features = [];
    my $ids = get_dataset_ids($self->{gid}, $dbh);
    foreach my $dsid (@$ids) {
        _add_features $self->{chromosome}, $self->{type}, $dsid, $features, $dbh;
    }
#	_write_index($file, $features);
	return $features;
}

sub _read_index {
	my $file = shift;
	my $features = [];
	open(my $fh, '<' . $file);
	binmode($fh);
	my $buf;
	read($fh, $buf, 4);
	my $size = 0 + $buf;
	for (my $i=0; $i<$size; $i++) {
		read($fh, $buf, 4);
		my @loc = unpack('LL', $buf);
		push(@$features, \@loc);
	}
	close($fh);
	return $features;
}

sub _write_index {
	my $file = shift;
	my $features = shift;
	open(my $fh, '>' . $file);
	binmode($fh);
	print $fh pack('L', scalar @$features);
	foreach my $loc (@$features) {
		print $fh pack('LL', $loc->[1], $loc->[2]);
	}
	close($fh);
}

1;
