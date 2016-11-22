package CoGe::Services::API::JBrowse::Search;

use Mojo::Base 'Mojolicious::Controller';

use CoGe::Core::Experiment;
use CoGe::Core::Storage qw($DATA_TYPE_QUANT $DATA_TYPE_POLY $DATA_TYPE_ALIGN $DATA_TYPE_MARKER);
use CoGe::Services::Auth;
use CoGeDBI qw( get_dataset_ids feature_type_names_to_id );

# returns a psuedo object with coderefs next() and line()
sub get_db_data {
	my ($gid, $type_names, $chr, $dbh) = @_;

	my $query = 'SELECT start,stop';
	$query .= ',chromosome' unless $chr;
	$query .= ' FROM feature WHERE dataset_id';
	my $ids = get_dataset_ids($gid, $dbh);
	$query .= (scalar @$ids == 1) ? '=' . $ids->[0] : ' IN(' . join(',', @$ids) . ')';
	$query .= " AND chromosome='" . $chr . "'" if $chr;
	if ($type_names && $type_names ne 'all') {
		my $type_ids = feature_type_names_to_id($type_names, $dbh);
		$query .= ' AND feature_type_id';
		$query .= (index($type_ids, ',') == -1) ? '=' . $type_ids : ' IN(' . $type_ids . ')';
	}
	$query .= ' ORDER BY ';
	$query .= 'chromosome,' unless $chr;
	$query .= 'start,stop';
	my $sth = $dbh->prepare($query);
	$sth->execute();
	return {
		sth => $sth,
		line => sub {
				my $self = shift;
				return $self->{row};
			},
		next => sub {
			my $self = shift;
			my $row = $self->{sth}->fetchrow_arrayref;
			return undef unless $row;
			$self->{row} = $row;
			return $row->[0], $row->[1];
		}
	};
}

sub get_data {
	my ($self, $num, $db) = @_;
	my $type = $self->param('type' . $num);
	if ($type eq 'experiment') {
		my $experiment = $db->resultset('Experiment')->find($self->param('eid' . $num));
		return get_experiment_data($experiment->id, $experiment->data_type, $self->param('chr'), 0);
	}
	if ($type eq 'features') {
		return get_db_data($self->param('gid' . $num), $self->param('features' . $num), $self->param('chr'), $db->storage->dbh);
	}
}

# returns a psuedo object with coderefs next() and line()
sub get_experiment_data {
	my ($eid, $data_type, $chr, $all) = @_;

	return {
		all => $all,
		data => CoGe::Core::Experiment::query_data(
				eid => $eid,
				data_type => $data_type,
				chr => $chr,
				all => $all
			),
		data_type => $data_type,
		index => -1,
		line => sub {
				my $self = shift;
				return $self->{data}[$self->{index}];
			},
		next => sub {
				my $self = shift;
				$self->{index}++;
				return undef if $self->{index} == scalar @{$self->{data}};
				my $data_point = $self->{data}[$self->{index}];
				if ($self->{data_type} == $DATA_TYPE_ALIGN && $self->{all}) { # sam format
					my (undef, undef, undef, $pos, undef, undef, undef, undef, undef, $seq, undef, undef) = split(/\t/, $data_point);
					$pos = int($pos);
					return ($pos, $pos + length($seq));
				}
				my $c1 = index($data_point, ',') + 1;
				my $c2 = index($data_point, ',', $c1) + 1;
				my $c3 = index($data_point, ',', $c2);
				return (int(substr($data_point, $c1, $c2 - $c1 - 1)), int(substr($data_point, $c2, $c3 - $c2)));
			}
	};
}

# pass in two psuedo objects. data points in the first one that overlap those in the second are returned in an array
sub _in {
	my $data1 = shift;
	my $data2 = shift;
	my $hits = [];
	my ($start1, $stop1) = $data1->{next}->($data1);
	my ($start2, $stop2) = $data2->{next}->($data2);
	while ($start1 && $start2) {
		while ($start1 && $stop1 < $start2) {
			($start1, $stop1) = $data1->{next}->($data1);
		}
		last unless $start1;
		while ($start2 && $stop2 < $start1) {
			($start2, $stop2) = $data2->{next}->($data2);
		}
		last unless $start2;
		if ($stop1 >= $start2) {
			push @$hits, $data1->{line}->($data1);
		}
		($start1, $stop1) = $data1->{next}->($data1);
	}
	return $hits;
}

# pass in two psuedo objects. data points in the first one that don't overlap those in the second are returned in an array
sub _not_in {
	my $data1 = shift;
	my $data2 = shift;
	my $hits = [];
	my ($start1, $stop1) = $data1->{next}->($data1);
	my ($start2, $stop2) = $data2->{next}->($data2);
	while ($start1 && $start2) {
		while ($start1 && $stop1 < $start2) {
			push @$hits, $data1->{line}->($data1);
			($start1, $stop1) = $data1->{next}->($data1);
		}
		last unless $start1;
		while ($start2 && $stop2 < $start1) {
			($start2, $stop2) = $data2->{next}->($data2);
		}
		last unless $start2;
		if ($stop1 < $start2) {
			push @$hits, $data1->{line}->($data1);
		}
		($start1, $stop1) = $data1->{next}->($data1);
	}
	return $hits;
}

sub overlaps {
	my $self = shift;
	my $not = $self->param('not');

	my ($db, $user, $conf) = CoGe::Services::Auth::init($self);
	my $data1 = $self->get_data(1, $db);
	my $data2 = $self->get_data(2, $db);
	my $hits = $not ? _not_in($data1, $data2) : _in($data1, $data2);
	if ($hits) {
		$self->render(json => $hits);
	}
}

1;
