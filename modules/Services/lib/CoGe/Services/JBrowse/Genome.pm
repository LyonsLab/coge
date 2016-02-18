package CoGe::Services::JBrowse::Genome;
use base 'CGI::Application';

use CoGeX;
use CoGe::Accessory::Web;
use JSON::XS;
use Data::Dumper;

sub setup {
    my $self = shift;
    $self->run_modes(
        'features' => 'features',
        'genes' => 'genes',
    );
    $self->mode_param('rm');
}

sub add_features {
    my ($name, $chr, $type_ids, $dsid, $hits, $dbh) = @_;
    my $query = 'SELECT name,chromosome,start,stop FROM feature JOIN feature_name on feature.feature_id=feature_name.feature_id WHERE dataset_id=' . $dsid;
    if ($chr) {
    	$query .= " AND chromosome='" . $chr . "'";
    }
    if ($type_ids) {
    	if (index($type_ids, ',') != -1) {
		    $query .= ' AND feature_type_id IN(' . $type_ids . ')';
    	}
    	else {
    		$query .= ' AND feature_type_id=' . $type_ids;
    	}
    }
    $query .= ' AND lower(name) ' . (index($name, '%') != -1 ? 'LIKE' : '=') . " '" . $name . "'";
    my $sth = $dbh->prepare($query);
    $sth->execute();
    while (my $row = $sth->fetch) {
        push @$hits, { name => $row->[0], location => { ref => $row->[1], start => $row->[2], end => $row->[3] } };
    }
}

sub features {
    my $self = shift;
    my $name = scalar $self->query->param('name');
    my $chr = $self->query->param('chr');

    my ( $db, $user ) = CoGe::Accessory::Web->init;
    my $dbh = $db->storage->dbh;

	my $types = $self->query->param('features');
	my $type_ids;
	if ($types ne 'all') {
		$type_ids = join(',', @{$dbh->selectcol_arrayref('SELECT feature_type_id FROM feature_type WHERE name IN(' . $types . ')')});
	}
	
    my $hits = [];
    my $ids = $dbh->selectcol_arrayref('SELECT dataset_id FROM dataset_connector WHERE genome_id=' . $self->param('gid'));
    foreach my $dsid (@$ids) {
        add_features '%' . lc($name) . '%', $chr, $type_ids, $dsid, $hits, $dbh;
    }
    my @sorted = sort { $a->{name} cmp $b->{name} } @{$hits};
    return encode_json(\@sorted);
}

sub genes {
    my $self = shift;
    my $name = scalar $self->query->param('equals');
    if (!$name) {
        $name = scalar $self->query->param('startswith');
        if ($name) {
            $name .= '%';
        }
    }
    
    if (index($name, ':') != -1 && index($name, '..') != -1) {
        return '[]';
    }

    my ( $db, $user ) = CoGe::Accessory::Web->init;
    my $dbh = $db->storage->dbh;

    my $hits = [];
    my $ids = $dbh->selectcol_arrayref('SELECT dataset_id FROM dataset_connector WHERE genome_id=' . $self->param('gid'));
    foreach my $dsid (@$ids) {
        add_features lc($name), undef, '1', $dsid, $hits, $dbh;
    }
    return encode_json($hits);
}

1;
