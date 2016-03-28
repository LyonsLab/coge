package CoGe::Services::API::JBrowse::Genome;

use Mojo::Base 'Mojolicious::Controller';
use CoGeX;
use CoGe::Services::Auth qw(init);
use Data::Dumper;

sub _add_features {
    my ($name, $chr, $type_ids, $dsid, $hits, $dbh) = @_;
    #TODO move this query to CoGeDBI.pm
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
    my $name = scalar $self->param('name');
    my $chr = $self->param('chr');

    my ($db, $user) = CoGe::Accessory::Web->init;
    my $dbh = $db->storage->dbh;

	my $types = $self->param('features');
	my $type_ids;
	if ($types ne 'all') {
		$type_ids = join(',', @{$dbh->selectcol_arrayref('SELECT feature_type_id FROM feature_type WHERE name IN(' . $types . ')')});
	}
	
    my $hits = [];
    #TODO move this query to CoGeDBI.pm
    my $ids = $dbh->selectcol_arrayref('SELECT dataset_id FROM dataset_connector WHERE genome_id=' . $self->stash('gid'));
    foreach my $dsid (@$ids) {
        _add_features '%' . lc($name) . '%', $chr, $type_ids, $dsid, $hits, $dbh;
    }
    my @sorted = sort { $a->{name} cmp $b->{name} } @{$hits};
    $self->render(json => \@sorted);
}

sub genes {
    my $self = shift;
    my $name = scalar $self->param('equals');
    if (!$name) {
        $name = scalar $self->param('startswith');
        if ($name) {
            $name .= '%';
        }
    }
    
    if (index($name, ':') != -1 && index($name, '..') != -1) {
        $self->render(json => []);
        return;
    }

    my ( $db, $user ) = CoGe::Accessory::Web->init;
    my $dbh = $db->storage->dbh;

    my $hits = [];
    #TODO move this query to CoGeDBI.pm
    my $ids = $dbh->selectcol_arrayref('SELECT dataset_id FROM dataset_connector WHERE genome_id=' . $self->stash('gid'));
    foreach my $dsid (@$ids) {
        _add_features lc($name), undef, '1', $dsid, $hits, $dbh;
    }
    $self->render(json => $hits);
}

1;
