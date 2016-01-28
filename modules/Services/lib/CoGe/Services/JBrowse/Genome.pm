package CoGe::Services::JBrowse::Genome;
use base 'CGI::Application';

use CoGeX;
use CoGe::Accessory::Web;
use JSON::XS;
use Data::Dumper;

sub setup {
    my $self = shift;
    $self->run_modes(
        'genes' => 'genes',
    );
    $self->mode_param('rm');
}

sub add_genes {
    my ($name, $dsid, $hits, $dbh) = @_;
    my $query = 'SELECT name,chromosome,start,stop FROM feature JOIN feature_name on feature.feature_id=feature_name.feature_id WHERE dataset_id=' . $dsid . " AND feature_type_id=1 AND lower(name) LIKE lower('" . $name . "')";
    my $sth = $dbh->prepare($query);
    $sth->execute();
    while (my $row = $sth->fetch) {
        push @$hits, { name => $row->[0], location => { ref => $row->[1], start => $row->[2], end => $row->[3] } };
    }
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
        add_genes $name, $dsid, $hits, $dbh;
    }
    return encode_json($hits);
}

1;
