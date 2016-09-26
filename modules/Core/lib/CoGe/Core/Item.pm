package CoGe::Core::Item;

use strict;
use warnings;

BEGIN {
    our ( @EXPORT, @ISA, $VERSION );
    require Exporter;

    $VERSION = 0.1;
    @ISA = qw( Exporter );
    @EXPORT = qw( get_item );
}

my $node_types = CoGeX::node_types();

sub get_item {
    my ($db, $id, $type) = @_;
    return unless ($db and $id and $type);
    
    my $item;
    if ($type eq 'genome') {
        $item = $db->resultset('Genome')->find($id);
    }
    elsif ($type eq 'experiment') {
        $item = $db->resultset('Experiment')->find($id);
    }
    elsif ($type eq 'notebook') {
        $item = $db->resultset('List')->find($id);
    }
    else {
        warn "CoGe::Core::Item: unknown type";
        return;    
    }
    
    return $item;
}

1;
