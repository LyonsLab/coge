package CoGe::Core::Annotations;

use v5.14;
use strict;
use warnings;

use JSON::XS qw(encode_json);

BEGIN {
    our (@ISA, $VERSION, @EXPORT);
    require Exporter;

    $VERSION = 0.0.1;
    @ISA = qw(Exporter);
    @EXPORT = qw( get_annotation );
}

sub get_annotation {
    my ($aid, $object_type, $db) = @_;
    return unless $aid && $object_type && $db;

    #TODO check user access here

    my $annotation = $db->resultset($object_type . 'Annotation')->find($aid);
    return unless $annotation;

    my $type       = '';
    my $type_group = '';
    if ( $annotation->type ) {
        $type = $annotation->type->name;
        $type_group = $annotation->type->group->name if ( $annotation->type->group );
    }
    return encode_json({
        annotation => $annotation->annotation,
        link       => $annotation->link,
        type       => $type,
        type_group => $type_group
    });
}

1;
