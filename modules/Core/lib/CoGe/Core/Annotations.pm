package CoGe::Core::Annotations;

use v5.14;
use strict;
use warnings;

BEGIN {
    our (@ISA, $VERSION, @EXPORT);
    require Exporter;

    $VERSION = 0.0.1;
    @ISA = qw(Exporter);
    @EXPORT = qw( get_annotation get_annotations get_type_groups );
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
    return {
        annotation => $annotation->annotation,
        link       => $annotation->link,
        type       => $type,
        type_group => $type_group
    };
}

sub get_annotations {
    my ($id, $object_type, $db) = @_;
    return unless $id && $object_type && $db;

    my $object = $db->resultset($object_type)->find($id);
    return unless $object;

    # Categorize annotations based on type group and type
    my %groups;
    my $num_annot = 0;
    foreach my $a ( $object->annotations ) {
        my $group = ( $a->type->group ? $a->type->group->name : '');
        my $type = $a->type->name;
        push @{ $groups{$group}{$type} }, $a if (defined $group and defined $type);
        $num_annot++;
    }

    my @groups;
    if ($num_annot) {
        foreach my $group ( sort keys %groups ) { # groups
            my %types;
            foreach my $type ( sort keys %{ $groups{$group} } ) { # types
                my @annotations;
                foreach my $a ( sort { $a->id <=> $b->id } @{ $groups{$group}{$type} } ) { # annotations
                    push @annotations, {
                        annotation => $a->info,
                        bisque_id => $a->bisque_id,
                        id => $a->id,
                        image_id => $a->image_id,
                        link => $a->link,
                        locked => $a->locked
                    };
                }
                $types{$type} = \@annotations;
            }
            push @groups, { $group => \%types };
        }
    }
    return \@groups;
}

sub get_type_groups {
    my $db = shift;
    my %unique;

    my $rs = $db->resultset('AnnotationTypeGroup');
    while ( my $atg = $rs->next ) {
        $unique{ $atg->name }++;
    }

    return [ sort keys %unique ];
}

sub _get_object {
    my ($id, $object_type, $own_or_edit, $db, $user) = @_;
    my $object = $db->resultset($object_type)->find($id);
    return 'Resource not found' unless $object;
    if ($own_or_edit) {
        return 'User not logged in' unless $user;
        return 'Access denied' unless $user->is_owner_editor(lc($object_type) => $id);
    }
    return 'Access denied' unless !$object->restricted || ($user && $user->has_access_to($object));
    return undef, $object;
}

1;
