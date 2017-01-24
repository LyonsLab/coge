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
    my ($id, $object_type, $db, $user) = @_;
    return unless $id && $object_type && $db;

    my $object = $db->resultset($object_type)->find($id);
    return unless $object;

    my $user_can_edit = $object->is_editable($user);

    # Categorize annotations based on type group and type
    my %groups;
    my $num_annot = 0;
    foreach my $a ( $exp->annotations ) {
        my $group = ( $a->type->group ? $a->type->group->name : '');
        my $type = $a->type->name;
        push @{ $groups{$group}{$type} }, $a if (defined $group and defined $type);
        $num_annot++;
    }

    # Build annotation table
    my $html;
    if ($num_annot) {
        $html .= '<table id="experiment_annotation_table" class="border-top border-bottom small" style="max-width:800px;overflow:hidden;word-wrap:break-word;border-spacing:0;"><thead style="display:none"></thead><tbody>';
        foreach my $group ( sort keys %groups ) { # groups
            my $first_group = 1;
            foreach my $type ( sort keys %{ $groups{$group} } ) { # types
                my $first_type = 1;
                foreach my $a ( sort { $a->id <=> $b->id } @{ $groups{$group}{$type} } ) { # annotations
                    my $header = ($group and $first_group-- > 0 ? "<b>$group</b> " : '') . ($first_type-- > 0 ? $type : '');
                    $html .= "<tr style='vertical-align:top;'>";
                    $html .= "<th class='title4' style='padding-right:10px;white-space:nowrap;font-weight:normal;background-color:white;text-align:left;'>$header</th><td class='data4'>";
                    if ($a->image) {
                        my $image_link = ( $a->image ? 'image.pl?id=' . $a->image->id : '' );
                        $html .= "<a href='$image_link' target='_blank' title='click for full-size image'><img height='40' width='40' src='$image_link' onmouseover='image_preview(this, 1);' onmouseout='image_preview(this, 0);' style='float:left;padding:1px;border:1px solid lightgray margin-right:5px;'></a>";
                    }
                    elsif ($a->bisque_id) {
                        $html .= "<a href='http://bisque.iplantc.org/client_service/view?resource=http://bisque.iplantc.org/data_service/";
                        $html .= $a->bisque_id;
                        $html .= "5' target='_blank' title='click to view in BisQue'><img src='http://bisque.iplantc.org/image_service/";
                        $html .= $a->bisque_id;
                        $html .= "5?thumbnail=200,200' onmouseover='image_preview(this, 1);' onmouseout='image_preview(this, 0);' style='float:left;padding:1px;border:1px solid lightgray;width:42px;margin-right:5px;'></a>";
                    }
                    warn $a->bisque_id;
                    $html .= $a->info;
                    $html .= '</td><td style="padding-left:5px;">';
                    $html .= linkify( $a->link, 'Link' ) if $a->link;
                    $html .= '</td>';
                    if ($user_can_edit && !$a->locked) {
                        my $aid = $a->id;
                        $html .=
                            '<td style="padding-left:20px;white-space:nowrap;">'
                          . "<span onClick=\"edit_annotation_dialog($aid);\" class='link ui-icon ui-icon-gear'></span>"
                          . "<span onClick=\"\$(this).fadeOut(); remove_annotation($aid);\" class='link ui-icon ui-icon-trash'></span>"
                          . '</td>';
                    }
                    $html .= '</tr>';
                }
            }
        }
        $html .= '</tbody></table>';
    }
    elsif ($user_can_edit) {
        $html .= '<table class="border-top border-bottom small padded note"><tr><td>There are no additional metadata items for this experiment.</tr></td></table>';
    }

    if ($user_can_edit) {
        $html .= qq{<div class="padded"><span onClick="add_annotation_dialog();" class='coge-button'>Add</span></div>};
    }
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
