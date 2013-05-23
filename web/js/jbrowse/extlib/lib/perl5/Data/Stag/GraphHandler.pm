package Data::Stag::GraphHandler;

=head1 NAME

  Data::Stag::GraphHandler

=head1 SYNOPSIS

  my $p = Data::Stag->parser;
  $p->handler('graph');
  $p->handler->mapping([vertex=>[

  $p->parse

=cut

=head1 DESCRIPTION


=head1 PUBLIC METHODS -

=cut

use strict;
use base qw(Data::Stag::BaseHandler);

use vars qw($VERSION);
$VERSION="0.11";
use Graph;

sub init {
    my $self = shift;
    $self->SUPER::init(@_);
    $self->graph(new Graph);
}

sub graph {
    my $self = shift;
    $self->{graph} = shift if @_;
    return $self->{graph};
}

sub mapping {
    my $self = shift;
    $self->{_mapping} = shift if @_;
    return $self->{_mapping};
}


sub graph_edges {
    my $self = shift;
    $self->{_graph_edges} = shift if @_;
    return $self->{_graph_edges};
}

sub catch_end {
    my $self = shift;
    my $ev = shift;
    my $node = shift;

    my $g = $self->graph;

    my $map = $self->mapping;
    if ($map) {
        my $uid_tag = $map->get_uid;
        my @vertexmaps = $map->getnode_vertex;
        my @edgemaps = $map->getnode_edge;
        my ($emap) = grep {$ev eq $_->get_tag} @edgemaps;
        if ($emap) {
            my $parent_tag = $emap->get_parent;
            my $child_tag = $emap->get_child;
            
            my $parent = $node->get($parent_tag);
            my $child = $node->get($child_tag);
            if (!@vertexmaps) {
                # graph defined purely by edges
                $g->add_vertex($parent);
                $g->add_vertex($child);
            }
            $g->add_edge($parent, $child);
            $self->_set_attributes($g, $node, [$parent, $child], $emap->get_att);
        }
        my ($vmap) = grep {$ev eq $_->get_tag} @vertexmaps;
        if ($vmap) {
            my $uid = $node->get($uid_tag);
            $g->add_vertex($uid);
            $self->_set_attributes($g, $node, [$uid], $vmap->get_att);
        }
    }
    return;
}

sub _set_attributes {
    my $self = shift;
    my ($g, $node, $idref, @atts) = @_;

    if (!scalar(@atts) ||
        grep {$_ eq '*'} @atts) {
        my %atth = map { $_->name => 1 } $node->kids;
        @atts = keys %atth;
    }
    foreach my $att (@atts) {
        my @v = $node->get($att);
        foreach my $v (@v) {
            $g->set_attribute($att, @$idref, $v);
        }
    }
    return;
}

1;
