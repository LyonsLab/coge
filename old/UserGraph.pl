#! /usr/bin/perl -w

use strict;
use CGI;
use JSON::XS;
use HTML::Template;
use CoGeX;
use CoGe::Accessory::Web;
use Benchmark;

use vars
  qw($P $PAGE_TITLE $USER $coge %FUNCTION $FORM %ITEM_TYPE $MAX_SEARCH_RESULTS
     $node_types);

$PAGE_TITLE = 'UserGraph';

$FORM = new CGI;
( $coge, $USER, $P ) = CoGe::Accessory::Web->init(
    cgi => $FORM,
    page_title => $PAGE_TITLE
);

$node_types = CoGeX::node_types();

%FUNCTION = (
    get_user_nodes => \&get_user_nodes,
    get_all_nodes  => \&get_all_nodes,
);

CoGe::Accessory::Web->dispatch( $FORM, \%FUNCTION, \&gen_html );

sub gen_html {
    my $template =
      HTML::Template->new( filename => $P->{TMPLDIR} . 'generic_page.tmpl' );
    $template->param( HELP => "/wiki/index.php?title=$PAGE_TITLE" );
    my $name = $USER->user_name;
    $name = $USER->first_name if $USER->first_name;
    $name .= " " . $USER->last_name if $USER->first_name && $USER->last_name;
    $template->param(
        USER       => $name,
        PAGE_TITLE => $PAGE_TITLE,
    );
    $template->param( LOGON => 1 ) unless $USER->user_name eq "public";
    $template->param( BODY => gen_body() );
    $template->param( ADMIN_ONLY => $USER->is_admin );
    $template->param( CAS_URL    => $P->{CAS_URL} || '' );

    return $template->output;
}

sub gen_body {
    return "Permission denied" unless $USER->is_admin;

    my $template =
      HTML::Template->new( filename => $P->{TMPLDIR} . "$PAGE_TITLE.tmpl" );
    $template->param(
        PAGE_NAME => "$PAGE_TITLE.pl",
        MAIN      => 1
    );

    return $template->output;
}

sub get_user_nodes {
    my @users;

    foreach ( $coge->resultset('User')->all ) {
        push @users, { name => 'user_' . $_->id, size => 100 };
    }

    return encode_json( { name => 'root', children => \@users } );
}

sub get_all_nodes {
    print STDERR "get_all_nodes\n";

    my %childrenByList;
    foreach my $conn ( $coge->resultset('ListConnector')->all ) {
        push @{ $childrenByList{ $conn->parent_id } },
          { name => $conn->child_id, type => $conn->child_type };
    }

    my %childrenByUser;
    foreach my $conn (
        $coge->resultset('UserConnector')->search(
            {
                parent_type => $node_types->{user},
                -or         => [
                    child_type => $node_types->{genome},
                    child_type => $node_types->{experiment},
                    child_type => $node_types->{list}
                ]
            }
        )
      )
    {
        if ( $conn->child_type == $node_types->{list} ) {
            push @{ $childrenByUser{ $conn->parent_id } },
              {
                name     => $conn->child_id,
                type     => $conn->child_type,
                children => $childrenByList{ $conn->child_id }
              };
        }
        else {
            push @{ $childrenByUser{ $conn->parent_id } },
              { name => $conn->child_id, type => $conn->child_type };
        }
    }

    my @users;
    foreach my $user ( $coge->resultset('User')->all ) {
        push @users,
          {
            name     => $user->id,
            type     => 5,
            info     => $user->info,
            children => $childrenByUser{ $user->id }
          };
    }

    my %childrenByGroup;
    foreach my $conn (
        $coge->resultset('UserConnector')->search(
            {
                parent_type => $node_types->{group},
                -or         => [
                    child_type => $node_types->{genome},
                    child_type => $node_types->{experiment},
                    child_type => $node_types->{list}
                ]
            }
        )
      )
    {
        if ( $conn->child_type == $node_types->{list} ) {
            push @{ $childrenByGroup{ $conn->parent_id } },
              {
                name     => $conn->child_id,
                type     => $conn->child_type,
                children => $childrenByList{ $conn->child_id }
              };
        }
        else {
            push @{ $childrenByGroup{ $conn->parent_id } },
              { name => $conn->child_id, type => $conn->child_type };
        }
    }

    my @groups;
    foreach my $group ( $coge->resultset('UserGroup')->all ) {
        my $num_users = @{ $group->users };
        push @groups,
          {
            name     => $group->id,
            type     => 6,
            info     => $group->info,
            size     => $num_users * 1000,
            children => $childrenByGroup{ $group->id }
          };
    }

    return encode_json(
        {
            name     => 'root',
            children => [
                { name => 'users',  children => \@users },
                { name => 'groups', children => \@groups }
            ]
        }
    );
}

sub get_type_name {
    my $type_id = shift;
    foreach my $type_name ( keys %{$node_types} ) {
        return $type_name if ( $type_id == $node_types->{$type_name} );
    }
}

sub get_node {
    my %opts    = @_;
    my $node_id = $opts{id};
    return unless $node_id;

    my ( $type, $id ) = split( '_', $node_id );

    my ( $url, @children );

    if ( $type eq 'user' ) {
        foreach my $conn (
            $coge->resultset('UserConnector')->search(
                { parent_id => $id, parent_type => $node_types->{user} }
            )
          )
        {
            my $child = $conn->child;
            my $type  = get_type_name( $conn->child_type );
            my $name  = ( $child->name ? $child->name : '<unnamed>' );
            my $info  = ( $child->info ? $child->info : '<no info>' );
            my $id    = $type . '_' . $conn->child_id;
            push @children,
              { id => $id, name => $name, info => $info, type => $type };
        }
    }
    elsif ( $type eq 'group' ) {
        $url = 'GroupView.pl?ugid=' . $id;
    }
    elsif ( $type eq 'list' ) {

        #$url = 'NotebookView.pl?nid='.$id;
        foreach my $conn (
            $coge->resultset('ListConnector')->search( { parent_id => $id } ) )
        {
            my $child = $conn->child;
            my $type  = get_type_name( $conn->child_type );
            my $name  = ( $child->name ? $child->name : '<unnamed>' );
            my $info  = ( $child->info ? $child->info : '<no info>' );
            my $id    = $type . '_' . $conn->child_id;
            push @children,
              { id => $id, name => $name, info => $info, type => $type };
        }
    }
    elsif ( $type eq 'genome' ) {
        $url = 'GenomeInfo.pl?gid=' . $id;
    }
    elsif ( $type eq 'experiment' ) {
        $url = 'ExperimentView.pl?eid=' . $id;
    }

    push @children, { info => '<empty>' } if ( not @children );

    return encode_json( { url => $url, children => \@children } );
}
