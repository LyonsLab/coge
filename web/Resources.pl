#! /usr/bin/perl -w

use strict;
use CGI;

use JSON::XS;
use HTML::Template;
use Digest::MD5 qw(md5_base64);
use Sort::Versions;
use List::Util qw(first);
use DBIxProfiler;
use URI::Escape::JavaScript qw(escape unescape);
use Data::Dumper;
use File::Path;
use File::stat;
use CoGeX;
use CoGe::Accessory::LogUser;
use CoGe::Accessory::Web;
use CoGeX::ResultSet::Experiment;
use CoGeX::ResultSet::Genome;
use CoGeX::ResultSet::Feature;
use Benchmark;
no warnings 'redefine';

use vars qw($P $DBNAME $DBHOST $DBPORT $DBUSER $DBPASS $connstr $PAGE_TITLE
  $TEMPDIR $USER $DATE $BASEFILE $coge $cogeweb %FUNCTION
  $COOKIE_NAME $FORM $URL $COGEDIR $TEMPDIR $TEMPURL %ITEM_TYPE
  $MAX_SEARCH_RESULTS);
$P = CoGe::Accessory::Web::get_defaults();

$DATE = sprintf(
    "%04d-%02d-%02d %02d:%02d:%02d",
    sub { ( $_[5] + 1900, $_[4] + 1, $_[3] ), $_[2], $_[1], $_[0] }
      ->(localtime)
);

$PAGE_TITLE = 'Resources';

$FORM = new CGI;

$DBNAME = $P->{DBNAME};
$DBHOST = $P->{DBHOST};
$DBPORT = $P->{DBPORT};
$DBUSER = $P->{DBUSER};
$DBPASS = $P->{DBPASS};
$connstr =
  "dbi:mysql:dbname=" . $DBNAME . ";host=" . $DBHOST . ";port=" . $DBPORT;
$coge = CoGeX->connect( $connstr, $DBUSER, $DBPASS );

#$coge->storage->debugobj(new DBIxProfiler());
#$coge->storage->debug(1);

$COOKIE_NAME = $P->{COOKIE_NAME};
$URL         = $P->{URL};
$COGEDIR     = $P->{COGEDIR};
$TEMPDIR     = $P->{TEMPDIR} . "$PAGE_TITLE/";
mkpath( $TEMPDIR, 0, 0777 ) unless -d $TEMPDIR;
$TEMPURL = $P->{TEMPURL} . "PAGE_TITLE/";

my ($cas_ticket) = $FORM->param('ticket');
$USER = undef;
($USER) = CoGe::Accessory::Web->login_cas(
    cookie_name => $COOKIE_NAME,
    ticket      => $cas_ticket,
    coge        => $coge,
    this_url    => $FORM->url()
) if ($cas_ticket);
($USER) = CoGe::Accessory::LogUser->get_user(
    cookie_name => $COOKIE_NAME,
    coge        => $coge
) unless $USER;
my $link = "http://" . $ENV{SERVER_NAME} . $ENV{REQUEST_URI};
$link = CoGe::Accessory::Web::get_tiny_link(
    db              => $coge,
    user_id         => $USER->id,
    page            => "$PAGE_TITLE.pl",
    url             => $link,
    disable_logging => 1
);

%FUNCTION = (
    gen_html => \&gen_html,
    get_node => \&get_node,
);

my $node_types = CoGeX::node_types();

dispatch();

sub dispatch {
    my %args  = $FORM->Vars;
    my $fname = $args{'fname'};
    if ($fname) {
        die if not defined $FUNCTION{$fname};

        #print STDERR Dumper \%args;
        if ( $args{args} ) {
            my @args_list = split( /,/, $args{args} );
            print $FORM->header, $FUNCTION{$fname}->(@args_list);
        }
        else {
            print $FORM->header, $FUNCTION{$fname}->(%args);
        }
    }
    else {
        print $FORM->header, gen_html();
    }
}

sub gen_html {
    my $template =
      HTML::Template->new( filename => $P->{TMPLDIR} . 'generic_page.tmpl' );
    $template->param( HELP => "/wiki/index.php?title=$PAGE_TITLE" );
    my $name = $USER->user_name;
    $name = $USER->first_name if $USER->first_name;
    $name .= " " . $USER->last_name if $USER->first_name && $USER->last_name;
    $template->param( USER => $name );

    #$template->param( TITLE      => 'User Profile' );
    $template->param( PAGE_TITLE => $PAGE_TITLE );
    $template->param( LOGO_PNG   => "$PAGE_TITLE-logo.png" );
    $template->param( LOGON      => 1 ) unless $USER->user_name eq "public";
    $template->param( DATE       => $DATE );
    $template->param( BODY       => gen_body() );
    $template->param( ADJUST_BOX => 1 );

    return $template->output;
}

sub gen_body {
    return "Top Secret!  Admin login required." unless $USER->is_admin;

    my $template =
      HTML::Template->new( filename => $P->{TMPLDIR} . "$PAGE_TITLE.tmpl" );
    $template->param(
        PAGE_NAME      => "$PAGE_TITLE.pl",
        MAIN           => 1,
        USER_NODE_LOOP => get_user_nodes()
    );

    return $template->output;
}

sub get_user_nodes {
    my @users;

    foreach ( sort { $a->info cmp $b->info } $coge->resultset('User')->all ) {
        push @users,
          { USER_NODE_ID => 'user_' . $_->id, USER_NODE_NAME => $_->info };
    }

    return \@users;
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
