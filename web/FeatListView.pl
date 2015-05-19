#! /usr/bin/perl -w
use strict;
use CGI;
use CGI::Carp 'fatalsToBrowser';
use CoGe::Accessory::LogUser;
use HTML::Template;
use Data::Dumper;
use CGI::Ajax;
use Digest::MD5 qw(md5_base64);
use CoGe::Graphics;
use File::Temp;
no warnings 'redefine';

$ENV{PATH} = "/opt/apache/CoGe/";

use vars
  qw( $DBNAME $DBHOST $DBPORT $DBUSER $DBPASS $connstr $coge $DATE $DEBUG $TEMPDIR $TEMPURL $USER $FORM $DB);

# set this to 1 to print verbose messages to logs
$DEBUG   = 0;
$TEMPDIR = "/opt/apache/CoGe/tmp";
$TEMPURL = "/CoGe/tmp";
$|       = 1;                        # turn off buffering
$DATE    = sprintf(
    "%04d-%02d-%02d %02d:%02d:%02d",
    sub { ( $_[5] + 1900, $_[4] + 1, $_[3] ), $_[2], $_[1], $_[0] }
      ->(localtime)
);

$FORM = new CGI;

$DBNAME = $P->{DBNAME};
$DBHOST = $P->{DBHOST};
$DBPORT = $P->{DBPORT};
$DBUSER = $P->{DBUSER};
$DBPASS = $P->{DBPASS};
$connstr =
  "dbi:$P->{DB}:dbname=" . $DBNAME . ";host=" . $DBHOST . ";port=" . $DBPORT;
$coge = CoGeX->connect( $connstr, $DBUSER, $DBPASS );

($USER) =
  CoGe::Accessory::LogUser->get_user( cookie_name => 'cogec', coge => $coge );

if ( $FORM->param('ticket') && $USER->user_name eq "public" ) {

    my @values = split( /'?'/, $FORM->url() );

    my ( $name, $fname, $lname, $email, $login_url ) =
      CoGe::Accessory::Web::login_cas( $FORM->param('ticket'), $values[0] );

    if ($name) {
        my ( $valid, $cookie, $urlx ) =
          login( name => $name, url => $login_url );

        if ( $valid eq 'true' ) {
            print STDERR 'valid';
        }
        else {

            my $new_row = $coge->resultset('User')->create(
                {
                    user_name  => $name,
                    first_name => $fname,
                    last_name  => $lname,
                    email      => $email
                }
            );
            $new_row->insert;
            print STDERR 'not valid';
            ( $valid, $cookie, $urlx ) =
              login( name => $name, url => $login_url );
        }

        print STDERR $cookie;
        print "Set-Cookie: $cookie\n";

    }
    $FORM->delete_all();

    ($USER) = CoGe::Accessory::LogUser->get_user(
        cookie_name => 'cogec',
        coge        => $coge
    );
    print 'Location:' . $FORM->redirect($login_url);
    print STDERR "***" . $USER->user_name;
}

my $pj = new CGI::Ajax(
    expand_list    => \&expand_list,
    close_list     => \&close_list,
    edit_list      => \&edit_list,
    delete_list    => \&delete_list,
    create_list    => \&create_list,
    delete_feature => \&delete_feature,
    add_feature    => \&add_feature,
    expand_group   => \&expand_group,
    close_group    => \&close_group,
    edit_group     => \&edit_group,
    delete_group   => \&delete_group,
    create_group   => \&create_group,
    group_view     => \&group_view,
    genomic_view   => \&genomic_view,
    add_image_box  => \&add_image_box,
    delete_image   => \&delete_image,
    edit_image     => \&edit_image,
);
$pj->JSDEBUG(0);
$pj->DEBUG(0);
print $pj->build_html( $FORM, \&gen_html );

#print "Content-Type: text/html\n\n";
#print gen_html();

sub gen_html {
    if ( $FORM->param('add_image') ) {
        add_image();
    }    #need to do this in order to deal with uploading an image file

    my ( $body, $seq_names, $seqs ) = gen_body();
    my $template =
      HTML::Template->new(
        filename => '/opt/apache/CoGe/tmpl/generic_page.tmpl' );

    $template->param( TITLE => 'Feature List Viewer' );
    my $name = $USER->user_name;
    $name = $USER->first_name if $USER->first_name;
    $name .= " " . $USER->last_name if $USER->first_name && $USER->last_name;
    $template->param( USER => $name );

    $template->param( LOGON => 1 ) unless $USER->user_name eq "public";
    $template->param( DATE => $DATE );
    $template->param( LOGO_PNG => "FeatListView-logo.png" );
    $template->param( BODY     => $body );
    $template->param( POSTBOX  => gen_foot() );
    $template->param( ADMIN_ONLY => $USER->is_admin );
    $template->param( CAS_URL    => $P->{CAS_URL} || '' );
    my $html;
    $html .= $template->output;
    return $html;
}

sub login {

    #$my $self= shift;

    my %opts = @_;
    my $name = $opts{name};
    my $url  = $opts{url};
    my ($u) = $coge->resultset('User')->search( { user_name => $name } );

    if ($u) {

        my $session = md5_base64( $name . $ENV{REMOTE_ADDR} );
        $session =~ s/\+/1/g;
        my $sid = $coge->log_user( user => $u, session => $session );

        my $c = CoGe::Accessory::LogUser->gen_cookie(
            session     => $session,
            cookie_name => 'cogec',
            url         => $url
        );

        return ( 'true', $c, $url );
    }
    else {
        my $c = CoGe::Accessory::LogUser->gen_cookie( session => "public" );
        return ( 'false', $c, $url );
    }

}

sub gen_foot {
    my $template =
      HTML::Template->new(
        filename => '/opt/apache/CoGe/tmpl/FeatListView.tmpl' );
    $template->param( FOOT => 1 );
    return $template->output;
}

sub gen_body {

    #    my $form = shift || $FORM;
    my $start = $FORM->param('start') || 1;
    my $limit = $FORM->param('ll')    || 200;
    my $html;
    if ( $FORM->param('flgid') ) {
        $html = group_view( $FORM->param('flgid') );
    }
    elsif ( $FORM->param('fl') ) {

        $html = feat_lists(
            start => $start,
            limit => $limit,
        );
    }
    elsif ( $FORM->param('flg') ) {
        $html = feat_list_groups();
    }
    else {
        $html = feat_list_groups();
    }

    #    print $html;
    return $html;
}

sub feat_lists {

    #need to make sure to do a user validation at some point
    my %opts  = @_;
    my $start = $opts{start};
    my $limit = $opts{limit};
    my $tmpl =
      HTML::Template->new(
        filename => '/opt/apache/CoGe/tmpl/FeatListView.tmpl' );

    my @lists;
    my $row_color = 'ltgreen';
    my $total     = $DB->get_feature_list_obj->count_all();
    $tmpl->param( "LIST_COUNT" => $total );
    my $stop = $start + $limit - 1;
    $stop = $total if $stop > $total;
    $tmpl->param( "LIST_START" => $start, "LIST_STOP" => $stop );
    my $pstart = $start - $limit;
    my $nstart = $start + $limit;
    my $lstart = $total - $limit + 1;
    $tmpl->param( "FIRST" =>
qq{<input type="submit" name="<<" value="<<" onClick="window.location='FeatListView.pl?fl=1&start=1'">}
    ) if ( $start > 1 );
    $tmpl->param( "PREVIOUS" =>
qq{<input type="submit" name="<" value="<" onClick="window.location='FeatListView.pl?fl=1&start=$pstart'">}
    ) if ( $start > 1 );
    $tmpl->param( "NEXT" =>
qq{<input type="submit" name=">" value=">" onClick="window.location='FeatListView.pl?fl=1&start=$nstart'">}
    ) if ( $stop < $total );
    $tmpl->param( "LAST" =>
qq{<input type="submit" name=">>" value=">>" onClick="window.location='FeatListView.pl?fl=1&start=$lstart'">}
    ) if ( $stop < $total );

    foreach my $fl (
        $DB->get_feature_list_obj->get_lists(
            start => $start,
            limit => $limit
        )
      )

      #    foreach my $fl ($DB->get_feature_list_obj->retrieve_all())
    {
        my $f_count    = scalar $fl->flc;
        my $readers    = join( ", ", map { $_->name } $fl->readers );
        my $editors    = join( ", ", map { $_->name } $fl->editors );
        my $group_name = $fl->group->name if $fl->group;
        push @lists,
          {
            LIST_NAME     => $fl->name,
            LIST_DESC     => $fl->desc,
            FEATURE_COUNT => $f_count,
            LIST_ID       => $fl->id,
            LIST_GROUP    => $group_name,
            LIST_EDITORS  => $readers,
            LIST_READERS  => $editors,
            LIST_NOTE     => $fl->notes,
            ROW_BG        => $row_color
          };
        $row_color = $row_color eq 'ltgreen' ? 'ltblue' : 'ltgreen';
    }
    $tmpl->param( FEAT_LIST  => 1 );
    $tmpl->param( FEAT_LISTS => \@lists );
    return $tmpl->output;
}

sub expand_list {
    my $lid = shift;
    my $fl  = $DB->get_feature_list_obj->retrieve($lid);
    my @feats;
    foreach my $f ( $fl->features ) {
        push @feats, {
            FEAT_TYPE => $f->type->name,
            FEAT_ORG  => $f->org->name,
            FEAT_NAME => join(
                ", ",
                sort map {
                    "<a href = FeatView.pl?accn="
                      . $_->name . ">"
                      . $_->name . "</a>"
                  } $f->names
            ),
            LOC        => gelo_link( feat => $f ),
            FEAT_PNAME => $fl->preferred_name($f),
            FEAT_DESC  => $fl->preferred_desc($f),
        };
    }
    my $tmpl =
      HTML::Template->new(
        filename => '/opt/apache/CoGe/tmpl/FeatListView.tmpl' );

    #    $tmpl->param(LIST_ID=>$fl->id);
    $tmpl->param( EXPAND_LIST  => 1 );
    $tmpl->param( FEATURE_INFO => \@feats );
    my $close = qq{
<input type="IMAGE" src="picts/open.png" onClick="close_list(['args__$lid'],['FLID$lid', 'showfeat$lid'])">
};
    return $tmpl->output, $close;
}

sub close_list {
    my $lid = shift;
    my $fl  = $DB->get_feature_list_obj->retrieve($lid);
    my $open =
qq{<input type="IMAGE" src="picts/close.png" onClick="expand_list(['args__$lid'],['FLID$lid', 'showfeat$lid'])">};
    return scalar $fl->flc, $open;
}

sub edit_list {
    my $lid = shift;
    my $fl  = $DB->get_feature_list_obj->retrieve($lid);
    my @feats;
    my $row_color = 'ltgreen';
    foreach my $f ( $fl->features ) {

        push @feats, {
            FEAT_ID   => $f->id,
            FEAT_ORG  => $f->org->name,
            FEAT_DATA => $f->dataset->name . "(v" . $f->dataset->version . ")",
            FEAT_TYPE => $f->type->name,
            LIST_ID   => $lid,
            ROW_BG    => $row_color,
            FEAT_NAME => join(
                ", ",
                sort map {
                    "<a href = FeatView.pl?accn="
                      . $_->name . ">"
                      . $_->name . "</a>"
                  } $f->names
            ),
            FEAT_PNAME => $fl->preferred_name($f),
            FEAT_LOC   => gelo_link( feat => $f ),

            #FEAT_LOC => "Chr: ".$f->chr." (".$f->start." - ".$f->stop.")",
        };
        $row_color = $row_color eq 'ltgreen' ? 'ltblue' : 'ltgreen';
    }
    @feats = sort {
        $a->{FEAT_PNAME} cmp $b->{FEAT_PNAME}
          || $a->{FEAT_NAME} cmp $b->{FEAT_NAME}
    } @feats;
    my $tmpl =
      HTML::Template->new(
        filename => '/opt/apache/CoGe/tmpl/FeatListView.tmpl' );
    $tmpl->param( EDIT_LIST    => 1 );
    $tmpl->param( LIST_NAME    => $fl->name );
    $tmpl->param( LIST_DESC    => $fl->desc );
    $tmpl->param( LIST_ID      => $fl->id );
    $tmpl->param( FEATURE_INFO => \@feats );

    return $tmpl->output;
}

sub create_list {

}

sub delete_list {
    my $lid = shift;
    my $fl  = $DB->get_feature_list_obj->retrieve($lid);
    $fl->delete if $fl;
}

sub delete_feature {
    my $fid = shift;
    my $lid = shift;
    foreach my $flc (
        $DB->get_feature_list_connector_obj->retrieve(
            feature_list_id => $lid,
            feature_id      => $fid
        )
      )
    {
        $flc->delete;
    }
    return $lid;
}

sub add_feature {
    my $id  = shift;
    my $lid = shift;
    return $lid unless $id                                        =~ /^\d+$/;
    return $lid unless ref( $DB->get_feature_obj->retrieve($id) ) =~ /Feature/;
    return $lid unless $lid                                       =~ /^\d+$/;
    $DB->get_feature_list_connector_obj->find_or_create(
        feature_list_id => $lid,
        feature_id      => $id
    );
    return $lid;
}

sub search_feature {

}

sub feat_list_groups {
    my $tmpl =
      HTML::Template->new(
        filename => '/opt/apache/CoGe/tmpl/FeatListView.tmpl' );
    my @lists;
    my $row_color = 'ltgreen';
    foreach my $flg ( $DB->get_feature_list_group_obj->retrieve_all() ) {
        my $count = scalar $flg->fl;
        push @lists,
          {
            GROUP_NAME  => $flg->name,
            GROUP_DESC  => $flg->desc,
            GROUP_COUNT => $count,
            GROUP_ID    => $flg->id,
            GROUP_NOTE  => $flg->notes,
            ROW_BG      => $row_color
          };
        $row_color = $row_color eq 'ltgreen' ? 'ltblue' : 'ltgreen';
    }
    $tmpl->param( LIST_GROUP  => 1 );
    $tmpl->param( LIST_GROUPS => \@lists );
    return $tmpl->output;
}

sub expand_group {
    my $gid = shift;
    my $flg = $DB->get_feature_list_group_obj->retrieve($gid);
    my @lists;
    foreach my $fl ( $flg->lists ) {
        my $count = scalar $fl->flc;
        push @lists,
          {
            LIST_ID    => $fl->id,
            LIST_NAME  => $fl->name,
            LIST_DESC  => $fl->desc,
            FEAT_COUNT => $count,
            LIST_ID    => $fl->id
          };
    }
    my $tmpl =
      HTML::Template->new(
        filename => '/opt/apache/CoGe/tmpl/FeatListView.tmpl' );

    #    $tmpl->param(LIST_ID=>$fl->id);
    $tmpl->param( EXPAND_GROUP => 1 );
    $tmpl->param( LIST_INFO    => \@lists );
    my $close = qq{
<input type="IMAGE" src="picts/open.png" onClick="close_group(['args__$gid'],['GID$gid', 'showlist$gid'])">
};
    return $tmpl->output, $close;
}

sub close_group {
    my $gid = shift;
    my $flg = $DB->get_feature_list_group_obj->retrieve($gid);
    my $open =
qq{<input type="IMAGE" src="picts/close.png" onClick="expand_group(['args__$gid'],['GID$gid', 'showlist$gid'])">};
    return scalar $flg->fl, $open;
}

sub edit_group {
}

sub create_group {
}

sub delete_group {
    my $gid = shift;
    my $flg = $DB->get_feature_list_group_obj->retrieve($gid);
    $flg->delete if $flg;
}

sub group_view {
    my $gid = shift;

    #    print "\n\n$gid\n\n";
    my $tmpl =
      HTML::Template->new(
        filename => '/opt/apache/CoGe/tmpl/FeatListView.tmpl' );
    $tmpl->param( GROUP_VIEW => 1 );
    my $flg = $DB->get_feature_list_group_obj->retrieve($gid);
    $tmpl->param( NAME  => $flg->name );
    $tmpl->param( DESC  => $flg->desc );
    $tmpl->param( NOTES => $flg->notes );
    my @lists;
    my $row_color = 'ltgreen';

    foreach my $fl ( $flg->lists ) {
        my $count = 0;
        foreach ( $fl->flc ) { $count++; }
        push @lists,
          {
            LIST_NAME  => $fl->name . "</a>",
            LIST_DESC  => $fl->desc,
            LIST_NOTE  => $fl->notes,
            FEAT_COUNT => $count,
            ROW_BG     => $row_color,
            LIST_ID    => $fl->id
          };
        $row_color = $row_color eq 'ltgreen' ? 'ltblue' : 'ltgreen';

    }
    $tmpl->param( LIST_LOOP => \@lists );
    my %seen;
    foreach my $feat ( $flg->features ) {
        $seen{ $feat->org->name }{ $feat->dataset->id }{ $feat->chr }{count}++;
        push @{ $seen{ $feat->org->name }{ $feat->dataset->id }{ $feat->chr }
              {feats} }, $feat;
    }
    my @orgs;
    foreach my $org ( keys %seen ) {
        foreach my $dsid ( keys %{ $seen{$org} } ) {
            foreach my $chr ( keys %{ $seen{$org}{$dsid} } ) {
                my ($start) =
                  sort { $a <=> $b }
                  map  { $_->start } @{ $seen{$org}{$dsid}{$chr}{feats} };
                my ($stop) =
                  sort { $b <=> $a }
                  map  { $_->stop } @{ $seen{$org}{$dsid}{$chr}{feats} };

                #		my $loc = "Chr: ".$chr." ($start - $stop)";
                my $loc = gelo_link(
                    chr   => $chr,
                    start => $start,
                    stop  => $stop,
                    di    => $dsid
                );
                my $action =
qq{<input type="submit" name="Genomic View" value="Genomic View" onClick="genomic_view(['args__flg', 'args__$gid', 'args__$chr', 'args__$dsid'],['genomic_view'])">};
                push @orgs,
                  {
                    ORG_NAME   => $org,
                    FEAT_COUNT => $seen{$org}{$dsid}{$chr}{count},
                    LOC        => $loc,
                    ROW_BG     => $row_color,
                    ACTIONS    => $action,
                  };
                $tmpl->param( ORG_LOOP => \@orgs );
                $row_color = $row_color eq 'ltgreen' ? 'ltblue' : 'ltgreen';
            }
        }
    }
    my @images;
    foreach my $image ( $flg->images ) {
        my $img_link = qq{<IMG SRC="image.pl?id=} . $image->id . qq{">};
        push @images,
          {
            CID    => $image->flgic->next->id,
            IMAGE  => $img_link,
            NAME   => $image->name,
            DESC   => $image->desc,
            TYPE   => 'flg',
            RELOAD => "FeatListView.pl?flgid=" . $flg->id,
          };
    }
    if (@images) {
        $tmpl->param( SHOW_IMAGE => 1 );
        $tmpl->param( IMAGE_LOOP => \@images );
    }
    $tmpl->param( ADD_IMAGE => 1 );
    $tmpl->param( TYPE      => 'flg' );
    $tmpl->param( ID        => $gid );
    my $html = $tmpl->output;

    #    print "!!!!!";
    #    print $html;
    #    print "!!!!\n";

    return $html;
}

sub add_image_box {
    my $type = shift;
    my $id   = shift;
    my $tmpl =
      HTML::Template->new(
        filename => '/opt/apache/CoGe/tmpl/FeatListView.tmpl' );
    $tmpl->param( ADD_IMG => 1 );
    $tmpl->param( TYPE    => $type );
    $tmpl->param( ID      => $id );
    my $html = $tmpl->output;
    return $html;
}

sub genomic_view {
    my $type = shift;
    my $id   = shift;
    return "No feature list or feature list group id provided" unless $id;
    my $chr  = shift;
    my $dsid = shift;
    my $list =
        $type eq "flg"
      ? $DB->get_feature_list_group_obj->retrieve($id)
      : $DB->get_feature_list_obj->retrieve($id);
    my @feats;
    foreach my $feat ( $list->features ) {
        push @feats, $feat if $feat->dataset->id == $dsid && $feat->chr eq $chr;
    }
    return generate_image(
        feats => \@feats,
        dsid  => $dsid,
        chr   => $chr,
        list  => $list
    );
}

sub add_image {
    my $type = $FORM->param('type');
    my $id   = $FORM->param('id');
    my $name = $FORM->param('name');
    my $desc = $FORM->param('desc');
    my $content;
    if ( $FORM->param('file') ) {
        my $fh = $FORM->upload('file');
        while (<$fh>) { $content .= $_ }
    }
    $name = $FORM->param('file') unless $name;
    if ( $type eq "flg" ) {
        $FORM->param( 'flgid' => $id );
        $DB->get_feature_list_group_obj->retrieve($id)
          ->insert_image( name => $name, desc => $desc, img => $content );
    }
}

sub delete_image {
    my $id   = shift;
    my $type = shift;

    #    print STDERR $id,"-",$type,"\n";
    my $obj;
    $obj = $DB->get_feature_list_group_image_connector_obj->retrieve($id)
      if $type eq "flg";
    my $img = $obj->image;
    $img->delete;
    $obj->delete;
}

sub edit_image {   #should let user change the name and decsription of the image
}

sub generate_image {
    my %opts  = @_;
    my $feats = $opts{feats};
    my $dsid  = $opts{dsid};
    my $chr   = $opts{chr};
    my $list  = $opts{list};
    my $iw    = $opts{iw} || 1600;
    my $ih    = $opts{ih} || 200;
    my $fh    = $opts{fh} || 25;
    my ($start) = sort { $a <=> $b } map { $_->start } @$feats;
    my ($stop)  = sort { $b <=> $a } map { $_->stop } @$feats;
    my $fnames = $list->preferred_names;
    my $fids   = map { "fid=" . $_->id } @$feats unless $fnames;
    my $file   = new File::Temp(
        TEMPLATE => 'FeatListView__XXXXX',
        DIR      => $TEMPDIR,
        SUFFIX   => '.png',
        UNLINK   => 0
    );
    my ($filename) = $file->filename =~ /([^\/]*$)/;
    my $mapname    = $filename . "map";
    my ($map)      = CoGe::Graphics->genomic_view(
        start   => $start,
        stop    => $stop,
        di      => $dsid,
        chr     => $chr,
        iw      => $iw,
        ih      => $ih,
        fsh     => $fh,
        fids    => $fids,
        fnames  => $fnames,
        img_map => $mapname,
        file    => $file->filename,

        #					     debug=>1,
    );
    my $img_link =
      qq{<img src="tmp/$filename" border=0 ismap usemap="#$mapname" >} . "\n"
      . $map;
    $img_link .= "<br>";
    $img_link .= qq{<FORM NAME="info"><DIV id="info"></DIV></FORM>};

    #    print STDERR $img_link,"\n";
    return $img_link;
}

sub gelo_link {
    my %opts  = @_;
    my $feat  = $opts{feat};
    my $chr   = $opts{chr};
    my $start = $opts{start};
    my $stop  = $opts{stop};
    my $di    = $opts{di};
    $chr   = $feat->chr         if $feat;
    $start = $feat->start       if $feat;
    $stop  = $feat->stop        if $feat;
    $di    = $feat->dataset->id if $feat;
    my $link =
        qq{<a href="GenomeView.pl?}
      . join( "&", "chr=" . $chr, "ds=" . $di, "z=10", "x=$start" )
      . qq{" target=_new>}
      . qq{Chr: $chr ($start - $stop)</a>};
    return $link;
}
