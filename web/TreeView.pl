#! /usr/bin/perl -w
use strict;
use CGI;
use CGI::Carp 'fatalsToBrowser';
use CGI::Ajax;
use CoGe::Accessory::LogUser;
use CoGe::Accessory::Web;
use HTML::Template;
use Data::Dumper;
use Bio::TreeIO;
use File::Temp;
use Digest::MD5 qw(md5_base64);
use File::Path;
no warnings 'redefine';

use vars
  qw($P $DBNAME $DBHOST $DBPORT $DBUSER $DBPASS $connstr $DATE $DEBUG $TEMPDIR $TEMPURL $USER $FORM $NEATO $DOT $COOKIE_NAME $coge);
$P = CoGe::Accessory::Web::get_defaults();
$ENV{PATH} = $P->{COGEDIR};

# set this to 1 to print verbose messages to logs
$DEBUG   = 0;
$TEMPDIR = $P->{TEMPDIR} . "TreeView";
$TEMPURL = $P->{TEMPURL} . "TreeView";
mkpath( $TEMPDIR, 0, 0777 ) unless -d $TEMPDIR;

$NEATO = $P->{NEATO};
$DOT   = $P->{DOT};
$|     = 1;             # turn off buffering
$DATE  = sprintf(
    "%04d-%02d-%02d %02d:%02d:%02d",
    sub { ( $_[5] + 1900, $_[4] + 1, $_[3] ), $_[2], $_[1], $_[0] }
      ->(localtime)
);

$FORM   = new CGI;
$DBNAME = $P->{DBNAME};
$DBHOST = $P->{DBHOST};
$DBPORT = $P->{DBPORT};
$DBUSER = $P->{DBUSER};
$DBPASS = $P->{DBPASS};
$connstr =
  "dbi:$P->{DB}:dbname=" . $DBNAME . ";host=" . $DBHOST . ";port=" . $DBPORT;
$coge = CoGeX->connect( $connstr, $DBUSER, $DBPASS );

$COOKIE_NAME = $P->{COOKIE_NAME};

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

my $pj = new CGI::Ajax(
    gen_html => \&gen_html,
    zoom     => \&zoom,
);
$pj->JSDEBUG(0);
$pj->DEBUG(0);

print $pj->build_html( $FORM, \&gen_html );

#print "Content-Type: text/html\n\n". gen_html();
#gen_html();

sub gen_html {
    my $html;
    my ( $body, $seq_names, $seqs ) = gen_body();
    my $template =
      HTML::Template->new( filename => $P->{TMPLDIR} . 'generic_page.tmpl' );

    $template->param( TITLE => 'Phylogenetic Tree Viewer' );

    my $name = $USER->user_name;
    $name = $USER->first_name if $USER->first_name;
    $name .= " " . $USER->last_name if $USER->first_name && $USER->last_name;
    $template->param( USER => $name );

    $template->param( DATE     => $DATE );
    $template->param( LOGO_PNG => "TreeView-logo.png" );
    $template->param( BODY     => $body );
    $template->param( LOGON    => 1 ) unless $USER->user_name eq "public";
    $template->param( ADMIN_ONLY => $USER->is_admin );
    $template->param( CAS_URL    => $P->{CAS_URL} || '' );
    $html .= $template->output;
    return $html;
}

sub gen_body {
    my $content;
    if ( $FORM->param('file') ) {
        my $fh = $FORM->upload('file');
        while (<$fh>) { $content .= $_ }
    }
    elsif ( $FORM->param('tree') ) {
        $content = $FORM->param('tree');
    }
    my $zoom = $FORM->param('zoom');

    #    return unless $content;
    $content =~ s/#.*?;//s;
    $content =~ s/\[.*?\]//sg;
    $content =~ s/end;//g;
    $content =~ s/\r/\n/g;
    $content =~ s/\n+/\n/g;
    my ($translate_table) = $content =~ /Translate(.*?);/s;
    $content =~ s/Translate.*?;//s;
    $content =
      translate_content( table => $translate_table, content => $content );

    #    print STDERR $content,"\n";

    #    $content =1;
    my $template =
      HTML::Template->new( filename => $P->{TMPLDIR} . 'TreeView.tmpl' );
    $template->param( FRONT_PAGE => 1 ) unless $content;
    if ($content) {
        my $cogeweb =
          CoGe::Accessory::Web::initialize_basefile( tempdir => $TEMPDIR );
        my $basefile = $cogeweb->basefile;
        my ( $tree, $dot, $img, $imgmap ) = generate_tree(
            content  => $content,
            basefile => $basefile,
            zoom     => $zoom
        );
        $img  =~ s/$TEMPDIR/$TEMPURL/;
        $tree =~ s/$TEMPDIR/$TEMPURL/;
        $dot  =~ s/$TEMPDIR/$TEMPURL/;
        $template->param( IMG         => $img );
        $template->param( IMGMAP      => $imgmap ) if $imgmap;
        $template->param( 'tree_file' => $tree );
        $template->param( 'dot_file'  => $dot );
        $template->param( 'tfile'     => $tree );
        my $mod_file = generate_mod_file(
            mod      => $FORM->param('nodes'),
            basefile => $basefile
        ) if $FORM->param('nodes');
        $mod_file =~ s/$TEMPDIR/$TEMPURL/;
        $template->param( 'mod' => $mod_file ) if $mod_file;
        $template->param( 'ROOTED' => $FORM->param('rooted') )
          ;    # if $FORM->param('rooted');
        $template->param( 'OVERLAP' => $FORM->param('overlap') )
          ;    # if $FORM->param('overlap');
        $template->param( 'LR' => $FORM->param('LR') ); # if $FORM->param('LR');
    }
    return $template->output;
}

sub zoom {
    my $zoom    = shift;
    my $tfile   = shift;
    my $mod     = shift;
    my $extra   = shift;
    my $rooted  = shift;
    my $overlap = shift;
    my $LR      = shift;
    $FORM->param( 'rooted',  $rooted );
    $FORM->param( 'overlap', $overlap );
    $FORM->param( 'LR',      $LR );

    $zoom += $extra if $extra;
    $mod =~ s/$TEMPURL/$TEMPDIR/;
    my $mod_content = load_mod($mod);
    $FORM->param( 'nodes', $mod_content );
    $tfile =~ s/$TEMPURL/$TEMPDIR/;
    warn "can't read $tfile" unless -f $tfile;
    my ($dot_file) = generate_dot_file( file => $tfile, zoom => $zoom );
    my ( $img, $map ) = generate_image( dot_file => $dot_file );
    my $zoomin = '<input type="hidden" id = "zoi" value=' . ( $zoom + 1 ) . '>';
    my $zoomout =
      '<input type="hidden" id = "zoo" value=' . ( $zoom - 1 ) . '>';
    return ( $img, $map, $zoomin, $zoomout );
}

sub generate_tree {
    my %opts     = @_;
    my $content  = $opts{content};
    my $basefile = $opts{basefile};
    my $zoom     = $opts{zoom};
    my ( $dot_file, $tree_file ) = generate_dot_file(
        content  => $content,
        basefile => $basefile,
        zoom     => $zoom
    );
    my ( $img_file, $map ) = generate_image( dot_file => $dot_file );
    return ( $tree_file, $dot_file, $img_file, $map );
}

sub generate_image {
    my %opts = @_;
    my $file = $opts{dot_file};
    my ($basefile) = $file =~ /^(.*)\.dot$/;
    $basefile .= rand(100);
    my $cmd = $FORM->param("rooted") ? $DOT : $NEATO;
    my $ofn = $basefile . ".png";
    system "$cmd -Tpng $file -o $ofn";
    open( IN, "$cmd -Tcmap $file |" );
    my $map;
    while (<IN>) { $map .= $_; }
    close IN;
    my $map_name = "tree";
    $ofn =~ s/$TEMPDIR/$TEMPURL/;
    return "<img src=\"$ofn\" USEMAP=\"#" . $map_name . "\" border=0>",
      "<MAP NAME=$map_name>$map</MAP>";
}

sub generate_dot_file {
    my %opts     = @_;
    my $content  = $opts{content};
    my $basefile = $opts{basefile};
    $content =~ s/end;//g;
    my $file = $opts{file} || $basefile . ".nexus_tree";
    ($basefile) = $file =~ /^(.*)\.nexus_tree/ unless $basefile;
    print STDERR $file, "\n";
    my $zoom = $opts{zoom};

    if ($content) {
        open( OUT, ">$file" );
        print OUT $content;
        close OUT;
    }
    return unless -r $file;

    my $treeio = new Bio::TreeIO(    #-format =>'newick',
        -file => $file,
    );

    #    print STDERR Dumper $treeio->next_tree;
    while ( my $tree = $treeio->next_tree() ) {
        print STDERR "Tree nodes ", $tree->number_nodes, "\n";
        next unless $tree->number_nodes;
        my $node_names = get_names( $file, $tree );
        my $root_node = $tree->get_root_node;
        next unless $root_node->internal_id;
        my $dot;
        $dot .= "graph tree {\n";
        $dot .= qq{
     rankdir=LR;
} if $FORM->param('LR');
        my $size = 8;
        $size += $zoom if $zoom;
        $size = 1 if $size < 1;
        $dot .= "\tsize=\"$size,$size\";\n";
        $dot .= "\tfontsize=10;\n";
        $dot .= "\tgraph [splines=true ";
        $dot .= "overlap=scale" unless $FORM->param('overlap');
        $dot .= "];\n";
        $dot .= process_node( $root_node, $node_names );
        $dot .= label_nodes( $tree, $node_names );

        if ( $FORM->param("rooted") ) {

            #let's make sure all named nodes are at the same rank
            my $rank = join( "; ", map { "\"" . $_ . "\"" } keys %$node_names );
            $dot .= qq{
     { rank = same $rank; }
};
        }
        $dot .= "}\n";
        my $dot_file = $basefile . ".dot";
        open( OUT, ">$dot_file" );
        print OUT $dot;
        close OUT;
        return ( $dot_file, $file );
    }
}

sub label_nodes {
    my $tree        = shift;
    my $names       = shift;
    my $node_labels = get_labels();
    my $content;
    foreach my $node ( $tree->get_nodes ) {

        #	print STDERR $node->id,"\n";
        if ( $node->id ) {
            my $name = $names->{ $node->id } || $node->id;
            my $label = $name;
            $label =~ s/_/\./;
            $label =~ s/_/\\n/g;
            my ($accn1) = $label =~ /^([^\\]+)/;
            my $url = "SynView.pl?";
            $url .= "accn1=$accn1&" if $accn1;

            #	    print STDERR $name,": ", Dumper $node;
            foreach my $sister ( $node->ancestor->get_all_Descendents ) {
                next if $sister->internal_id == $node->internal_id;
                next unless $sister->id;
                my $tmp = $names->{ $sister->id };
                next unless $tmp;
                $tmp =~ s/_/ /g;
                my ($accn2) = $tmp =~ /^(\S+)/;
                $url .= "accn2=$accn2" if $accn2;
                last if $accn2;
            }
            while ( my ( $pat, $mod ) = each %$node_labels ) {
                if ( $name =~ /$pat/ ) {
                    $content .=
"\t\t\"$name\"\t[label = \"$label\", URL=\"$url\", $mod ];\n";
                }
            }
        }
        else {
            my $id = $node->internal_id;
            $content .= "\t\t\"$id\"\t[label=\"$id\", shape=point ];\n";
        }
    }
    return $content;
}

sub get_labels {
    my $text = $FORM->param('nodes');
    my %mods;
    return \%mods unless $text;
    foreach my $line ( split /\n/, $text ) {
        my ( $key, $mod ) = split /\s+/, $line, 2;
        $mod =~ s/,\s+/,/g;
        $mod =~ s/\s+/,/g;
        next unless $key;
        next unless $mod;
        $mods{$key} = $mod;
    }
    return \%mods;
}

sub process_node {
    my $node  = shift;
    my $names = shift;
    my $name  = $node->id ? $names->{ $node->id } : $node->internal_id;

    #    my $len =  ($node->branch_length) if $node->branch_length;
    my $content;
    foreach my $cnode ( $node->each_Descendent ) {
        my $cname = $cnode->id ? $names->{ $cnode->id } : $cnode->internal_id;
        $content .= "\t\"$name\" -- \"$cname\" ";
        my $len = $cnode->branch_length if $cnode->branch_length;
        if ($len) {
            $content .= "[label=\"" . sprintf( "%.2f", $len ) . "\"";

            #	$content .= ", len=$len" if defined $len;
            $content .= "];";
        }
        $content .= "\n";
        my $tmp = process_node( $cnode, $names );
        $content .= $tmp if $tmp;
    }
    return $content;
}

sub get_names {
    my $file = shift;
    my $tree = shift;

    #    open (IN, $file);
    my %names;

    #     my $flag = 0;
    #     while (<IN>)
    #       {
    # 	chomp;
    # 	if ($flag)
    # 	  {
    # 	    s/^\s+//g;
    # 	    my ($num, $name) = split /\s+/, $_, 2;
    # 	    next unless ($num =~ /^\d+$/);
    # 	    $name =~ s/,$// if $name;
    # 	    $name =~ s/\|/_/g if $name;
    # 	    $name =~ s/jgi_//g if $name;
    # 	    $name =~ s/(Chlre3.\d+)_.*/$1/ if $name;
    # 	    $name =~ s/'//g;
    # 	    $name =~ s/-/_/g;
    # 	    $name =~ s/\./_/g;
    # 	    last if $_ =~ ";";
    # 	    $names{$num} = $name;
    # 	    $names{$name} = $num;

# 	  }
# 	$flag = 1 if /Translate/;
#       }
#     close IN;
#     #check to see if we got the node names, otherwise, grab the ids from the tree
#     print STDERR Dumper \%names;
#    unless (keys %names)
#      {
    foreach my $node ( $tree->get_nodes ) {
        $names{ $node->id } = $node->id if $node->id;
    }

    #      }
    return \%names;
}

sub generate_mod_file {
    my %opts     = @_;
    my $stuff    = $opts{mod};
    my $basefile = $opts{basefile} . ".mod";
    open( OUT, ">$basefile" );
    print OUT $stuff;
    close OUT;
    return $basefile;
}

sub load_mod {
    my $file = shift;
    open( IN, $file );
    my $data;
    while (<IN>) { $data .= $_; }
    close IN;
    return $data;
}

sub translate_content {
    my %opts    = @_;
    my $content = $opts{content};
    my $table   = $opts{table};
    foreach ( split /\n/, $table ) {
        s/^\s+//g;
        my ( $num, $name ) = split /\s+/, $_, 2;
        next unless ( $num =~ /^\d+$/ );
        $name =~ s/,$//                if $name;
        $name =~ s/\|/_/g              if $name;
        $name =~ s/jgi_//g             if $name;
        $name =~ s/(Chlre3.\d+)_.*/$1/ if $name;
        $name =~ s/'//g;
        $name =~ s/-/_/g;
        $name =~ s/\./_/g;

        #	    $name =~ s/_/\./g;
        #	    $name =~ s/_/\\n/ if $name;
        last if $_ =~ ";";
        $content =~ s/([\(\):,])$num([\(\):,])/$1$name$2/g;
    }
    return $content;
}
