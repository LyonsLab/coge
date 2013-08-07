#!/usr/bin/perl -wT

use strict;
use File::Spec::Functions;
use File::Path;
use LWP::Simple;
use Data::Dumper;
use CoGe::Accessory::Web;
use CGI;
no warnings 'redefine';

use vars qw($P $IMGURL $BASEDIR);

$ENV{'PATH'} = '';

my $cgi = CGI->new;

$P = CoGe::Accessory::Web::get_defaults();

$IMGURL = $P->{SERVER} . 'GenomePNG.pl?';

# where to start the caching
$BASEDIR = $P->{IMAGE_CACHE};
( my $x, $BASEDIR ) = CoGe::Accessory::Web::check_taint($BASEDIR);
if ( !-e $BASEDIR ) { mkdir($BASEDIR); }
if ( !-e $BASEDIR ) { warn "unable to find and create $BASEDIR"; exit; }
print "Content-type: image/png; mode=24bit\n\n";

my $basedir = [ split( /[\\\/]/, $BASEDIR ) ];

my ($dir) = get_dir_array( cgi => $cgi );
my $basepath = catfile(@$basedir);
my $fn = catfile( $basepath, @$dir ) . '.png';

my $data;
if ( !-e $fn ) {
    pop @$dir;    # get rid of the file name
    umask(0);
    mkpath( $basepath . '/' . join( "/", @$dir ) );
    LWP::Simple::getstore( $IMGURL . $ENV{QUERY_STRING}, $fn );

    #   chmod (0777, $fn);
}

# DUMP THE IMAGE
{
    local ( *IMG, $/ );
    my $mesg = "cant open $fn\nmake sure _cache_ dir is web writeable\n";
    open( IMG, "<", $fn ) or warn "$mesg: $!";
    binmode IMG;
    print <IMG>;
    close(IMG);
}

##################################################
# Use %ENV to find where in directory structure
# we are and make an array of directory names
##################################################
sub get_dir_array {
    my %opts   = @_;
    my $cgi    = $opts{cgi};
    my @dir    = ();
    my $xmin   = $cgi->param('xmin');     #delete $query_pairs{xmin};
    my $xmax   = $cgi->param('xmax');     #delete $query_pairs{xmax};
    my $ds     = $cgi->param('ds');       #delete $query_pairs{ds};
    my $dsg    = $cgi->param('dsg');      #delete $query_pairs{dsg};
    my $tilew  = $cgi->param('width');    #$query_pairs{width};
    my $gstid  = $cgi->param('gstid');    #$query_pairs{gstid};
    my $expid  = $cgi->param('expid');
    my $chr    = $cgi->param('chr');
    my $layers = $cgi->param('layers');
    my $MAX = int( log( 1000 * abs( $xmax - $xmin ) / $tilew ) / log(10) );
    push @dir, ( "dsg", $dsg ) if $dsg;
    push @dir, ( "ds",  $ds )  if $ds;
    push @dir, ("chr__$chr");
    push @dir, ("layers__$layers");
    push @dir, ("gstid__$gstid");
    push @dir, ("expid__$expid") if $expid;
    push @dir, ("width__$tilew");

    my @vals;
    my $val = $xmin;
    $val =~ s/-/n/g;
    $val = scalar reverse($val);
    while ( $val =~ /(\d{1,$MAX}n?)/g ) {
        unshift( @vals, scalar reverse($1) );
    }
    my $first = shift @vals;
    push( @dir, 'xmin__' . $first );
    map { push( @dir, $_ ) } @vals;

    @vals = undef;
    $val  = $xmax;
    $val =~ s/-/n/g;
    $val = scalar reverse($val);
    while ( $val =~ /(\d{1,$MAX}n?)/g ) {
        unshift( @vals, scalar reverse($1) );
    }
    $first = shift @vals;
    push( @dir, 'xmax__' . $first );
    map { push( @dir, $_ ) } @vals;

    my $x;
    my @tested;
    foreach my $item (@dir) {
        ( $x, $item ) = CoGe::Accessory::Web::check_taint($item);
        push @tested, $item if $x && $item;
    }
    return \@tested;
}
