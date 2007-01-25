#! /usr/bin/perl -w

use strict;
use File::Spec::Functions;
require LWP::Simple;
use vars qw($IMGURL $ZOOM_PAR $SPATIAL_PARS $BUPT);

$IMGURL = 'http://biocon.berkeley.edu/CoGe/GenomePNG.pl?layers=all,gc&';
$ZOOM_PAR = 'z';
$SPATIAL_PARS = ['x'];
$BUPT = 10;

###DEBUGGING testing to make sure that 'x' extrapolates to 0 based on
#zoom and $BUPT
my %param = map { split('=', $_) } split(/&/,$ENV{QUERY_STRING});
if ($param{x}%($BUPT*(2**$param{z})))
  {
   my $mod = $param{x}%($BUPT*(2**$param{z}));
    print STDERR "x location error in tiler.pl: x=".$param{x}."\tz=".$param{z}."\tmod=$mod\n ";
  }

print "Content-type: image/png\n\n";

my ($basedir,$dir) = get_dir_array();
my $reldir = catfile(@$dir) . '.png';
my $basepath = catfile(@$basedir);

my $fn = catfile($basepath,@$dir) . '.png' ;
($fn) = $fn =~ /(.*)/;
my $data;
if(! -e $fn){
   make_dir($basepath,$dir);
   LWP::Simple::getstore($IMGURL . $ENV{QUERY_STRING},$fn);
}


# DUMP THE IMAGE
{
local( *IMG,$/ ); 
open( IMG, $fn) or die "cant open $fn\n"; binmode IMG;
print <IMG>; close(IMG);
}


#############################################
# Make the directory structure 
#############################################
sub make_dir {
    my ($dirstr,$dir) = @_;
    pop @$dir;
    $dirstr = catfile($dirstr,shift @$dir);
    umask 0;
    foreach my $subdir (@$dir){
        ($dirstr) = catdir($dirstr,$subdir) =~ /(.*)/;
        my $status = mkdir($dirstr,0777) if !-e $dirstr;
    }
}

#############################################
# Find the number of numbers so that there
# are never too many files in a single dir
#############################################
sub get_max {
    my $z = shift;
    my $MAX_FILES_PER_DIR = 1000;
    my $upt = $BUPT * 2 ** $z;
    return int(log($MAX_FILES_PER_DIR * $upt)/log(10));
}

# find see if a value is in an array;
sub grep_match {
    my ($needle,$haystack) = @_;
    return grep { $needle eq $_ } @$haystack;
}

##################################################
# Use %ENV to find where in directory structure
# we are and make an array of directory names  
##################################################
sub get_dir_array {
    my $qstr = $ENV{QUERY_STRING};
    my $base_dir = "/opt/apache";#$ENV{DOCUMENT_ROOT};
    my @basedir = split(/[\/\\]/,$base_dir);
    my @dir = split(/[\/\\]/,$ENV{SCRIPT_NAME});
    shift @dir;
    # get rid of file name;
    $dir[$#dir] = '_cache_';

    # parse the url;
    # keep order of entries but put the spatial pars last becuase
    # it just works out best that way.
    my @keyvals = map { split('=', $_) } split(/&/,$qstr);
    my %query_pairs = @keyvals;
    #@keyvals = map { $_ % 2 == 0 && !grep_match($keyvals[$_],$SPATIAL_PARS) ? $keyvals[$_] : '' } 0.. $#keyvals;
    #@keyvals = sort @keyvals;
    #push(@keyvals,@$SPATIAL_PARS);
    @keyvals = qw(ds chr iw z x);
    my $MAX = get_max($query_pairs{$ZOOM_PAR});
    foreach my $key( @keyvals ){
        next if !$key;
        my $val = $query_pairs{$key};
        my @vals;
        # 123456 becomes /x__123/456/  if $MAX == 3
        # 123456 becomes /x__123456/   if $MAX  > 5
        if(grep_match($key,$SPATIAL_PARS)){
            if($val !~ /(\D|-)+/){
                $val =~ s/-/n/g;
                $val = scalar reverse($val);
                while($val =~ /(\d{1,$MAX}n?)/g){
                    unshift(@vals,scalar reverse($1));
                }
            }
            # &layer=fred becomes layer__fred
        }
        @vals = ($val) unless @vals;
        my $first = shift @vals;
        push(@dir,$key . '__' . $first);
        map { push(@dir,$_) } @vals;
    }
    return (\@basedir,\@dir); 
}
