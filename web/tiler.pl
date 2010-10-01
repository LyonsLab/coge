#!/usr/bin/perl -wT

use strict;
use File::Spec::Functions;
use File::Path;
use LWP::Simple;
use Data::Dumper;
use CoGe::Accessory::Web;
no warnings 'redefine';

$ENV{'PATH'} = '';
use vars qw($P $IMGURL $BASEDIR);
$P = CoGe::Accessory::Web::get_defaults();

$IMGURL = $P->{SERVER}.'GenomePNG.pl?';
# where to start the caching
$BASEDIR = $P->{IMAGE_CACHE};
if (! -e $BASEDIR ){ mkdir($BASEDIR);  }
if (! -e $BASEDIR ){ warn "unable to find and create $BASEDIR";exit;}
print "Content-type: image/png; mode=24bit\n\n";

my $basedir = [split(/[\\\/]/,$BASEDIR)];

my ($dir) = get_dir_array();
my $reldir = catfile(@$dir) . '.png';
my $basepath = catfile(@$basedir);


my $fn = catfile($basepath,@$dir) . '.png' ;
my $data;
if(!-e $fn){
   pop @$dir; # get rid of the file name
   umask (0);
   mkpath($basepath . '/' . join("/",@$dir));
   LWP::Simple::getstore($IMGURL . $ENV{QUERY_STRING},$fn);
#   chmod (0777, $fn);
#   print STDERR $IMGURL.$ENV{QUERY_STRING},"\n";
}


# DUMP THE IMAGE
{
local( *IMG,$/ ); 
my $mesg = "cant open $fn\nmake sure _cache_ dir is web writeable\n";
open( IMG,"<", $fn) or warn "$mesg: $!"; binmode IMG;
print <IMG>; close(IMG);
}


##################################################
# Use %ENV to find where in directory structure
# we are and make an array of directory names  
##################################################
sub get_dir_array {
    my ($qstr) = $ENV{QUERY_STRING} =~ /(.*)/;
    #print STDERR $IMGURL.$ENV{QUERY_STRING},"\n";
    $qstr =~ s/[\s-]//;
    my @dir = ();
#    print STDERR $qstr,"\n";
    # parse the url;
    my %query_pairs = map { split('=', $_) } split(/&/,$qstr);
#    print STDERR Dumper \%query_pairs;
    # keep order of entries but put the spatial pars last becuase
    # it just works out best that way.
    my $xmin = delete $query_pairs{xmin};
    my $xmax = delete $query_pairs{xmax};
    my $ds = delete $query_pairs{ds};
    my $tilew = $query_pairs{width};
    my $gstid = $query_pairs{gstid};
    my $MAX = int(log(1000*abs($xmax - $xmin)/$tilew)/log(10));
    my @keyvals = ('ds',sort keys %query_pairs);
    $query_pairs{xmin} = $xmin;
    $query_pairs{xmax} = $xmax;
    $query_pairs{ds} = $ds;
    push(@keyvals,"xmin");
    push(@keyvals,"xmax");
    
    # &layer=fred becomes layer__fred
    foreach my $key( @keyvals ){
        next unless $key;
        my $val = $query_pairs{$key};
        my @vals;
        # 123456 becomes /x__123/456/  if $MAX == 3
        # 123456 becomes /x__123456/   if $MAX  > 5
        if($key eq "xmin" || $key eq "xmax"){
            if($val !~ /(\D|-)+/){
                $val =~ s/-/n/g;
                $val = scalar reverse($val);
                while($val =~ /(\d{1,$MAX}n?)/g){
                    unshift(@vals,scalar reverse($1));
                }
            }
        }
        @vals = ($val) unless @vals;
        my $first = shift @vals;
        push(@dir,$key . '__' . $first);
        map { push(@dir,$_) } @vals;
    }
    return \@dir; 
}
