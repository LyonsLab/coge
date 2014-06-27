package CoGe::Accessory::Tile::Cache;

require LWP::Simple;
use base qw(Class::Accessor);
use strict;
use Data::Dumper;
use BerkeleyDB;
use Storable qw(freeze thaw);
use Benchmark;

$ENV{'PATH'} = '/bin:/usr/bin:/usr/local/bin';

my $db;
#my $lock;
BEGIN {
  my $VERSION = 0.01;
  my @ACCESS = qw(
    bypass_cache cache_file  img_generator_url DEBUG  force_retile
    normalize config
               );

  __PACKAGE__->mk_accessors( @ACCESS );
}

=for AccessorMethods:
    bypass_cache   # if set to 1 will request a new image without caching it
    cache_file     # the path to the cache file
    DEBUG   # turn on (1) off (0) debugging
    force_retile # if (1), fetch the remote tile and overwrite it in the cache
    normalize  # if(1) changes the pixel value to the value at the left corner to allow cacing
    zoom_power  # power change between zoom z(n) = upp_base * zoom_power ^ (n-1)

=cut

sub new {
    my $class = shift;
    my $this = {};
    $this->{db_has_changed} = 0;
    bless $this,$class;
    $this->initialize(@_);
    return $this;
}

sub delete_config {
  my $this = shift;
  $db->db_del('Tile_Cache_Config');
  if($this->DEBUG()){
    print STDERR 'DELETED KEY' . "\n\n";
  }
}

# only called if the config info is not stored in the cache.
sub parse_js_config {
  my $JS;
  eval {
    my $d = $ENV{DOCUMENT_ROOT};
    my $scriptdir = $ENV{SCRIPT_NAME};
    $scriptdir =~ s/(.*)([\\|\/])\w+\.pl/$1$2/;
    my $file = $d . $scriptdir . "js/tilerConfig.js";
    $file =~ s/htdocs\/// unless (-r $file);
    $file = "/opt/apache/coge/web/js/tilerConfig.js" unless -f $file; #last ditch
    my $open = open ($JS, "<" , $file );
    $open = open ($JS, "<" , $file ) or die "cant find tilerConfig.js: file tried $file";
  };

  #SLURPEE  (switch to sysread)??
  my $javascript = join("\n",<$JS>);
  close($JS);
  return 'cant open (or find) config.js' if $@;

  # remove /* comments */
  # from the perlop manpage;
  $javascript =~ s {
      /\*     # Match the opening delimiter.
      .*?     # Match a minimal number of characters.
      \*/     # Match the closing delimiter.
  } []gsxm;

  # use JSON::Syck qw(Load);

  use JSON;
  local $JSON::BareKey = 1;

  my $str = "{";
  my @J = split("\n",$javascript);
  foreach my $line (@J){
    # get rid of the inital assignment statement
    next if $line =~ /var\s+config\s*=\s*{/;
    chomp $line;
    $line =~ s/\s+//;
    $line =~ s/http:\/\//http_ss/g;
#    next if $line =~ s/\/\//$1/; I dunno what Brent is doing here.  switched to line below
    next if $line =~ /^\/\//;
    $line =~ s/\An(.*)/$1/;
    $str .= $line;
  }
  chomp $str;
  $str =~ s/,\}/}/ms;
  $str =~ s/;//ms;
  $str =~ s/'/"/gms;
  $str =~ s/http_ss/http:\/\//ms;
  # return Load( "$str");
  return  jsonToObj( "$str" );

}

sub initialize {
    my $this = shift;
    my ($cache_file) = @_;
    my $config;

    # defaults, can all be set with accessorized methods
     my $env = new BerkeleyDB::Env
#            -Cachesize  => 131072,
            # this one works best.. but not great
            -Flags      => DB_CREATE |  DB_INIT_LOCK | DB_INIT_MPOOL;
#            -Flags      => DB_CREATE |   DB_INIT_MPOOL | DB_INIT_CDB;

     $db = new BerkeleyDB::Hash
                -Filename   => $cache_file,
                -Pagesize   => 65536,
                -Flags      => DB_CREATE,
                -Env        => $env
      or die "Cannot open $cache_file: $!\n check that directory is web-writable" ;
     my $error = $db->db_get('Tile_Cache_Config',$config);
     #$lock = $db->cds_lock();
     if($error){
       $config = __PACKAGE__->parse_js_config();
       $db->db_put('Tile_Cache_Config',freeze($config));
       #$this->{db_has_changed} = 1;
     }else{
       print STDERR 'TILE_CONFIG_PARS in cache' . "\n" if $this->DEBUG;
       $config = thaw($config);
     }

     # uggghhhh need to go from "['x','y']" to ['x','y']
     # JSON::Syck ???
     $config->{SPATIAL_PARS_NAME} = eval($config->{SPATIAL_PARS_NAME});
     $this->config($config);

    $this->{zoom_power} = $this->config()->{ZOOM_POWER} || 2;
    $this->{upp_x_base}= $this->config()->{BASE_UNITS_PER_PIXEL_X};

    $this->{upp_y_base} =  $this->config()->{BASE_UNITS_PER_PIXEL_Y}
        || $this->{upp_x_base};

    $this->{spatial_pars} = $this->config()->{SPATIAL_PARS};

    # accessor methods
    $this->cache_file($this->config()->{PATH_TO_CACHE_FILE});
    $this->img_generator_url($this->config()->{IMG_GENERATOR_URL});
    $this->normalize(0);
    $this->bypass_cache(0);
    $this->force_retile(0);
}

sub DESTROY {
  my $this = shift;
  if($this->{db_has_changed}){
    #$lock->cds_unlock();
    $db->db_sync();

  }

}

sub get_tile
  {
    my $t0 = new Benchmark;
    my $this = shift;
    my $base_url = $this->img_generator_url();
    my $opts_href = shift;
    my $IMG = $opts_href->{img};  # the return value;
    delete $opts_href->{img};
    $this->{tile_size} = $opts_href->{tile_size};

    # cache unless specifically told not to.  ie. bypass_cache => 1
    my $bypass_cache = $this->bypass_cache();

    # if force_cache overwrite image even if it's already in cache
    my $force_retile = $this->force_retile();

    # dont need this in the hash any longer
    delete $opts_href->{bypass_cache} if defined $opts_href->{bypass_cache};
    # get the url that will be the key in the db and check to
    my $param_url = $this->default_url_generator($opts_href);
    if($this->DEBUG()){
      print STDERR 'param_url: ' . $param_url . "\n\n";
    }
    my $full_url = $base_url . $param_url;
    chomp $full_url;

    if($force_retile){
      $db->db_del($param_url);
    }

    # cache  BerkeleyDB returns 0 on success (for C compatiblity)
    my $in_cache = ! $db->db_get($param_url, $IMG) unless $IMG;
    return $IMG if $in_cache;
    print STDERR "Getting image from web $full_url\n" if $this->DEBUG;
    my $tmp = $IMG ? $IMG : LWP::Simple::get($full_url);
    return $tmp if $bypass_cache;

    # if we made it this far, it's not in the cache;
    # but it wants to be.
    #my $cursor = $db->db_cursor(DB_WRITECURSOR);
    #if(!$cursor){ die "Could not get cursor\n" }
    my $t1 = new Benchmark;
    my $success = ! $db->db_put($param_url, $tmp);
    my $t2 = new Benchmark;
    $db->db_sync();
#    print STDERR Dumper $db->db_stat;
    $this->{db_has_changed} = 1;
    if($this->DEBUG()){
        print STDERR 'status of put: ', $success,"\n";
    }
    #$cursor->c_close();
    my $t3 = new Benchmark;
    my $it_time = timestr(timediff($t1, $t0));
    my $st_time = timestr(timediff($t2, $t1));
    my $ft_time = timestr(timediff($t3, $t2));
    print STDERR qq{
Initialization: $it_time
Storage:        $st_time
Finalization:   $ft_time
} if 0;
    return $tmp;
}

sub get_tile_from_db
  {
    my $this = shift;
    my $base_url = $this->img_generator_url();
    my $opts_href = shift;
    $this->{tile_size} = $opts_href->{tile_size};
    my $param_url = $this->default_url_generator($opts_href);
    my $full_url = $base_url . $param_url;
    chomp $full_url;
    my $TMP;
    my $in_cache = ! $db->db_get($param_url, $TMP);
    return $TMP;
  }

sub default_url_generator {
  my $this = shift;
  my $spatials = $this->{spatial_pars};

  my $h_ref = shift;

  my $purl = "";
  # here, it's best with the BTREE algorithm if tiles that will be
  # accessed togehter are stored together. this will sort of be the
  # case if x, and y are hte last. but also have zoomlevel...
  # when using tile, just use ztilex, ztiley

  # even if the request is by tile, we convert to pixels and store
  # that way
#  my $do_normalize = $this->normalize();
#  if($this->DEBUG()){
#      print STDERR "normalize: $do_normalize\n";
#  }
#  if($do_normalize eq 'by_tile' ){
#    ($h_ref->{'x'},$h_ref->{'y'})
#      = $this->xy_from_tilez( $h_ref->{'x'} || 0, $h_ref->{'y'}
#      || 0,$h_ref->{'z'});
#  }
#  elsif( $do_normalize ){
#    ($h_ref->{'x'},$h_ref->{'y'})
#      = $this->normalize_xyz( $h_ref->{'x'} || 0, $h_ref->{'y'}
#      || 0,$h_ref->{'z'});
#  }
  map { $purl .= '&' . $_ . '=' . $h_ref->{$_} } sort keys %{$h_ref};
  return $purl;
}

# at the lowest zoom level (z == 0) the units/pixel is upp;
# so if the ZOOM_POWER == 2, then that is upp*2^0;
# at zoom_level 1, it is upp is upp*2^1. etc.
# so we need:
#   the upp at level 0:     $units_per_pixel_base
#   the power of each zoom  $ZOOM_POWER
#
# then its no problemo
sub normalize_xyz {
  my $this = shift;
  my ($x,$y,$z) = @_;
  my ($tile_x,$tile_y);
  my $units_per_pixel = $this->{upp_x_base}  * $this->{zoom_power} ** $z;
  my $units_per_tile  = $units_per_pixel * $this->{tile_size};
 # might want to add a units per tile y
  $tile_x = int($x / $units_per_tile);
  $tile_y = int($y / $units_per_tile);
  return ($tile_x*$units_per_tile,$tile_y*$units_per_tile);

}

sub tile_from_xyz {
  my $this = shift;
  my ($x,$y,$z) = @_;
  my ($tile_x,$tile_y);
  my $units_per_pixel = $this->{upp_x_base} * $this->{zoom_power} ** $z;
  my $units_per_tile = $units_per_pixel * $this->{tile_size};
  if($this->DEBUG()){
    print STDERR "normalizing x and y ...\n";
    print STDERR "x original : $x\n";
    print STDERR "units [per pixel, per tile]: $units_per_pixel , $units_per_tile\n";
    print STDERR "tiles [x, y]: $tile_x, $tile_y\n";
  }
  # TODO: always add 1? maybe some systems start with 0;
  # ++$x; ++$y;
  $tile_x = int($x / $units_per_tile);
  $tile_y = int($y / $units_per_tile);
  return ($tile_x,$tile_y);
}

sub xy_from_tilez {
  my $this = shift;
  my ($tile_x,$tile_y,$z) = @_;
  # this line works for z == 0
  #my $units_per_pixel = $this->upp_base * $this->zoom_power ** $z;
  my $units_per_pixel = $this->{upp_x_base} * $this->{zoom_power} ** $z;
  my $units_per_tile = $units_per_pixel * $this->{tile_size};
  # TODO: always add 1? maybe some systems start with 0;

  my $x = $units_per_tile * $tile_x; # + 1;
  my $y = $units_per_tile * $tile_y; # + 1;
  if($this->DEBUG()){
    print STDERR "getting and x and y from tiles ...\n";
    print STDERR "x normalized : $x\n";
    print STDERR "units [per pixel, per tile]: $units_per_pixel , $units_per_tile\n";
    print STDERR "tiles [x, y]: $tile_x, $tile_y\n";
  }
  return ($x,$y);
}

1;

=head1 NAME

Cache

=head1 SYNOPSIS

use Cache

=head1 DESCRIPTION

=head1 USAGE

=head1 BUGS

=head1 SUPPORT

=head1 AUTHOR

Eric Lyons

=head1 COPYRIGHT

This program is free software; you can redistribute
it and/or modify it under the same terms as Perl itself.

The full text of the license can be found in the
LICENSE file included with this module.

=head1 SEE ALSO

=cut
