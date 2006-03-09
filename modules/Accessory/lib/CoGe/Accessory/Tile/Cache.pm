package CoGe::Accessory::Tile::Cache;

require LWP::Simple;
use base qw(Class::Accessor);
use strict; 
use Data::Dumper;
use BerkeleyDB;

$ENV{'PATH'} = '/bin:/usr/bin:/usr/local/bin';

my $db;
BEGIN {
  my $VERSION = 0.01;
  my @METHODS = qw(
    bypass_cache cache_file  img_generator_url DEBUG  force_retile
    normalize upp_base zoom_power);

  __PACKAGE__->mk_accessors( @METHODS );
}
=for AccessorMethods:
    bypass_cache   # if set to 1 will request a new image without caching it
    cache_file     # the path to the cache file
    img_generator_url  # the web address of the script that generates images
    DEBUG   # turn on (1) off (0) debugging
    force_retile # if (1), fetch the remote tile and overwrite it in the cache
    normalize  # if(1) changes the pixel value to the value at the left corner to allow cacing
    upp_base  # pixels per unit at the lowest zoom level
    zoom_power  # power change between zoom z(n) = upp_base * zoom_power ^ (n-1)
=cut


sub new {
    my $class = shift;
    my $this = {};
    bless $this,$class;
    $this->initialize(@_);
    return $this;
}

sub initialize {
    my $this = shift;
    my ($img_generator_url,$upp_base,$cache_file) = @_;

    # defaults, can all be set with accessorized methods
    $this->zoom_power(2);
    $this->img_generator_url($img_generator_url);
    $this->upp_base($upp_base);
    $this->cache_file($cache_file ||'tile.cache');
    $this->normalize(1);
    $this->bypass_cache(0);
    $this->force_retile(0);
     my $env = new BerkeleyDB::Env
            -Cachesize  => 65536, 
            -Flags      => DB_INIT_LOCK,
            -LockDetect => DB_LOCK_YOUNGEST;
     $db = new BerkeleyDB::Hash
                -Filename   => $cache_file, 
                -Pagesize   => 65536,
                -Flags      => DB_CREATE,
                -Env        => $env
      or die "Cannot open $cache_file: $!\n check that directory is web-writable" ;
#		tie %internal_cache => 'DB_File' , $this->cache_file, O_RDWR|O_CREAT, 0666, $DB_HASH or die "cant access $cache_file for the Cache file\n";
}

sub DESTROY {
}

sub get_tile {
    my $this = shift;
    my $base_url = $this->img_generator_url();
    my $opts_href = shift;

    my $IMG;  # the return value;
    $this->{tile_size} = $opts_href->{tile_size};

    # cache unless specifically told not to.  ie. bypass_cache => 1
    my $bypass_cache = $this->bypass_cache();

    # if force_cache overwrite image even if it's already in cache 
    my $force_retile = $this->force_retile();

    # dont need this in the hash any longer
    delete $opts_href->{bypass_cache} if defined $opts_href->{bypass_cache};

    # get the url that will be the key in the db and check to 
    my $param_url = $this->default_url_generator($opts_href);
    my $full_url = $base_url . $param_url; 
    chomp $full_url;

    if($this->DEBUG()){
      print STDERR "bypass_cache: ", $bypass_cache, "\n";
      print STDERR "full_url: ", $full_url ,"|\n";
    }
    # cache  BerkeleyDB returns 0 on success (for C compatiblity)
    my $in_cache = ! $db->db_get($param_url, $IMG);
    if(! $bypass_cache ) {
        
        if($this->DEBUG()){
          if( $in_cache ) {
            print STDERR "status: fetching from cache\n";
          }
          else{
            print STDERR "status: not in cache, fetching and caching \n";
          }
        }

        if(! $force_retile ){
          # TODO: fix, this is a temporary hack to fix non-caching
            if(! $in_cache ) {
                my $tmp ='';
                my $maxtries = 4;
                while(length($tmp)<100 && $maxtries){
                  $tmp = LWP::Simple::get($full_url);
                  $maxtries--;
                }
                $maxtries = 4; 
                my $success = 0;
                while(!$success && $maxtries){
                  $success = !$db->db_put($param_url, $tmp);
                  print STDERR "in sucess : $success \t\t";
                  $maxtries--;
                }
                if($this->DEBUG()){
                    print STDERR 'status of put: ', $success,"\n";
                }
            }
        }
        else {  # force retile
           $db->db_del($param_url);
           my $success = !$db->db_put($param_url, LWP::Simple::get($full_url));
           if($this->DEBUG()){
              print STDERR 'success of forced retile: ', $success,"\n";
           }
        }

        if($this->DEBUG()){
          print STDERR "status: got image, returning\n\n";
        }
       $db->db_sync();
       $db->db_get($param_url, $IMG);
       return $IMG;
    }

    # dont cache just return it  # if bypass_cache => 1
    if($this->DEBUG()){
       print STDERR "status: fetching but not caching\n\n";
    }
    return LWP::Simple::get($full_url);
}

sub default_url_generator {
  my $this = shift;
  my $h_ref = shift;

  my $purl = "";
  # here, it's best with the BTREE algorithm if tiles that will be
  # accessed togehter are stored together. this will sort of be the
  # case if x, and y are hte last. but also have zoomlevel...
  # when using tile, just use ztilex, ztiley
  
  # even if the request is by tile, we convert to pixels and store
  # that way
  my $do_normalize = $this->normalize();
  if($this->DEBUG()){
      print STDERR "normalize: $do_normalize\n";
  }
  if($do_normalize eq 'by_tile' ){
    ($h_ref->{'x'},$h_ref->{'y'}) 
      = $this->xy_from_tilez( $h_ref->{'x'} || 0, $h_ref->{'y'} 
      || 0,$h_ref->{'z'});
  }
  elsif( $do_normalize ){ 
    ($h_ref->{'x'},$h_ref->{'y'}) 
      = $this->normalize_xyz( $h_ref->{'x'} || 0, $h_ref->{'y'} 
      || 0,$h_ref->{'z'});
  }
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
  my $units_per_pixel = $this->upp_base  * $this->zoom_power ** $z;
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
  my $units_per_pixel = $this->upp_base * $this->zoom_power ** $z;
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
  my $units_per_pixel = $this->upp_base * $this->zoom_power ** $z;
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

