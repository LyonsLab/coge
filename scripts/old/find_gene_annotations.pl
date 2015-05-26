#!/usr/bin/perl -w

use strict;
use CoGeX;
use DBIxProfiler;
use Getopt::Long;

my ($help, $org, $version, $type, $chr, $DEBUG, $name, $anno);

GetOptions ("h|help" =>  \$help,
            "o|org=s" => \$org,
            "v|version=s" => \$version,
            "t|type=s"    => \$type,
            "d|debug"     => \$DEBUG,
            "c|chr=s"     => \$chr,
	    "n|name|accn=s" =>\$name,
	    "a|anno=s"=>\$anno,
            );

my $connstr = 'dbi:mysql:dbname=DB;host=HOST;port=PORT';
my $coge = CoGeX->connect($connstr, 'USER', 'PASSWORD' );
#$coge->storage->debugobj(new DBIxProfiler());
#$coge->storage->debug(1);

cogesearch (accn=>$name, anno=>$anno,org=>$org, ver=>$version, type=>$type);

sub cogesearch
  {
    my %opts = @_;
    my $accn = $opts{accn};
    my $anno = $opts{anno};
    my $type = $opts{type};
    my $org = $opts{org};
    my $version = $opts{ver};
    my $feat_accn_wild = $opts{feat_name_wild} || "both";
    my $feat_anno_wild = $opts{feat_anno_wild} || "both";
#    print STDERR Dumper \%opts;
    my $blank = qq{<input type="hidden" id="accn_select">};
#    print STDERR "cogesearch: $accn\n";
#    print STDERR Dumper @_;
    my $weak_query = "Query needs to be better defined.";
    return $weak_query.$blank unless length($accn) > 2 || $type || $org || length($anno) > 5;
    if (!$accn && !$anno)
      {
	return $weak_query.$blank unless $org && $type;
      }
    my $html;
    my %seen;
    my @opts;
    $accn = "%".$accn if $accn && ($feat_accn_wild eq "both" || $feat_accn_wild eq "left");
    $accn = $accn."%" if $accn && ($feat_accn_wild eq "both" || $feat_accn_wild eq "right");
    $anno = "%".$anno if $anno && ($feat_anno_wild eq "both" || $feat_anno_wild eq "left");
    $anno = $anno."%" if $anno && ($feat_anno_wild eq "both" || $feat_anno_wild eq "right");
    my $search = {'me.name'=>{like=>$accn}} if $accn;
    $search->{annotation}={like=>$anno} if $anno;
    $search->{feature_type_id}=$type if $type;
    $search->{organism_id}=$org if $org;
    $search->{version} = $version if $version;
    my $join = {'feature'=>'dataset'};
    $join->{'feature'} = ['dataset','annotations'] if $anno;
    foreach my $name ($coge->resultset('FeatureName')->search(
							      $search,
							      {join=>$join,
							       order_by=>'name ASC',
							       prefetch=>'feature',
							      },
							     ))
      {
	my $item = $name->feature->id;
	next if $seen{uc($item)};
	next unless $name->feature->type->name =~ /gene/;
	$seen{uc($item)}++;
	print join ("\t", $name->name, $name->name." ". join (", ", map {$_->annotation} $name->feature->annotations)),"\n";
      }
  }
