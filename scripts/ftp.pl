#!/usr/bin/env perl
use v5.14;
use strict;
use warnings;

use Getopt::Long qw(GetOptions);
use CoGe::Accessory::Web qw(ftp_get_file);

our ($url, $username, $password, $dest_path);

GetOptions(
    "url=s"       => \$url,
    "username=s"  => \$username,
    "password=s"  => \$password,
    "dest_path=s" => \$dest_path
);

$| = 1;
print STDOUT "Starting $0 (pid $$)\n", qx/ps -o args $$/;

unless ($url && $dest_path) {
    print STDERR "Missing required argument\n";
    exit(-1);
}

my $result = CoGe::Accessory::Web::ftp_get_file(
    url => $url, 
    username => $username, 
    password => $password, 
    dest_path => $dest_path
);

if ( !$result || $result->{error} ) {
    print STDERR "Failed: ", $result->{error}, "\n";
    exit(-1);
}

exit;
