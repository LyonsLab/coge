#! /usr/bin/perl -w

use strict;
use CGI;
use Data::Dumper;
use LWP::UserAgent;

use CoGe::Accessory::Web qw(get_defaults);

my $request = CGI->new;
my $id = $request->param('id');

my $conf = get_defaults();

print CGI::header();

my $ua = LWP::UserAgent->new();
my $req = HTTP::Request->new(GET => 'https://bisque.cyverse.org/client_service/embedded?view=2d&resource=https://bisque.cyverse.org/data_service/' . $id);
$req->authorization_basic('coge', $conf->{BISQUE_PASS});
my $res = $ua->request($req,
sub {
    my($chunk, $res) = @_;
    print $chunk;
});
