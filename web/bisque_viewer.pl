#! /usr/bin/perl -w

use strict;
use CGI;
use Data::Dumper;
use LWP::UserAgent;

my $request = CGI->new;
my $id = $request->param('id');

binmode(STDOUT);

my $ua = LWP::UserAgent->new();
my $req = HTTP::Request->new(GET => 'https://bisque.cyverse.org/client_service/embedded?view=2d&resource=https://bisque.cyverse.org/data_service/' . $id);
$req->authorization_basic('coge', 'C0G3Rules!');
my $header_sent = 0;
my $res = $ua->request($req,
sub {
    my($chunk, $res) = @_;
    if (!$header_sent) {
        print CGI::header($res->header('Content-Type'));
        $header_sent = 1;
    }
    print $chunk;
    warn $chunk;
});
