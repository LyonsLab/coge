use v5.14;
use strict;
use warnings;

use CoGe::Accessory::Web qw(set_defaults url_for);
use Test::Exception;
use Test::More;
use URI;

my $config = {
    SERVER =>  "host.com",
    URL    =>  "test"
};

throws_ok { url_for("") } qr/CONFIG was not found/, "bail no config";

set_defaults({});
throws_ok { url_for("") } qr/SERVER option not found in CONFIG/, "bail SERVER not set";

set_defaults($config);
my $u1 = URI->new(url_for(""));
is $u1->scheme, "http",     "has correct scheme (default http)";
is $u1->host,   "host.com", "has correct host";
is $u1->path,   "/test/",   "has correct path";
is $u1->query,  undef,      "has correct query string";
is url_for(""), "http://host.com/test/", "has correct url";

# Set secure
$config->{SECURE} = 1;
set_defaults($config);
my $u2 = URI->new(url_for(""));
is $u2->scheme, "https",    "has correct scheme (default https)";
is $u2->host,   "host.com", "has correct host";
is $u2->path,   "/test/",   "has correct path";
is $u2->query,  undef,      "has correct query string";
is url_for(""), "https://host.com/test/", "has correct url";

# Force default
$config->{SERVER} = "http://localhost";
set_defaults($config);
my $u3 = URI->new(url_for(""));
is $u3->scheme, "http",      "has correct scheme (overrides default scheme)";
is $u3->host,   "localhost", "has correct host";
is $u3->path,   "/test/",    "has correct path";
is $u3->query,  undef,       "has correct query string";
is url_for("Test.pl"), "http://localhost/test/Test.pl", "has correct url";

# VERIFY QUERY STRING
my $url1 = url_for("Test.pl", a005 => 1, a05 => 2, a5 => 3);
my $u4 = URI->new($url1);
is $u4->scheme, "http",              "has correct scheme (overrides default scheme)";
is $u4->host,   "localhost",         "has correct host";
is $u4->path,   "/test/Test.pl",     "has correct path";
is $u4->query,  "a005=1&a05=2&a5=3", "has correct query string";
is $url1, "http://localhost/test/Test.pl?a005=1&a05=2&a5=3", "has correct url";

my $url2 = url_for("Test.pl", a05 => 2, a5 => 3, a005 => 1);
my $u5 = URI->new($url2);
is $u5->scheme, "http",              "has correct scheme (overrides default scheme)";
is $u5->host,   "localhost",         "has correct host";
is $u5->path,   "/test/Test.pl",     "has correct path";
is $u5->query,  "a005=1&a05=2&a5=3", "has correct query string regardless of order";
is $url2, "http://localhost/test/Test.pl?a005=1&a05=2&a5=3", "has correct url";

# Base url set in server
$config->{SERVER} = "http://localhost/test/";
$config->{URL} = "/test/";
set_defaults($config);
my $url3 = url_for("Test.pl", a05 => 2, a5 => 3, a005 => 1);
my $u6 = URI->new($url3);
is $u6->scheme, "http",              "has correct scheme (overrides default scheme)";
is $u6->host,   "localhost",         "has correct host";
is $u6->path,   "/test/Test.pl",     "has correct path";
is $u6->query,  "a005=1&a05=2&a5=3", "has correct query string regardless of order";
is $url3, "http://localhost/test/Test.pl?a005=1&a05=2&a5=3", "has correct url";

# Base url set in server in different case
$config->{SERVER} = "http://localhost/TEST/";
$config->{URL} = "/test/";
set_defaults($config);
my $url4 = url_for("Test.pl", a05 => 2, a5 => 3, a005 => 1);
my $u7 = URI->new($url4);
is $u7->scheme, "http",              "has correct scheme (overrides default scheme)";
is $u7->host,   "localhost",         "has correct host";
is $u7->path,   "/test/Test.pl",     "has correct path";
is $u7->query,  "a005=1&a05=2&a5=3", "has correct query string regardless of order";
is $url4, "http://localhost/test/Test.pl?a005=1&a05=2&a5=3", "has correct url regardless of case";

done_testing();
