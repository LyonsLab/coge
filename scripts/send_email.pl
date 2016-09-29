#!/usr/bin/env perl
use v5.14;
use strict;
use warnings;

use File::Path qw(mkpath);
use Getopt::Long qw(GetOptions);

use CoGe::Accessory::Web qw(send_email);
use CoGe::Accessory::Utils qw(to_pathname);
use URI::Escape::JavaScript qw(unescape);

our ($from, $to, $subject, $body, $done_file);

GetOptions(
    # Required params
    "from=s"      => \$from,
    "to=s"        => \$to,
    "subject=s"   => \$subject,
    "body=s"      => \$body,
    "done_file=s" => \$done_file
);

$| = 1;
print STDOUT "Starting $0 (pid $$)\n", qx/ps -o args $$/;

# Check required parameters
die "ERROR: from not specified" unless $from;
die "ERROR: to not specified" unless $to;
die "ERROR: subject not specified" unless $subject;
die "ERROR: body not specified" unless $body;

$from = unescape($from);
$to = unescape($to);
$subject = unescape($subject);
$body = unescape($body);

#send email
send_email(from => $from, to => $to, subject => $subject, body => $body);
    
# Save something in done file -- signals task completion to JEX
if ($done_file) {
    my $done_path = to_pathname($done_file);
    mkpath($done_path);
    open(my $fh, ">>", $done_file);
    say $fh "email sent: from:" . $from . " to:" . $to . " subject:" . $subject;
    close($fh);
}

exit;
