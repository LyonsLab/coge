#!/usr/bin/env perl
#-------------------------------------------------------------------------------
# Usage:
#    ./alert_on_filesystem_full.pl -names <fs1[,fs2,...]> -limit <percentage> -email <email_address>
#
#	 ./alert_on_filesystem_full.pl -names /,/scratch,/mysql -limit 85% -email coge@gmail.com
#
# Check given file-system partitions usage and send email notification if over limit
#
#-------------------------------------------------------------------------------

use strict;
use warnings;

use Getopt::Long qw(GetOptions);
use Sys::Hostname;
use CoGe::Accessory::Web qw(send_email);

our ($names, $limit, $email);

GetOptions(
    "names=s" => \$names,
    "limit=s" => \$limit,
    "email=s" => \$email,
);

die "Missing required arguments" unless ($names && $limit && $email);

$limit =~ s/\%$//; # remove trailing %

foreach my $name (split(',', $names)) {
    my $usage = `df $name | awk 'NR>1 {print \$5}'`;
    chomp($usage);
    $usage =~ s/\%$//; # remove trailing %

    if ($usage > $limit) {
        print "Alert: filesystem $name is $usage% full, sending notification to $email\n";
        send_email(
            from    => 'root@'.hostname,
            to      => $email,
            subject => "Alert: filesystem $name usage $usage% is above $limit% limit",
            body    => "This automated message was sent by the script $0");
    }
}

exit;