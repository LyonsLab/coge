#!/usr/bin/perl -w

use v5.10;
use strict;
use warnings;

use Getopt::Long;
use JSON qw(decode_json);
use CoGe::Core::Storage qw(add_workflow_result);

our ($USER_NAME, $WORKFLOW_ID, $RESULT_JSON);

GetOptions(
    "user_name=s" => \$USER_NAME,
    "wid=s"       => \$WORKFLOW_ID,
    "result=s"    => \$RESULT_JSON
);

print STDERR $RESULT_JSON;

my $result = decode_json($RESULT_JSON);
unless (add_workflow_result($USER_NAME, $WORKFLOW_ID, $result))
{
    print STDOUT "log: error: could not add workflow result\n";
    exit(-1);
}

exit(0);