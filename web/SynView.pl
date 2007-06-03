#!/usr/bin/perl -w

use strict;
use CGI;
my $form = new CGI;
my $query = $form->query_string;
#print STDERR $query;

print "Content-Type: text/html\n\n";

print qq{
<head>
<META HTTP-EQUIV="Refresh" CONTENT="0; URL=GEvo.pl?$query">
</head>
};
