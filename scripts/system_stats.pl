#!/usr/bin/env perl
#-------------------------------------------------------------------------------
# Purpose:  Capture system stats for Admin page
# Author:   Matt Bomhoff
# Created:  7/22/15
#------------------------------------------------------------------------------

use warnings;
use strict;

use Sys::Load qw(getload); 
use POSIX qw(strftime); 

# Print timestamp
print strftime("\%H:\%M:\%S \%m/\%d/\%Y", localtime()); 

# Print load averages
my @loads = getload(); 
print "\t", join("\t", @loads);

# Print used memory (in MB)
print "\t", `free -m | awk 'FNR == 3 {printf "%s", \$3}'`; #`free -m | awk 'FNR == 2 {printf "%s", \$3}'`; # mdb changed 11/19/15

# Print disk load average
print "\t", `iostat -c | awk 'FNR == 4 {print \$4}'`;

exit;

#-------------------------------------------------------------------------------
