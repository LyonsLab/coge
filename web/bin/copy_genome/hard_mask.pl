#!/usr/bin/perl -w

use strict;

while (<>)
  {
    if (/^>/)
      {
	s/lcl\|//g;
	s/gi\|//g;
	print $_;
      }
    else
      {
	s/[a-z]/X/g;
	print $_;
      }
  }
