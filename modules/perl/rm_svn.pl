#!/usr/bin/perl -w

use strict;
my ($dir) = @ARGV;
process_dir($dir);

sub process_dir
  {
    my $dir = shift;
    my @files;
    opendir (DIR, $dir) || return;
    while (my $f = readdir (DIR))
      {
	next if $f =~ /^\.\.?$/;
	system "rm -rf $dir/$f" if $f =~ /^\.svn$/;
	next if $f =~ /^\.svn$/;
	push @files,"$dir/$f";
      }
    closedir (DIR);
    foreach (@files)
      {
	if (-d $_)
	  {
	    process_dir($_);
	  }
      }
  }
