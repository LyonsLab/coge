#!/usr/bin/perl -w

use strict;
my $dir = shift;
$dir = "/opt/apache/CoGe/data/diags" unless $dir;
process_dir($dir);

sub process_dir
  {
    my $dir = shift;
    my @dir;
    opendir (DIR, $dir) || die "Can't opendir $dir: $!";
    while (my $item =readdir(DIR))
      {
	next if $item =~ /^\.\.?$/;
	if (-d "$dir/$item")
	  {
	    push @dir, $item;
	  }
	elsif ($item =~ /blast$/)
	  {
	    print "Blast file: $item\n";
	  }
	else
	  {
	    print "File to delete: $item\n";
	    unlink "$dir/$item";
	  }
      }
    closedir DIR;
    foreach my $item (@dir)
      {
	process_dir ("$dir/$item");
      }
  }
