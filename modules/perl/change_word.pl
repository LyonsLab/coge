#!/usr/bin/perl -w

use strict;
my ($dir, $old_word, $new_word) = @ARGV;
my $TMPDIR = "/tmp";
process_dir($dir);

sub process_dir
  {
    my $dir = shift;
    my @files;
    opendir (DIR, $dir) || return;
    while (my $f = readdir (DIR))
      {
	next if $f =~ /^\.\.?$/;
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
	else
	  {
	    process_file($_);
	  }
      }
  }

sub process_file
  {
    my $file = shift;
    open (OUT, ">$TMPDIR/tmp") || die "Can't open $TMPDIR/tmp for writing: $!";
    open (IN, $file) || warn "Can't open $file for reading: $!";
    while (<IN>)
      {
	s/$old_word/$new_word/g;
	print OUT $_;
      }
    close IN;
    close OUT;
    system "mv $TMPDIR/tmp $file";
  }
