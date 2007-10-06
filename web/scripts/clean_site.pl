#!/usr/bin/perl -w

use strict;

my $DIR = shift || ".";
my $new_dir = shift || "new";

mkdir $new_dir unless -d $new_dir;
mkdir "$new_dir/$DIR" unless -d "$new_dir/$DIR";
process_dir($DIR);

sub process_dir
  {
    my $dir = shift;
    my @files;
    opendir (DIR, $dir) || die;
    while (my $file = readdir(DIR))
      {
	next if $file =~ /^\.\.?$/;
	push @files, "$dir/$file";
      }
    closedir DIR;
    foreach my $file (@files)
      {
	if (-d $file)
	  {
	    mkdir "$new_dir/$file";
	    process_dir($file);
	  }
	else
	  {
	    process_file($file);
	  }
      }
  }

sub process_file
  {
    my $file = shift;
    if ($file =~ /\.html$/)
      {
	open (IN, $file) || die "can't open $file for reading";
	open (OUT, ">$new_dir/$file") || die;
	while (<IN>)
	  {
	    s/Â//g;
	    s/â€™/'/g;
	    s/â€œ/"/g;
	    s/â€/"/g;
	    s/>\xa0</><br></g;
	    s/\xa0/ /g;
#	    s/å//g;
#	    s/overflow: hidden/overflow: visible/g;
	    s/http:\/\/\/?CoGe/\/CoGe/ig;
#	    s/http:\/\/\/?coge/\/CoGe/ig;
	    print OUT $_;
	  }
	close IN;
	close OUT;
      }
    else
      {
	`cp "$file" "$new_dir/$file"`;
      }
  }
