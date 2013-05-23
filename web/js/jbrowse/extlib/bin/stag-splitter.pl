#!/usr/bin/perl -w

eval 'exec /usr/bin/perl -w -S $0 ${1+"$@"}'
    if 0; # not running under some shell

use strict;
use Data::Stag qw(:all);
use Data::Stag::XMLParser;
use Getopt::Long;

my $split;
my $name;
my $dir;
my $fmt = '';
my $outfmt = '';
my $help;
GetOptions("split|s=s"=>\$split,
	   "name|n=s"=>\$name,
	   "dir|d=s"=>\$dir,
	   "format|parser|p=s"=>\$fmt,
	   "outformat|writer|w=s"=>\$outfmt,
	   "help|h"=>\$help,
	  );
if ($help) {
    system("perldoc $0");
    exit 0;
}
my $h = Splitter->new;
if ($dir) {
    `mkdir $dir` unless -d $dir;
    $h->dir($dir);
}
if (!@ARGV) {
    die "no files passed!";
}
$h->split_on_element($split);
$h->name_by_element($name);
$h->fmt($outfmt || $fmt);
foreach my $file (@ARGV) {
    my @pargs = (-file=>$file, -format=>$fmt, -handler=>$h);
    if ($file eq '-') {
	if (!$fmt) {
	    $fmt = 'xml';
	}
	@pargs = (-format=>$fmt, -fh=>\*STDIN, -handler=>$h);
    }
#    $split = $split || shift @ARGV;
    Data::Stag->parse(@pargs);
#    $p->handler($h);
#    $p->parse($file);
#    print stag_xml($h->tree);
}

sub usage {
    return <<EOM
EOM
}

package Splitter;
use base qw(Data::Stag::BaseHandler);
use Data::Stag qw(:all);

sub dir {
    my $self = shift;
    $self->{_dir} = shift if @_;
    return $self->{_dir};
}

sub split_on_element {
    my $self = shift;
    $self->{_split_on_element} = shift if @_;
    return $self->{_split_on_element};
}

sub fmt {
    my $self = shift;
    $self->{_fmt} = shift if @_;
    return $self->{_fmt};
}


sub name_by_element {
    my $self = shift;
    $self->{_name_by_element} = shift if @_;
    return $self->{_name_by_element};
}

sub i {
    my $self = shift;
    $self->{_i} = shift if @_;
    return $self->{_i} || 0;
}

sub end_event {
    my $self = shift;
    my $ev = shift;
    if ($ev eq $self->split_on_element) {
	my $topnode = $self->popnode;
	my $name_elt = $self->name_by_element;
	my $name;
	if ($name_elt) {
	    $name = stag_get($topnode, $name_elt);
	}
	if (!$name) {
	    $self->i($self->i()+1);
	    $name = $ev."_".$self->i;
	}
	my $dir = $self->dir || '.';
	my $fmt = $self->fmt;
	$fmt = $fmt || 'xml';
	$name = safe($name);
	my $fh = FileHandle->new(">$dir/$name.$fmt") || die("can't open >$dir/$name.$fmt");
	if ($fmt eq 'xml') {
	    my $out = stag_xml($topnode);
	    print $fh $out;
	}
	else {
	    stag_generate($topnode, -fh=>$fh, -fmt=>$fmt);
	}
	$fh->close;
	return [];
    }
    else {
	return $self->SUPER::end_event($ev, @_);
    }
}

sub safe {
    my $n = shift;
    $n =~ s/\///g;
    $n;
}

1;

__END__

=head1 NAME 

stag-splitter.pl - splits a stag file into multiple files

=head1 SYNOPSIS

  stag-splitter.pl -split person -name social_security_no file.xml

=head1 DESCRIPTION

Splits a file using a user specified parser (default xml) around a
specified split node, naming each file according to the name argument

the files will be named anonymously, unless the '-name' switch is specified; this will use the value of the specified element as the filename

eg; if we have


  <top>
    <a>
      <b>foo</b>
      <c>yah</c>
      <d>
        <e>xxx</e>
      </d>
    </a>
    <a>
      <b>bar</b>
      <d>
        <e>wibble</e>
      </d>
    </a>
  </top>

if we run

  stag-splitter.pl -split a -name b

it will generate two files, "foo.xml" and "bar.xml"

input format can be 'xml', 'sxpr' or 'itext' - if this is left blank
the format will be guessed from the file suffix

the output format defaults to the same as the input format, but
another can be chosen.

files go in the current directory, but this can be overridden with the
'-dir' switch

=head1 USAGE

   stag-splitter.pl [-split <ELEMENT-NAME>] [-name <ELEMENT-NAME>] [-dir <DIR>] [-format <INPUT-FORMAT>] [-outformat <OUTPUT-FORMAT>] <FILENAMES>

=over

=item -p|parser FORMAT

FORMAT is one of xml, sxpr or itext, or the name of a perl module

xml assumed as default

=item -w|writer FORMAT

FORMAT is one of xml, sxpr or itext, or the name of a perl module

=item -split|s NODE

node to split on

=item -name|n NODE

field/element to use when naming files

will use surrogate IDs if this argument not specified

=item -dir|d DIR

write files to this directory

=back


=cut
