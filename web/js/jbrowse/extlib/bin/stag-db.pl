#!/usr/bin/perl -w

eval 'exec /usr/bin/perl -w -S $0 ${1+"$@"}'
    if 0; # not running under some shell

# POD docs at bottom of file

use strict;
use Data::Stag qw(:all);
use Data::Stag::StagDB;
use Getopt::Long;

my $record_type;
my $unique_key;
my $dir;
my $fmt = '';
my $outfmt;
my $help;
my $top;
my $indexfile;
my $qf;
my @query = ();
my $keys;
my $reset;
GetOptions("record_type|r=s"=>\$record_type,
	   "unique_key|unique|u|k=s"=>\$unique_key,
	   "parser|format|p=s"=>\$fmt,
	   "handler|writer|w=s"=>\$outfmt,
	   "indexfile|index|i=s"=>\$indexfile,
	   "top=s"=>\$top,
	   "query|q=s@"=>\@query,
	   "qf=s"=>\$qf,
	   "help|h"=>\$help,
	   "keys"=>\$keys,
	   "reset"=>\$reset,
	  );
if ($help) {
    system("perldoc $0");
    exit 0;
}

my $sdb = Data::Stag::StagDB->new;

$sdb->record_type($record_type) if $record_type;
$sdb->unique_key($unique_key) if $unique_key;
$sdb->indexfile($indexfile) if $indexfile;
$sdb->reset() if $reset;

foreach my $file (@ARGV) {
    my $p;
    if ($file eq '-') {
	$fmt ||= 'xml';
	$p = Data::Stag->parser(-format=>$fmt, -fh=>\*STDIN);
	$p->handler($sdb);
	$p->parse(-fh=>\*STDIN);
    }
    else {
	if (!-f $file) {
	    print "the file \"$file\" does not exist\n";
	}
	$p = Data::Stag->parser($file, $fmt);
	$p->handler($sdb);
	$p->parse($file);
    }
}

if ($qf) {
    open(F, $qf) || die "cannot open queryfile: $qf";
    @query = map {chomp;$_} <F>;
    close(F);
}

if ($keys) {
    my $idx = $sdb->index_hash;
    printf "$_\n", $_ foreach (keys %$idx);
}

if (@query) {
    
    my $w;
    if ($outfmt) {
	$w = Data::Stag->getformathandler($outfmt);
    }
    else {
	$w = Data::Stag->makehandler;
    }
    if ($top) {
	$w->start_event($top);
    }
    my $idx = $sdb->index_hash;
    my $n_found = 0;
    foreach my $q (@query) {
	my $nodes = $idx->{$q} || [];
	if (!@$nodes) {
	    print STDERR "Could not find a record indexed by key: \"$q\"\n";
	    next;
	}
	foreach my $node (@$nodes) {
	    $n_found++;
	    if ($w) {
		$node->sax($w);
	    }
	    else {
		print $node->xml;
	    }
	}
    }
    if ($top) {
	$w->end_event($top);
    }
    if (!$n_found && !$top) {
	print STDERR "NONE FOUND!\n";
    }
    else {
	if (!$outfmt) {
	    print $w->stag->xml;
	}
    }
}


=head1 NAME 

stag-db.pl - persistent storage and retrieval for stag data (xml, sxpr, itext)

=head1 SYNOPSIS

  stag-db.pl -r person -k social_security_no -i ./person-idx myrecords.xml
  stag-db.pl -i ./person-idx -q 999-9999-9999 -q 888-8888-8888

=head1 DESCRIPTION

Builds a simple file-based database for persistent storage and
retrieval of nodes from a stag compatible document.

Imagine you have a very large file of data, in a stag compatible
format such as XML. You want to index all the elements of type
B<person>; each person can be uniquely identified by
B<social_security_no>, which is a direct subnode of B<person>

The first thing to do is to build an index file, which will be stored
in your current directory:

  stag-db.pl -r person -k social_security_no -i ./person-idx myrecords.xml

You can then use the index "person-idx" to retrieve B<person> nodes by
their social security number

  stag-db.pl -i ./person-idx -q 999-9999-9999 > some-person.xml

You can export using different stag formats

  stag-db.pl -i ./person-idx -q 999-9999-9999 -w sxpr > some-person.xml

You can retrieve multiple nodes (although these need to be rooted to
make a valid file)

  stag-db.pl -i ./person-idx -q 999-9999-9999 -q 888-8888-8888 -top personset

Or you can use a list of IDs from a file (newline delimited)

  stag-db.pl -i ./person-idx -qf my_ss_nmbrs.txt -top personset

=head2 ARGUMENTS

=head3 -i INDEXFILE

This file will be used as the persistent index for storage/retrieval

=head3 -r RELATION-NAME

This is the name of the stag node (XML element) that will be stored in
the index; for example, with the XML below you may want to use the
node name B<person> and the unique key B<id>

  <person_set>
    <person>
      <id>...</id>
    </person>
    <person>
      <id>...</id>
    </person>
    ...
  </person_set>

This flag should only be used when you want to store data

=head3 -k UNIQUE-KEY

This node will be used as the unique/primary key for the data

This node should be nested directly below the node that is being
stored in the index - if it is more that one below, specify a path

This flag should only be used when you want to store data

=head3 -u UNIQUE-KEY

Synonym for B<-k>

=head3 -p PARSER

This can be the name of a stag supported format (xml, sxpr, itext) -
XML is assumed by default

It can also be a module name - this module is used to parse the input
file into a stag stream; see L<Data::Stag::BaseGenerator> for details
on writing your own parsers/event generators

This flag should only be used when you want to store data

=head3 -q QUERY-ID

Fetches the relation/node with unique key value equal to query-id

Multiple arguments can be passed by specifying -q multple times

This flag should only be used when you want to query data

=head3 -top NODE-NAME

If this is specified in conjunction with B<-q> or B<-qf> then all the
query result nodes will be nested inside a node with this name (ie
this provides a root for the resulting document tree)

=head3 -qf QUERY-FILE

This is a file of newline-seperated IDs; this is useful for querying
the index in batch

=head3 -keys

This will write a list of all primary keys in the index

=head3 -w WRITER

This format will be used to write the data; can be any stag format
(xml, sxpr, itext) - default XML.

Can also be a module that catches the incoming stag event stream and
does something with it (for example, this could be a module you write
yourself that transforms the stag events into HTML)


=head1 SEE ALSO

L<Data::Stag>

For more complex stag to database mapping, see L<DBIx::DBStag> and the
scripts

L<stag-storenode.pl>

L<selectall_xml>

=cut

