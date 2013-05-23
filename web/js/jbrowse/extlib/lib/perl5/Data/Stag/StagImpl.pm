# $Id: StagImpl.pm,v 1.66 2008/06/03 17:31:15 cmungall Exp $
#
# Author: Chris Mungall <cjm@fruitfly.org>
#
# You may distribute this module under the same terms as perl itself

package Data::Stag::StagImpl;

=head1 NAME

  Data::Stag::StagImpl

=head1 SYNOPSIS

  use Data::Stag qw(:all);

=head1 DESCRIPTION

This is the default implementation for Data::Stag - please see L<Data::Stag>

=cut

use FileHandle;
use IO::String;
use Carp;
use strict;
use vars qw($AUTOLOAD $DEBUG);
use Data::Stag::Base;
use Data::Stag::Util qw(rearrange);
use base qw(Data::Stag::StagI);

use vars qw($VERSION);
$VERSION="0.11";


sub new {
    my $proto = shift; 
    my $class = ref($proto) || $proto;
    if (@_ == 2 || @_ == 0) {
        return bless [@_], $class;
    }
    else {
        confess @_;
    }
}

sub unflatten {
    my $proto = shift; 
    my $class = ref($proto) || $proto;
    my ($name, $flist) = @_;
    my @uflist = ();
    if (!ref($flist)) {
	return $class->new($name=>$flist);
    }
    if (ref($flist) eq 'HASH') {
	# unpack hash into array
	$flist = [%$flist];
    }
    if (ref($flist) ne 'ARRAY') {
	confess("$name => $flist not array");
    }
    while (@$flist) {
	my $k = shift @$flist;
        my $v = shift @$flist;
	if (ref($v)) {
	    push(@uflist,
		 $class->unflatten($k=>$v));
	}
	else {
	    push(@uflist,
		 [$k=>$v]);
	}
    }
    return $class->new($name=>[@uflist]);
}

sub unhash {
    my $proto = shift; 
    my $class = ref($proto) || $proto;
    my %hash = @_;

    my @tags = ();
    foreach my $k (keys %hash) {
	my $v = $hash{$k};
	if (ref($v)) {
	    if (ref($v) eq 'ARRAY') {
		push(@tags, [$k=>$_]) foreach @$v;
	    }
	    elsif (ref($v) eq 'HASH') {
		my $stag = unhash($class, %$v);
		push(@tags, [$k=>$stag->data]);
	    }
	    else {
		confess("cannot unhash $v");
	    }
	}
	else {
	    push(@tags, [$k=>$v]);
	}
    }
    return $class->new(stag=>[@tags]);
}

sub unstone {
    my $tree = shift; 
    my $stone = shift;
    my $xml = $stone->asXML;
    from($tree, xmlstr=>$xml);
}

sub stone {
    my $tree = shift;
    my %h = hash($tree);
    load_module("Stone");
    Stone->new(%h);
}

sub load_module {

    my $classname = shift;
    confess("no class") unless $classname;
    my $mod = $classname;
    $mod =~ s/::/\//g;

    if (!$mod) {
        confess("must supply as mod as argument");
    }

    if ($main::{"_<$mod.pm"}) {
    }
    else {
        require "$mod.pm";
    }
}

# -------------------------------------
# ----   INPUT/OUTPUT FUNCTIONS   -----
# -------------------------------------


sub parser {
    my $tree = shift;
    my ($fn, $fmt, $h, $eh, $str, $fh) = 
      rearrange([qw(file format handler errhandler str fh)], @_);

    if ($fn && $fn eq '-') {
	$fn = '';
	$fh = \*STDIN;
    }
    # GUESS FORMAT BASED ON FILENAME
    if (!$fmt && $fn) {
	if ($fn =~ /\.xml$/) {
            $fmt = "xml";
        }
        elsif ($fn =~ /\.indent$/) {
            $fmt = "indent";
        }
        elsif ($fn =~ /\.ind$/) {
            $fmt = "indent";
        }
        elsif ($fn =~ /\.xtc$/) {
            $fmt = "indent";
        }
        elsif ($fn =~ /\.ite?xt$/) {
            $fmt = "itext";
        }
        elsif ($fn =~ /\.se?xpr$/) {
            $fmt = "sxpr";
        }
        elsif ($fn =~ /\.el$/) {
            $fmt = "sxpr";
        }
        elsif ($fn =~ /\.pl$/) {
            $fmt = "perl";
        }
        elsif ($fn =~ /\.perl$/) {
            $fmt = "perl";
        }
        else {
            # default to xml
	     if (!$str && $fn) {
		 my $fh = FileHandle->new($fn) || 
                   confess("no such file $fn");
		 # get the first line
		 $str = <$fh>;
		 chomp $str;
		 $fh->close;
	     }
	     else {
		 $fmt = "xml";
	     }
        }
    }

    # GUESS FORMAT BASED ON STR
    if (!$fmt && $str) {
        if ($str =~ /^\s*\'/) {
            $fmt = "sxpr";
        }
        elsif ($str =~ /^\s*\(/) {
            $fmt = "sxpr";
        }
        elsif ($str =~ /^\s*\;/) {
            $fmt = "sxpr";
        }
        elsif ($str =~ /^\s*\</) {
            $fmt = "xml";
        }
        elsif ($str =~ /^\s*\w+\:\s/) {
            $fmt = "itext";
        }
        elsif ($str =~ /^\s*\#/) {
            $fmt = "indent";
        }
        elsif ($str =~ /^\w+/) {
            $fmt = "indent";
        }
        else {
        }
    }

    my $parser;
    if (!ref($fmt) && $fmt =~ /::/) {
        load_module($fmt);
        $fmt = $fmt->new;
    }

    if (ref($fmt)) {
        $parser = $fmt;
    }
    elsif ($fmt eq "xml") {
        $parser = "Data::Stag::XMLParser";
    }
    elsif ($fmt eq "indent") {
        $parser = "Data::Stag::IndentParser";
    }
    elsif ($fmt eq "itext") {
        $parser = "Data::Stag::ITextParser";
    }
    elsif ($fmt eq "perl") {
        $parser = "Data::Stag::PerlParser";
    }
    elsif ($fmt eq "sxpr") {
        $parser = "Data::Stag::SxprParser";
    }
    else {
	confess("cannot guess parser from fmt=\"$fmt\" @_") unless $parser;
    }
    unless (ref($parser)) {
	load_module($parser);
	$parser = $parser->new;
    }
    $parser->file($fn) if $fn;
#    $parser->fh($fh) if $fh;
    return $parser;
}

sub parse {
    my $tree = shift;
    my ($fn, $fmt, $h, $eh, $str, $fh) = 
      rearrange([qw(file format handler errhandler str fh)], @_);

    if (!$tree || !ref($tree)) {
        $tree = [];
    }

    my $p = parser($tree, @_);
    $h = Data::Stag::Base->new unless $h;
    if (!$eh) {
	$eh = getformathandler($tree, 'xml');
	$eh->fh(\*STDERR);
    }
    $p->handler($h);
    $p->errhandler($eh);
    $p->parse(
              -file=>$fn,
              -str=>$str,
              -fh=>$fh,
             );
    my $htree = $h->can("tree") ? $h->tree : [];
    Nodify($tree);
    @$tree = @{$htree || []};
    return $tree;
}
*parseFile = \&parse;

sub parsestr {
    my $tree = shift;
    my ($str, $fmt, $h, $eh) = 
      rearrange([qw(str format handler errhandler)], @_);
    return 
      $tree->parse(-str=>$str,
		   -format=>$fmt,
		   -handler=>$h,
		   -errhandler=>$eh);
		 
}
*parseStr = \&parsestr;

sub from {
    my $class = shift;
    my ($fmt, $file, $str) = 
      rearrange([qw(fmt file str)], @_);
    if ($fmt eq 'xmlstr') {
        return xmlstr2tree($file);
    }
    elsif ($fmt =~ /(.*)str/) {
        $fmt = $1;
        return parse([], 
                     -str=>$file,
                     -format=>$fmt);
        
    }
    elsif ($fmt eq 'xml') {
        return xml2tree($file);
    }
    else {
        return parse([], 
                     -file=>$file,
                     -format=>$fmt);
    }
}

sub _gethandlerobj {
    my $tree = shift || [];
    my ($fn, $fmt, $fh) = 
      rearrange([qw(file fmt fh)], @_);
    if (!$fmt) {
	if (!$fn) {
	    $fmt = "xml";
	}
        elsif ($fn =~ /\.xml$/) {
            $fmt = "xml";
        }
        elsif ($fn =~ /\.ite?xt$/) {
            $fmt = "itext";
        }
        elsif ($fn =~ /\.ind/) {
            $fmt = "indent";
        }
        else {
        }
    }
    my $writer;
    if (ref($fmt)) {
        return $fmt;
    }
    elsif ($fmt =~ /xml/i) {
        $writer = "Data::Stag::XMLWriter";
    }
    elsif ($fmt =~ /itext/i) {
        $writer = "Data::Stag::ITextWriter";
    }
    elsif ($fmt =~ /indent/i) {
        $writer = "Data::Stag::IndentWriter";
    }
    elsif ($fmt =~ /sxpr/i) {
        $writer = "Data::Stag::SxprWriter";
    }
    elsif ($fmt =~ /yaml/i) {
        $writer = "Data::Stag::YAMLWriter";
    }
    elsif ($fmt =~ /perl/i) {
        $writer = "Data::Stag::PerlWriter";
    }
    elsif ($fmt =~ /dtd/i) {
        $writer = "Data::Stag::DTDWriter";
    }
    elsif ($fmt =~ /simple/i) {
        $writer = "Data::Stag::Simple";
    }
    elsif ($fmt =~ /xslt\/(.+)/i) {
	my $xslt_file = $1;
        $writer = "Data::Stag::XSLTHandler";
	load_module($writer);
	my $w = $writer->new(-file=>$fn, -fh=>$fh);
	$w->xslt_file($xslt_file);
	return $w;
    }
    elsif ($fmt =~ /::/) {
        $writer = $fmt;
    }
    elsif (!$fmt) {
	confess("no format/writer!");
    }
    else {
	confess("unrecognised:$fmt");
    }
    load_module($writer);

    my $w = $writer->new(-file=>$fn, -fh=>$fh);
    return $w;
}
*findhandler = \&_gethandlerobj;

sub getformathandler {
    my $tree = shift;
    my $fmt = shift;
    return findhandler($tree, -fmt=>$fmt);
}

sub generate {
    my $tree = shift || [];
    my $w = _gethandlerobj($tree, @_);
    $w->is_buffered(1);
    $w->event(@$tree);
    return $w->popbuffer || '';
}
*gen = \&generate;

sub write {
    my $tree = shift || [];
    my $w = _gethandlerobj($tree, @_);

    $w->is_buffered(0);
    $w->event(@$tree);
    $w->close_fh;
    return;
}
    
sub makexslhandler {
    my $tree = shift;
    my $xslt_file = shift;
    load_module("Data::Stag::XSLTHandler");
    my $handler = Data::Stag::XSLTHandler->new;
    $handler->xslt_file($xslt_file);
    return $handler;
}

sub makehandler {
    my $tree = shift;
    my $handler;
    if (@_ == 1) {
	my $module = shift;
        if ($module =~ /\.xsl/) {
            $handler = makexslhandler($tree, $module);
        }
        else {
            load_module($module);
            $handler = $module->new;
        }
    }
    else {
	my %trap_h = @_;
	my %opt_h = ();
	%trap_h =
	  map {
	      if ($_ =~ /^-(.*)/) {
		  $opt_h{lc($1)} = $trap_h{$_};
		  ();
	      }
	      else {
		  ($_ => $trap_h{$_})
	      }
	  } keys %trap_h;
	load_module("Data::Stag::BaseHandler");
	$handler = Data::Stag::BaseHandler->new;
	$handler->trap_h(\%trap_h);
	
	if ($opt_h{notree}) {
	    load_module("Data::Stag::null");
	    my $null = Data::Stag::null->new;
	    my $ch = Data::Stag->chainhandlers([keys %trap_h],
					       $handler,
					       $null);
	    return $ch;
	}
    }
    return $handler;
}
*mh = \&makehandler;

sub chainhandlers {
    my $tree = shift;
    my $block = shift;
    my @sh = @_;

    load_module("Data::Stag::ChainHandler");
    my $handler = Data::Stag::ChainHandler->new;
    $handler->subhandlers([
                           map {
                               if (ref($_)) {
                                   if (ref($_) eq 'HASH') {
                                       # make a new handler
                                       makehandler($tree, %$_);
                                   }
                                   else {
                                       # assume it is an object
                                       $_;
                                   }
                               }
			       elsif (!$_) {
				   ()
			       }
                               else {
                                   # assume it is string specifying format
                                   _gethandlerobj($tree, -fmt=>$_)
                               }
                             } @sh
                          ]);
    $handler->blocked_event($block);

    # if no explicit blocked events set, then introspect
    # the subhandlers to see if they declare what they emit
    if (ref($block) && !@$block) {
        my @emits = map {$_->CONSUMES} @{$handler->subhandlers};
        $handler->blocked_event(\@emits);
    }
    return $handler;
}

sub transform {
    my $tree = shift;
    my @T = @_;
    my %trap_h = 
      map {
	  my ($from, $to) = @$_;
	  $from=> sub {
	      my $self = shift;
	      my $stag = shift;
#	      print STDERR "Transforming $from => $to\n";
#	      print STDERR $stag->sxpr;
	      my $data = $stag->data;
	      my @path = splitpath($to);
	      my $node = [];
	      my $p = $node;
	      while (@path) {
		  my $elt = shift @path;
		  $p->[0] = $elt;
		  if (@path) {
		      my $newpath = [];
		      $p->[1] = [$newpath];
		      $p = $newpath;
		  }
		  else {
		      $p->[1] = $data;
		  }
	      }
#	      @$stag = @$node;
#	      print STDERR $stag->sxpr;
	      return $node;
#	      return 0;
	  }
      } @T;
    load_module("Data::Stag::BaseHandler");
    my $handler = Data::Stag::BaseHandler->new;
    $handler->trap_h(\%trap_h);
    $tree->events($handler);
    my $nu = $handler->stag;
    @$tree = @$nu;
    return;
}
*t = \&transform;

# transform stag into hash datastruct;
# stag keys become hash keys (unordered)
# single valued keys map to single value (itself a hash or primitive)
# multivalued map to arrayrefs
sub hash {
    my $tree = shift;
    my ($ev, $subtree) = @$tree;

    # make sure we have non-terminal
    if (ref($subtree)) {
	# make hash using stag keys
        my %h = ();
	foreach my $subnode (@$subtree) {
	    my $k = $subnode->[0];
	    my $v;
	    
	    # terminals map to data value; non-terminals
	    # get recursively mapped to hashes
	    if (isterminal($subnode)) {
		$v = $subnode->[1];
	    }
	    else {
		$v = {hash($subnode)};
	    }

	    # determine if it is single-valued or multi-valued hash
	    my $curr = $h{$k};
	    if ($curr) {
		if (ref($curr) && ref($curr) eq 'ARRAY') {
		    push(@$curr, $v);
		}
		else {
		    $h{$k} = [$curr, $v];
		}
	    }
	    else {
		# if there is only one value, don't use array -- yet
		$h{$k} = $v;
	    }
        }
	return %h;
    }
    else {
	warn("can't make a hash from a terminal");
        return ();
    }
}
*tree2hash = \&hash;

# this FLATTENS everything
# any non terminal node is flattened and lost
# does not check for duplicates!
sub pairs {
    return @{(_pairs(@_))[1]};
}
*tree2pairs = \&pairs;
sub _pairs {
    my $tree = shift;
    my ($ev, $subtree) = @$tree;
    if (ref($subtree)) {
        my @pairs = map { _pairs($_) } @$subtree;
        return $ev=>[@pairs];
    }
    else {
        return $ev=>$subtree;
    }
}

# PRIVATE
sub tab {
    my $t=shift;
    return "  " x $t;
}

sub dlist {
    my $tree = shift;
    my $del = shift || '/';
    _dlist($tree, $del, $del, 1);
}

sub _dlist {
    my $tree = shift;
    my $del = shift || '/';
    my $root = shift;
    my $nextid = shift;
    my @kids = $tree->kids;
    if (isterminal($tree)) {
        my $data = $tree->data;
        $data =~ s/$del/\\$del/g;
        #        my $stem = $root . $tree->element . "[$nextid]";
        #        $root = $stem . $del;
        #        return ($root, $root . $del . $data);

        return ($root . $tree->element . ': ' . $data . "[$nextid]");
    }
    else {
        my $id = 1;
        my $stem = $root . $tree->element . "[$nextid]";
        $root = $stem . $del;
        my @dlist =
          map {
              _dlist($_, $del, $root, $id++);
          } @kids;
        return ($stem, @dlist);
    }
#    return (@dlist);
}

sub xml {
    my $tree = shift;
    my $fn = shift;
    generate($tree, $fn, 'xml', @_);
}
*tree2xml = \&xml;

sub sxpr {
    my $tree = shift;
    my $fn = shift;
    generate($tree, $fn, 'sxpr', @_);
}

sub itext {
    my $tree = shift;
    my $fn = shift;
    generate($tree, $fn, 'itext', @_);
}

sub indent {
    my $tree = shift;
    my $fn = shift;
    generate($tree, $fn, 'indent', @_);
}

sub as {
    my $tree = shift;
    my $fmt = shift;
    my $fn = shift;
    generate($tree, $fn, $fmt, @_);
}

sub perldump {
    my $tree = shift;
    return 
      _tree2perldump($tree, 1) . ";\n";
}
*tree2perldump = \&perldump;

sub _tree2perldump {
    my $tree = shift;
    my $indent = shift || 0;
    my ($ev, $subtree) = @$tree;
    if (ref($subtree)) {
	$indent++;
        return 
          sprintf("%s[ '$ev' => [\n%s%s]]",
#                  tab($indent++),
		  "",
                  join(",\n", map {
		      tab($indent) .
			_tree2perldump($_, $indent)
		  } @$subtree),
#                  tab($indent-1),
		  "",
                 );
    }
    else {
        return
          sprintf("%s[ '$ev' => %s ]",
#                  tab($indent),
		  "",
                  perlesc($subtree) || "");
    }
}

sub perlesc {
    my $val = shift;
    return "undef" if !defined $val;
    $val =~ s/\'/\\\'/g;
    return "'$val'";
}

sub sax {
    my $tree = shift;
    my $saxhandler = shift;
    $saxhandler->start_document;
    _tree2sax($tree, $saxhandler);
    $saxhandler->end_document;
}
*tree2sax = \&sax;

sub _tree2sax {
    my $tree = shift;
    my $saxhandler = shift;
    my ($ev, $subtree) = @$tree;
    $saxhandler->start_element({Name => $ev});
    if (ref($subtree)) {
        map { _tree2sax($_, $saxhandler) } @$subtree;
    }
    else {
        $saxhandler->characters({Data => $subtree});
        
    }
    $saxhandler->end_element({Name => $ev});
}

sub xslt {
    my $tree = shift;
    my $xsltstr = xsltstr($tree,@_);
    return parsestr($tree,
                    -str=>$xsltstr,
                    -format=>'xml');
}

sub xsltstr {
    my $stag = shift;
    my $xslt_file = shift;
    
    load_module("XML::LibXML");
    load_module("XML::LibXSLT");
    my $parser = XML::LibXML->new();
    my $source = $parser->parse_string($stag->xml);
    
    my $xslt = XML::LibXSLT->new();
    my $styledoc = $parser->parse_file($xslt_file);
    my $stylesheet = $xslt->parse_stylesheet($styledoc);

    my $results = $stylesheet->transform($source);
    return $results->toString;
}


sub events {
    my $tree = shift;
    my $handler = shift;
    $handler->event(@$tree);
}

sub xmlesc {
    my $word = shift;
    return '' unless defined $word;
    $word =~ s/\&/\&amp;/g;
    $word =~ s/\</\&lt;/g;
    $word =~ s/\>/\&gt;/g;
    return $word;
    
}

sub xml2tree {
    my $file = shift;
    my $handler = Data::Stag::Base->new;
    load_module("XML::Parser::PerlSAX");
    my $parser = XML::Parser::PerlSAX->new();
    my %parser_args = (Source => {SystemId => $file},
                       Handler => $handler);
    $parser->parse(%parser_args);
    return Node(@{$handler->tree});
}

sub xmlstr2tree {
    my $str = shift;
    my $handler = Data::Stag::Base->new;
    load_module("XML::Parser::PerlSAX");
    my $parser = XML::Parser::PerlSAX->new();
    my %parser_args = (Source => {String => $str},
                       Handler => $handler);
    $parser->parse(%parser_args);
    return Node(@{$handler->tree});
}

# WANTARRAY
sub sxpr2tree {
    my $sxpr = shift;
    my $indent = shift || 0;
    my $i=0;
    my @args = ();
    while ($i < length($sxpr)) {
        my $c = substr($sxpr, $i, 1);
#################################        print STDERR "c;$c i=$i\n";
        $i++;
        if ($c eq ')') {
            my $funcnode = shift @args;
#            print STDERR "f = $funcnode->[1]\n";
            map {print xml($_)} @args;
            return [$funcnode->[1] =>[@args]], $i;
        }
        if ($c =~ /\s/) {
            next;
        }
        if ($c eq '(') {
            my ($tree, $extra) = sxpr2tree(substr($sxpr, $i), $indent+1);
            push(@args, $tree);
            $i += $extra;
#            printf STDERR "tail: %s\n", substr($sxpr, $i);
        }
        else {
            # look ahead
            my $v = "$c";
            my $p=0;
            while ($i+$p < length($sxpr)) {
                my $c = substr($sxpr, $i+$p, 1);
                if ($c =~ /\s/) {
                    last;
                }
                if ($c eq '(') {
                    last;
                }
                if ($c eq ')') {
                    last;
                }
                $p++;
                $v.= $c;
            }
            $i+=$p;
            push(@args, [arg=>$v]);
        }
    }
#    map {print xml($_)} @args;
    map {Nodify($_)} @args;
    if (wantarray) {
        return @args;
    }
    
    return $args[0];
}

sub addkid {
    my $tree = shift;
    my $newtree = shift;
    confess("problem: $tree not arr") unless ref($tree) && ref($tree) eq "ARRAY" || isastag($tree);
    my ($ev, $subtree) = @$tree;
    if (ref($newtree) && $newtree->[0] eq '@') {
        unshift(@$subtree, $newtree);
    }
    else {
        push(@$subtree, $newtree);
    }
    $newtree;
}
*addChildTree = \&addkid;
*addchild = \&addkid;
*ak = \&addkid;

#sub findChildren {
#    my $tree = shift;
#    confess("problem: $tree not arr") unless ref($tree) && ref($tree) eq "ARRAY" || isaNode($tree);
#    my ($ev, $subtree) = @$tree;
#    return @$subtree;
#}

sub findnode {
    my $tree = shift;
    my ($node, @path) = splitpath(shift);

    my $replace = shift;
    confess("problem: $tree not arr") unless ref($tree) && ref($tree) eq "ARRAY" || isastag($tree);

    if (@path) {
        my @r = map { $_->findnode(\@path, $replace) } findnode($tree, $node);        
        return @r;
    }

    my ($ev, $subtree) = @$tree;
    my @r = ();
    if ($DEBUG) {
        print STDERR "$ev, $subtree;; replace = $replace\n";
    }
    if (test_eq($ev, $node)) {
        if ($DEBUG) {
            print STDERR "  MATCH\n";
        }
        if (defined $replace) {
            my @old = @$tree;
            @$tree = @$replace;
            return Nodify([@old]);
        }
#        return [$ev=>$subtree] ;
#        return Nodify($tree);
	@r = (Nodify($tree));
    }
    my @nextlevel = ();
    if ( ref($subtree)) {
	@nextlevel =
	  map { 
	      map {
		  Nodify($_)
	      } findnode($_, $node, $replace);
	      
	  } @$subtree;
	# get rid of empty nodes
	# (can be caused by replacing)
	@$subtree = map { ref($_) && !scalar(@$_) ? () : $_ } @$subtree;
    }
#    if (wantarray) {
#        return $nextlevel[0];
#    }
    return (@r, @nextlevel);
}
*fn = \&findnode;
*findSubTree = \&findnode;
*fst = \&findnode;

sub remove {
    my $tree = shift;
    my ($node, @path) = splitpath(shift);

    my $replace = shift;
    confess("problem: $tree not arr") unless ref($tree) && ref($tree) eq "ARRAY" || isastag($tree);

    if (@path) {
        $_->remove(\@path) foreach findnode($tree, $node);        
	return;
    }

    my ($ev, $subtree) = @$tree;
    return unless ref $subtree;
    my @subnodes = @$subtree;
    @subnodes =
      grep {
	  $_->[0] ne $node;
      } @subnodes;
    @$subtree = @subnodes;
    remove($_, $node) foreach @subnodes;
    return;
}

sub set {
    my $tree = shift || confess;
    my ($node, @path) = splitpath(shift);
    my @replace = @_;

    if (@path) {
        my $last = pop @path;
        my @nodes = getnode($tree, [$node, @path]);
        foreach (@nodes) {
            set($_, $last, @replace);
        }
        return @replace;
    }

    confess("problem: $tree not arr") unless ref($tree) && ref($tree) eq "ARRAY" || isastag($tree);
    my ($ev, $subtree) = @$tree;
    confess("$subtree not arr [$ev IS A TERMINAL NODE!!]") unless ref($subtree);
    my $is_set;
    my @nu = ();
    foreach my $st (@$subtree) {
        push(@nu, $st);
        my ($ev, $subtree) = @$st;
        if (test_eq($ev, $node)) {
            if (!$is_set) {
                pop @nu;
                push(@nu,
                     map {
                         [$node => $_]
                     } @replace);
                $is_set = 1;
            }
            else {
                pop @nu;
            }
        }
    }
    
    # place at the end if not already present
    if (!$is_set) {
        map {
            addkid($tree, [$node=>$_]);
        } @replace;
    }
    else {
        @$tree = ($ev, \@nu);
    }
    return @replace;
}
*s = \&set;
*setSubTreeVal = \&set;

sub setl {
    my $tree = shift || confess;
    my @args = @_;
    while (@args) {
	set($tree, splice(@args, 0, 2));
    }
    return;
}
*sl = \&setl;
*setlist = \&setl;

sub setnode {
    my $tree = shift;
    my $elt = shift;
    my $nunode = shift;
    set($tree, $elt, data($nunode)); 
}
*sn = \&setnode;
*setn = \&setnode;
*settree = \&setnode;

# EXPERIMENTAL
sub free {
    my $tree = shift;
    @$tree = ();
}

sub add {
    my $tree = shift || confess;
    my $node = shift;
    # usage1: stag_add($tree, $name, @nodes)
    my @v = @_; # nodes to be added
    if (ref($node)) {
        # usage2: stag_add($tree, $node)
	if ($node->isnull) {
	    confess("cannot add null node");
	}
        # split node into name and data
        ($node, @v) = ($node->[0], [$node->[1]]);
    }
    if (ref($v[0]) && !ref($v[0]->[0])) {
	@v = map { $_->[1] } @v;
    }
    confess("problem: $tree not arr") unless ref($tree) && ref($tree) eq "ARRAY" || isastag($tree);
    my ($ev, $subtree) = @$tree;

    my @nu_subtree = ();
    my $has_been_set = 0;
    for (my $i=0; $i<@$subtree; $i++) {
	my $st = $subtree->[$i];
	my $next_st = $subtree->[$i+1];
        my ($ev, $subtree) = @$st;
	push(@nu_subtree, $st);
        if (!$has_been_set &&
            test_eq($ev, $node) &&
	    (!$next_st ||
	     $next_st->[0] ne $ev)) {
	    push(@nu_subtree, 
		 map { [$ev=>$_] } @v);
	    $has_been_set = 1;
	}
    }
    if (!$has_been_set) {
	addkid($tree, [$node=>$_]) foreach @v; 
    }
    else {
	@$subtree = @nu_subtree;
    }
    return;
}
*a = \&add;
*addSubTreeVal = \&add;

sub addnode {
    my $tree = shift;
    my $elt = shift;
    my $node = shift;
    my $nodename = $node->name;
    if ($nodename ne $elt) {
        confess("$nodename ne $elt");
    }
    add($tree, $elt, $node->data);
}
*an = \&addnode;
*addn = \&addnode;

sub unset {
    my $tree = shift || confess;
    my ($node, @path) = splitpath(shift);
    confess("problem: $tree not arr") unless ref($tree) && ref($tree) eq "ARRAY" || isastag($tree);

    if (@path) {
	$_->unset(\@path) foreach findnode($tree, $node);
	return;
    }

    my ($ev, $subtree) = @$tree;
    return unless ref $subtree;
    my @nu_tree = ();
    foreach my $st (@$subtree) {
        my ($ev, $subtree) = @$st;
        if ($ev ne $node) {
            push(@nu_tree, $st);
        }
    }
    $tree->[1] = \@nu_tree;
    return;
}
*unsetSubTreeVal = \&unset;
*u = \&unset;


# WANTARRAY
sub find {
    my $tree = shift || confess;
    my ($node, @path) = splitpath(shift);
    my $replace = shift;

    confess("problem: $tree not arr") unless (ref($tree) && ref($tree) eq "ARRAY") || isastag($tree);

    my @r = ();
    if (@path) {
        @r = map { $_->find(\@path, $replace) } findnode($tree, $node)
    }
    else {
        my ($ev, $subtree) = @$tree;
        if (test_eq($ev, $node)) {
            my $is_nt = ref($subtree);
            if (defined $replace)  {
                if ($is_nt) {
                    confess("use findval") unless ref($replace);
                    @$tree = @$replace;
                }
                else {
                    $tree->[1] = $replace;
                }
            }
            return $is_nt ? Nodify($tree) : $subtree;
#	    @r = ($is_nt ? Nodify($tree) : $subtree);
        }
        return unless ref($subtree);
        push(@r, map { find($_, $node, $replace) } @$subtree);
    }
    if (wantarray) {
        return @r;
    }
    $r[0];
}
*f = \&find;

sub findval {
    my $tree = shift || confess;
    my ($node, @path) = splitpath(shift);
    my $replace = shift;

    confess("problem: $tree not arr") unless ref($tree) && ref($tree) eq "ARRAY" || isastag($tree);

    my @r = ();
    if (@path) {
        @r = map { $_->findval(\@path, $replace) } findnode($tree, $node)
    }
    else {
        my ($ev, $subtree) = @$tree;
        if (test_eq($ev, $node)) {
	    my $dataref = \$tree->[1];
	    if (ref($subtree)) {
		# check if it is the data node of
		# an element with attributes
		my @kids = grep {$_->[0] eq '.'} @$subtree;
		if (@kids == 1) {
		    $dataref = \$kids[0]->[1];
		}
	    }
            if (defined $replace)  {
                $$dataref = $replace;
            }
            return $$dataref;
        }
        return unless ref($subtree);
        @r = map { findval($_, $node, $replace) } @$subtree;
    }
    if (wantarray) {
        return @r;
    }
    $r[0];
}
*fv = \&findval;
*finddata = \&findval;
*fd = \&findval;
*findSubTreeVal = \&findval;

# WANTARRAY
sub getdata {
    my $tree = shift || confess;
    my ($node, @path) = splitpath(shift);
    my $replace = shift;

    confess("problem: $tree not arr") unless ref($tree) && ref($tree) eq "ARRAY" || isastag($tree);

    my @v = ();
    if (@path) {
        @v = map { $_->getdata(\@path, $replace) } getnode($tree, $node)
    }
    else {

        my ($top_ev, $children) = @$tree;
        @v = ();
        foreach my $child (@$children) {
            confess unless ref $child;
            my ($ev, $subtree) = @$child;
            if (test_eq($ev, $node)) {
                if (defined $replace)  {
                    $child->[1] = $replace;
                }
                push(@v, $subtree);
            }
        }
    }
    if (wantarray) {
        return @v;
    }
    $v[0];
}
*gd = \&getdata;

sub get {
    my $tree = shift || confess;
    my ($node, @path) = splitpath(shift);
    my $replace = shift;
#    confess("problem: $tree not arr") unless ref($tree) && (ref($tree) eq "ARRAY" || isastag($tree));
    if (!ref($tree->[1])) {
        # terminal node - always returns undef
        return;
    }
    my @v = ();
    if (@path) {
        @v = map { $_->get(\@path, $replace) } getnode($tree, $node)
    }
    else {

        my ($top_ev, $children) = @$tree;
        @v = ();
	if (!ref($children)) {
	    confess("problem with $node/$top_ev => $children; cannot call get on terminal node");
	}
        foreach my $child (@$children) {
            confess unless ref $child;
            my ($ev, $subtree) = @$child;
            if (test_eq($ev, $node)) {
                my $is_nt = 0;
                if (ref($subtree)) {
                    $is_nt = 1;
                }
                if (defined $replace)  {
#                    $tree->[1] = $replace;
                    if ($is_nt) {
                        if (!ref($replace)) {
                            confess("use getdata instead");
                        }
                        @$child = @$replace;
                    }
                    else {
                        $child->[1] = $replace;
                    }
                }
                if ($is_nt) {
                    push(@v, Nodify($child));
                }
                else {
                    push(@v, $subtree);
                }
            }
        }
    }
    if (wantarray) {
        return @v;
    }
    $v[0];
}
*g = \&get;

# WANTARRAY
sub getnode {
    my $tree = shift || confess;
    my ($node, @path) = splitpath(shift);
    my $replace = shift;

    confess("problem: $tree not arr") unless ref($tree) && ref($tree) eq "ARRAY" || isastag($tree);

    my @v = ();
    if (@path) {
        @v = map { $_->getnode(\@path, $replace) } getnode($tree, $node)
    }
    else {

        my ($top_ev, $children) = @$tree;
        if (!ref($children)) {
            confess("problem: $top_ev => \"$children\" not a tree");
        }
        foreach my $child (@$children) {
            my ($ev, $subtree) = @$child;
            if (test_eq($ev, $node)) {
                if (defined $replace)  {
                    $tree->[1] = $replace;
                }
                push(@v, Nodify($child));
            }
        }
    }

    if (wantarray) {
        return @v;
    }
    $v[0];
}
*getn = \&getnode;
*gn = \&getnode;
*gettree = \&getnode;

sub sgetnode {
    my $tree = shift;
    my @v = getnode($tree, @_);
    return $v[0];
}
*sgetn = \&sgetnode;
*sgn = \&sgetnode;
*sgettree = \&sgetnode;


sub getl {
    my $tree = shift || confess;
    my @elts = @_;
    my %elth = map{$_=>1} @elts;
    my %valh = ();
    confess("problem: $tree not arr") unless ref($tree) && ref($tree) eq "ARRAY" || isastag($tree);
    my ($top_ev, $children) = @$tree;
    my @v = ();
    foreach my $child (@$children) {
        my ($ev, $subtree) = @$child;
        my $is_nt = ref($subtree);
        if ($elth{$ev}) {
            # warn if dupl?
            $valh{$ev} = $is_nt ? $child : $subtree;
        }
    }
    return map {$valh{$_}} @elts;
}
*getlist = \&getl;
*gl = \&getl;

sub sget {
    my $tree = shift;
    my @v = get($tree, @_);
    # warn if multivalued?
    return $v[0];
}
*sg = \&sget;

sub sgetdata {
    my $tree = shift;
    my @v = getdata($tree, @_);
    # warn if multivalued?
    return $v[0];
}
*sgd = \&sgetdata;

sub mapv {
    my $tree = shift;
    my %maph = @_;
    foreach my $oldkey (keys %maph) {
	my $newkey = $maph{$oldkey};
	my @currv = get($tree, $newkey);
	next if @currv;
	my @v = get($tree, $oldkey);
	set($tree, $newkey, @v);
	unset($tree, $oldkey);
    }
    return;
}

sub sgetmap {
    my $tree = shift;
    my %maph = @_;
    my %vh = ();
    foreach my $oldkey (keys %maph) {
	# if 0 is supplied, use oldkey
	my $newkey = $maph{$oldkey} || $oldkey;
	my @v = get($tree, $oldkey);
	if (@v > 1) {
	    $tree->throw("multivalued key $oldkey");
	}
	$vh{$newkey} = $v[0];
    }
    return %vh;
}
*sgm = \&sgetmap;

sub sfindval {
    my $tree = shift;
    my @v = findval($tree, @_);
    # warn if multivalued?
    return $v[0];
}
*sfv = \&sfindval;
*singlevalFindSubTreeVal = \&sfindval;

# private
sub indexOn {
    my $tree = shift;
    my $key = shift;

    my %h = ();
    my ($evParent, $stParent) = @$tree;
    foreach my $subtree (@$stParent) {
	my @vl = get($subtree, $key);
	foreach my $v (@vl) {
	    if (!$h{$v}) { $h{$v} = [] }
	    push(@{$h{$v}}, $subtree);
	}
    }
    return \%h;
}

# does a relational style join
sub ijoin {
    my $tree = shift;
    my $element = shift;      # name of element to join
    my $key = shift;          # name of join element
    my $searchstruct = shift; # structure
    my @elts = $tree->fst($element);
    paste($_, $key, $searchstruct)
      foreach @elts;
    
    return;
}
*ij = \&ijoin;
*j = \&ijoin;
*nj = \&ijoin;
*njoin = \&ijoin;

sub paste {
    my $tree = shift;
    my $key = shift;
    my $searchstruct = shift;
    # use indexing?
    my ($key1, $key2) = ($key, $key);
    if ($key =~ /(.*)=(.*)/) {
	($key1, $key2) = ($1, $2);
    }
    my $ssidx = indexOn($searchstruct, $key2);

    my ($evParent, $stParent) = @$tree;
    my @children = ();
    foreach my $subtree (@$stParent) {
	my @nu = ($subtree);
	my ($ev, $st) = @$subtree;
	if ($ev eq $key1) {
	    $tree->throw("can't join on $ev - $st is not primitive")
	      if ref $st;
	    my $replace = $ssidx->{$st} || [];
	    @nu = @$replace;
	}
	push(@children, @nu);
    }
    @$tree = ($evParent, \@children);
    return;
}

# 
sub maptree {
    my $tree = shift;
    my $code = shift;
    my $parent = shift;
    my $next = $code->($tree, $parent);
    if ($next) {
        if (ref($next) && ref($next) eq 'ARRAY') {
            return @$next;
        }
        my $elt = $next->element;
        my @subnodes = subnodes($next);
        if (!@subnodes) {
            return $next;
        }
        my @new_subnodes = ();
        foreach (@subnodes) {
            my @nu = maptree($_, $code, $tree);
            push(@new_subnodes,@nu);
        }
        return Data::Stag->new($elt=>[@new_subnodes]);
    }
    else {
        return ();
    }
}

# iterate depth first through tree executing code
sub iterate {
    my $tree = shift;
    my $code = shift;
    my $parent = shift;
    $code->($tree, $parent);
    my @subnodes = subnodes($tree);
    foreach (@subnodes) {
        iterate($_, $code, $tree);
    }
    return;
}
*i = \&iterate;

# takes a denormalized flat table of rows/columns
# and turns it back into its original tree structure;
# useful for querying databases
# DEPRECATED
sub normalize {
    my $tree = shift || [];
    my ($schema,
        $rows,
        $top,
        $cols,
        $constraints,
        $path) =
          rearrange([qw(schema
                        rows
                        top
                        cols
                        constraints
                        path)], @_);
    if (!$schema) {
        $schema = $tree->new(schema=>[]);
    }
    if (!isastag($schema)) {
        if (!ref($schema)) {
            # it's a string - parse it
            # (assume sxpr)
        }
        $schema = $tree->from('sxprstr', $schema);
    }
    # TOP - this is the element name
    # to group the structs under.
    # [override if specified explicitly]
    if ($top) {
        $schema->set_top($top);
    }
    $top = $schema->get_top || "set";
    my $topstruct = $tree->new($top, []);

    # COLS - this is the columns (attribute names)
    # in the order they appear
    # [override if specified explicitly]
    if ($cols) {
        my @ncols =
          map {
              if (ref($_)) {
                  $_
              }
              else {
                  # presume it's a string
                  # format = GROUP.ATTRIBUTENAME
                  if (/(\w+)\.(\w+)/) {
                      $tree->new(col=>[
                                       [group=>$1],
                                       [name=>$2]]);
                  }
                  else {
                      confess $_;
                  }
              }
          } @$cols;
        $schema->set_cols([@ncols]);
    }


    # PATH - this is the tree structure in
    # which the groups are structured
    # [override if specified explicitly]
    if ($path) {
        if (ref($path)) {
        }
        else {
            $path = $tree->from('sxprstr', $path);
        }
        $schema->set_path([$path]);
    }
    $path = $schema->sgetnode_path;
    if (!$path) {
        confess("no path!");
    }
    
    # column headings
    my @cols = $schema->sgetnode_cols->getnode_col();

    # set the primary key for each group;
    # the default is all the columns in that group
    my %pkey_by_groupname = ();
    my %cols_by_groupname = ();
    foreach my $col (@cols) {
        my $groupname = $col->get_group;
        my $colname = $col->get_name;
        $pkey_by_groupname{$groupname} = []
          unless $pkey_by_groupname{$groupname};
        push(@{$pkey_by_groupname{$groupname}},
             $colname);
        $cols_by_groupname{$groupname} = []
          unless $cols_by_groupname{$groupname};
        push(@{$cols_by_groupname{$groupname}},
             $colname);
    }
    my @groupnames = keys %pkey_by_groupname;

    # override PK if set as a constraint
    my @pks = $schema->findnode("primarykey");
    foreach my $pk (@pks) {
        my $groupname = $pk->get_group;
        my @cols = $pk->get_col;
        $pkey_by_groupname{$groupname} = [@cols];
    }

    # ------------------
    #
    # loop through denormalised rows,
    # grouping the columns into their
    # respecive groups
    #
    # eg
    #
    #  <----- a ----->   <-- b -->
    #  a.1   a.2   a.3   b.1   b.2
    #
    # algorithm:
    #  use path/tree to walk through
    #
    # ------------------

    # keep a hash of all groups by their primary key vals
    #  outer key = groupname
    #  inner key = pkval
    #  hash val  = group structure
    my %all_group_hh = ();
    foreach my $groupname (@groupnames) {
        $all_group_hh{$groupname} = {};
    }

    # keep an array of all groups
    #  outer key = groupname
    #  inner array = ordered list of groups
#    my %all_group_ah = ();
#    foreach my $groupname (keys %pkey_by_groupname) {
#        $all_group_ah{$groupname} = [];
#    }

    my ($first_in_path) = $path->subnodes;
    my $top_record_h = 
      {
       child_h=>{
                 $first_in_path->name=>{}
                },
       struct=>$topstruct
      };
    # loop through rows
    foreach my $row (@$rows) {
        my @colvals = @$row;

        # keep a record of all groups in
        # this row
        my %current_group_h = ();
        for (my $i=0; $i<@cols; $i++) {
            my $colval = $colvals[$i];
            my $col = $cols[$i];
            my $groupname = $col->get_group;
            my $colname = $col->get_name;
            my $group = $current_group_h{$groupname};
            if (!$group) {
                $group = {};
                $current_group_h{$groupname} = $group;
            }
            $group->{$colname} = $colval;
        }

        # we now have a hash of hashes -
        #  outer keyed by group id
        #  inner keyed by group attribute name
        
        # traverse depth first down path;
        # add new nodes as children of the parent
        sub make_a_tree {
            my $class = shift;
            my $parent_rec_h = shift;
            my $node = shift;
            my %current_group_h= %{shift ||{}};
            my %pkey_by_groupname = %{shift ||{}};
            my %cols_by_groupname = %{shift ||{}};
            my $groupname = $node->name;
            my $grouprec = $current_group_h{$groupname};
            my $pkcols = $pkey_by_groupname{$groupname};
            my $pkval = 
              CORE::join("\t",
                         map {
                             esctab($grouprec->{$_} || '')
                         } @$pkcols);
            my $rec = $parent_rec_h->{child_h}->{$groupname}->{$pkval};
            if (!$rec) {
                my $groupcols = $cols_by_groupname{$groupname};
                my $groupstruct =
                  $class->new($groupname=>[
                                          map {
                                              [$_ => $grouprec->{$_}]
                                          } @$groupcols
                                         ]);
                my $parent_groupstruct = $parent_rec_h->{struct};
                if (!$parent_groupstruct) {
                    confess("no parent for $groupname");
                }
                add($parent_groupstruct,
                    $groupstruct->name,
                    $groupstruct->data);
                $rec =
                  {struct=>$groupstruct,
                   child_h=>{}};
                foreach ($node->subnodes) {
                    # keep index of children by PK
                    $rec->{child_h}->{$_->name} = {};
                }
                $parent_rec_h->{child_h}->{$groupname}->{$pkval} = $rec;
            }
            foreach ($node->subnodes) {
                make_a_tree($class,
                            $rec, $_, \%current_group_h,
                            \%pkey_by_groupname, \%cols_by_groupname);
            }
        }
        make_a_tree($tree,
                    $top_record_h, $first_in_path, \%current_group_h,
                    \%pkey_by_groupname, \%cols_by_groupname);
    }
    return $topstruct;
}
*norm = \&normalize;
*normalise = \&normalize;


sub esctab {
    my $w=shift;
    $w =~ s/\t/__MAGICTAB__/g;
    $w;
}

sub findvallist {
    my $tree = shift || confess;
    my @nodes = @_;
    my @vals =
      map {
          my @v = findval($tree, $_);
          if (@v > 1) {
              confess(">1 val for $_: @v");
          }
          $v[0];
      } @nodes;
    return @vals;
}
*findSubTreeValList = \&findvallist;
*fvl = \&findvallist;

sub findChildVal {
    confess("deprecated - use get");
    my $tree = shift;
    my $node = shift;
    my $replace = shift;
    confess unless ref($tree);
    my ($ev, $subtree) = @$tree;
    return $subtree if test_eq($ev, $node);
    return unless ref($subtree);
    my @children = grep { $_->[0] eq $node } @$subtree;
#    print "@children\n";
    if (defined $replace) {
        return 
          map { $_->[1] = $replace } @children;
    }
    my @r = map { $_->[1] } @children;
    if (wantarray) {
        return @r;
    }
    $r[0];
}

sub qmatch {
    my $tree = shift;
    my $elt = shift;
    my $matchkey = shift;
    my $matchval = shift;
    my $replace = shift;

    my @st = findnode($tree, $elt);
    my @match =
      grep {
#          tmatch($_, $matchkey, $matchval);
	  grep {$_ eq $matchval} $_->get($matchkey);
      } @st;
    if ($replace) {
	map {
	    @$_ = @$replace;
	} @match;
    }
    return @match;
}
*qm = \&qmatch;
*findSubTreeMatch = \&qmatch;


sub tmatch {
    my $tree = shift;
    my $elt = shift;
    my $matchval = shift;
    my @vals = findval($tree, $elt);
    return grep {$_ eq $matchval} @vals;
}
*testSubTreeMatch = \&tmatch;
*tm = \&tmatch;


sub tmatchhash {
    my $tree = shift;
    my $match = shift;
    my @mkeys = keys %$match;
    my @mvals = map {$match->{$_}} @mkeys;

    my @rvals = findvallist($tree, @mkeys);
    my $pass =  1;
    for (my $i=0; $i<@mvals; $i++) {
        $pass = 0 if $mvals[$i] ne $rvals[$i];
    }
#    print "CHECK @mvals eq @rvals [$pass]\n";
    return $pass;
}
*tmh = \&tmatchhash;
*testSubTreeMatchHash = \&tmatchhash;

sub tmatchnode {
    my $tree = shift;
    my $matchtree = shift;
    my ($node, $subtree) = @$tree;
    confess unless ref $matchtree;
    my ($mnode, $msubtree) = @$matchtree;
    if ($node ne $mnode) {
        return unless ref $subtree;
        return 
          grep {
              testSubTreeMatchTree($_,
                                   $matchtree)
          } @$subtree;
    }
    if (!ref($subtree) && !ref($msubtree)) {
        return $subtree eq $msubtree;
    }
    if (ref($subtree) && ref($msubtree)) {
        my @got = ();
        for (my $i=0; $i<@$msubtree; $i++) {
            my $n = $msubtree->[$i]->[0];
            my ($x) = grep {$_->[0] eq $n} @$subtree;
            return unless $x;
            my $ok =
              testSubTreeMatchTree($x,
                                   $msubtree->[$i]);
            return unless $ok;
                                   
        }
        return 1;
    }

}
*tmn = \&tmatchnode;
*testSubTreeMatchTree = \&tmatchnode;

sub cmatch {
    my $tree = shift;
    my $node = shift;
    my $matchval = shift;
    my @vals = findval($tree, $node);
    return scalar(grep {$_ eq $matchval} @vals);
}
*cm = \&cmatch;
*countSubTreeMatch = \&cmatch;


sub OLDwhere {
    my $tree = shift;
    my $node = shift;
    my $testcode = shift;
    my $replace = shift;
    my @subtrees = findnode($tree, $node);
    my @match =
      grep {
          $testcode->($_);
      } @subtrees;
    if (defined $replace) {
        map {
            @$_ = @$replace;
        } @match;
    }
    return @match;
}
sub where {
    my $tree = shift;
    my $node = shift;
    my $testcode = shift;
    my $replace = shift;

    my @match = ();
    iterate($tree,
	    sub {
		my $stag = shift;
		if (name($stag) eq $node) {
		    if ($testcode->($stag)) {
			push(@match, $stag);
		    }
		}
	    });

    if (defined $replace) {
        map {
            @$_ = @$replace;
        } @match;
    }
    return @match;
}
*w = \&where;
*findSubTreeWhere = \&where;


#sub findSubTreeWhere {
#    my $tree = shift;
#    my $node = shift;
#    my $testcode = shift;
#    my $replace = shift;
##    use Data::Dumper;
##    print Dumper($expr);
#    my $call = $expr->name;
#    my @p = $expr->findChildVal('arg');
##    print Dumper(\@p);
#    my @subtrees = findSubTree($tree, $node);
#    no strict 'refs'; 
#    my @match =
#      grep {
#          &$call($_, @p);
#      } @subtrees;
#    if (defined $replace) {
#        map {
#            print "WAS:";
#            print xml($_);
#            print "NOW:";
#            print xml($replace);
#            @$_ = @$replace;
#        } @match;
#    }
#    return @match;
#}

sub run {
    my $tree = shift;
    my $code = shift;
    my $func = $code->name;
    my @p = $code->children();
    my @args = ();
    foreach my $p (@p) {
        if ($p->name eq 'arg') {
            if (isastag($p->children)) {
                die;
                push(@args,
                     evalTree($tree,
                              $p->children));
                              
            }
            else {
                push(@args, $p->children);
            }
        }
        else {
#            print "rcall $tree $p\n";
            push(@args,
                 evalTree($tree,
                          $p));
        }
    }
    no strict 'refs'; 
    my @r = &$func($tree, @p);
    return @r if wantarray;
    return shift @r;
}
*evalTree = \&run;


# ------------------------
# collapseElement($tree, $elt)
#
# eg if we have
#
# +personset
#    +person
#       +name jim
#       +job  sparky
#    +person
#       +name jim
#       +hair green
#    +person
#       +name eck
#       +hair blue
#
# then execute
# collapseElement($tree, 'person', 'name')
#
# we end up with 
#
# +personset
#    +person
#       +name jim
#       +job  sparky
#       +hair green
#    +person
#       +name eck
#       +hair blue
#
# OR if we have
#
# +personset
#    +person
#       +name jim
#       +job  sparky
#       +petname  bubbles
#       +pettype  chimp
#    +person
#       +name jim
#       +job  sparky
#       +petname  flossie
#       +pettype  sheep
#    +person
#       +name eck
#       +job  bus driver
#       +petname  gnasher
#       +pettype  dug
#
#
# then execute
# collapseElement($tree, 'name', 'job')
#
# we end up with 
#
# +personset
#    +person
#       +name jim
#       +job  sparky
#       +pet
#          +petname  bubbles
#          +pettype  chimp
#       +pet
#          +petname  flossie
#          +pettype  sheep
#    +person
#       +name eck
#       +job  bus driver
#       +pet
#          +petname  gnasher
#          +pettype  dug
#
#
# warning: element should be unique
# todo: allow person/name
#
# ------------------------
sub collapse {
    my $tree = shift;
    my $elt = shift;
    my $merge_elt = shift;
    my @subtrees = findnode($tree, $elt);
    my %treeh = ();
    my @elt_order = (); # preserve ordering
    map {
        my @v = findval($_, $merge_elt);
        die unless scalar(@v) == 1;
        my $val = shift @v;
        push(@elt_order, $val) unless $treeh{$val};
        $treeh{$val} = [] unless $treeh{$val};
        push(@{$treeh{$val}}, $_);
    } @subtrees;
    my @new_subtrees =
      map {
          my $trees = $treeh{$_};
          my ($v) = findval($trees->[0], $merge_elt);
          [
           $elt=>[
                  [$merge_elt=>$v],
                  map {
                      my ($ev, $subtree) = @$_;
                      grep {
                          $_->[0] ne $merge_elt
                      } @$subtree;
                  } @$trees,
                 ]
          ]
      } @elt_order;

    findnode($tree, $elt, []);
    push(@{$tree->[1]}, @new_subtrees);
    return $tree;
}
*collapseElement = \&collapse;

sub merge {
    my $tree = shift;
    my @elts = @{shift || []};
    my $merge_key = shift;

    return unless @elts;
    my $e1 = $elts[0];
    my ($type, $v) = @$e1;
    my @cur_elts = findnode($tree, $type);
    foreach my $elt (@elts) {
        my ($v) = findval($elt, $merge_key);
        foreach my $cur_elt (@cur_elts) {
            my ($cv) = findval($cur_elt, $merge_key);
            if ($cv eq $v) {
                # merge
                my $cur_children = $cur_elt->[1];
                my $children = $elt->[1];
                push(@$cur_children,
                     grep {
                         $_->[0] ne $merge_key
                     } @$children);
            }
        }0
    }
    return $tree;
}
*mergeElements = \&merge;

sub makeattrsnodes {
    my $tree = shift;
    return if isterminal($tree);
    my @attrs = get($tree,'@');
    if (@attrs) {
	my @nu = ();
	foreach (@attrs) {
	    push(@nu, kids($_));
	}
	unset($tree, '@');
	unshift(@{$tree->[1]},@nu);
    }
    my @subnodes = subnodes($tree);
    makeattrsnodes($_) foreach @subnodes;
    return;
}

sub duplicate {
    my $tree = shift;
    load_module('Data::Dumper');
    use Data::Dumper;
    my $nu;
    my $d = Data::Dumper->new( [$tree], [qw($nu)] );
    my $dump = $d->Dump;
    eval $dump;
    if ($@) {
	confess $@;
    }
    return stagify($nu);
}
*d = \&duplicate;
*clone = \&duplicate;

sub isastag {
    my $node = shift;
    return UNIVERSAL::isa($node, "Data::Stag::StagI") ||
      UNIVERSAL::isa($node, "Node");
}
*isanode = \&isastag;
*isa_node = \&isanode;
*isaNode = \&isanode;

sub isnull {
    my $node = shift;
    if (@$node) {
	return 0;
    }
    return 1;
}

sub node {
    return Data::Stag::StagImpl->new(@_);
}
*Node = \&node;

sub stagify {
    my $tree = shift;
    if (!ref($tree)) {
        # allow static or nonstatic usage
        $tree = shift;
    }
    confess unless ref $tree;
    my $class = "Data::Stag::StagImpl";
    bless $tree, $class
}
*nodify = \&stagify;
*Nodify = \&stagify;

sub xpath {
    my $tree = shift;
    load_module("XML::XPath");
    my $xp = XML::XPath->new(xml=>xml($tree));
    return $xp;
}
*xp = \&xpath;
*getxpath = \&xpath;
*tree2xpath = \&xpath;


sub xpquery {
    my $tree = shift;
    my @args = @_;
    load_module("XML::XPath::XMLParser");
    my $xp = $tree->getXPath;
    my $nodeset = $xp->find(@args);
    my @nodes =
      map {
	  xmlstr2tree(XML::XPath::XMLParser::as_string($_));
      } $nodeset->get_nodelist;
    return @nodes;
}
*xpq = \&xpquery;
*xpathquery = \&xpquery;
*xpfind = \&xpquery;
*xpFind = \&xpquery;

sub grammarparser {
    my $tree = shift;
    my $grammar = shift;
    load_module("Parse::RecDescent");
    $::RD_AUTOACTION = q 
      { use Data::Stag;
	$#item == 1 && !ref $item[1] ? $item[1] : Data::Stag->new(shift @item, [map {if(ref($_)) {$_} else {[arg=>$_]}} @item ]); 
      };
    my $parser = Parse::RecDescent->new($grammar) or confess "Bad grammar!\n";
    return $parser;
}

#use overload
#  '.' => sub {my @r=findnodeVal($_[0],$_[1]);$r[0]},
#  '-' => sub {my @r=findnodeVal($_[0],$_[1]);$r[0]},
#  '+' => sub {my @r=$_[0]->findnode($_[1]);return $r[0]},
#  '/' => sub {[findnodeVal($_[0],$_[1])]},
#  '*' => sub {[$_[0]->findnode($_[1])]},
#  '<' => sub {Data::Stag::testSubTreeMatchTree($_[0], $_[1])},
#  qw("" stringify);

#sub stringify {
#   $_[0];
#}


sub kids {
    my $self = shift;
    if (@_) {
        @$self = ($self->[0], [map {Nodify($_)} @_]);
    }
    my ($name, $kids) = @$self;
    if (!ref($kids)) {
        return $kids;
    }
    return map {Nodify($_)} @$kids;
}
*k = \&kids;
*children = \&kids;
*getChildren = \&kids;
*getKids = \&kids;

sub subnodes {
    my $self = shift;
    return grep {ref($_)} kids($self);
}

# non-terminal nodes
sub ntnodes {
    my $self = shift;
    my @subnodes = $self->subnodes;
    return grep {!$_->isterminal} @subnodes;
}
*nonterminalnodes = \&ntnodes;
*nonterminals = \&ntnodes;

# terminal nodes
sub tnodes {
    my $self = shift;
    my @subnodes = $self->subnodes;
    return grep {$_->isterminal} @subnodes;
}
*terminalnodes = \&tnodes;
*terminals = \&tnodes;

sub element {
    my $self = shift;
    if (@_) {
        $self->[0] = shift;
    }
    return $self->[0];
}
*e = \&element;
*name = \&element;
*tagname = \&element;

sub data {
    my $self = shift;
    if (@_) {
        $self->[1] = shift;
    }
    return $self->[1];
}

sub rename {
    my $tree = shift;
    my $from = shift;
    my $to = shift;
    foreach (kids($tree)) {
	if ($_->[0] eq $from) {
	    $_->[0] = $to;
	}
    }
    return;
}

sub isterminal {
    my $self = shift;
    return !ref($self->[1]);
}

sub _min {
    my $x = shift;
    my $y = shift;
    return $x if !defined($y);
    return $y if !defined($x);
    $x < $y ? $x : $y;
}
sub _max {
    my $x = shift || 0;
    my $y = shift || 0;
    $x > $y ? $x : $y;
}

# automatically deduces schema
# FORMAT MAY CHANGE!!!!!
sub autoschema {
    my $tree = shift;
    my $schema = _autoschema($tree);
#    use Data::Dumper;
#    print Dumper $schema;
    my $nu = genschema($tree, undef, $tree->element, $schema);
    $nu->iterate(sub{
		     my $node = shift;
		     if ($node->name =~ /(\S+)\.(\S+)/) {
			 $node->name($2);		     }
		 });
    return $nu;
}

sub dtd {
    my $tree = shift;
    my ($name, $subtree) = @$tree;
    my $done_h = shift || {};
    $name =~ s/[\+\?\*]$//;
    return '' if $done_h->{$name};
    my $is_nt = ref($subtree);
    my $s;
    if ($is_nt) {
        my $s2 = join('', map {dtd($_,$done_h)} @$subtree);
        my @subnames = map {$_->[0]} @$subtree;
        my $S = join('|',@subnames);
        if (@subnames < 2) {
            $S = "@subnames";
            if (!@subnames) {
                $S = 'EMPTY';
            }
        }

        $s = "\n<!-- $name: (node) -->\n<!ELEMENT $name ($S)+>\n$s2";
    }
    else {
        $s = "\n<!-- $name: ($subtree) -->\n<!ELEMENT $name (#PCDATA)>\n";
    }
    $done_h->{$name} = 1;
    $s;
}

sub genschema {
    my $tree = shift;
    my $parent = shift;
    my $root = shift;
    my $schema = shift;
    my $path = shift || [];

    my $cycle = 0;
    if (grep {$_ eq $root} @$path) {
        $cycle = 1;
        warn "cycle detected: @$path $root";
        return (Data::Stag->new($root=>[]));
    }
    push(@$path, $root);

    my $data = $schema->{data};
    my $childh = $schema->{childh};

    my $card = "";
    if ($parent) {
        my $link = "$parent $root";
        my $min = $schema->{mincard}->{$link};
        my $max = $schema->{maxcard}->{$link};
        die if !$max;
        if ($min == 0) {
            if ($max == 1) {
                $card = '?';
            }
            else {
                $card = '*';
            }
        }
        else {
            # $min >= 1
            if ($max == 1) {
                $card = '';
            }
            else {
                $card = '+';
            }
        }
    }
    my $ss = Data::Stag->new($root.$card=>[]);
    if ($data->{$root}) {
        $ss->data($data->{$root});
    }
    else {
        my $c = $childh->{$root};
        my @sn =
          map {
              genschema($tree, $root, $_, $schema, $path);
          } @$c;
        $ss->data([@sn]);
    }
    return $ss;
}

# automatically deduces schema
sub _autoschema {
    my $tree = shift;
    my $schema = shift || {data=>{}, 
                           childh=>{},
                           mincard=>{},
                           maxcard=>{},
                          };
    my $data = $schema->{data};
    my $childh = $schema->{childh};
    my $mincard = $schema->{mincard};
    my $maxcard = $schema->{maxcard};

    my $elt = $tree->element;
    my @sn = $tree->subnodes;
    # CARD:   default(blank) : 1
    #                      + : 1 or more
    #                      * : 0 or more
    #                      ? : 0 or one
    my %lcard = ();  # local cardinality
    foreach (@sn) {
	# nonterminal nodes are uniquely defined
	# by the node name
        my $se = element($_);
	if (isterminal($_)) {
	    # a terminal node is uniquely defined
	    # by parent.node
	    $se = "$elt.$se";
	}
        push(@{$childh->{$elt}}, $se)
          unless $childh->{$elt} &&
            grep { $_ eq $se } @{$childh->{$elt}};
        $lcard{$se} = 0 unless $lcard{$se};
        $lcard{$se}++;
    }
#    foreach (keys %lcard) {
    foreach (@{$childh->{$elt}}) {
        my $link = "$elt $_";
#        print "$link :: $lcard{$_}\n";
        $mincard->{$link} =
          _min($lcard{$_} || 0,
               $mincard->{$link});
        $maxcard->{$link} =
          _max($lcard{$_},
               $maxcard->{$link});
    }
    foreach (grep {isterminal($_)} @sn) {
        my $elt = $elt.'.'.element($_);
#        push(@{$data->{$elt}}, $_->data);
        my $in = $_->data;
        my $d = $data->{$elt} || 'INT';
        if (!$in) {
        }
        elsif ($in =~ /^\-?\d+$/) {  # LOOKS LIKE INT
	    # CHANGE SIZE IF IT IS
            if ($d eq 'INT') {
		$d = 'INT';
		my $lin = length($in);
		if ($lin > 10) {
		    # too big for an int
		    # TODO: largeint?
		    $d = _mkvarchar($lin);
		}
	    }
        }
        elsif ($in =~ /^\-?\d*\.\d+$/) { # LOOKS LIKE FLOAT
	    # PROMOTE TO FLOAT IF INT
            if ($d eq 'INT') {
		$d = 'FLOAT';
	    }
        }
        else {
            my $lin = length($in) || 0;
            if ($d =~ /VARCHAR\((\d+)\)/) {
                if ($lin > $1) {
		    $d = _mkvarchar($lin);
                }
                else {
                }
            }
            elsif ($d =~ /TEXT/) {
            }
            else {
                $d = _mkvarchar($lin);
            }
        }
        $data->{$elt} = $d;
    }
    foreach (grep {!isterminal($_)} @sn) {
        _autoschema($_, $schema);
    }
    return $schema;
}

sub _mkvarchar {
    my $size = shift || 1;
    # round up to log2
    my $s2 = 2**(int(log($size) / log(2))+2) -1;
    if ($s2 > 255) {
        return 'TEXT';
    }
    return "VARCHAR($s2)";
}

sub AUTOLOAD {
    my $self = shift;
    my @args = @_;

    my $name = $AUTOLOAD;
    $name =~ s/.*://;   # strip fully-qualified portion

    if ($name eq "DESTROY") {
	# we dont want to propagate this!!
	return;
    }

    if ($name =~ /^([a-zA-Z]+)_(\w+)/) {
        if ($self->can($1)) {
            return $self->$1($2, @args);
        }
    }
    confess("no such method:$name");
}

# --MISC--

sub splitpath {
    my $node = shift;
    if (ref($node) && ref($node) eq 'ARRAY') {
        @$node;
    }
    elsif ($node =~ /\//) {
        (split(/\//, $node));
    }
    else {
        ($node);
    }
}

sub test_eq {
    my ($ev, $node) = @_;
    $ev = '' unless defined $ev;
    $node = '' unless defined $node;
    return $ev eq $node || $node eq '*';
}

1;

