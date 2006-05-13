#! /usr/bin/perl -w
use strict;
use CGI;
use CGI::Carp 'fatalsToBrowser';
use CGI::Ajax;
use CoGe::Accessory::LogUser;
use HTML::Template;
use Data::Dumper;
use Bio::TreeIO;
use File::Temp;

$ENV{PATH} = "/opt/apache2/CoGe/";

use vars qw( $DATE $DEBUG $TEMPDIR $TEMPURL $USER $FORM $NEATO $DOT);

# set this to 1 to print verbose messages to logs
$DEBUG = 0;
$TEMPDIR = "/opt/apache/CoGe/tmp";
$TEMPURL = "/CoGe/tmp";
$NEATO = "/opt/apache/CoGe/bin/neato";
$DOT = "/opt/apache/CoGe/bin/dot";
$| = 1; # turn off buffering
$DATE = sprintf( "%04d-%02d-%02d %02d:%02d:%02d",
		 sub { ($_[5]+1900, $_[4]+1, $_[3]),$_[2],$_[1],$_[0] }->(localtime));

$FORM = new CGI;
$USER = CoGe::Accessory::LogUser->get_user();
my $pj = new CGI::Ajax(
		       gen_html=>\&gen_html,
		       zoom=>\&zoom,
		      );
$pj->JSDEBUG(0);
$pj->DEBUG(0);

print $pj->build_html($FORM, \&gen_html);
#print "Content-Type: text/html\n\n". gen_html();
#gen_html();

sub gen_html
  {
    my ($body, $seq_names, $seqs) = gen_body();
    my $template = HTML::Template->new(filename=>'/opt/apache/CoGe/tmpl/generic_page.tmpl');

    $template->param(TITLE=>'CoGe: Phylogenetic Tree Viewer');
    $template->param(USER=>$USER);
    $template->param(DATE=>$DATE);
    $template->param(LOGO_PNG=>"TreeView-logo.png");
    $template->param(BODY=>$body);
    my $html;# =  "Content-Type: text/html\n\n";
    $html .= $template->output;
    return $html;
  }

sub gen_body
  {
    my $content;
    if ($FORM->param('file') )
      {
	my $fh = $FORM->upload('file');
	while (<$fh>) {$content .= $_};
      }
    elsif ($FORM->param('tree') )
      {
	$content = $FORM->param('tree');
      }
    $content =~ s/\r/\n/g;
    $content =~ s/\n+/\n/g;

#    $content =1;
    my $template = HTML::Template->new(filename=>'/opt/apache/CoGe/tmpl/TreeView.tmpl');
    $template->param(FRONT_PAGE => 1) unless $content;
    if ($content)
      {
	my ($tree, $dot, $img, $imgmap) = generate_tree($content);
	$img =~ s/$TEMPDIR/$TEMPURL/;
	$tree =~ s/$TEMPDIR/$TEMPURL/;
	$dot =~ s/$TEMPDIR/$TEMPURL/;
	$template->param(IMG=>$img);
	$template->param(IMGMAP=>$imgmap) if $imgmap;
	$template->param('tree_file'=>$tree);
	$template->param('dot_file'=>$dot);
	$template->param('tfile'=>$tree);
	my $mod_file = generate_mod_file($FORM->param('nodes')) if $FORM->param('nodes');
	$mod_file =~ s/$TEMPDIR/$TEMPURL/;
	$template->param('mod'=>$mod_file) if $mod_file;
	
      }
    return $template->output;
  }



sub zoom
  {
    my $zoom = shift;
    my $tfile = shift;
    my $mod = shift;
    my $extra = shift;
    $zoom += $extra if $extra;
    print STDERR $zoom,"\n";
    $mod =~ s/$TEMPURL/$TEMPDIR/;
    my $mod_content = load_mod($mod);
    $FORM->param('nodes', $mod_content);
    $tfile =~ s/$TEMPURL/$TEMPDIR/;
    warn "can't read $tfile" unless -f $tfile;
    my ($dot_file) = generate_dot_file(file=>$tfile, zoom=>$zoom);
    my ($img, $map) = generate_image($dot_file);
    my $zoomin = '<input type="hidden" id = "zoi" value='.($zoom+1).'>';
    my $zoomout = '<input type="hidden" id = "zoo" value='.($zoom-1).'>';
    return ($img, $map, $zoomin, $zoomout);
  }

sub generate_tree
  {
    my $content = shift;
    return unless $content;
    
    my ($dot_file, $tree_file) = generate_dot_file(content=>$content);
    my ($img_file, $map) = generate_image($dot_file);
    return ($tree_file, $dot_file, $img_file, $map);
  }
sub generate_image
  {
    my $file = shift;
    my $cmd = $FORM->param("rooted") ? $DOT : $NEATO;
    my $of = new File::Temp ( TEMPLATE=>'TV__XXXXX',
			      DIR=>$TEMPDIR,
			      SUFFIX=>'.png',
			      UNLINK=>0);
    my $ofn = $of->filename;
    system "$cmd -Tpng $file -o $ofn";
    close $of;
    open (IN,  "$cmd -Tcmap $file |");
    my $map;
    while (<IN>) {$map .= $_;}
    close IN;
    my $map_name = "tree";
    $ofn =~ s/$TEMPDIR/$TEMPURL/;
    return "<img src=\"$ofn\" USEMAP=\"#".$map_name."\" border=0>", "<MAP NAME=$map_name>$map</MAP>";
  }

sub generate_dot_file
  {
    my %opts = @_;
    my $content = $opts{content};
    my $file = $opts{file};
    my $zoom = $opts{zoom};
    if ($content)
      {
	my $tmp_file = new File::Temp ( TEMPLATE=>'TV__XXXXX',
					DIR=>$TEMPDIR,
					SUFFIX=>'.nexus_tree',
					UNLINK=>0);
	print $tmp_file $content;
	$file = $tmp_file->filename;
	close $tmp_file;
      }
    return unless -r $file;
    my $node_names = get_names($file);
    my $treeio = new Bio::TreeIO (#-format =>'newick',
				  -file=>$file,
				 );
    $treeio->next_tree;
    while (my $tree = $treeio->next_tree())
      {
	my $root_node = $tree->get_root_node;
	next unless $root_node->internal_id;
	my $dot;
	$dot .= "graph tree {\n";
	my $size = 8;
	$size += $zoom if $zoom;
	$size = 1 if $size < 1;
	$dot .= "\tsize=\"$size,$size\";\n";
	$dot .= "\tfontsize=10;\n";
	$dot .= "\tgraph [splines=true ";
	$dot .= "overlap=scale" unless $FORM->param('overlap');
	$dot .= "];\n";
	$dot .= process_node($root_node, $node_names);
	$dot .= label_nodes($tree, $node_names);
	$dot .= "}\n";
	my $dot_file = new File::Temp ( TEMPLATE=>'TV__XXXXX',
					DIR=>$TEMPDIR,
					SUFFIX=>'.dot',
					UNLINK=>0);
	print $dot_file $dot;
	my $dot_fn = $dot_file->filename;
	close $dot_file;
	return ($dot_fn, $file);
      }
  }
sub label_nodes
  {
    my $tree = shift;
    my $names = shift;
    my $node_labels = get_labels();
    my $content;
    foreach my $node ($tree->get_nodes)
      {
	if ($node->id)
	  {
	    my $name = $names->{$node->id};
	    my $label = $name;
	    $label =~ s/_/\./;
	    $label =~ s/_/\\n/g;
	    my ($accn1) = $label =~ /^([^\\]+)/;
	    my $url = "B2V.pl?";
	    $url .= "accn1=$accn1&" if $accn1;
	    foreach my $sister ($node->ancestor->get_all_Descendents)
	      {
		next if $sister->internal_id == $node->internal_id;
		my $tmp = $names->{$sister->id};
		next unless $tmp;
		$tmp =~ s/_/ /g;
		my ($accn2) = $tmp =~ /^(\S+)/;
		$url .= "accn2=$accn2" if $accn2;
		last if $accn2;
	      }
	    while (my ($pat, $mod) = each %$node_labels)
	      {
		if ( $name =~ /$pat/i)
		  {
		    $content .= "\t\t$name\t[label = \"$label\", URL=\"$url\", $mod ];\n";
		  }
	      }
	  }
	else
	  {
	    my $id = $node->internal_id;
	    $content .= "\t\t$id\t[label=\"$id\", shape=point ];\n";
	  }
      }
    return $content;
  }

sub get_labels
  {
    my $text = $FORM->param('nodes');
    my %mods;
    return \%mods unless $text;
    foreach my $line (split /\n/,$text)
      {
	my ($key, $mod) = split /\s+/, $line, 2;
	$mod =~ s/,\s+/,/g;
	$mod =~ s/\s+/,/g;
	next unless $key;
	next unless $mod;
	$mods{$key} = $mod;
      }
    return \%mods;
  }

sub process_node
  {
    my $node = shift;
    my $names = shift;
    my $name = $node->id ? $names->{$node->id} : $node->internal_id;
    my $len =  ($node->branch_length) if $node->branch_length;
    my $content;
    foreach my $cnode ($node->each_Descendent)
      {
	my $cname = $cnode->id ? $names->{$cnode->id} : $cnode->internal_id;
	$content .= "\t$name -- $cname ";
#	print "[label=\"".sprintf("%.2f",$len)."\",len=$len]" if defined $len;
	$content .= "[label=\"".sprintf("%.2f",$len)."\"]" if defined $len;
	$content .= ";\n";
	my $tmp = process_node($cnode, $names);
	$content .= $tmp if $tmp;
      }
    return $content;
  }

sub get_names
  {
    my $file = shift;
    open (IN, $file);
    my %names;
    my $flag = 0;
    while (<IN>)
      {
	chomp;
	if ($flag)
	  {
	    s/^\s+//g;
	    my ($num, $name) = split /\s+/, $_, 2;
	    $name =~ s/,$// if $name;
	    $name =~ s/\|/_/g if $name;
	    $name =~ s/jgi_//g if $name;
	    $name =~ s/(Chlre3.\d+)_.*/$1/ if $name;
#	    $name =~ s/_/\\n/ if $name;
	    last if $_ =~ ";";
	    next unless ($num =~ /^\d+$/);
	    $names{$num} = $name;
	    $names{$name} = $num;

	  }
	$flag = 1 if /Translate/;
      }
    close IN;
    return \%names;
  }

sub generate_mod_file
  {
    my $stuff = shift;
    my $file = new File::Temp ( TEMPLATE=>'TV__XXXXX',
				    DIR=>$TEMPDIR,
				    SUFFIX=>'.mod',
				    UNLINK=>0);
    print $file $stuff;
    my $fn = $file->filename;
    close $file;
    return $fn;
  }

sub load_mod
  {
    my $file = shift;
    open (IN, $file);
    my $data;
    while (<IN>) {$data .= $_;}
    return $data;
  }
