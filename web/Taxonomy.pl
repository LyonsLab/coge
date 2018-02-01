#! /usr/bin/perl -w

use strict;

use CGI;
use Data::Dumper;
use JSON qw(encode_json);

use CoGe::Accessory::Web;

no warnings 'redefine';

$|=1;

use vars qw(
    $config $PAGE_TITLE $PAGE_NAME $user $BASEFILE $db %FUNCTION
    $FORM $MAX_SEARCH_RESULTS %ITEM_TYPE $JEX $link
);

$PAGE_TITLE = 'Taxonomy';
$PAGE_NAME  = "$PAGE_TITLE.pl";

$FORM = new CGI;
( $db, $user, $config, $link ) = CoGe::Accessory::Web->init( cgi => $FORM, page_title => $PAGE_TITLE );

%FUNCTION = (
	user_is_admin					=> \&user_is_admin,
    gen_tree_json					=> \&gen_tree_json
);

CoGe::Accessory::Web->dispatch( $FORM, \%FUNCTION, \&gen_html );


sub gen_html {
	my $html;
	my $template =
	  HTML::Template->new( filename => $config->{TMPLDIR} . 'generic_page.tmpl' );
	$template->param( USER       => $user->display_name || '',
	                  PAGE_TITLE => $PAGE_TITLE,
	                  PAGE_LINK  => $link,
	                  SUPPORT_EMAIL => $config->{SUPPORT_EMAIL},
	                  TITLE      => "TAXONOMY",
	                  HOME       => $config->{SERVER},
                      HELP       => '',
                      WIKI_URL   => $config->{WIKI_URL} || '',
                      CAS_URL    => $config->{CAS_URL} || '',
                      ADMIN_ONLY => $user->is_admin,
                      COOKIE_NAME => $config->{COOKIE_NAME} || '');
	$template->param( LOGON      => 1 ) unless $user->user_name eq "public";
	$template->param( BODY       => gen_body() );
	$html .= $template->output;
}

sub gen_body {
	my $template =
	  HTML::Template->new( filename => $config->{TMPLDIR} . 'Taxonomy.tmpl' );
	$template->param( 	API_BASE_URL  	=> $config->{SERVER} . 'api/v1/',  
						USER			=> $user->user_name,	
					);

	return $template->output;
}

sub gen_tree_json {
	my %taxonomic_tree  = (
		name => "root", 
		children => [],
	);
	
	### 
	# Helper functions:
	###

	*check_tree = sub {
		my $search_term = shift;
		if (scalar(@_) == 0) {
			return undef;
		}
		
		my @new_roots;
		foreach my $root (@_) {
			foreach my $child (@{$root->{children}}) {
				push (@new_roots, $child);
				if(lc $search_term eq lc $child->{name}) {
					return $child;
				}
			}
		}
		return check_tree($search_term, @new_roots);
	};
	
	*gen_subtree = sub {
		my @array = @{$_[0]};
		if (scalar(@array) > 0) {
			my $hash = {
				name => $array[0],
				children => [],
			};
			shift @array;
			if (scalar(@array) > 0) {
				push ($hash->{children}, gen_subtree(\@array));
			}
			return $hash;
		}
	};

	*add_fix = sub {
		my $add_tree = $_[0];
		my $move_tree;
		my $i;
		for ($i = scalar(@{$taxonomic_tree{children}} - 1); $i > -1; $i--) {
			if ((lc $add_tree->{name}) eq (lc @{$taxonomic_tree{children}}[$i]->{name})) {
				$move_tree = splice(@{$taxonomic_tree{children}}, $i, 1);
				
				foreach my $child (@{$move_tree->{children}}) {
					add_to_tree($add_tree, $child);
				}
			}
		}
		#check for further inconsistencies
		foreach my $child (@{$add_tree->{children}}) {
			add_fix($child);
		}
	};
	
	*add_to_tree = sub {
		my $root = $_[0];
		my $sub_tree = $_[1];
	    my $root_children = \@{$root->{children}}; #array reference
		my $top_value = $sub_tree->{name};
		my $find = check_tree($top_value, $root);
		if (!$find) {
			#check for inconsistencies
			add_fix($sub_tree);
			
			#add the subtree as a child of the root
			push ($root_children, $sub_tree);
		} else {
			#recurse to the next level of the tree
			foreach my $child (@{$sub_tree->{children}}) {
				add_to_tree($find, $child);
			}
		}
	};
	
	###
	# Main code
	###

	my @organism_descriptions = CoGeDBI::get_table($db->storage->dbh, 'organism', undef, {description => "\"%;%\""}, {description => " like "});
	my $organism_descriptions = $organism_descriptions[0];
	
	foreach my $organism (values %$organism_descriptions) {
		my $taxon_string = $organism->{description};
		my @taxons = split(/\s*;\s*/, $taxon_string);
		
		if($taxons[0]) {
			$taxons[0] =~ s/^\s+|\s+$//g;
			my $size = scalar @taxons;
			$taxons[$size - 1] =~ s/^\s+|\s+$//g;
			my $sub_tree = gen_subtree(\@taxons);
			add_to_tree(\%taxonomic_tree, $sub_tree);
		}
	}
    
	return encode_json(\%taxonomic_tree);
}

sub user_is_admin {
	return $user->is_admin;
}
