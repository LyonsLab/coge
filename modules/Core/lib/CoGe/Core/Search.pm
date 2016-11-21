package CoGe::Core::Search;

use strict;
use warnings;

use CoGeX;
use CoGeDBI;
use CoGe::Core::Favorites;
use CoGe::Core::Feature qw( search_features );
use Data::Dumper;
use Text::ParseWords;

BEGIN {
    our ( @EXPORT, @ISA, $VERSION );
    require Exporter;

    $VERSION = 0.1;
    @ISA = qw( Exporter );
    @EXPORT = qw( search );
}

sub add_join {
	my $attributes = shift || {};
	my $new_join = shift;

	my $join = $attributes->{'join'};
	if ($join) {
		if (ref($join) eq 'ARRAY') {
			push @$join, $new_join;
		}
		else {
			$attributes->{'join'} = [$join, $new_join];
		}
	}
	else {
		$attributes->{'join'} = $new_join;
	}
	$attributes->{'distinct'} = 1;
	return $attributes;
}

sub build_search {
	my ($columns, $query_terms) = @_;
	return undef unless @$query_terms;

	my @search;
	for ( my $i = 0 ; $i < @$query_terms ; $i++ ) {
		my @items;
		for (@$columns) {
			push @items, $_;
			push @items, $query_terms->[$i];
		}
		push @search, { ($query_terms->[$i]{'not_like'} ? '-and' : '-or') => \@items };
	}
	return \@search;
}

sub do_search {
    my ($db, $user, $type, $search, $attributes, $deleted, $restricted, $metadata_key, $metadata_value, $role, $favorite, $certified) = @_;
    
    my $and = $search || [];
    push @$and, { 'me.certified' => $certified } if defined $certified;
    push @$and, { 'me.deleted' => $deleted } if defined $deleted;
    push @$and, { 'me.restricted' => $restricted } if defined $restricted;
	if ($metadata_key || $metadata_value) {
		$attributes = add_join($attributes, lc($type) . '_annotations');
		if ($metadata_key) {
			my $dbh = $db->storage->dbh;
			my @row = $dbh->selectrow_array('SELECT annotation_type_id FROM annotation_type WHERE name=' . $dbh->quote($metadata_key));
			push @$and, { lc($type) . '_annotations.annotation_type_id' => $row[0] };
		}
		if ($metadata_value) {
			push @$and, { 'annotation' => $metadata_value };
		}
	}
	if ($role && $user) {
		$attributes = add_join($attributes, 'user_connectors');
		push @$and, { 'user_connectors.child_type' => $type eq 'Genome' ? 2 : $type eq 'Experiment' ? 3 : 1 };
		push @$and, { 'user_connectors.parent_type' => 5 };
		push @$and, { 'user_connectors.parent_id' => $user->id };		
		push @$and, { 'user_connectors.role_id' => $role eq 'owner' ? 2 : $role eq 'editor' ? 3 : 4 };		
	}
	if ($favorite) { # only allow searching for favorites, specifying favorite::0 is not allowed
		$attributes = add_join($attributes, 'favorite_connectors');
		push @$and, { 'favorite_connectors.child_type' => $type eq 'Genome' ? 2 : $type eq 'Experiment' ? 3 : 1 };
		push @$and, { 'favorite_connectors.user_id' => $user->id };
	}
	# my ($sql, @b) = $db->resultset($type)->search_rs({ -and => $and }, $attributes)->as_query();
	# warn Dumper $sql;
	# warn Dumper \@b;
	return $db->resultset($type)->search_rs({ -and => $and }, $attributes);
}

sub info_cmp {
	my $info_a = lc($a->info);
	my $info_b = lc($b->info);
	$info_a = substr($info_a, 6) if substr($info_a, 0, 2) eq '🔒 ';
	$info_b = substr($info_b, 6) if substr($info_b, 0, 2) eq '🔒 ';
	$info_a cmp $info_b;
}

sub push_results {
	my ($results, $objects, $type, $user, $access_method, $favorites) = @_;
	foreach ( @$objects ) {
		if (!$access_method || !$user || $user->$access_method($_)) {
			my $result = {
				'type'          => $type,
				'name'          => $_->info(hideRestrictedSymbol=>1),
				'id'            => int $_->id
			};
			$result->{'certified'} = $_->certified ? Mojo::JSON->true : Mojo::JSON->false if $_->can('certified');
			$result->{'deleted'} = $_->deleted ? Mojo::JSON->true : Mojo::JSON->false if $_->can('deleted');
			$result->{'favorite'} = ($favorites->is_favorite($_) ? Mojo::JSON->true : Mojo::JSON->false) if defined($favorites) && ($type eq 'genome' || $type eq 'experiment' || $type eq 'notebook');
			$result->{'restricted'} = $_->restricted ? Mojo::JSON->true : Mojo::JSON->false if $_->can('restricted');
			
			push @$results, $result;
		}
	}
}

sub search {
	my %opts = @_;

    my $search_term = $opts{search_term};
	my $db 			= $opts{db};
	my $user		= $opts{user};
	my $show_users	= $opts{show_users};

    my @results;
	my $certified;
	my $deleted;
	my $favorite;
    my $favorites;
	$favorites = CoGe::Core::Favorites->new(user => $user) if ($user && !$user->is_public);
	my $metadata_key;
    my $metadata_value;
	my $restricted;
	my $role;
	my $type = "none";

	my @search_terms;
    my @query_terms = parse_line('\s+', 0, $search_term);
	for ( my $i = 0 ; $i < @query_terms ; $i++ ) {
    	my $handled = 0;
        if ( index( $query_terms[$i], "::" ) != -1 ) {
            my @splitTerm = split( '::', $query_terms[$i] );
			if ($splitTerm[0] eq 'certified') {
				$certified = $splitTerm[1];
				$handled = 1;
			}
			elsif ($splitTerm[0] eq 'deleted') {
				$deleted = $splitTerm[1];
				$handled = 1;
			}
			elsif ($splitTerm[0] eq 'favorite') {
				$favorite = $splitTerm[1];
				$handled = 1;
			}
			elsif ($splitTerm[0] eq 'metadata_key') {
				$metadata_key = $splitTerm[1];
				$handled = 1;
			}
			elsif ($splitTerm[0] eq 'metadata_value') {
				$metadata_value = $splitTerm[1];
				$handled = 1;
			}
			elsif ($splitTerm[0] eq 'restricted') {
				$restricted = $splitTerm[1];
				$handled = 1;
			}
			elsif ($splitTerm[0] eq 'role') {
				$role = $splitTerm[1];
				$handled = 1;
			}
			elsif ($splitTerm[0] eq 'type') {
				$type = $splitTerm[1];
				$handled = 1;
			}
        }
        if ($handled) {
            splice( @query_terms, $i, 1 );
            $i--;
			next;
        }
		push @search_terms, $query_terms[$i];
		if (index( $query_terms[$i], '!' ) == -1) {
            $query_terms[$i] = { 'like', '%' . $query_terms[$i] . '%' };
        } else {
        	my $bang_index = index($query_terms[$i], '!');
        	my $new_term = substr($query_terms[$i], 0, $bang_index) . substr($query_terms[$i], $bang_index + 1);
        	$query_terms[$i] = { 'not_like', '%' . $new_term . '%' };
        }
    }

    # organisms
	if (($type eq 'none' || $type eq 'organism') && @query_terms && !defined($certified) && !defined($deleted) && !defined($favorite) && !defined($restricted) && !defined($metadata_key) && !defined($metadata_value) && !defined($role)) {
		my $search = build_search(['me.name', 'me.description', 'me.organism_id'],  \@query_terms);
		my @organisms = $db->resultset("Organism")->search( { -and => $search } );
		foreach ( sort { lc($a->name) cmp lc($b->name) } @organisms ) {
			push @results, {
				'type' => 'organism',
				'name' => $_->name,
				'id' => int $_->id,
				'description' => $_->description
			};
		}
	}

    # genomes
	if ($type eq 'none' || $type eq 'genome') {
		my $search = build_search(['me.name', 'me.description', 'me.genome_id', 'organism.name', 'organism.description'],  \@query_terms);
		my $rs = do_search($db, $user, 'Genome', $search, $search ? { join => 'organism' } : undef, $deleted || 0, $user ? $restricted : 0, $metadata_key, $metadata_value, $role, $favorite, $certified);
		if ($rs) {
			my @genomes = sort info_cmp $rs->all();
			push_results(\@results, \@genomes, 'genome', $user, 'has_access_to_genome', $favorites);
		}
	}

    # experiments
	if (($type eq 'none' || $type eq 'experiment') && !defined($certified)) {
		my $search = build_search(['me.name', 'me.description', 'me.experiment_id', 'genome.name', 'genome.description', 'organism.name', 'organism.description'],  \@query_terms);
		my $rs = do_search($db, $user, 'Experiment', $search, $search ? { join => { 'genome' => 'organism' } } : undef, $deleted || 0, $user ? $restricted : 0, $metadata_key, $metadata_value, $role, $favorite);
		if ($rs) {
			my @experiments = sort info_cmp $rs->all();
			push_results(\@results, \@experiments, 'experiment', $user, 'has_access_to_experiment', $favorites);
		}
	}

    # notebooks
	if (($type eq 'none' || $type eq 'notebook') && !defined($certified)) {
		my $search = build_search(['name', 'description', 'list_id'],  \@query_terms);
		my $rs = do_search($db, $user, 'List', $search, undef, $deleted || 0, $user ? $restricted : 0, $metadata_key, $metadata_value, $role, $favorite);
		if ($rs) {
			my @notebooks = sort { lc($a->name) cmp lc($b->name) } $rs->all();
			push_results(\@results, \@notebooks, 'notebook', $user, 'has_access_to_list', $favorites);
		}
	}

    # user groups
	if ($show_users && ($type eq 'none' || $type eq 'usergroup') && !defined($certified) && !defined($favorite) && !defined($restricted) && !defined($metadata_key) && !defined($metadata_value) && !defined($role)) {
		my $search = build_search(['name', 'description', 'user_group_id'],  \@query_terms);
		my $rs = do_search($db, $user, 'UserGroup', $search, undef, $deleted || 0);
		if ($rs) {
			my @user_groups = sort info_cmp $rs->all();
			push_results(\@results, \@user_groups, 'user_group', $user, undef, $favorites);
		}
	}

    # features
	if (($type eq 'none' || $type eq 'feature') && !defined($certified) && !defined($deleted) && !defined($favorite) && !defined($restricted) && !defined($metadata_key) && !defined($metadata_value) && !defined($role)) {
		my $sql = 'SELECT feature_name.name,feature.feature_id,feature_type.name,organism.name,data_source.name,genome.version,genomic_sequence_type.name ' .
			'FROM feature_name ' .
				'JOIN feature USING(feature_id) ' .
				'JOIN feature_type USING(feature_type_id) ' .
				'JOIN dataset USING(dataset_id) ' .
				'JOIN data_source USING(data_source_id) ' .
				'JOIN dataset_connector USING(dataset_id) ' .
				'JOIN genome ON dataset_connector.genome_id=genome.genome_id AND !genome.deleted ' .
				'JOIN organism USING(organism_id) ' .
				'JOIN genomic_sequence_type USING(genomic_sequence_type_id) ' .
			'WHERE MATCH(feature_name.name) AGAINST (\'' . (join ',', @search_terms) . '\') ' .
			'GROUP BY feature_name.name,feature.feature_id';
			warn $sql;
		my $rows = $db->storage->dbh->selectall_arrayref($sql);
		foreach (@$rows) {
			push @results, {
				type         => 'feature',
				name         => $_->[0],
				id           => int $_->[1],
				feature_type => $_->[2],
				genome       => $_->[3] ? $_->[3] . ' (' . $_->[4] . ', ' . $_->[5] . ', ' .  $_->[6] . ')' : ''
			}
		}
	}

	return @results;
}

1;
