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
    @EXPORT = qw( parse_query search );
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

sub build_conditions {
	my ($columns, $terms) = @_;
	return [] unless $terms;

	my @conditions;
	foreach (@$terms) {
		my %expr;
        my $index = index($_, '!');
		if ($index == -1) {
            $expr{'like'} = '%' . $_ . '%';
        } else {
        	$expr{'not_like'} = '%' . substr($_, 0, $index) . substr($_, $index + 1) . '%';
        }
		my @items;
		foreach (@$columns) {
			push @items, $_;
			push @items, \%expr;
		}
		push @conditions, { ($expr{'not_like'} ? '-and' : '-or') => \@items };
	}
	return \@conditions;
}

sub contains_none_of {
	my ($query, $keys) = @_;
	foreach (@$keys) {
		return 0 if exists $query->{$_}; 
	}
	return 1;
}

sub do_search {
    my ($type, $conditions, $attributes, $query, $db, $user) = @_;

	push @$conditions, { 'me.certified' => $query->{'certified'} } if exists $query->{'certified'};
    push @$conditions, { 'me.deleted' => exists $query->{'deleted'} ? $query->{'deleted'} : 0 }; # don't search for deleted items by default
    push @$conditions, { 'me.restricted' => $query->{'restricted'} } if exists $query->{'restricted'};
	my $metadata_key = $query->{'metadata_key'};
	my $metadata_value = $query->{'metadata_value'};
	if ($metadata_key || $metadata_value) {
		$attributes = add_join($attributes, lc($type) . '_annotations');
		if ($metadata_key) {
			my $dbh = $db->storage->dbh;
			my @row = $dbh->selectrow_array('SELECT annotation_type_id FROM annotation_type WHERE name=' . $dbh->quote($metadata_key));
			push @$conditions, { lc($type) . '_annotations.annotation_type_id' => $row[0] };
		}
		if ($metadata_value) {
			push @$conditions, { 'annotation' => $metadata_value };
		}
	}
	if ($query->{'tag'}) {
		$attributes = add_join($attributes, 'experiment_type_connectors');
		my $dbh = $db->storage->dbh;
		my @row = $dbh->selectrow_array('SELECT experiment_type_id FROM experiment_type WHERE name=' . $dbh->quote($query->{'tag'}));
		push @$conditions, { 'experiment_type_connectors.experiment_type_id' => $row[0] };
	}
	my $role = $query->{'role'};
	if ($role && $user) {
		$attributes = add_join($attributes, 'user_connectors');
		push @$conditions, { 'user_connectors.child_type' => $type eq 'Genome' ? 2 : $type eq 'Experiment' ? 3 : 1 };
		push @$conditions, { 'user_connectors.parent_type' => 5 };
		push @$conditions, { 'user_connectors.parent_id' => $user->id };		
		push @$conditions, { 'user_connectors.role_id' => $role eq 'owner' ? 2 : $role eq 'editor' ? 3 : 4 };		
	}
	if (exists $query->{'favorite'} && $user && !$user->is_public) { # only allow searching for favorites, specifying favorite::0 is not allowed
		$attributes = add_join($attributes, 'favorite_connectors');
		push @$conditions, { 'favorite_connectors.child_type' => $type eq 'Genome' ? 2 : $type eq 'Experiment' ? 3 : 1 };
		push @$conditions, { 'favorite_connectors.user_id' => $user->id };
	}
	# my ($sql, @b) = $db->resultset($type)->search_rs({ -and => $conditions }, $attributes)->as_query();
	# warn Dumper $sql;
	# warn Dumper \@b;
	return $db->resultset($type)->search_rs({ -and => $conditions }, $attributes);
}

sub info_cmp {
	my $info_a = lc($a->info);
	my $info_b = lc($b->info);
	$info_a = substr($info_a, 6) if substr($info_a, 0, 2) eq '🔒 ';
	$info_b = substr($info_b, 6) if substr($info_b, 0, 2) eq '🔒 ';
	$info_a cmp $info_b;
}

sub parse_query {
	my $search_text = shift;
	my $query = {};
	my @search_terms;
	foreach (parse_line('\s+', 0, $search_text)) {
		next if !$_;
		my $index = index($_, '::');
		if ($index != -1) {
			$query->{substr($_, 0, $index)} = substr($_, $index + 2);
		}
		else {
			push @search_terms, $_;
		}
    }
	$query->{'search_terms'} = \@search_terms if @search_terms;
	return $query;
}

sub push_results {
	my ($results, $objects, $type, $user, $access_method, $favorites) = @_;
	foreach (@$objects) {
		if (!$access_method || !$user || $user->$access_method($_)) {
			my $result = {
				'type'          => $type,
				'name'          => $_->info(hideRestrictedSymbol=>1),
				'id'            => int $_->id
			};
			$result->{'certified'}  = $_->certified ? Mojo::JSON->true : Mojo::JSON->false if $_->can('certified');
			$result->{'deleted'}    = $_->deleted ? Mojo::JSON->true : Mojo::JSON->false if $_->can('deleted');
			$result->{'favorite'}   = ($favorites->is_favorite($_) ? Mojo::JSON->true : Mojo::JSON->false) if defined $favorites;
			if ($_->can('restricted')) {
				$result->{'restricted'} = $_->restricted ? Mojo::JSON->true : Mojo::JSON->false;
				next if (!$user && $result->{'restricted'}); # mdb added 1/3/17 COGE-809 -- remove from results for public user
			}
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

	my $query = parse_query($search_term, $user);
	my $type = $query->{'type'};

    my @results;
    my $favorites;
	$favorites = CoGe::Core::Favorites->new(user => $user) if ($user && !$user->is_public && (!$type || $type eq 'genome' || $type eq 'experiment' || $type eq 'notebook'));

    # Genomes
	if ((!$type || $type eq 'genome') && contains_none_of($query, ['feature_type', 'tag'])) {
		my $conditions = build_conditions(['me.name', 'me.description', 'me.genome_id', 'organism.name', 'organism.description'],  $query->{'search_terms'});
		my $rs = do_search('Genome', $conditions, $conditions ? { join => 'organism' } : undef, $query, $db, $user);
		if ($rs) {
			my @genomes = sort info_cmp $rs->all();
			push_results(\@results, \@genomes, 'genome', $user, 'has_access_to_genome', $favorites);
		}
	}

    # Organisms
	if ((!$type || $type eq 'organism') && $query->{'search_terms'} && contains_none_of($query, ['certified', 'deleted', 'favorite', 'feature_type', 'restricted', 'metadata_key', 'metadata_value', 'role', 'tag'])) {
		my $conditions = build_conditions(['me.name', 'me.description', 'me.organism_id'],  $query->{'search_terms'});
		my @organisms = $db->resultset("Organism")->search( { -and => $conditions } );
		foreach ( sort { lc($a->name) cmp lc($b->name) } @organisms ) {
			push @results, {
				'type' => 'organism',
				'name' => $_->name,
				'id' => int $_->id,
				'description' => $_->description
			};
		}
	}

    # Experiments
	if ((!$type || $type eq 'experiment') && contains_none_of($query, ['certified', 'feature_type'])) {
		my $conditions = build_conditions(['me.name', 'me.description', 'me.experiment_id', 'genome.name', 'genome.description', 'organism.name', 'organism.description'],  $query->{'search_terms'});
		my $rs = do_search('Experiment', $conditions, $conditions ? { join => { 'genome' => 'organism' } } : undef, $query, $db, $user);
		if ($rs) {
			my @experiments = sort info_cmp $rs->all();
			push_results(\@results, \@experiments, 'experiment', $user, 'has_access_to_experiment', $favorites);
		}
	}

    # Notebooks
	if ((!$type || $type eq 'notebook') && contains_none_of($query, ['certified', 'feature_type', 'tag'])) {
		my $conditions = build_conditions(['name', 'description', 'list_id'],  $query->{'search_terms'});
		my $rs = do_search('List', $conditions, undef, $query, $db, $user);
		if ($rs) {
			my @notebooks = sort { lc($a->name) cmp lc($b->name) } $rs->all();
			push_results(\@results, \@notebooks, 'notebook', $user, 'has_access_to_notebook', $favorites);
		}
	}

    # User groups
	if ($show_users && (!$type || $type eq 'usergroup') && contains_none_of($query, ['certified', 'favorite', 'feature_type', 'restricted', 'metadata_key', 'metadata_value', 'role', 'tag'])) {
		my $conditions = build_conditions(['name', 'description', 'user_group_id'],  $query->{'search_terms'});
		my $rs = do_search('UserGroup', $conditions, undef, $query, $db, $user);
		if ($rs) {
			my @user_groups = sort info_cmp $rs->all();
			push_results(\@results, \@user_groups, 'group', $user, undef, $favorites);
		}
	}

    # Features
	if ((!$type || $type eq 'feature') && $query->{'search_terms'} && contains_none_of($query, ['certified', 'deleted', 'favorite', 'restricted', 'metadata_key', 'metadata_value', 'role', 'tag'])) {
		my $dbh = $db->storage->dbh;
		my $feature_type = $query->{'feature_type'};

		#FIXME move literal query into CoGeDBI, mdb 12/22/16
		my $sql = 'SELECT feature_name.name,feature.feature_id,' . ($feature_type ? "'" . $feature_type . "'" : 'feature_type.name') . ',organism.name,data_source.name,genome.version,genomic_sequence_type.name ' .
			'FROM feature_name ' .
				'JOIN feature USING(feature_id) ';
		$sql .= 'JOIN feature_type USING(feature_type_id) ' unless $feature_type;
		$sql .= 'JOIN dataset USING(dataset_id) ' .
				'JOIN data_source USING(data_source_id) ' .
				'JOIN dataset_connector USING(dataset_id) ' .
				'JOIN genome ON dataset_connector.genome_id=genome.genome_id AND !genome.deleted ';
		$sql .= 'AND !genome.restricted ' if !$user || $user->is_public;
		$sql .= 'JOIN organism USING(organism_id) ' .
				'JOIN genomic_sequence_type USING(genomic_sequence_type_id) ' .
			'WHERE MATCH(feature_name.name) AGAINST (\'' . (join ',', @{$query->{'search_terms'}}) . '\') ';
		if ($feature_type) {
			my @row = $dbh->selectrow_array('SELECT feature_type_id FROM feature_type WHERE name=\'' . $feature_type . '\'');
			$sql .= 'AND feature.feature_type_id=' . $row[0] . ' ';
		}
		$sql .= 'GROUP BY feature_name.name,feature.feature_id';
		my $rows = $dbh->selectall_arrayref($sql);

		foreach (@$rows) { #TODO use fetchall_hashref and map for performance improvement, mdb 12/22/16
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
