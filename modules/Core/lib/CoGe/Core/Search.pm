package CoGe::Core::Search;

use strict;
use warnings;

use CoGeX;
use CoGeDBI;
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
		push @search, { ($query_terms->[$i]{"not_like"} ? '-and' : '-or') => \@items };
	}
	warn Dumper \@search;
	return \@search;
}

sub do_search {
    my ($db, $user, $type, $search, $attributes, $deleted, $restricted, $metadata_key, $metadata_value, $role) = @_;
    
    my $and = [];
	push @$and, $search if $search;
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
		push @$and, { 'user_connectors.parent_id' => $user->id };		
		push @$and, { 'user_connectors.role_id' => $role eq 'owner' ? 2 : $role eq 'editor' ? 3 : 4 };		
	}
	return undef if !@$and;
	return $db->resultset($type)->search_rs({ -and => $and }, $attributes);
}

sub info_cmp {
	my $info_a = lc($a->info);
	my $info_b = lc($b->info);
	$info_a = substr($info_a, 6) if substr($info_a, 0, 6) eq '&reg; ';
	$info_b = substr($info_b, 6) if substr($info_b, 0, 6) eq '&reg; ';
	$info_a cmp $info_b;
}

sub name_cmp {
	my $name_a = lc($a->{'name'});
	my $name_b = lc($b->{'name'});
	$name_a = substr($name_a, 6) if substr($name_a, 0, 6) eq '&reg; ';
	$name_b = substr($name_b, 6) if substr($name_b, 0, 6) eq '&reg; ';
	$name_a cmp $name_b;
}

sub push_results {
	my ($results, $objects, $type, $user, $access_method) = @_;
	foreach ( @$objects ) {
		if (!$user || $user->$access_method($_)) {
			push @$results, {
				'type'          => $type,
				'name'          => $_->info,
				'id'            => int $_->id,
				'deleted'       => $_->deleted ? Mojo::JSON->true : Mojo::JSON->false,
				'restricted'    => $_->restricted ? Mojo::JSON->true : Mojo::JSON->false,
			};
		}
	}
}

sub search {
	my %opts = @_;

    my $search_term = $opts{search_term};
	my $db 			= $opts{db};
	my $user		= $opts{user};
	my $show_users	= $opts{show_users};

    my @organism_ids;
    my @results;
	my $deleted;
    my $metadata_key;
    my $metadata_value;
	my $restricted;
	my $role;
	my $type = "none";

    my @query_terms = parse_line('\s+', 0, $search_term);
	my @query_ids;
    for ( my $i = 0 ; $i < @query_terms ; $i++ ) {
    	my $handled = 0;
        if ( index( $query_terms[$i], "::" ) != -1 ) {
            my @splitTerm = split( '::', $query_terms[$i] );
			if ($splitTerm[0] eq 'deleted') {
				$deleted = $splitTerm[1];
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
		if ($query_terms[$i] =~ /^\d+$/) {
			push @query_ids, $query_terms[$i];
		}
        if (index( $query_terms[$i], '!' ) == -1) {
            $query_terms[$i] = { 'like', '%' . $query_terms[$i] . '%' };
        } else {
        	my $bang_index = index($query_terms[$i], '!');
        	my $new_term = substr($query_terms[$i], 0, $bang_index) . substr($query_terms[$i], $bang_index + 1);
        	$query_terms[$i] = { 'not_like', '%' . $new_term . '%' };
        }
    }

	#only show public data if $user is undefined
	$restricted = 0 if !$user;

    # organisms
	if (($type eq 'none' || $type eq 'organism' || $type eq 'genome') && scalar @query_terms) {
        my @orgArray;
        for ( my $i = 0 ; $i < @query_terms ; $i++ ) {
            my $and_or;
            if ($query_terms[$i]{"not_like"}) {
                $and_or = "-and";
			} else {
				$and_or = "-or";
			}
			push @orgArray, {
				$and_or => [
					name        => $query_terms[$i],
					description => $query_terms[$i],
					organism_id => $query_terms[$i]
				]
			};
		}
		my @organisms = $db->resultset("Organism")->search( { -and => \@orgArray } );
		if ($type ne 'genome') {
			foreach ( sort { lc($a->name) cmp lc($b->name) } @organisms ) {
				push @results, {
					'type' => "organism",
					'name' => $_->name,
					'id' => int $_->id,
					'description' => $_->description
				};
			}
		}
		@organism_ids = map { $_->id } @organisms unless $type eq 'organism';
	}

    # genomes
	if ($type eq 'none' || $type eq 'genome') {
		my @genome_results;

		# get genomes linked to found organisms
		if (scalar @organism_ids) {
			my $rs = do_search($db, $user, 'Genome', @organism_ids ? { organism_id => { -in => \@organism_ids }} : undef, undef, $deleted, $restricted, $metadata_key, $metadata_value, $role);
			if ($rs) {
				my @genomes = $rs->all();
				push_results(\@genome_results, \@genomes, 'genome', $user, 'has_access_to_genome');
			}
		}

		# get genomes by id
		if (scalar @query_ids) {
			my $rs = do_search($db, $user, 'Genome', { genome_id => { -in => @query_ids }}, undef, $deleted, $restricted, $metadata_key, $metadata_value, $role);
			if ($rs) {
				my @genomes = $rs->all();
				push_results(\@genome_results, \@genomes, 'genome', $user, 'has_access_to_genome');
			}
		}
		push @results, sort name_cmp @genome_results;
	}

    # experiments
	if ($type eq 'none' || $type eq 'experiment') {
		my $search = build_search(['me.name', 'me.description', 'experiment_id', 'genome.name', 'genome.description', 'organism.name', 'organism.description'],  \@query_terms);
		my $rs = do_search($db, $user, 'Experiment', $search, { join => { 'genome' => 'organism' } }, $deleted, $restricted, $metadata_key, $metadata_value, $role);
		if ($rs) {
			my @experiments = sort info_cmp $rs->all();
			push_results(\@results, \@experiments, 'experiment', $user, 'has_access_to_experiment');
		}
	}

    # notebooks
	if ($type eq 'none' || $type eq 'notebook') {
		my @noteArray;
		for ( my $i = 0 ; $i < @query_terms ; $i++ ) {
			my $and_or;
            if ($query_terms[$i]{"not_like"}) {
                $and_or = "-and";
            } else {
                $and_or = "-or";
            }
			push @noteArray,
			  [
				$and_or => [
					name        => $query_terms[$i],
					description => $query_terms[$i],
					list_id     => $query_terms[$i]
				]
			  ];
		}
		my $rs = do_search($db, $user, 'List', scalar @noteArray ? \@noteArray : undef, undef, $deleted, $restricted, $metadata_key, $metadata_value, $role);
		if ($rs) {
			my @notebooks = sort { lc($a->name) cmp lc($b->name) } $rs->all();
			push_results(\@results, \@notebooks, 'notebook', $user, 'has_access_to_list');
		}
	}

    # user groups
	if ($show_users) {
		if ( $type eq 'none' || $type eq 'usergroup' ) {
			my @usrGArray;
			for ( my $i = 0 ; $i < @query_terms ; $i++ ) {
				my $and_or;
	            if ($query_terms[$i]{"not_like"}) {
	                $and_or = "-and";
	            } else {
	                $and_or = "-or";
	            }
				push @usrGArray,
				  [
					$and_or => [
						name          => $query_terms[$i],
						description   => $query_terms[$i],
						user_group_id => $query_terms[$i]
					]
				  ];
			}
			my $rs = do_search($db, $user, 'UserGroup', scalar @usrGArray ? \@usrGArray : undef, undef, $deleted);
			if ($rs) {
				foreach ( sort { lc($a->name) cmp lc($b->name) } $rs->all() ) {
					push @results, {
						'type'    => "user_group",
						'name'    => $_->name,
						'id'      => int  $_->id,
						'deleted' => $_->deleted ? Mojo::JSON->true : Mojo::JSON->false,
					};
				}
			}
		}
	}
	return @results;
}

1;
