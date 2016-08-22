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
		$attributes->{'join'} = { join => $new_join };
	}
	$attributes->{'distinct'} = 1;
	return $attributes;
}

sub do_search {
    my ($db, $user, $type, $search, $attributes, $deleted, $restricted, $metadata_key, $metadata_value, $role) = @_;
    
    my $and = [];
	push @$and,  $search if $search;
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
	warn Dumper $attributes;
	return $db->resultset($type)->search({ -and => $and }, $attributes);
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

    my $search_term = $opts{search_term};	# Takes in the entire string, to be processed later
	my $db 			= $opts{db};
	my $user		= $opts{user};
	my $show_users	= $opts{show_users};

    my @idList;
    my @results;
	my $deleted;
    my $metadata_key;
    my $metadata_value;
	my $restricted;
	my $role;
	my $type = "none";

    my @searchArray = parse_line('\s+', 0, $search_term);
    for ( my $i = 0 ; $i < @searchArray ; $i++ ) {
    	my $handled = 0;
        if ( index( $searchArray[$i], "::" ) != -1 ) {
            my @splitTerm = split( '::', $searchArray[$i] );
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
            splice( @searchArray, $i, 1 );
            $i--;
        }
        elsif ( index( $searchArray[$i], '!' ) == -1 ) {
            $searchArray[$i] = { 'like', '%' . $searchArray[$i] . '%' };
        } else {
        	my $bang_index = index($searchArray[$i], '!');
        	my $new_term = substr($searchArray[$i], 0, $bang_index) . substr($searchArray[$i], $bang_index + 1);
        	$searchArray[$i] = { 'not_like', '%' . $new_term . '%' };
        }
    }

	#only show public data if $user is undefined
	$restricted = 1 if !$user;

    # Perform organism search
    if (($type eq 'none' || $type eq 'organism' || $type eq 'genome') && scalar @searchArray) {
        my @orgArray;
        for ( my $i = 0 ; $i < @searchArray ; $i++ ) {
            my $and_or;
            if ($searchArray[$i]{"not_like"}) {
                $and_or = "-and";
			} else {
				$and_or = "-or";
			}
			push @orgArray, {
				$and_or => [
					name        => $searchArray[$i],
					description => $searchArray[$i],
					organism_id => $searchArray[$i]
				]
			};
		}
		my @organisms = $db->resultset("Organism")->search( { -and => @orgArray } );
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
		@idList = map { $_->id } @organisms;
	}

	# Perform user search
	# if ($show_users) {
	# 	if ( $type eq 'none' || $type eq 'user' ) {
	# 		my @usrArray;
	# 		for ( my $i = 0 ; $i < @searchArray ; $i++ ) {
	# 			my $and_or;
	#             if ($searchArray[$i]{"not_like"}) {
	#                 $and_or = "-and";
	#             } else {
	#                 $and_or = "-or";
	#             }
	# 			push @usrArray,
	# 			  [
	# 				$and_or => [
	# 					user_name  => $searchArray[$i],
	# 					first_name => $searchArray[$i],
	# 					user_id    => $searchArray[$i],
	# 					last_name => $searchArray[$i]
	# 				]
	# 			  ];
	# 		}
	# 		my @users = $db->resultset("User")->search( { -and => [ @usrArray, ], } );
    #
	# 		foreach ( sort { lc($a->user_name) cmp lc($b->user_name) } @users ) {
	# 			push @results, {
	# 				'type' => "user",
	# 				'first' => $_->first_name,
	# 				'last' => $_->last_name,
	# 				'username' => $_->user_name,
	# 				'id' => int $_->id,
	# 				'email' => $_->email
	# 			};
	# 		}
	# 	}
	# }

	# Perform genome search (corresponding to Organism results)
	if ($type eq 'none' || $type eq 'genome') {
		my @genome_results;
		my @genomes = do_search($db, $user, 'Genome', @idList ? { organism_id => { -in => \@idList }} : undef, undef, $deleted, $restricted, $metadata_key, $metadata_value, $role);
		push_results(\@genome_results, \@genomes, 'genome', $user, 'has_access_to_genome');

		# Perform direct genome search (by genome ID)
		if (scalar @searchArray) {
			my @genomes = do_search($db, $user, 'Genome', { genome_id => { -in => grep(/d+/, @searchArray) }}, undef, $deleted, $restricted, $metadata_key, $metadata_value, $role);
			push_results(\@genome_results, \@genomes, 'genome', $user, 'has_access_to_genome');
		}
		push @results, sort name_cmp @genome_results;
	}

	# Perform experiment search
	if ($type eq 'none' || $type eq 'experiment') {
		my @expArray;
		for ( my $i = 0 ; $i < @searchArray ; $i++ ) {
			my $and_or;
            if ($searchArray[$i]{"not_like"}) {
                $and_or = "-and";
            } else {
                $and_or = "-or";
            }
			push @expArray,
			  [
				$and_or => [
					'me.name'        => $searchArray[$i],
					'me.description' => $searchArray[$i],
					'experiment_id'  => $searchArray[$i],
					'genome.name'  => $searchArray[$i],
                    'genome.description' => $searchArray[$i],
					'organism.name'  => $searchArray[$i],
					'organism.description' => $searchArray[$i]
				]
			  ];
		}
		my @experiments = do_search($db, $user, 'Experiment', scalar @expArray ? \@expArray : undef, { join => { 'genome' => 'organism' } }, $deleted, $restricted, $metadata_key, $metadata_value, $role);
		@experiments = sort info_cmp @experiments;
		push_results(\@results, \@experiments, 'experiment', $user, 'has_access_to_experiment');
	}

	# Perform notebook search
	if ($type eq 'none' || $type eq 'notebook') {
		my @noteArray;
		for ( my $i = 0 ; $i < @searchArray ; $i++ ) {
			my $and_or;
            if ($searchArray[$i]{"not_like"}) {
                $and_or = "-and";
            } else {
                $and_or = "-or";
            }
			push @noteArray,
			  [
				$and_or => [
					name        => $searchArray[$i],
					description => $searchArray[$i],
					list_id     => $searchArray[$i]
				]
			  ];
		}
		my @notebooks = do_search($db, $user, 'List', scalar @noteArray ? \@noteArray : undef, undef, $deleted, $restricted, $metadata_key, $metadata_value, $role);
		@notebooks = sort { lc($a->name) cmp lc($b->name) } @notebooks;
		push_results(\@results, \@notebooks, 'notebook', $user, 'has_access_to_list');
	}

	# Perform user group search
	if ($show_users) {
		if ( $type eq 'none' || $type eq 'usergroup' ) {
			my @usrGArray;
			for ( my $i = 0 ; $i < @searchArray ; $i++ ) {
				my $and_or;
	            if ($searchArray[$i]{"not_like"}) {
	                $and_or = "-and";
	            } else {
	                $and_or = "-or";
	            }
				push @usrGArray,
				  [
					$and_or => [
						name          => $searchArray[$i],
						description   => $searchArray[$i],
						user_group_id => $searchArray[$i]
					]
				  ];
			}
			my @userGroup = $db->resultset("UserGroup")->search( build_search(@usrGArray, $deleted) );

			foreach ( sort { lc($a->name) cmp lc($b->name) } @userGroup ) {
				push @results, {
					'type'    => "user_group",
					'name'    => $_->name,
					'id'      => int  $_->id,
					'deleted' => $_->deleted ? Mojo::JSON->true : Mojo::JSON->false,
				};
			}
		}
	}
	return @results;
}

1;
