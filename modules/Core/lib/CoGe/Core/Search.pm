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

sub search {
	my %opts = @_;
	
    my $search_term = $opts{search_term};	# Takes in the entire string, to be processed later
	my $db 			= $opts{db};
	my $user		= $opts{user};
	my $show_users	= $opts{show_users};

    my @idList;
    my @results;
    my $metadata_key;

	my $type = "none";
	my @restricted =
	  [ -or => [ 'me.restricted' => 0, 'me.restricted' => 1 ] ]; 
	  #Specifies either restricted(1) OR unrestriced(0) results. Default is all.
	my @deleted =
	  [ -or => [ 'me.deleted' => 0, 'me.deleted' => 1 ] ]; 
	  #Specifies either deleted(1) OR non-deleted(0) results. Default is all.

    my @searchArray = parse_line('\s+', 0, $search_term);
    for ( my $i = 0 ; $i < @searchArray ; $i++ ) {
    	my $handled = 0;
        if ( index( $searchArray[$i], "::" ) != -1 ) {
            my @splitTerm = split( '::', $searchArray[$i] );
			if ($splitTerm[0] eq 'deleted') {
				@deleted = [ deleted => $splitTerm[1] ];
				if ( $type eq "none" ) {
					$type = 'deleted'; 
					#Sets the "type" so that only fields relevant to deletion are shown. (i.e. there are no deleted users).
				}
				$handled = 1;
			}
			elsif (index($splitTerm[0], 'metadata_key') != -1) {
				$type = $splitTerm[0];
				$metadata_key = $splitTerm[1];
				$handled = 1;
			}
			elsif ($splitTerm[0] eq 'restricted') {
				@restricted = [ 'me.restricted' => $splitTerm[1] ];
				if ( $type eq "none" ) {
					$type = 'restricted'; 
					#Sets the "type" so that only fields relevant to restriction are shown. (i.e. there are no restricted users.)
				}
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
	if (!$user) {
		$type = "restricted";
		@restricted = [ 'me.restricted' => 0 ];
	}

    # Perform organism search
    if (   $type eq 'none'
        || $type eq 'organism'
        || $type eq 'genome'
        || $type eq 'restricted'
        || $type eq 'deleted' )
    {
        my @orgArray;
        for ( my $i = 0 ; $i < @searchArray ; $i++ ) {
            my $and_or;
            if ($searchArray[$i]{"not_like"}) {
                $and_or = "-and";
			} else {
				$and_or = "-or";
			}
			push @orgArray,
			  [
				$and_or => [
					name        => $searchArray[$i],
					description => $searchArray[$i],
					organism_id => $searchArray[$i]
				]
			  ];
		}
		my @organisms = $db->resultset("Organism")->search( { -and => [ @orgArray, ], } );

		#if ( $type eq 'none' || $type eq 'organism' ) { # mdb removed 2/18/16 -- wasn't showing organism results unless logged in
			foreach ( sort { lc($a->name) cmp lc($b->name) } @organisms ) {
				push @results, { 
					'type' => "organism", 
				  	'name' => $_->name, 
				  	'id' => int $_->id, 
				  	'description' => $_->description 
				};
			}
		#}
		@idList = map { $_->id } @organisms;
	}

	# Perform user search
	if ($show_users) {
		if ( $type eq 'none' || $type eq 'user' ) {
			my @usrArray;
			for ( my $i = 0 ; $i < @searchArray ; $i++ ) {
				my $and_or;
	            if ($searchArray[$i]{"not_like"}) {
	                $and_or = "-and";
	            } else {
	                $and_or = "-or";
	            }
				push @usrArray,
				  [
					$and_or => [
						user_name  => $searchArray[$i],
						first_name => $searchArray[$i],
						user_id    => $searchArray[$i],
						last_name => $searchArray[$i]
					]
				  ];
			}
			my @users = $db->resultset("User")->search( { -and => [ @usrArray, ], } );
	
			foreach ( sort { lc($a->user_name) cmp lc($b->user_name) } @users ) {
				push @results, { 
					'type' => "user", 
					'first' => $_->first_name, 
					'last' => $_->last_name, 
					'username' => $_->user_name, 
					'id' => int $_->id, 
					'email' => $_->email 
				};
			}
		}
	}
	

	# Perform genome search (corresponding to Organism results)
	if (   $type eq 'none'
		|| $type eq 'genome'
		|| $type eq 'genome_metadata_key'
		|| $type eq 'restricted'
		|| $type eq 'deleted' )
	{
		my $attributes;
		my $search;
		if ($type eq 'genome_metadata_key') {
			$attributes = { distinct => 1, join => 'genome_annotations' };
			my $dbh = $db->storage->dbh;
			my @row = $dbh->selectrow_array('SELECT annotation_type_id FROM annotation_type WHERE name=' . $dbh->quote($metadata_key));
			$search = { 'genome_annotations.annotation_type_id' => $row[0] };
		}
		else {
			$search = {
				-and => [
					organism_id => { -in => \@idList },
					@restricted,
					@deleted,
				]
			};
		}
		my @genome_results;
		my @genomes = $db->resultset("Genome")->search( $search, $attributes );

		foreach ( @genomes ) {
			if (!$user || $user->has_access_to_genome($_)) {
				push @genome_results, {
					'type'          => "genome",
					'name'          => $_->info,
					'id'            => int $_->id,
					'deleted'       => $_->deleted ? Mojo::JSON->true : Mojo::JSON->false,
					'restricted' 	=> $_->restricted ? Mojo::JSON->true : Mojo::JSON->false,
				};
			}
		}

		# Perform direct genome search (by genome ID)
		if (scalar @searchArray) {
			my @genIDArray;
			for ( my $i = 0 ; $i < @searchArray ; $i++ ) {
				push @genIDArray, [ -or => [ genome_id => $searchArray[$i] ] ];
			}
			my @genomeIDs = $db->resultset("Genome")->search( { -and => [ @genIDArray, @restricted, @deleted, ], } );

			foreach ( @genomeIDs ) {
				if (!$user || $user->has_access_to_genome($_)) {
					push @genome_results, {
						'type'          => "genome",
						'name'          => $_->info,
						'id'            => int $_->id,
						'deleted'       => $_->deleted ? Mojo::JSON->true : Mojo::JSON->false,
						'restricted'    => $_->restricted ? Mojo::JSON->true : Mojo::JSON->false,
					};
				}
			}
		}
		push @results, sort name_cmp @genome_results;
	}

	# Perform experiment search
	if (   $type eq 'none'
		|| $type eq 'experiment'
		|| $type eq 'experiment_metadata_key'
		|| $type eq 'restricted'
		|| $type eq 'deleted' )
	{
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
		my $attributes;
		my $search;
		if ($type eq 'experiment_metadata_key') {
			$attributes = { distinct => 1, join => [ 'experiment_annotations' ] };
			my $dbh = $db->storage->dbh;
			my @row = $dbh->selectrow_array('SELECT annotation_type_id FROM annotation_type WHERE name=' . $dbh->quote($metadata_key));
			$search = { 'experiment_annotations.annotation_type_id' => $row[0] };
		}
		else {
		    $attributes = { join => { 'genome' => 'organism' } };
			$search = { -and => [ @expArray, @restricted, @deleted] };
		}
		my @experiments = $db->resultset("Experiment")->search( $search, $attributes );

		foreach ( sort info_cmp @experiments ) {
			if (!$user || $user->has_access_to_experiment($_)) {
				push @results, {
					'type'          => "experiment",
					'name'          => $_->info,
					'id'            => int $_->id,
					'deleted'       => $_->deleted ? Mojo::JSON->true : Mojo::JSON->false,
					'restricted'    => $_->restricted ? Mojo::JSON->true : Mojo::JSON->false,
				};
			}
		}
	}

	# Perform notebook search
	if (   $type eq 'none'
		|| $type eq 'notebook'
		|| $type eq 'list_metadata_key'
		|| $type eq 'restricted'
		|| $type eq 'deleted' )
	{
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
		my $attributes;
		my $search;
		if ($type eq 'list_metadata_key') {
			$attributes = { distinct => 1, join => 'list_annotations' };
			my $dbh = $db->storage->dbh;
			my @row = $dbh->selectrow_array('SELECT annotation_type_id FROM annotation_type WHERE name=' . $dbh->quote($metadata_key));
			$search = { 'list_annotations.annotation_type_id' => $row[0] };
		}
		else {
			$search = { -and => [ @noteArray, @restricted, @deleted, ] };
		}
		my @notebooks = $db->resultset("List")->search( $search, $attributes );

		foreach ( sort { lc($a->name) cmp lc($b->name) } @notebooks ) {
			if (!$user || $user->has_access_to_list($_)) {
				push @results, {
					'type'          => "notebook",
					'name'          => $_->info,
					'id'            => int $_->id,
					'deleted'       => $_->deleted ? Mojo::JSON->true : Mojo::JSON->false,
					'restricted'    => $_->restricted ? Mojo::JSON->true : Mojo::JSON->false,
				};
			}
		}
	}

	# Perform user group search
	if ($show_users) {
		if ( $type eq 'none' || $type eq 'usergroup' || $type eq 'deleted' ) {
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
			my @userGroup = $db->resultset("UserGroup")->search( { -and => [ @usrGArray, @deleted, ], } );
	
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
	#return encode_json( { timestamp => $timestamp, items => \@results } );
	return @results;
}

1;
