package CoGe::Services::Data::Experiment;

use Mojo::Base 'Mojolicious::Controller';
use Data::Dumper;
#use IO::Compress::Gzip 'gzip';
use CoGeX;
use CoGe::Accessory::Utils;
use CoGe::Core::Experiment;
use CoGe::Services::Auth;
use CoGe::Services::Data::Job;

sub search {
    my $self = shift;
    my $search_term = $self->stash('term');

    # Validate input
    if (!$search_term or length($search_term) < 3) {
        $self->render(status => 400, json => { error => { Error => 'Search term is shorter than 3 characters' } });
        return;
    }

    # Authenticate user and connect to the database
    my ($db, $user) = CoGe::Services::Auth::init($self);

    # Search experiments
    my $search_term2 = '%' . $search_term . '%';
    my @experiments = $db->resultset("Experiment")->search(
        \[
            'experiment_id = ? OR name LIKE ? OR description LIKE ?',
            [ 'experiment_id', $search_term  ],
            [ 'name',        $search_term2 ],
            [ 'description', $search_term2 ]
        ]
    );

    # Filter response
    my @filtered = grep {
         !$_->restricted || (defined $user && $user->has_access_to_experiment($_))
    } @experiments;

    # Format response
    my @result = map {
      {
        id => int($_->id),
        name => $_->name,
        description => $_->description,
      }
    } @filtered;

    $self->render(json => { experiments => \@result });
}

sub fetch {
    my $self = shift;
    my $id = int($self->stash('id'));
    
    # Validate input
    unless ($id) {
        $self->render(status => 400, json => {
            error => { Error => "Invalid input"}
        });
        return;
    }

    # Authenticate user and connect to the database
    my ($db, $user) = CoGe::Services::Auth::init($self);

    # Get experiment
    my $experiment = $db->resultset("Experiment")->find($id);
    unless (defined $experiment) {
        $self->render(status => 404, json => {
            error => { Error => "Resource not found" }
        });
        return;
    }

    # Check permissions
    unless ( !$experiment->restricted || (defined $user && $user->has_access_to_experiment($experiment)) ) {
        $self->render(json => {
            error => { Auth => "Access denied" }
        }, status => 401);
        return;
    }

    # Format metadata
    my @metadata = map {
        {
            text => $_->annotation,
            link => $_->link,
            type => $_->type->name,
            type_group => $_->type->group
        }
    } $experiment->annotations;

    # Format types
    my @types = map {
        {
            name => $_->name,
            description => $_->description
        }
    } $experiment->types;

    $self->render(json => {
        id => int($experiment->id),
        name => $experiment->name,
        description => $experiment->description,
        version => $experiment->version,
        genome_id  => int($experiment->genome->id),
        source => {
            name => $experiment->source->name,
            description => $experiment->source->description,
            link => $experiment->source->link
        },
        types => \@types,
        additional_metadata => \@metadata,
        restricted => $experiment->restricted ? Mojo::JSON->true : Mojo::JSON->false,
    });
}

sub add {
    my $self = shift;
    my $data = $self->req->json;

# mdb removed 9/17/15 -- auth is handled by Job::add below, redundant token validation breaks CAS proxyValidate
#    # Authenticate user and connect to the database
#    my ($db, $user, $conf) = CoGe::Services::Auth::init($self);
#
#    # User authentication is required to add experiment
#    unless (defined $user) {
#        $self->render(json => {
#            error => { Auth => "Access denied" }
#        });
#        return;
#    }

    # Valid data items
    unless ($data->{source_data} && @{$data->{source_data}}) {
        $self->render(status => 400, json => {
            error => { Error => "No data items specified" }
        });
        return;
    }
    
    # Marshall incoming payload into format expected by Job Submit.
    # Note: This is kind of a kludge -- is there a better way to do this using
    # Mojolicious routing?
    my $request = {
        type => 'load_experiment',
        parameters => $data
    };
    
    return CoGe::Services::Data::Job::add($self, $request);
}

sub update {
	my $self = shift;
    my $id = int($self->stash('id'));
    
    # Validate input
    unless ($id) {
        $self->render(status => 400, json => {
            error => { Error => "Invalid input"}
        });
        return;
    }

    # Authenticate user and connect to the database
    my ($db, $user) = CoGe::Services::Auth::init($self);

    # Get experiment
    my $experiment = $db->resultset("Experiment")->find($id);
    unless (defined $experiment) {
        $self->render(status => 404, json => {
            error => { Error => "Resource not found" }
        });
        return;
    }

    # Check permissions
    unless ($user->is_owner_editor(experiment => $id)) {
        $self->render(json => {
            error => { Auth => "Access denied" }
        }, status => 401);
        return;
    }

    my $data = $self->req->json;
    if (exists($data->{metadata}->{id})) {
	    delete $data->{metadata}->{id};
    }
	$experiment->update($data->{metadata});
	$self->render(json => {
		success => Mojo::JSON->true
	});
}

sub debug {
	my $data = shift;
	my $new_file = shift;
	my $OUTFILE;
	open $OUTFILE, ($new_file ? ">/tmp/sean" : ">>/tmp/sean");
	print {$OUTFILE} Dumper $data;
	print {$OUTFILE} "\n";
	close $OUTFILE;
}

sub data {
	my $self = shift;
    my $id = int($self->stash('id'));
    my $chr = $self->param('chr');
    my $data_type = $self->param('data_type');
    my $type = $self->param('type');
    my $gte = $self->param('gte');
    my $lte = $self->param('lte');
    my $transform = $self->param('transform');

    # Authenticate user and connect to the database
    my ($db, $user) = CoGe::Services::Auth::init($self);
    unless ($user) {
        $self->render(json => {
            error => { Error => "User not logged in" }
        });
        return;
    }    	

    # Get experiment
    my $experiment = $db->resultset("Experiment")->find($id);
    unless (defined $experiment) {
        $self->render(json => {
            error => { Error => "Experiment not found" }
        });
        return;
    }

    # Check permissions
    unless ($user->is_admin() || $user->is_owner_editor(experiment => $id)) {
        $self->render(json => {
            error => { Auth => "Access denied" }
        }, status => 401);
        return;
    }

	$self->res->headers->content_disposition('attachment; filename=experiment.csv;');
	$self->write('# experiment: ' . $experiment->name . "\n");
	$self->write("# chromosome: $chr\n") if $chr;
	if ($type) {
		$self->write("# search: type = $type");
		$self->write(", gte = $gte") if $gte;
		$self->write(", lte = $lte") if $lte;
		$self->write("\n");
	}
	$self->write("# transform: $transform\n") if $transform;
	my $cols = CoGe::Core::Experiment::get_fastbit_format()->{columns};
	my @columns = map { $_->{name} } @{$cols};
	$self->write('# columns: ');
	for (my $i=0; $i<scalar @columns; $i++) {
		$self->write(',') if $i;
		$self->write($columns[$i]);
	}
	$self->write("\n");

	my $lines = CoGe::Core::Experiment::query_data(
		eid => $id,
		data_type => $data_type,
		col => join(',', @columns),
		chr => $chr,
		type => $type,
		gte => $gte,
		lte => $lte,
	);
	my $score_column = CoGe::Core::Experiment::get_fastbit_score_column($data_type);
	my $log10 = log(10);
	foreach my $line (@{$lines}) {
		if ($transform) {
			my @tokens = split ',', $line;
			if ($transform eq 'Inflate') {
				$tokens[$score_column] = 1;
			}
			elsif ($transform eq 'Log2') {
				$tokens[$score_column] = log(1 + $tokens[$score_column]);
			}
			elsif ($transform eq 'Log10') {
				$tokens[$score_column] = log(1 + $tokens[$score_column]) / $log10;
			}
			$line = join ',', @tokens;
		}
		$self->write($line);
		$self->write("\n");
	}
	$self->finish();
}

1;
