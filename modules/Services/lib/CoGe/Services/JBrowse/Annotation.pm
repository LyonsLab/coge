package CoGe::Services::JBrowse::Annotation;
use base 'CGI::Application';

use CoGeX;
use CoGe::Accessory::Web;
use JSON qq{encode_json};
use Data::Dumper;

sub setup {
	my $self = shift;
	$self->run_modes(
		'stats_global' 	=> 'stats_global',
		'features' 		=> 'features',
	);
	$self->start_mode('stats_global');
	$self->mode_param('rm');
}

sub stats_global {
	print STDERR "Annotation::stats_global\n";
	return qq{{
	}};
}

sub features {
	my $self = shift;
	my $dsid = $self->param('ds');
	my $chr = $self->param('chr');
	my $start = $self->query->param('start');
	my $end = $self->query->param('end');
	print STDERR "Annotation::features $chr:$start:$end\n";

	# Load config file
#	unless ($coge_conf) {
#		$coge_conf = abs_path($0);
#		$coge_conf =~ s/services\/features.pl//;
#		$coge_conf .= 'coge.conf';
#	}
	my $coge_conf = '/home/mbomhoff/public/CoGe/coge.conf';
	my $P = CoGe::Accessory::Web::get_defaults($coge_conf);
	my $DBNAME = $P->{DBNAME};
	my $DBHOST = $P->{DBHOST};
	my $DBPORT = $P->{DBPORT};
	my $DBUSER = $P->{DBUSER};
	my $DBPASS = $P->{DBPASS};

	# Connect to the database
	my $connstr = "dbi:mysql:dbname=$DBNAME;host=$DBHOST;port=$DBPORT";
	my $coge = CoGeX->connect($connstr, $DBUSER, $DBPASS);
	#$coge->storage->debugobj(new DBIxProfiler());
	#$coge->storage->debug(1);

	# Get features
	my @raw_feats = $coge->resultset('Feature')->search(
		{	"me.chromosome" => $chr,
			"me.dataset_id" => $dsid,
			-and => [
				"me.start" => {"<=" => $end},
				"me.stop"  => {">=" => $start},
			],
		},
		{	prefetch => [ "locations", 
						  "feature_names", 
						  "feature_type" 
						],
#			order_by=>"me.start",
		}
	);

	# Output JSON
	my @feats;
	my %loc;
	foreach my $f (@raw_feats) {
		my ($n) = sort {$a->name cmp $b->name} $f->feature_names({primary_name=>1});
		my $name = ($n ? $n->name : '');
		my $type = $f->feature_type->name;
		#print STDERR "feat: " . $f->id . " " . $type . " " . $name . " start=" . $f->start . " end=" . $f->stop . "\n";
		push @{$loc{$name}{$type}}, $f->locations;
	}
	
	foreach my $name (sort keys %loc) {
		my ($rootName) = $name =~ /^(\w+)\.?/;
		next unless $rootName;
		next if ($rootName eq $name);
		my $gene = $loc{$rootName}{'gene'}->[0];
		next unless $gene;
		
		my @subfeats;
		foreach my $f (@{$loc{$name}{'CDS'}}) {
			push @subfeats,
				{ start  => $f->start,
				  end 	 => $f->stop,
				  strand => $f->strand,
				  id 	 => $f->id,
				  type 	 => 'CDS',
				};
		}
		
		push @feats, { 
			start 	=> $gene->start,
			end 	=> $gene->stop,
			strand 	=> $gene->strand,
			name 	=> $name,
			type 	=> 'mRNA',
			id 		=> $gene->id,
			subfeatures => \@subfeats
		};
		print STDERR "transcript: " . $gene->id . " " . $name . " start=" .$gene->start . " end=" . $gene->stop . "\n";
	}

#	print STDERR encode_json({"features" => \@feats});
	return encode_json({"features" => \@feats});
}

sub features_old {
	my $self = shift;
	my $dsid = $self->param('ds');
	my $chr = $self->param('chr');
	my $start = $self->query->param('start');
	my $end = $self->query->param('end');
	print STDERR "Annotation::features $chr:$start:$end\n";

	# Load config file
#	unless ($coge_conf) {
#		$coge_conf = abs_path($0);
#		$coge_conf =~ s/services\/features.pl//;
#		$coge_conf .= 'coge.conf';
#	}
	my $coge_conf = '/home/mbomhoff/public/CoGe/coge.conf';
	my $P = CoGe::Accessory::Web::get_defaults($coge_conf);
	my $DBNAME = $P->{DBNAME};
	my $DBHOST = $P->{DBHOST};
	my $DBPORT = $P->{DBPORT};
	my $DBUSER = $P->{DBUSER};
	my $DBPASS = $P->{DBPASS};

	# Connect to the database
	my $connstr = "dbi:mysql:dbname=$DBNAME;host=$DBHOST;port=$DBPORT";
	my $coge = CoGeX->connect($connstr, $DBUSER, $DBPASS);
	#$coge->storage->debugobj(new DBIxProfiler());
	#$coge->storage->debug(1);

	# Get features
#	my @raw_feats = $coge->resultset('Feature')->search({
#		dataset_id => $dsid,
#		chromosome => $chr,
#		{ or => [
#			{start => {between => [$start,$end]}},
#			{stop => {between => [$start,$end]}}
#		]}
#	});
	my @raw_feats = $coge->get_features_in_region(
		dataset => $dsid,
		chr => $chr, 
		start => $start, 
		end => $end, 
		#ftid => 1#ftid => $ftid
	);

	my %type_to_name = map { $_->id => $_->name } $coge->resultset('FeatureType')->all;

#	foreach (@raw_feats) {
#		print STDERR "feat: " . $type_to_name{$_->feature_type_id} . " start=" . $_->start . " end=" . $_->stop . "\n";
#	}
	
	# Output JSON
	my @feats;
	my %featsByType;
	map { push @{$featsByType{$type_to_name{$_->feature_type_id}}}, $_ } @raw_feats;
	foreach my $gene (@{$featsByType{'mRNA'}}) {
		print STDERR "mRNA " . $gene->start . ":" . $gene->stop . "\n";
		
		my @subfeats;
		foreach my $type ('CDS') {
			foreach (@{$featsByType{$type}}) {
				next if (!isOverlapping($_->start, $_->stop, $start, $end));
				push @subfeats, 
					{ start => $_->start,
					  end => $_->stop,
					  strand => $_->strand,
					  type => $type,
					};
				print STDERR "   CDS " . $_->start . ":" . $_->stop . "\n";
			}
		}
		
		my $name = $gene->id;#$_->primary_name->name;
		push @feats, { 
			start 	=> $gene->start,
			end 	=> $gene->stop,
			strand 	=> $gene->strand,
			name 	=> $name,
			#label 	=> $name,
			type 	=> $type_to_name{$gene->feature_type_id},
			subfeatures => \@subfeats
		};
	}
	
	print STDERR encode_json({"features" => \@feats});
	return encode_json({"features" => \@feats});
}

sub isOverlapping {
	my ($s1, $e1, $s2, $e2) = @_;
	return ($s1 <= $e2 && $s2 <= $e1);
}

1;