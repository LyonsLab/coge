package CoGe::Services::JBrowse::Annotation;
use base 'CGI::Application';

use CoGeX;
use CoGe::Accessory::Web;
use JSON qq{encode_json};
use Data::Dumper;
use Cwd 'abs_path';

my $coge_conf;

sub setup {
	my $self = shift;
	
	#FIXME - move this into service.pl
	$coge_conf = abs_path($0);
	$coge_conf =~ s/services\/service\.pl/coge\.conf/;	
	
	$self->run_modes(
		'stats_global' 	=> 'stats_global',
		'features' 		=> 'features',
	);
	$self->start_mode('stats_global');
	$self->mode_param('rm');
}

sub stats_global {
	print STDERR "Annotation::stats_global\n";
	return qq{{}};
}

sub features {
	my $self = shift;
	my $dsid = $self->param('ds');
	my $chr = $self->param('chr');
	my $start = $self->query->param('start');
	my $end = $self->query->param('end');
	print STDERR "Annotation::features $chr:$start:$end\n";

	# Load config file
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
				{ 
				  start		=> $f->start,
				  end		=> $f->stop,
				  strand 	=> $f->strand,
				  uniqueID	=> $f->id,
				  type		=> 'CDS',
				};
		}
		
		push @feats, { 
			chr 		=> $chr,
			start 		=> $gene->start,
			end 		=> $gene->stop,
			strand 		=> '0'.$gene->strand,
			name 		=> $name,
			type 		=> 'mRNA',
			uniqueID 	=> $gene->id,
			subfeatures => \@subfeats
		};
		print STDERR "transcript: " . $gene->id . " " . $name . " start=" .$gene->start . " end=" . $gene->stop . "\n";
	}

#	print STDERR encode_json({"features" => \@feats});
	return encode_json({"features" => \@feats});
}

sub isOverlapping {
	my ($s1, $e1, $s2, $e2) = @_;
	return ($s1 <= $e2 && $s2 <= $e1);
}

1;
