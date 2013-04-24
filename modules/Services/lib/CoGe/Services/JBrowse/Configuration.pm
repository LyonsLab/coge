package CoGe::Services::JBrowse::Configuration;
use base 'CGI::Application';

use CoGeX;
use CoGe::Accessory::Web;
use JSON;
use Data::Dumper;
use Sort::Versions;

sub setup {
	my $self = shift;
	$self->run_modes(
		'refseq_config'	=> 'refseq_config',
		'track_config' 	=> 'track_config',
	);
	$self->start_mode('track_config');
	$self->mode_param('rm');
}

sub refseq_config {
	my $self = shift;
	my $gid = $self->query->param('gid');
	my $SEQ_CHUNK_SIZE = 20000;
	print STDERR "Configuration::refseq_config gid=$gid\n";
	#return qq{[{"length":4686137,"name":"1","seqChunkSize":20000,"end":4686137,"start":0}]};
	
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
	
	my $genome = $coge->resultset('Genome')->find($gid);
	return unless $genome;	

	my @chromosomes;
	foreach my $chr ($genome->genomic_sequences) {
		push @chromosomes, { 
			name => $chr->chromosome,
			length => $chr->sequence_length,
			seqChunkSize => $SEQ_CHUNK_SIZE,
			start => 0,
			end => $chr->sequence_length - 1
		};
	}
	
	return encode_json(\@chromosomes);
}

sub track_config {
	my $self = shift;
	my $gid = $self->query->param('gid');
	print STDERR "Configuration::track_config gid=$gid\n";
	
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
	
	my $genome = $coge->resultset('Genome')->find($gid);
	return unless $genome;
	
	my @tracks;
	
	# Add reference sequence track
	push @tracks,
	{
		chunkSize => 20000,
		baseUrl => "http://geco.iplantcollaborative.org/mbomhoff/CoGe/services/service.pl/sequence/$gid/",
		type => "SequenceTrack",
		storeClass => "JBrowse/Store/SeqFeature/REST",
		label => "sequence",
		key => "Sequence",
		formatVersion => 1
	};

	# Add gene annotation tracks
	my $num;
	foreach $ds ($genome->datasets) {
		next if (not $ds->has_gene_annotation);
		my $dsid = $ds->id;
		push @tracks, 
		{
			baseUrl => "http://geco.iplantcollaborative.org/mbomhoff/CoGe/services/service.pl/annotation/$dsid/",
            autocomplete => "all",
            track => "genes",
            label => "genes",
            key => "Genes" . ($num++ ? " $num" : ''),
         	type => "FeatureTrack",
            storeClass => "JBrowse/Store/SeqFeature/REST",
         	style => {
            	className => "cds",
         	},
        }
	}
	
	# Add experiment tracks
	foreach $e (sort experimentcmp $genome->experiments) {
		next if ($e->deleted);
		my $eid = $e->id;
		push @tracks, 
		{ 
			baseUrl => "http://geco.iplantcollaborative.org/mbomhoff/CoGe/services/service.pl/experiment/$eid/",
		    autocomplete => "all",
		    track => "exp$eid",
		    label => "exp$eid",
		    key => $e->name,
		    type => "JBrowse/View/Track/Wiggle/XYPlot",
		    storeClass => "JBrowse/Store/SeqFeature/REST"				
		}
	}
	
	return encode_json({ tracks => \@tracks, dataset_id => 'coge', formatVersion => 1 });
}

# FIXME this comparison routine is duplicated elsewhere
sub experimentcmp {
	no warnings 'uninitialized'; # disable warnings for undef values in sort
	versioncmp($b->version, $a->version) || $a->name cmp $b->name || $b->id cmp $a->id
}

1;
