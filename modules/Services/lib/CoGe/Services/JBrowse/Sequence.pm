package CoGe::Services::JBrowse::Sequence;
use base 'CGI::Application';

use CoGeX;
use CoGe::Accessory::Web;
use JSON qq{encode_json};
use Data::Dumper;
use Cwd 'abs_path';
use File::Basename 'fileparse';

my $coge_conf;

sub setup {
	my $self = shift;
	
	#FIXME - move this into service.pl
	$coge_conf = abs_path($0);
	$coge_conf =~ s/services\/JBrowse\/service\.pl/coge\.conf/;	
	
	$self->run_modes(
		'stats_global' 	=> 'stats_global',
		'features'		=> 'features',
	);
	#$self->start_mode('stats_global');
	$self->mode_param('rm');
}

sub stats_global {
	print STDERR "Sequence::stats_global\n";
	return qq{{}};	
}

sub features {
	my $self = shift;
	my $gid = $self->param('gid');
	my $chr = $self->param('chr');
	my $size = $self->query->param('seqChunkSize');
	my $start = $self->query->param('start');
	my $end = $self->query->param('end');
	my $len = $end-$start;
	print STDERR "Sequence::features gid=$gid chr=$chr size=$size start=$start end=$end\n";
	
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
	
	my $genome = $coge->resultset('Genome')->find($gid);
	return unless $genome;

	# Extract requested piece of sequence file
	my (undef, $storagepath) = fileparse($genome->file_path);
	my $seqfile = $storagepath.'/chr/'.$chr;
	open(my $fh, $seqfile) or die;
	seek($fh, $start, 0);
	read($fh, my $seq, $len);
	close($fh);
	#print STDERR "$seqfile $len $seq\n";
	
	return qq{
		{ "features" : [
			{ "start": $start, "end": $end, "seq": "$seq" }
			]
		}
	};
}

1;
