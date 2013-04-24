package CoGe::Services::JBrowse::Experiment;
use base 'CGI::Application';

use Cwd 'abs_path';

my $NUM_COL = 6;
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
	print STDERR "experiment stats_global\n";
	return qq{{
		"scoreMin" : -1,
		"scoreMax" : 1
	}};
}

sub features {
	my $self = shift;
	my $exp = $self->param('exp');
	my $chr = $self->param('chr');
	my $start = $self->query->param('start');
	my $end = $self->query->param('end');
	print STDERR "experiment features $chr:$start:$end\n";
	
	my $exp_storage_path = '/storage/coge/data/experiments/0/0/0/'.$exp;
	
	# Call FastBit to do query
	# Issue 61: query string must contain a "." for fastbit to use consistent
	# format (see jira thread for explanation) so "0.0=0.0" was added along with
	# -v option.  Output parsing was modified accordingly for new output format.
	#my $cmd = "$CMDPATH -d $exp_storage_path -q \"select chr,start,stop,strand,value1,value2 where chr='$chr' and start > $start and stop < $stop\" 2>&1"; # mdb removed 3/27/13 issue 61
	my $cmd = "/usr/local/bin/ibis -v 1 -d $exp_storage_path -q \"select chr,start,stop,strand,value1,value2 where 0.0=0.0 and chr='$chr' and start > $start and stop < $end order by start limit 999999999\" 2>&1"; # mdb added 3/27/13 issue 61
#	print STDERR "$cmd\n";
	my @cmdOut = qx{$cmd};
#	print STDERR @cmdOut;
	my $cmdStatus = $?;
	die "Error executing command $CMDPATH ($cmdStatus)" if ($cmdStatus != 0);
	
	# Convert FastBit output into JSON
	my $results = '';
	foreach (@cmdOut) { # mdb rewritten 3/27/13 issue 61
		chomp;
		if (/^\"/) { #if (/^\"$chr\"/) { # potential result line
			s/"//g;
			my @items = split(/,\s*/);
			next if (@items != $NUM_COL);# || $items[0] !~ /^\"?$chr/); # make sure it's a row output line
			for (my $i =0; $i<@items; $i++) {
				$items[$i] = 1 if $items[$i] !~ /\w/; # what's this for?
			}
			#$results .= '[' . join(',', map {"\"$_\""} @items) . ']';
			
			my ($chr, $start, $end, $strand, $value1, $value2) = @items;
			$end = $start+1 if ($end == $start); #FIXME revisit this
			$value1 = $strand*$value1;
			$results .= ($results ? ',' : '') . qq{{ "start": $start, "end": $end, "score": $value1 }};
		}
	}
	
	#print STDERR "{ 'features' : [ $results ] }\n";
	return qq{{ "features" : [ $results ] }};
}

1;