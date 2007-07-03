package CoGe::Accessory::anchors;

use strict;
use warnings;
#use File::Temp;
use CoGe::Accessory::chaos_report;
use CoGe::Accessory::dialign_report;
use base qw(Class::Accessor);

#use Data::Dumper;

BEGIN
	{
		use vars qw($VERSION);
		$VERSION = "0.01";
	}
__PACKAGE__->mk_accessors qw(file1 file2 run_chaos run_dialign base_name output_dir chaos_file anchor_file fasta_file dialign_file run_chaos_opts run_dialign_opts chaos_report_opts dialign_report_opts chaos_report dialign_report);

###############################################################################
# anchors -- Josh Kane  UC Berkeley
###############################################################################
sub new {
	my $proto = shift;
	my $opts = shift;
	$opts = {} unless $opts;
	my $class = ref($proto) || $proto;
	my $self = bless ({%$opts}, $class);
	$self->run_chaos("/opt/apache/CoGe/bin/lagan/chaos") unless $self->run_chaos;
	$self->run_dialign("/home/jkane/src/dialign_package/src/dialign2-2") unless $self->run_dialign;
	$self->output_dir("/opt/apache/CoGe/tmp") unless $self->output_dir;
	$self->run_chaos_opts("-v") unless $self->run_chaos_opts;
	$self->run_dialign_opts("-anc -fn") unless $self->run_dialign_opts;
	$self->run_program();
	return $self;
}

sub run_program
{
	my $self = shift;
	$self->generate_anchors;
	$self->run_dialign_with_anchors;
	return $self;
}
	

sub generate_anchors
{
	my $self = shift;
	my $file1 = shift || $self->file1;
	my $file2 = shift || $self->file2;
	my $output_dir = $self->output_dir;
	my $run_chaos = $self->run_chaos;
	my $base_name = $self->base_name;
	my $chaos_opts = $self->run_chaos_opts;
	my $parser_opts = $self->chaos_report_opts;
	$parser_opts = {} unless ref($parser_opts) =~ /hash/i;
	my ($query_start,$query_stop,$subject_start,$length,$score) = (0,0,0,0,0);
	my $anchor_file = "";
	
	print "file1 is $file1, file2 is $file2\n";
	
	`$run_chaos $file1 $file2 $chaos_opts > $output_dir/$base_name.chaos`;
	
	open(IN,"< $output_dir/$base_name.chaos");
	my $data = join ("", <IN>);
	close IN;
	print "data is $data\n";
	$self->chaos_file($data);
	
	`cat $file1 > $output_dir/$base_name.fasta`;
	`cat $file2 >> $output_dir/$base_name.fasta`;
	
	open(IN,"< $output_dir/$base_name.fasta");
	$data = join ("", <IN>);
	close IN;
	$self->fasta_file($data);
	
	my $chaos_report = new CoGe::Accessory::chaos_report({file=>"$output_dir/$base_name.chaos", %$parser_opts});
	$self->chaos_report($chaos_report);
	
	foreach my $hsp (@{$chaos_report->hsps})
	{
	  $query_start = $hsp->query_start;
	  $query_stop = $hsp->query_stop;
	  $subject_start = $hsp->subject_start;
	  $score = $hsp->score;
	  $length = ($query_stop - $query_start) + 1;
	  $anchor_file .= "1 2 $query_start $subject_start $length $score\n";
	}
	$self->anchor_file($anchor_file);
	return $self;
}
	
sub run_dialign_with_anchors
{
	my $self = shift;
	my $run_dialign = $self->run_dialign;
	my $dialign_opts = $self->run_dialign_opts;
	my $parser_opts = $self->dialign_report_opts;
	$parser_opts = {} unless ref($parser_opts) =~ /hash/i;
	my $base_name = $self->base_name;
	my $output_dir = $self->output_dir;
	
	`$run_dialign $dialign_opts $output_dir/$base_name.dialign $output_dir/$base_name.fasta`;
	#`mv $basename.ali $output_dir/$basename.ali`;
	#some move command to output_dir directory
	
	open(IN,"< $output_dir/$base_name.dialign");
	my $data = join ("", <IN>);
	close IN;
	$self->dialign_file($data);
	
	
	my $dialign_report = new CoGe::Accessory::dialign_report({file=>"$output_dir/$base_name.ali", %$parser_opts});
	
	$self->dialign_report($dialign_report);
	return $self;
}
	
1;

__END__

	