package CoGe::Accessory::anchors;

use strict;
use warnings;
#use File::Temp;
use CoGe::Accessory::bl2seq_report;
use CoGe::Accessory::chaos_report;
use CoGe::Accessory::dialign_report;
use base qw(Class::Accessor);

#use Data::Dumper;

BEGIN
	{
		use vars qw($VERSION $DEBUG);
		$VERSION = "0.01";
	}
__PACKAGE__->mk_accessors qw(file1 file2 run_anchor run_dialign base_name extension output_dir anchor_file anchor_file fasta_file dialign_file run_anchor_opts run_dialign_opts anchor_report_opts dialign_report_opts anchor_report dialign_report);

###############################################################################
# anchors -- Josh Kane  UC Berkeley
###############################################################################
sub new {
	my $proto = shift;
	my $opts = shift;
	$opts = {} unless $opts;
	my $class = ref($proto) || $proto;
	my $self = bless ({%$opts}, $class);
	$self->run_anchor("/opt/apache/CoGe/bin/lagan/chaos_coge") unless $self->run_anchor;
	$self->run_dialign("/opt/apache/CoGe/bin/dialign2_dir/dialign2-2_coge") unless $self->run_dialign;
	$self->output_dir("/opt/apache/CoGe/tmp") unless $self->output_dir;
	$self->run_anchor_opts("-v") unless $self->run_anchor_opts;
	$self->run_dialign_opts("-n -anc -fn") unless $self->run_dialign_opts;
	$self->extension("chaos") unless $self->extension;
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
	my $run_anchor = $self->run_anchor;
	my $base_name = $self->base_name;
	my $anchor_opts = $self->run_anchor_opts;
	my $parser_opts = $self->anchor_report_opts;
	my $extension = $self->extension;
	$parser_opts = {} unless ref($parser_opts) =~ /hash/i;
	my ($query_start,$query_stop,$subject_start,$length,$score) = (0,0,0,0,0);
	my $anchor_file = "";
	my $anchor_report = [];
	
	print "file1 is $file1, file2 is $file2\n" if $DEBUG;
	
	#print "call to $extension: $run_anchor $file1 $file2 $chaos_opts > $output_dir/$base_name.chaos\n" if $DEBUG;
	`$run_anchor $file1 $file2 $anchor_opts > $output_dir/$base_name.$extension` if $extension=~/chaos/i;
	`$run_anchor -i $file1 -j $file2 $anchor_opts > $output_dir/$base_name.$extension` if $extension=~/bl2/i;
	
	open(IN,"< $output_dir/$base_name.$extension") || die "can't open $base_name.$extension for reading: $!";
	my $data = join ("", <IN>);
	close IN;
	
	$self->anchor_file($data);
	
	`cat $file1 > $output_dir/$base_name.fasta`;
	`cat $file2 >> $output_dir/$base_name.fasta`;
	
	open(IN,"< $output_dir/$base_name.fasta") || die "can't open $base_name.fasta for reading: $!";
	$data = join ("", <IN>);
	close IN;
	$self->fasta_file($data);
	
	$anchor_report = new CoGe::Accessory::chaos_report({file=>"$output_dir/$base_name.chaos", %$parser_opts}) if $extension=~/chaos/i;
	
	$anchor_report = new CoGe::Accessory::bl2seq_report({file=>"$output_dir/$base_name.chaos", %$parser_opts}) if $extension=~/bl2/i;
	
	$self->anchor_report($anchor_report);
	
	foreach my $hsp (@{$anchor_report->hsps})
	{
	  $query_start = $hsp->query_start;
	  $query_stop = $hsp->query_stop;
	  $subject_start = $hsp->subject_start;
	  $score = $hsp->score;
	  $length = ($query_stop - $query_start) + 1;
	  $anchor_file .= "1 2 $query_start $subject_start $length $score\n";
	}
	$self->anchor_file($anchor_file);
	
	open(NEW,"> $output_dir/$base_name.anc");
	print NEW $anchor_file;
	close NEW;
	
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
	
	print "call to dialign: $run_dialign $dialign_opts $output_dir/$base_name.dialign $output_dir/$base_name.fasta\n" if $DEBUG;
	`$run_dialign $dialign_opts $output_dir/$base_name.dialign $output_dir/$base_name.fasta`;
	#`mv $basename.ali $output_dir/$basename.ali`;
	#some move command to output_dir directory
	
	open(IN,"< $output_dir/$base_name.dialign") || die "can't open $base_name.dialign for reading: $!";
	my $data = join ("", <IN>);
	close IN;
	$self->dialign_file($data);
	
	
	my $dialign_report = new CoGe::Accessory::dialign_report({file=>"$output_dir/$base_name.dialign", %$parser_opts});
	
	$self->dialign_report($dialign_report);
	return $self;
}
	
1;

__END__

	