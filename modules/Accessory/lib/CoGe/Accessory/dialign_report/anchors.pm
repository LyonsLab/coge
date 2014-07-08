package CoGe::Accessory::dialign_report::anchors;

use strict;
use warnings;
#use File::Temp;
use CoGe::Accessory::bl2seq_report;
use CoGe::Accessory::chaos_report;
use CoGe::Accessory::blastz_report;
use CoGe::Accessory::dialign_report;
use base qw(Class::Accessor);

#use Data::Dumper;

BEGIN
	{
		use vars qw($VERSION $DEBUG);
		$VERSION = "0.01";
	}
__PACKAGE__->mk_accessors(qw(file1 file2 run_anchor run_dialign base_name extension output_dir anchor_output anchor_file fasta_file dialign_file run_anchor_opts run_dialign_opts anchor_report_opts dialign_report_opts anchor_report dialign_report log_file DEBUG));

###############################################################################
# anchors -- Josh Kane  UC Berkeley
###############################################################################
sub new {
	my $proto = shift;
	my $opts = shift;
	$opts = {} unless $opts;
	my $class = ref($proto) || $proto;
	my $self = bless ({%$opts}, $class);
	$self->run_anchor("/opt/apache/coge/web/bin/lagan/chaos_coge") unless $self->run_anchor;
	$self->run_dialign("/opt/apache/coge/web/bin/dialign2_dir/dialign2-2_coge") unless $self->run_dialign;
	$self->output_dir("/opt/apache/coge/web/tmp") unless $self->output_dir;
	$self->run_anchor_opts("-v") unless $self->run_anchor_opts;
	$self->run_dialign_opts("-n") unless $self->run_dialign_opts;
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

	print STDERR  "file1 is $file1, file2 is $file2\n" if $self->DEBUG;
	my $command;
	if ($extension=~/bl2/i)
	  {$command = "$run_anchor -p blastn -i $file1 -j $file2 $anchor_opts > $output_dir/$base_name.$extension";}
	else{
	$command = "$run_anchor $file1 $file2 $anchor_opts > $output_dir/$base_name.$extension";}

	print STDERR  "call to $extension: $command\n" if $self->DEBUG;
	$self->write_log("running $command");
	`$command`;

	$self->anchor_output("$output_dir/$base_name.$extension");
	print STDERR "Anchor output: ",$self->anchor_output,"\n" if $self->DEBUG;

	`cat $file1 > $output_dir/$base_name.fasta`;
	`cat $file2 >> $output_dir/$base_name.fasta`;

	$self->fasta_file("$output_dir/$base_name.fasta");

	$anchor_report = new CoGe::Accessory::chaos_report({file=>"$output_dir/$base_name.$extension", %$parser_opts}) if $extension=~/chaos/i;

	$anchor_report = new CoGe::Accessory::bl2seq_report({file=>"$output_dir/$base_name.$extension", %$parser_opts}) if $extension=~/bl2/i;

	$anchor_report = new CoGe::Accessory::blastz_report({file=>"$output_dir/$base_name.$extension", %$parser_opts}) if $extension=~/blastz/i;

	$self->anchor_report($anchor_report);

	foreach my $hsp (@{$anchor_report->hsps})
	{
	  next if $hsp->strand =~ /-/;
	  $query_start = $hsp->query_start;
	  $query_stop = $hsp->query_stop;
	  $subject_start = $hsp->subject_start;
	  $score = $hsp->score;
	  $length = ($query_stop - $query_start) + 1;
	  $anchor_file .= "1 2 $query_start $subject_start $length $score\n";
	}

	print STDERR "Anchor file is: \n",$anchor_file,"\n" if $self->DEBUG;
	$self->write_log("creating dialign anchor file");
	open(NEW,"> $output_dir/$base_name.anc");
	print  NEW $anchor_file;
	close NEW;

	$self->anchor_file("$output_dir/$base_name.anc");
	print STDERR "Anchor file: ",$self->anchor_file,"\n" if $self->DEBUG;
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
	my $command = "$run_dialign $dialign_opts -anc -fn $output_dir/$base_name.dialign $output_dir/$base_name.fasta";
	print STDERR  "call to dialign: $command\n" if $self->DEBUG;
	$self->write_log("running $command");
	`$command`;
	#`mv $basename.ali $output_dir/$basename.ali`;
	#some move command to output_dir directory

	$self->dialign_file("$output_dir/$base_name.dialign");

	my $dialign_report = new CoGe::Accessory::dialign_report({file=>"$output_dir/$base_name.dialign", %$parser_opts});

	$self->dialign_report($dialign_report);
	return $self;
}

sub write_log
  {
    my $self = shift;
    my $message = shift;
    return unless $self->log_file;
    open (OUT, ">>".$self->log_file) || return;
    print OUT $message,"\n";
    close OUT;
  }

1;

__END__
