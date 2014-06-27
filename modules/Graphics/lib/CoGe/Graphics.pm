package CoGe::Graphics;
use strict;
use base qw(Class::Accessor);
use CoGe::Graphics::Chromosome;
use CoGe::Graphics::Feature;
use CoGe::Graphics::Feature::Gene;
use CoGe::Graphics::Feature::NucTide;
use CoGe::Graphics::Feature::GAGA;
use CoGe::Graphics::Feature::GBox;
use CoGe::Graphics::Feature::Sigma54;
use CoGe::Graphics::Feature::Exon_motifs;
use CoGe::Graphics::Feature::AminoAcid;
use CoGe::Graphics::Feature::Domain;
use CoGe::Graphics::Feature::Block;
use CoGe::Graphics::Feature::Outline;
use CoGe::Graphics::Feature::Link;
use CoGe::Graphics::Feature::Quant;
use CoGeX;
use DBIxProfiler;
use Data::Dumper;
use Benchmark;
use JSON;
#use LWP::Simple;
use Carp;

use vars qw($VERSION $MAX_FEATURES $MAX_NT $DEBUG $BENCHMARK);
$VERSION = '0.1';
__PACKAGE__->mk_accessors( "MAX_FEATURES", "MAX_NT", "DEBUG", );
$BENCHMARK = 0;

#################### subroutine header begin ####################

=head2 sample_function

 Usage     : How to use this function/method
 Purpose   : What it does
 Returns   : What it returns
 Argument  : What it wants to know
 Throws    : Exceptions and other anomolies
 Comment   : This is a sample subroutine header.
           : It is polite to include more pod and fewer comments.

See Also   :

=cut

#################### subroutine header end ####################

sub new
{
	my ( $class, %parameters ) = @_;

	my $self = bless( {}, ref($class) || $class );

	return $self;
}

#################### main pod documentation begin ###################
## Below is the stub of documentation for your module.
## You better edit it!

=head1 NAME

CoGe::Graphics - CoGe::Graphics

=head1 SYNOPSIS

  use CoGe::Graphics;
  blah blah blah

=head1 DESCRIPTION

=head1 USAGE

=head1 BUGS

=head1 SUPPORT

=head1 AUTHOR

	Eric Lyons
	elyons@nature.berkeley.edu

=head1 COPYRIGHT

Permission to use, copy, modify, and distribute this software and its documentation for educational, research, and not-for-profit purposes, without fee and without a signed licensing agreement, is hereby granted, provided that the above copyright notice, this paragraph and the following two paragraphs appear in all copies, modifications, and distributions. Contact The Office of Technology Licensing, UC Berkeley, 2150 Shattuck Avenue, Suite 510, Berkeley, CA 94720-1620, (510) 643-7201, for commercial licensing opportunities.

The full text of the license can be found in the
LICENSE file included with this module.

=head1 SEE ALSO

perl(1).

=cut

#################### main pod documentation end ###################

#################### subroutine header begin ####################

=head2 genomic_view

 Usage     : CoGe::Graphics->genomic_view(
				   start=>$start,
				   stop=>$stop,
				   di=>$di,
				   version=>$version,
				   org=>$org_id,
				   chr=>$chr,
				   iw=>$iw,
				   ih=>$ih,
				   file=>$file,
				   simple=>$simple,
				   ch=>$chr_height,
				   fh=>$feat_height,
                                   gstid=>$gstid, #genoic sequence type id
				   fids=>\@fids,
				   fns=>\@fnames,
                                   layers=>\@layers,
				   );
 Purpose   : mama of a function.  Docs to come.
 Returns   :
 Argument  :
 Throws    :
 Comment   :
           :

See Also   :

=cut

#################### subroutine header end ####################

sub genomic_view
{
	my $self = shift;
	$self = new("CoGe::Graphics") unless ref($self) =~ /Graphics/;
	my %opts  = @_;
	my $start = $opts{'start'} || $opts{'x'} || 0;    #28520458;
	my $stop  = $opts{'stop'};                        # || 6948000;#6949600;#190000;
	$start = 1 unless defined $start;

	#$stop = $start unless $stop;
	my $ds      = $opts{'ds'}  || $opts{'di'};
	my $version = $opts{'v'}   || $opts{'version'};
	my $org     = $opts{'o'}   || $opts{'org'} || $opts{'organism'} || $opts{'org_id'} || $opts{'oid'};
	my $chr     = $opts{'chr'} || $opts{'chromosome'};
	my $iw      = $opts{'iw'}  || $opts{'width'} || $opts{'tile size'} || $opts{'tile_size'} || 256;
	my $ih          = $opts{'ih'};
	my $file        = $opts{'file'};                  # || "./tmp/pict.png";
	my $simple      = $opts{'simple'};
	my $chr_height  = $opts{'ch'} || 200;
	my $feat_height = $opts{'fh'} || 10;
	$BENCHMARK = $opts{'bm'} || 0;
	my $fids = $opts{'fid'} || $opts{'fids'};                       #used to highlight special features by their database id
	my $fnames = $opts{'fn'} || $opts{'fns'} || $opts{'fnames'};    #used to highlight special features by their name
	my $img_map_name      = $opts{'img_map'};
	my $layers            = process_layers( $opts{layers} );
	my $major_tick_labels = $opts{major_tick_labels} || 1;
	my $minor_tick_labels = $opts{minor_tick_labels} || -1;
	my $gstid             = $opts{gstid};                           #genomic_sequence_type_id for finding the right genomic sequence
	my $dsgid             = $opts{dsgid};                           #dataset_group for finding the right genomic sequence
	my $coge              = $opts{coge};                            #database connector object
	my $ftid              = $opts{ftid};
	my $expids = $opts{expid};
	my $color             = $opts{color};

	unless ($coge)
	{
		print STDERR "Need to pass in a coge database object!\n";
		return;
	}
	$DEBUG = $opts{debug} || $opts{DEBUG} || 0;
	$self->DEBUG($DEBUG);
	print STDERR "Options: " . Dumper \%opts if $self->DEBUG;

	$fids   = [$fids]   unless ref($fids)   =~ /array/i;
	$fnames = [$fnames] unless ref($fnames) =~ /array/i;
	##some limits to keep from blowing our stack
	$self->MAX_FEATURES( $iw * 10 );    #10 features per pixel
	$self->MAX_NT( $iw * 100000 );      #100,000 nts per pixel
	my $t0 = new Benchmark if $BENCHMARK;

	#CoGe objects that we will need
	my $c = new CoGe::Graphics::Chromosome;

	($org) = $coge->resultset('Organism')->resolve($org) if $org;
	($ds)  = $coge->resultset('Dataset')->resolve($ds)   if $ds;
	my $dsg = $coge->resultset('Genome')->find($dsgid) if $dsgid;

	#    print STDERR"!", $dsgid, "\n";#, $dsg->id,"\n";
	if ( $ds && !$dsg )
	{
		foreach my $item ( sort { $a->genomic_sequence_type_id <=> $b->genomic_sequence_type_id } $ds->dataset_groups )
		{
			last if $dsg;
			$dsg = $item if $gstid && $item->genomic_sequence_type_id == $gstid;
			$dsg = $item unless $gstid;
		}
	}
	##we need a dataset group to make sure we get all the appropriate features that may be spread over multiple datasets

	if ($ds)    #there can be additional information about a chromosome for a particular version of the organism that is not in the same data_information.  Let's go find the organism_id and version for the specified data information item.
	{
		($org) = $ds->organism unless $org;
		$version = $ds->version unless $version;
		($chr) = $ds->get_chromosomes() unless $chr;
	}

	if ( $org && !$ds )
	{
		($ds) = $org->current_datasets();
		$version = $ds->version;
		($chr) = $ds->get_chromosomes() unless $chr;
	}

	unless ( defined $start && $chr && ($ds || $dsg) )
	{
		my $dsid = $ds ? $ds->id : "none";
		carp "missing needed parameters: Start: $start, Dataset_id: " . $dsid . ", Genome_id: ".$dsgid.", Chr: $chr\n";
		return 0;
	}
	if (!$org && $dsg)
	{
	  $org = $dsg->organism;
	  $version = $dsg->version;
	}

	print STDERR "generating image for ds: " . $ds->name . " (" . $ds->id . ")\n" if $self->DEBUG;
	my $ta = new Benchmark if $BENCHMARK;
	my $tb = new Benchmark if $BENCHMARK;
	my $finddid_time = timestr( timediff( $tb, $ta ) ) if $BENCHMARK;
	my $last_position;
	if ($ds) {$last_position  = $ds->last_chromosome_position($chr);}
	if ($dsg) {$last_position = $dsg->last_chromosome_position($chr);}
	$self->initialize_c(
		 chr                => $chr,
		 iw                 => $iw,
		 ih                 => $ih,
		 start              => $start,
		 stop               => $stop,
		 c                  => $c,
		 ch                 => $chr_height,
		 fh                 => $feat_height,
		 feature_labels     => 1,
		 major_tick_labels  => $major_tick_labels,
		 minor_tick_labels  => $minor_tick_labels,
		 max_track          => 2,
		 overlap_adjustment => $layers->{features}{overlap_check},
	);

	if ( $last_position < $stop )
	{
		$c->chr_end($last_position);
		#	$c->stop($last_position);
		#	print STDERR join ("\t", $c->start, $c->stop, $c->chr_end, $c->chr_length),"\n" ;
	}
	my $tc = new Benchmark if $BENCHMARK;
	print STDERR "processing nucleotides\n" if $self->DEBUG;
	$self->process_nucleotides( start => $start, stop => $stop, chr => $chr, ds => $ds, dsg=>$dsg, c => $c, layers => $layers, gstid => $gstid ) if $layers->{get_nt_seq};
	my $td = new Benchmark if $BENCHMARK;
	my $init_c_time = timestr( timediff( $td, $tc ) ) if $BENCHMARK;

	unless ( $c->chr_length )
	{
		warn "error initializing the chromosome object.  Failed for valid chr_length\n";
		return (0);
	}

	my $t1  = new Benchmark if $BENCHMARK;
	my $taa = new Benchmark if $BENCHMARK;
	print STDERR "processing features\n" if $self->DEBUG;
	foreach my $item ( $dsg->datasets( chr => $chr ) )
	{
		$self->process_features( start => $start, stop => $c->stop, chr => $chr, ds => $item, dsg=>$dsg, coge => $coge, c => $c, fids => $fids, fnames => $fnames, layers => $layers, gstid => $gstid, ftid => $ftid, color => $color ) unless $simple;
		my $tab = new Benchmark if $BENCHMARK;
		my $feat_time = timestr( timediff( $tab, $taa ) ) if $BENCHMARK;
		print STDERR " processing features for dsid " . $item->name, " ", Dumper($layers), " (", $item->id, "):   $feat_time\n" if $BENCHMARK;
	}

	if ($expids && $layers->{quant} && ref ($expids) =~ /array/i)
	{
		foreach my $expid (@$expids)
		{
			$self->process_experiment(start => $start, stop => $c->stop, chr => $chr, coge => $coge, c => $c, color => $color, expid=>$expid ) unless $simple;
		}
	}
	unless ( $layers->{ruler} || $layers->{chromosome} || $layers->{all} )
	{
		$c->draw_ruler(0);
		$c->draw_chromosome(0);
	}

	my $t2 = new Benchmark if $BENCHMARK;
	print STDERR "generating image: $file\n" if $self->DEBUG;
	$self->generate_output( file => $file, c => $c );
	print STDERR "generating image map : $img_map_name\n" if $self->DEBUG && $img_map_name;
	my $img_map = $c->generate_imagemap( name => $img_map_name ) if $img_map_name;
	my $t3 = new Benchmark if $BENCHMARK;
	my $init_time    = timestr( timediff( $t1, $t0 ) ) if $BENCHMARK;
	my $feature_time = timestr( timediff( $t2, $t1 ) ) if $BENCHMARK;
	my $draw_time    = timestr( timediff( $t3, $t2 ) ) if $BENCHMARK;
	print STDERR "Image height: ", $c->ih, "\n" if $self->DEBUG;

	print STDERR qq{
BENCHMARKING GenomePNG.pl

  Finding DID:          $finddid_time
  Init chr obj:         $init_c_time
Total Initialization:   $init_time
Processing features:    $feature_time
Image generation:       $draw_time

} if $BENCHMARK;

	return $img_map if $img_map;
}

sub initialize_c
{
	my $self = shift;

	#    unless ref($self) =~ /Graphics/ unshift @_, $self;
	my %opts               = @_;
	my $chr                = $opts{chr};
	my $iw                 = $opts{iw};
	my $ih                 = $opts{ih};
	my $start              = $opts{start};
	my $stop               = $opts{stop};
	my $c                  = $opts{c};
	my $csh                = $opts{csh} || $opts{chr_height};
	my $fh                 = $opts{fh} || $opts{feature_height};
	my $draw_chr           = $opts{draw_chr};
	my $draw_ruler         = $opts{draw_ruler};
	my $draw_chr_end       = $opts{draw_chr_end};
	my $feature_labels     = $opts{feature_labels};
	my $fill_labels        = $opts{fill_labels};
	my $debug              = $opts{debug} || 0;
	my $invert_chromosome  = $opts{invert_chromosome};
	my $major_tick_labels  = $opts{major_tick_labels};
	my $minor_tick_labels  = $opts{minor_tick_labels};
	my $draw_hi_qual       = $opts{draw_hi_qual};
	my $padding            = $opts{padding};
	my $max_track          = $opts{max_track};
	my $overlap_adjustment = $opts{overlap_adjustment};
	$debug = 1 if $c->DEBUG;
	$draw_ruler = 1 unless defined $opts{draw_ruler};
	$c->iw($iw);
	$c->ih($ih) if $ih;
	$c->DEBUG($debug);
	$c->feature_labels($feature_labels)         if defined $feature_labels;
	$c->fill_labels($fill_labels)               if defined $fill_labels;
	$c->draw_chromosome($draw_chr)              if defined $draw_chr;
	$c->draw_ruler($draw_ruler)                 if defined $draw_ruler;
	$c->draw_chr_end($draw_chr_end)             if defined $draw_chr_end;
	$c->chr_height($csh)                        if defined $csh;
	$c->feature_height($fh)                     if defined $fh;
	$c->major_tick_labels($major_tick_labels)   if defined $major_tick_labels;
	$c->minor_tick_labels($minor_tick_labels)   if defined $minor_tick_labels;
	$c->overlap_adjustment($overlap_adjustment) if defined $overlap_adjustment;
	$c->invert_chromosome($invert_chromosome)   if defined $invert_chromosome;
	$c->draw_hi_qual($draw_hi_qual)             if defined $draw_hi_qual;
	$c->padding($padding)                       if defined $padding;
	$c->_max_track($max_track)                  if $max_track;
	$c->set_region( start => $start, stop => $stop );
	$start = $c->region_start;
	$stop  = $c->region_stop;
	return ( $start, $stop );
}

sub process_nucleotides
{

	#layers that require nt sequence:
	# gc
	# nt
	# all
	# gbox
	# gaga
	my $self     = shift;
	my %opts     = @_;
	my $start    = $opts{start};
	my $stop     = $opts{stop};
	my $seq      = $opts{seq};
	my $gstid    = $opts{gstid};
	my $at_color = $opts{at_color};
	my $gc_color = $opts{gc_color};
	my $n_color  = $opts{n_color};
	my $x_color  = $opts{x_color};
	$start = 1 unless $start;
	$stop = length $seq if $seq && !$stop;

	if ( $self->MAX_NT && abs( $stop - $start ) > $self->MAX_NT() )
	{
		warn "exceeded nucleotide limit of ", $self->MAX_NT(), " (requested: " . abs( $stop - $start ) . ")\n";
		return;
	}
	my $chr    = $opts{chr};
	my $ds     = $opts{ds};
	my $dsg     = $opts{dsg};
	my $layers = $opts{layers};
	my $c      = $opts{c};

	#process nucleotides
	my $t7 = new Benchmark if $BENCHMARK;
	unless ($seq)
	{
		#return unless $ds;
		if ($ds) {$seq = uc( $ds->get_genomic_sequence( start => $start, end => $stop, chr => $chr, gstid => $gstid ) );}
		if ($dsg) {$seq = uc( $dsg->get_genomic_sequence( start => $start, end => $stop, chr => $chr) );}
	}
	my $t8      = new Benchmark if $BENCHMARK;
	my $seq_len = length $seq;

	#    my $rstop = $c->chr_end ? $c->chr_end : $c->stop;
	my $chrs = int( ( $c->region_length ) / $c->iw );
	$chrs = 1 if $chrs < 1;
	my $pos = 0;
	$start = 1 if $start < 1;
	if ( $layers->{gc} || $layers->{nt} || $layers->{all} )
	{
		while ( $pos < $seq_len )
		{
			my $subseq = substr( $seq, $pos, $chrs );
			my $rcseq  = substr( $seq, $pos, $chrs );
			$rcseq =~ tr/ATCG/TAGC/;
			next unless $subseq && $rcseq;
			if ( !$layers->{gc} && !$layers->{gaga} && !$layers->{gbox} && !$layers->{all} && $subseq !~ /N/i && $subseq !~ /X/i )
			{
				$pos += $chrs;
				next;
			}

			my $options = $layers->{gc} ? "gc" : "nt";
			my $f1 = CoGe::Graphics::Feature::NucTide->new( { nt => $subseq, strand => 1,  start => $pos + $start, options => $options } );
			my $f2 = CoGe::Graphics::Feature::NucTide->new( { nt => $rcseq,  strand => -1, start => $pos + $start, options => $options } );
			$f1->at_color($at_color) if $at_color && ref($at_color) =~ /array/i && scalar( @{$at_color} ) eq 3;
			$f1->gc_color($gc_color) if $gc_color && ref($gc_color) =~ /array/i && scalar( @{$gc_color} ) eq 3;
			$f1->n_color($n_color)   if $n_color  && ref($n_color)  =~ /array/i && scalar( @{$n_color} )  eq 3;
			$f1->x_color($x_color)   if $x_color  && ref($x_color)  =~ /array/i && scalar( @{$x_color} )  eq 3;
			$f2->at_color($at_color) if $at_color && ref($at_color) =~ /array/i && scalar( @{$at_color} ) eq 3;
			$f2->gc_color($gc_color) if $gc_color && ref($gc_color) =~ /array/i && scalar( @{$gc_color} ) eq 3;
			$f2->n_color($n_color)   if $n_color  && ref($n_color)  =~ /array/i && scalar( @{$n_color} )  eq 3;
			$f2->x_color($x_color)   if $x_color  && ref($x_color)  =~ /array/i && scalar( @{$x_color} )  eq 3;
			$f1->show_label(1);
			$f2->show_label(1);
			$f1->font_size(13);
			$f2->font_size(13);
			$f1->use_external_image(1);
			$f2->use_external_image(1);
			$c->add_feature($f1) if $f1;
			$c->add_feature($f2) if $f2;
			$pos += $chrs;
		}
	}
	my $step = $chrs;
	$step = 6 if $step < 6;
	$pos = 0;
	if ( $layers->{gbox} )
	{
		my $startseq;
		if ($ds) {$startseq = uc( $ds->get_genomic_sequence( start => $start - 5, end => $start - 1, chr => $chr, gstid => $gstid ) );}
		if ($dsg) {$startseq = uc( $dsg->get_genomic_sequence( start => $start - 5, end => $start - 1, chr => $chr, gstid => $gstid ) );}
		my $stopseq;
		if ($ds) {$stopseq = uc( $ds->get_genomic_sequence( start => $stop + 1,  end => $stop + 5,  chr => $chr, gstid => $gstid ) );}
		if ($dsg) {$stopseq = uc( $dsg->get_genomic_sequence( start => $stop + 1,  end => $stop + 5,  chr => $chr, gstid => $gstid ) );}
		my $tmpseq   = $seq;
		$tmpseq = $startseq . $tmpseq if $startseq;
		$tmpseq .= $stopseq if $stopseq;

		while ( $pos < length($tmpseq) )
		{

			#	    print STDERR $pos,":",":",$chrs,"\n";
			my $subseq = substr( $tmpseq, $pos, $step );
			my $rcseq = $subseq;
			$rcseq =~ tr/ATCG/TAGC/;
			next unless $subseq && $rcseq;
			if ( $subseq =~ /CACGTG/ || $subseq =~ /GTGCAC/ )
			{
				my $f1 = CoGe::Graphics::Feature::GBox->new( { nt => $subseq, strand => 1,  start => $pos + $start - length($startseq), stop => $pos + $start + $step - length($startseq) - 1 } );
				my $f2 = CoGe::Graphics::Feature::GBox->new( { nt => $rcseq,  strand => -1, start => $pos + $start - length($startseq), stop => $pos + $start + $step - length($startseq) - 1 } );
				$f1->show_label(1);
				$f2->show_label(1);

				#	    $f1->use_external_image(1);
				#	    $f2->use_external_image(1);
				$c->add_feature($f1) if $f1;
				$c->add_feature($f2) if $f2;
			}
			if ( $chrs < 6 )
			{
				$pos += 1;
			}
			else
			{
				$pos += $chrs - 5;
			}
		}
	}
	$step = $chrs;
	$step = 4 if $step < 4;
	$pos  = 0;
	if ( $layers->{gaga} || $layers->{all} )
	{
		while ( $pos < $seq_len )
		{
			my $subseq = substr( $seq, $pos, $step );
			my $rcseq  = substr( $seq, $pos, $step );
			$rcseq =~ tr/ATCG/TAGC/;
			next unless $subseq && $rcseq;
			my $f1 = CoGe::Graphics::Feature::GAGA->new( { nt => $subseq, strand => 1,  start => $pos + $start, stop => $pos + $start + $step - 1 } );
			my $f2 = CoGe::Graphics::Feature::GAGA->new( { nt => $rcseq,  strand => -1, start => $pos + $start, stop => $pos + $start + $step - 1 } );
			$f1->gd;
			$f2->gd;
			$c->add_feature($f1) if $f1 && $f1->gaga_count;
			$c->add_feature($f2) if $f2 && $f2->gaga_count;
			$pos += $chrs;
		}
	}
	my $t9 = new Benchmark if $BENCHMARK;
	my $time1 = timestr( timediff( $t8, $t7 ) ) if $BENCHMARK;
	my $time2 = timestr( timediff( $t9, $t8 ) ) if $BENCHMARK;
	print STDERR qq{
 Time to get nt seq              :  $time1
 Time to process nt for graphics :  $time2
} if $BENCHMARK;

	return $chrs;
}

sub process_features
{
	my $self = shift;

	#process features
	my %opts        = @_;
	my $start       = $opts{start};
	my $stop        = $opts{stop};
	my $chr         = $opts{chr};
	my $ds          = $opts{ds};
	my $dsg          = $opts{dsg};
	my $dsid        = $opts{dsid};
	my $coge        = $opts{coge};
	my $c           = $opts{c};
	my $accn        = $opts{accn};
	my $print_names = $opts{print_names};
	my $fids        = $opts{fids};
	my $fnames      = $opts{fnames};
	my $layers      = $opts{layers};
	my $gstid       = $opts{gstid};
	my $ftid        = $opts{ftid};          #used in quant only
	my $color       = $opts{color};         #used in quant only
	return unless $layers->{all} || $layers->{features};    #no need to get features if no layers requiring them are requested
	$start = $c->region_start unless $start;
	$stop  = $c->region_stop  unless $stop;
	$start = 0 if $start < 0;
	my $tf1 = new Benchmark if $BENCHMARK;
	my $tf2 = new Benchmark if $BENCHMARK;
	my $tf3 = new Benchmark if $BENCHMARK;
	my @cds_feats;
	$dsid = $ds->id if $ds;
	my @dsids;
	push @dsids, $dsid if $dsid;
	push @dsids, map {$_->id} $dsg->datasets if $dsg;
	my @feats;
	foreach my $id (@dsids)
	  {
	 	push @feats, $coge->get_features_in_region( start => $start, end => $stop, dataset => $id, chr => $chr, ftid => $ftid );
	  }
	my @tmp1;
	my @tmp2;
	shift @feats while ( scalar @feats && !$feats[0] );
	return unless scalar @feats;
	shift @feats while ( scalar @feats && $feats[0]->type->name =~ /(contig|chr|source)/ );
	return unless scalar @feats;
	my $research = 0;
	my ($tmpfeat) = sort { $a->start <=> $b->start } @feats;

	#    print STDERR "$start - $stop\n";
	if ( $c->overlap_adjustment && $tmpfeat->start < $start )
	{
		$start    = $tmpfeat->start;
		$research = 1;
	}
	($tmpfeat) = sort { $b->stop <=> $a->stop } @feats;
	if ( $c->overlap_adjustment && $tmpfeat->stop > $stop )
	{
		$stop     = $tmpfeat->stop;
		$research = 1;
	}
	if ($research)
	  {
            @feats = undef;
            foreach my $id (@dsids)
             {
                push @feats, $coge->get_features_in_region( start => $start, end => $stop, dataset => $id, chr => $chr, ftid => $ftid );
             }
	  }

	#    print STDERR "$start - $stop\n";
	my %feats = map { $_->id, $_ } @feats;
	my $tf4 = new Benchmark if $BENCHMARK;

	#let's find and color local duplications
	my ($tandem_type_group) = $coge->resultset('AnnotationTypeGroup')->search( { name => "Tandem duplicates" } );

	#    my ($parent_type) = $coge->resultset('AnnotationType')->search({name=>"Parent", annotation_type_group_id=>$anno_type_group->id});
	#    my ($daughter_type) = $coge->resultset('AnnotationType')->search({name=>"Daughter", annotation_type_group_id=>$anno_type_group->id});
	my ($gevo_link_group) = $coge->resultset('AnnotationTypeGroup')->search( { name => "GEvo link" } );
feat: foreach my $feat ( values %feats )
	{
		my $tf4a = new Benchmark if $BENCHMARK;
		my @f;

		#	print STDERR "Feat info: Name: ",$feat->type->name,", Type: ",$feat->type->name, ", Loc: ", $feat->start,"-",$feat->stop,"\n" if $feat->type->name =~/gene/i;# if $self->DEBUG;
		if ( ( $layers->{features}{pseudogene} || $layers->{all} ) && $feat->type->name =~ /pseudogene/i )
		{
			my $f = CoGe::Graphics::Feature::Gene->new();
			$f->color( [ 255, 33, 0, 50 ] );
			foreach my $loc ( $feat->locs )
			{
				$f->add_segment( start => $loc->start, stop => $loc->stop );
				$f->strand( $loc->strand );
			}
			$f->order(1);
			$f->overlay(1);
			$f->mag(0.5);
			$f->no_3D(1) if $layers->{features}{flat};
			push @f, $f;
		}
		elsif ( ( $layers->{features}{gene} || $layers->{all} ) && $feat->type->name =~ /^Gene$/i )
		{
			my $f = CoGe::Graphics::Feature::Gene->new();
			$f->color( [ 219, 219, 219, 50 ] );
			foreach my $loc ( $feat->locs )
			{
				$f->add_segment( start => $loc->start, stop => $loc->stop );
				$f->strand( $loc->strand );
			}
			$f->order(1);
			$f->overlay(1);
			$f->mag(0.5);
			$f->no_3D(1) if $layers->{features}{flat};
			push @f, $f;
		}
		elsif ( ( $layers->{features}{transposable} || $layers->{all} ) && $feat->type->name =~ /transposable/i )
		{
			my $f = CoGe::Graphics::Feature::Gene->new();
			$f->color( [ 200, 0, 100 ] );
			foreach my $loc ( $feat->locs )
			{
				$f->add_segment( start => $loc->start, stop => $loc->stop );
				$f->strand( $loc->strand );
			}
			$f->order(1);
			$f->overlay(1);
			$f->mag(0.75);
			$f->no_3D(1) if $layers->{features}{flat};
			push @f, $f;
		}
		elsif ( ( $layers->{features}{quant} ) && $feat->type->name =~ /quantitation/i )
		{
			my $f = CoGe::Graphics::Feature::Quant->new();
			$color = [ 200, 0, 200 ] unless defined $color;
			$f->color($color);
			foreach my $loc ( $feat->locs )
			{
				$f->strand( $loc->strand );
			}
			foreach my $quant ( $feat->quantitations )    #this is seriously flawed.  Only the last quant value will be shown.  Need to write in a whole support structure for selecting sets of quant data, and having a special way to visualize them.  Perhaps a combo of selecting the right features and sending an array of values to the Quant.pm object.  But for a future time.
			{
				$f->fill_height( $quant->value );
			}
			push @f, $f;
		}
		elsif ( ( $layers->{features}{gevo_link} || $layers->{all} ) )
		{

			#next unless $coge->resultset('Annotation')->count({feature_id=>$feat->id, annotation_type_id=>$gevo_link_type->id});
			next unless $coge->resultset('FeatureAnnotation')->count( { feature_id => $feat->id, annotation_type_group_id => $gevo_link_group->id }, { join => 'annotation_type' } );
			my $f = CoGe::Graphics::Feature::Link->new();
			$f->color( [ 50, 200, 200, 50 ] );
			foreach my $loc ( $feat->locs )
			{
				$f->add_segment( start => $loc->start, stop => $loc->stop );
				$f->strand( $loc->strand );
			}
			$f->order(2);
			$f->overlay(1);
			$f->mag(.6);
			push @f, $f;
		}
		elsif ( ( $layers->{features}{cds} || $layers->{all} ) && $feat->type->name =~ /CDS/i )
		{
			my $f = CoGe::Graphics::Feature::Gene->new();
			my $color = [ 0, 255, 0, 50 ];
			$f->color($color);
			if ( $layers->{features}{cbc} || $layers->{features}{cbc50} )
			{
				my $seq = $feat->genomic_sequence( gstid => $gstid );
				$f->sequence($seq);
				$f->color_by_codon(1);
				$f->codon_limit(50) if $layers->{features}{cbc50};
			}
			$f->order(1);
			$f->overlay(3);
			$f->no_3D(1) if $layers->{features}{flat};
			push @f,         $f;
			push @cds_feats, $feat;
		}
		elsif ( ( $layers->{features}{mrna} || $layers->{all} ) && $feat->type->name =~ /mrna/i || $feat->type->name =~ /exon/i )
		{
			my $f = CoGe::Graphics::Feature::Gene->new();
			$f->color( [ 0, 0, 255, 50 ] );
			$f->order(1);
			$f->overlay(2);
			$f->mag(0.75);
			$f->no_3D(1) if $layers->{features}{flat};
			push @f, $f;
		}
		elsif ( ( $layers->{features}{rna} || $layers->{all} ) && $feat->type->name =~ /[^m]rna/i && !$feat->type->name =~ /quantitation/i )
		{
			my $f = CoGe::Graphics::Feature::Gene->new();
			$f->color( [ 200, 200, 200, 50 ] );
			$f->order(1);
			$f->overlay(2);
			$f->no_3D(1) if $layers->{features}{flat};
			push @f, $f;
		}
		elsif ( ( $layers->{features}{domain} || $layers->{all} ) && $feat->type->name =~ /functional domains/i )
		{
			my $f = CoGe::Graphics::Feature::Domain->new();
			$f->order(2);
			$f->overlay(2);
			push @f, $f;
		}
		elsif ( ( $layers->{features}{repeat} || $layers->{all} ) && $feat->type->name =~ /repeat/i )
		{
			my $f = CoGe::Graphics::Feature::Outline->new();
			$f->order(1);
			$f->overlay(1);
			my $color = [ 0, 0, 255 ];
			$f->color($color);
			push @f, $f;
		}
		elsif ( ( $layers->{features}{cns} || $layers->{all} ) && $feat->type->name =~ /CNS/i )
		{
			my $f = CoGe::Graphics::Feature::Block->new();
			$f->order(1);

			#$f->overlay(0);
			my $color = [ 255, 102, 0 ];
			$f->color($color);
			push @f, $f;
		}
		elsif ( ( $layers->{features}{gene_space} || $layers->{all} ) && $feat->type->name =~ /gene_space/i )
		{
			my $f = CoGe::Graphics::Feature::Domain->new();
			$f->order(1);
			$f->overlay(0);
			$f->transparency(.5);
			my $color = [ 255, 255, 200 ];
			$f->color($color);
			push @f, $f;
		}
		elsif ( $feat->type->name =~ /source/i || $feat->type->name =~ /chromosome/i )
		{
			next;
		}
		elsif ( $layers->{features}{other} || $layers->{all} )
		{
			my @not_other = qw(rna gene cds functional repeat CNS gene_space peptide protein);
			foreach my $item (@not_other)
			{
				next feat if $feat->type->name =~ /$item/i;
			}

			my $f = CoGe::Graphics::Feature::Block->new();
			$f->order(2);
			$f->overlay(3);
			my $color = [ 255, 100, 0 ];
			$f->color($color);
			push @f, $f;
		}
		if ( ( $layers->{features}{protein} || $layers->{all} ) && $feat->type->name =~ /CDS/i )
		{

			$self->process_proteins( genomic_feat => $feat, c => $c );

		}
		if ( ( $layers->{features}{local_dup} || $layers->{all} ) && $feat->type->name =~ /CDS/i )
		{
			my $color;
			if ( $coge->resultset('FeatureAnnotation')->count( { feature_id => $feat->id, annotation_type_group_id => $tandem_type_group->id }, { join => 'annotation_type' } ) )
			{
				$color = [ 0, 225, 255, 50 ];
			}

			#	    elsif ($coge->resultset('Annotation')->count({feature_id=>$feat->id, annotation_type_id=>$daughter_type->id}))
			#	      {
			#		$color = [0,255,255, 50];
			#	      }
			if ($color)
			{
				my $f = CoGe::Graphics::Feature::Gene->new();
				$f->color($color);
				$f->order(1);
				$f->overlay(4);
				$f->no_3D(1) if $layers->{features}{flat};
				push @f, $f;
			}
		}

		my $tf4b = new Benchmark if $BENCHMARK;
		next unless scalar @f;
		foreach my $id (@$fids)
		{
			next unless $id;
			foreach my $f (@f)
			{
				$f->color( [ 255, 255, 0 ] ) if $feat->id eq $id;
			}
		}
		my $tf4c = new Benchmark if $BENCHMARK;
		foreach my $name (@$fnames)
		{
			next unless $name;
			foreach my $n ( map { $_->name } $feat->names )
			{
				if ( $n =~ /^$name$/i )
				{
					foreach my $f (@f)
					{
						$f->color( [ 255, 255, 0 ] );
						$f->label($name);
					}
				}
			}
		}
		my $tf4d = new Benchmark if $BENCHMARK;
		foreach my $loc ( $feat->locs )
		{
			print STDERR "\tadding feature location: ", $loc->start, "-", $loc->stop, ". . ." if $self->DEBUG;
			foreach my $f (@f)
			{
				$f->add_segment( start => $loc->start, stop => $loc->stop );
				$f->strand( $loc->strand );
			}
			print "done!\n" if $self->DEBUG;
		}
		my $tf4e = new Benchmark if $BENCHMARK;
		my ($name) = $feat->names;

		#	$f->description($feat->annotation_pretty_print);
		foreach my $f (@f)
		{
			$f->link( "GeLo.pl?" . join( "&", "chr=" . $feat->chr, "ds=" . $feat->dataset_id, "z=10", "INITIAL_CENTER=" . $feat->start . ",0" ) );
			$f->label($name) if $print_names;
			$f->type( $feat->type->name );
			$f->skip_overlap_search(1) unless $c->overlap_adjustment;

			#	print STDERR Dumper $f;

			$c->add_feature($f);
		}
		my $tf4f = new Benchmark if $BENCHMARK;
		print STDERR "\tfinished process feature.\n" if $self->DEBUG;
	}
	my $tf5 = new Benchmark if $BENCHMARK;

	my $time1 = timestr( timediff( $tf2, $tf1 ) ) if $BENCHMARK;
	my $time2 = timestr( timediff( $tf3, $tf2 ) ) if $BENCHMARK;
	my $time3 = timestr( timediff( $tf4, $tf3 ) ) if $BENCHMARK;
	my $time4 = timestr( timediff( $tf5, $tf4 ) ) if $BENCHMARK;
	print STDERR qq{
sub process_features:
 Time to get chr length            :  $time1
 Time to count features            :  $time2
 Time to get features              :  $time3
 Time to process features for image:  $time4
} if $BENCHMARK;
	return \@cds_feats;
}

sub process_experiment
{
	my $self = shift;
	my %opts = @_;
	my $start = $opts{start};
	my $stop = $opts{stop};
	my $chr = $opts{chr};
	my $expid = $opts{expid};
	my $c = $opts{c};
	my $color=$opts{color};

	# We will need to abstract out the call to the data engine in order for
	# this to be compatible on systems with multiple intsallations of coge
	# using different databases
	my $url = "http://genomevolution.org/CoGe/bin/fastbit_query.pl?exp_id=$expid;chr=$chr;start=$start;stop=$stop"; #FIXME hardcoded server

	# mdb removed 3/26/13
#	my $cmd = "curl '$url'";
#	my $result = `$cmd`;
#	print STDERR "$cmd\n$result\n";

	# mdb added 3/26/13
	#my $url = "http://geco.iplantcollaborative.org/mbomhoff/CoGe/bin/fastbit_query.pl?exp_id=$expid;chr=$chr;start=$start;stop=$stop";
	my $result = LWP::Simple::get($url);
	#print STDERR "$url\n$result\n";
	return unless $result;

	my $data = decode_json($result);
	foreach my $quant (@{$data->{results}})
	{
		my $f = CoGe::Graphics::Feature::Quant->new();
		$color = [ 200, 0, 200 ] unless defined $color;
		$f->color($color);
		$f->fill_height( $quant->[4] );
		print STDERR "\tadding quant location: ", $quant->[1], "-", $quant->[2], ". . ." if $self->DEBUG;
		$f->add_segment( start => $quant->[1], stop => $quant->[2] );
		$f->strand( $quant->[3] );
		$f->skip_overlap_search(1);
		$c->add_feature($f);
	}
}

sub generate_output
{
	my $self = shift;
	my %opts = @_;
	my $file = $opts{file};
	my $c    = $opts{c};
	if ($file)
	{
		$c->generate_png( file => $file );
	}
	else
	{
		print "Content-type: image/png\n\n";
		$c->generate_png();
	}
}

sub process_proteins
{
	my $self = shift;
	return $self->draw_prots(@_);
}

sub draw_prots
{
	my $self = shift;
	my %opts = @_;
	my $feat = $opts{genomic_feat};
	my $c    = $opts{c};

	#Do we have any protein sequence we can use?
	my $chrs;
	foreach my $seq ( $feat->sequences )
	{
		next unless $seq->sequence_type->name =~ /prot/i;
		my ($pseq) = $seq->sequence_data;
		my $lastaa = 0;
		$chrs = int( ( ( $c->region_length ) / $c->iw ) / 3 + .5 );
		$chrs = 1 if $chrs < 1;
		my $pos = 0;
		while ( $pos < length $pseq )
		{
			my $aseq = substr( $pseq, $pos, $chrs );
			if ( $pos + $chrs == ( length $pseq ) ) { $lastaa = 1; }
			foreach my $loc ( $seq->get_genomic_locations( start => $pos + 1, stop => $pos + $chrs ) )
			{
				next if $loc->{stop} < $c->start;
				next if $loc->{start} > $c->stop;
				print STDERR "Adding protein feature: $aseq at position ", $loc->{start}, "-", $loc->{stop}, ": region: ", $c->start, "-", $c->stop, "\n" if $self->DEBUG;
				my $ao = CoGe::Graphics::Feature::AminoAcid->new( { aa => $aseq, start => $loc->{start}, stop => $loc->{stop}, strand => $loc->{strand}, order => 2, lastaa => $lastaa } );
				$ao->skip_overlap_search(1);
				$ao->mag(0.75);
				$ao->overlay(2);
				$c->add_feature($ao);
				delete $loc->{__Changed};    #silence the warning from Class::DBI when location is destroyed
			}

			$pos += $chrs;
		}
	}
	return $chrs;
}

sub process_layers
{
	my $layers = shift;
	my %valid_layers = (
		ruler              => "ruler",
		chromosome         => "chromosome",
		chr                => "chromosome",
		cds                => "cds",
		coding             => "cds",
		exons              => "cds",
		exon               => "cds",
		rna                => "rna",
		mrna               => "mrna",
		gene               => "gene",
		proteins           => "protein",
		protein            => "protein",
		pro                => "protein",
		funcitonal_domains => "domain",
		domains            => "domain",
		domain             => "domain",
		other              => "other",
		cns                => "cns",
		nt                 => "nt",               #requres nt sequence
		nucleotides        => "nt",               #requres nt sequence
		nucleotide         => "nt",               #requres nt sequence
		nuc                => "nt",               #requres nt sequence
		gc                 => "gc",               #requres nt sequence
		background         => "background",
		all                => "all",              #requres nt sequence
		pseudogene         => "pseudogene",
		gene_space         => "gene_space",
		"local_dup"        => "local_dup",
		"local_dups"       => "local_dup",
		"tandem"           => "local_dup",
		"gaga"             => "gaga",             #requres nt sequence
		"gbox"             => "gbox",             #requres nt sequence
		"cbc"              => "cbc",              #color CDS by codon
		"cbc50"            => "cbc50",            #color CDS by codon
		"flat"             => "flat",             #are gene models draw 'flat' or pseudo-3D?
		"olc"              => "overlap_check",    #are features checked for overlap when drawing image?
		"overlap_check"    => "overlap_check",
		"repeats"          => "repeat",
		"repeat"           => "repeat",
		"repeat_region"    => "repeat",
		"gevo_link"        => "gevo_link",
		"TE"               => "transposable",
		"transposable"     => "transposable",
		"quant"            => "quant",            #quantitation
	);
	my %features = (
		cds             => 1,
		mrna            => 1,
		protein         => 1,
		gene            => 1,
		pseudogene      => 1,
		domain          => 1,
		other           => 1,
		cns             => 1,
		rna             => 1,
		local_dup       => 1,
		cbc             => 1,
		cbc50           => 1,
		gene_space      => 1,
		"flat"          => 1,
		"overlap_check" => 1,
		"repeat"        => 1,
		"gevo_link"     => 1,
		"transposable"  => 1,
		"quant"         => 1,
	);

	#determine of nt sequence is needed
	my %nt = (
		all  => 1,
		nt   => 1,
		gc   => 1,
		gaga => 1,
		gboc => 1,
	);
	my %layers;
	foreach my $layer (@$layers)
	{

		#	my $strand;
		$layer = lc($layer);

		next unless $valid_layers{$layer};

		#	$layer = "rna" if $layer =~ /rna/;
		if ( $features{ $valid_layers{$layer} } )
		{
			$layers{features}{ $valid_layers{$layer} } = 1;
		}
		else
		{
			$layers{ $valid_layers{$layer} } = 1;
		}
		$layers{quant} =1 if $layer eq "quant"; #FIXME this is temporary while "quant" is a feature option.  to remove when demo is removed from CoGe server.  Note to Matt and Eric.  Remove this.
		$layers{get_nt_seq} = 1 if ( $nt{ $valid_layers{$layer} } );    #flag for whether nt sequence needs to be retrieved
	}
	$layers{all} = 1 unless keys %layers;
	return \%layers;
}

1;
