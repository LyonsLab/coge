package CoGeX::DatasetGroup;

# Created by DBIx::Class::Schema::Loader v0.03009 @ 2006-12-01 18:13:38

use strict;
use warnings;
use base 'DBIx::Class';
use File::Spec::Functions;
use Data::Dumper;

__PACKAGE__->load_components("PK::Auto", "ResultSetManager", "Core");
__PACKAGE__->table("dataset_group");
__PACKAGE__->add_columns(
  "dataset_group_id",{ data_type => "INT", default_value => undef, is_nullable => 0, size => 11 },
  "name",{ data_type => "VARCHAR", default_value => "", is_nullable => 0, size => 200 },
  "description",
  {
    data_type => "VARCHAR",
    default_value => undef,
    is_nullable => 1,
    size => 255,
  },
  "version",
  {
    data_type => "VARCHAR",
    default_value => undef,
    is_nullable => 0,
    size => 50,
  },
  "organism_id",{ data_type => "INT", default_value => 0, is_nullable => 0, size => 11 },
  "genomic_sequence_type_id",{ data_type => "INT", default_value => 0, is_nullable => 0, size => 11 },
  "file_path",
  {
    data_type => "VARCHAR",
    default_value => undef,
    is_nullable => 0,
    size => 255,
  },
);

__PACKAGE__->set_primary_key("dataset_group_id");
__PACKAGE__->has_many("dataset_connectors" => "CoGeX::DatasetConnector", 'dataset_group_id');
__PACKAGE__->has_many("genomic_sequences" => "CoGeX::GenomicSequence", 'dataset_group_id');
__PACKAGE__->belongs_to("organism" => "CoGeX::Organism", 'organism_id');
__PACKAGE__->belongs_to("genomic_sequence_type" => "CoGeX::GenomicSequenceType", 'genomic_sequence_type_id');


sub datasets
  {
    my $self = shift;
    my %opts = @_;
    my $chr = $opts{chr} || $opts{chromosome};
    my @ds;
    foreach my $dsc ($self->dataset_connectors())
      {
	if ($chr)
	  {
	    foreach my $ds_chr ($dsc->dataset->chromosomes)
	      {
		return $dsc->dataset if $ds_chr eq $chr;
	      }
	  }
	else
	  {
	    push @ds, $dsc->dataset;
	  }
      }
    return wantarray ? @ds : \@ds;
  }

sub get_genomic_sequence {
  my $self = shift;
  my %opts = @_;
  my $start = $opts{start} || $opts{begin};
  my $stop = $opts{stop} || $opts{end};
  my $chr = $opts{chr} || $opts{chromosome};
  my $strand = $opts{strand};
  my $debug = $opts{debug};
  my $str = "";
  $start = 1 unless $start;
  my $last = $self->sequence_length($chr);
  return if $start > $last && $stop > $last; #outside of chromosome
  return if $start < 1 && $stop < 1;
  $stop = $last unless $stop;
  $stop = $last if $stop > $last;
  $start = $last if $start > $last;
  if (defined $start && defined $stop)
    {
      $chr = "1" unless defined $chr;
      $start = 1 if $start < 1;
      $stop = $start unless $stop;
      return undef unless ($start =~ /^\d+$/ and  $stop =~ /^\d+$/);
      ($start, $stop) = ($stop, $start) if $stop < $start;
    }
  else
    {
      warn "missing parameters in sub get_genomic_sequence\n";
      warn Dumper \%opts;
      return;
    }
#  my ($seq) = $self->genomic_sequences({chromosome=>$chr});
  my $file = $self->file_path();
  return $self->get_seq(db=>$file, chr=>$chr, start=>$start, stop=>$stop, strand=>$strand, debug=>$debug);
}

sub get_seq #using fastacmd to get the sequence
  { 
    my $self = shift;
    my %opts = @_;
#    use Data::Dumper;
#    print STDERR Dumper \%opts;
    my $fastacmd = $opts{fastacmd} || '/usr/bin/fastacmd';
    my $blastdb = $opts{blastdb} || $opts{db};
    my $seqid   = $opts{seqid} || $opts{chr} || $opts{chromosome}; # chr 
    my $debug = $opts{debug};
    my $start = $opts{start};
    ($start) = $start=~/(\d+)/ if $start;
    my $stop = $opts{stop} || $opts{end};
    ($stop) = $stop=~/(\d+)/ if $stop;
    my $cmd = "$fastacmd -d $blastdb -s \"$seqid\" ";
    if($start && $stop)
      {
        $cmd .= "-L " . $start . "," . $stop. " ";
      }
    if($opts{reverse_complement} || $opts{rc} || ($opts{strand} && $opts{strand} =~/-/) )
      {
        $cmd .= "-S 2";
      }
    print STDERR $cmd,"\n" if $debug;

    open(FASTA, $cmd . "|") || die "can't run $cmd";
    # get rid of the header line...
    <FASTA>;
    my $seq = join ("",<FASTA>);
    close FASTA;
    $seq =~ s/\n//g;
    return $seq;
}

sub get_genome_sequence
  {
    return shift->get_genomic_sequence(@_);
  }
sub genomic_sequence
  {
    return shift->get_genomic_sequence(@_);
  }


################################################## subroutine header start ##

=head2 last_chromsome_position

 Usage     : my $last = $genome_seq_obj->last_chromosome_position($chr);
 Purpose   : gets the last genomic sequence position for a dataset given a chromosome
 Returns   : an integer that refers to the last position in the genomic sequence refered
             to by a dataset given a chromosome
 Argument  : string => chromsome for which the last position is sought
 Throws    : 
 Comments  : 

See Also   : 

=cut

################################################## subroutine header end ##


 sub sequence_length
   {
     my $self = shift;
     my $chr = shift;
     return unless $chr;
     my ($item) =  $self->genomic_sequences(
					    {
					     chromosome=>"$chr",
					    },
					   );
     my $stop = $item->sequence_length;
     unless ($stop)
      {
        warn "No genomic sequence for ",$self->name," for chr $chr\n";
        return;
      }
     return $stop;
   }

sub sequence_type
  {
    shift->genomic_sequence_type(@_);
  }

sub type
  {
    shift->genomic_sequence_type(@_);
  }

sub resolve : ResultSet {
    my $self = shift;
    my $info = shift;
    return $info if ref($info) =~ /Dataset/i;
    return $self->find($info) if $info =~ /^\d+$/;
    return $self->search({ 'name' => { '-like' => '%' . $info . '%'}},
			 ,{});
  }

sub get_chromosomes
  {
    my $self = shift;
    my @data =  map {$_->chromosome} $self->genomic_sequences(
					{},
					{
					 select=>{distinct=>"chromosome"},
					 as=>"chromosome",
					},
				       );
    return wantarray ? @data : \@data;
  }


sub chromosomes
  {
    my $self = shift;
    $self->get_chromosomes(@_);
  }
    

sub percent_gc
  {
    my $self = shift;
    my %opts = @_;
    my $count = $opts{count};
    my $seq = $self->genomic_sequence(%opts);
    my $length = length $seq;
    return unless $length;
    my ($gc) = $seq =~ tr/GCgc/GCgc/;
    my ($at) = $seq =~ tr/ATat/ATat/;
    my ($n) = $seq =~ tr/nNxX/nNxX/;
    return ($gc,$at, $n) if $count;
    return sprintf("%.4f", $gc/$length),sprintf("%.4f", $at/$length),sprintf("%.4f", $n/$length);
  }

sub fasta
  {
    my $self = shift;
    my %opts = @_;
    my $col = $opts{col};
    #$col can be set to zero so we want to test for defined variable
    $col = $opts{column} unless defined $col;
    $col = $opts{wrap} unless defined $col;
    $col = 100 unless defined $col;
    my $chr = $opts{chr};
    ($chr) = $self->get_chromosomes unless defined $chr;
    my $strand = $opts{strand} || 1;
    my $start = $opts{start} || 1;
    $start =1 if $start < 1;
    my $stop = $opts{stop} || $self->sequence_length($chr);
    my $prot = $opts{prot};
    my $rc = $opts{rc};
    $strand = -1 if $rc;
    my $seq = $self->genomic_sequence(start=>$start, stop=>$stop, chr=>$chr);
    $stop = $start + length($seq)-1 if $stop > $start+length($seq)-1;
    my $head = ">".$self->organism->name." (".$self->name;
    $head .= ", ".$self->description if $self->description;
    $head .= ", v".$self->version.")".", Location: ".$start."-".$stop." (length: ".($stop-$start+1)."), Chromosome: ".$chr.", Strand: ".$strand;

    $Text::Wrap::columns=$col;
    my $fasta;


    $seq = $self->reverse_complement($seq) if $rc;
    if ($prot)
      {
	my $trans_type = $self->trans_type;
	my $feat = new CoGeX::Feature;
	my ($seqs, $type) = $feat->frame6_trans(seq=>$seq, trans_type=>$trans_type);
	foreach my $frame (sort {length($a) <=> length($b) || $a cmp $b} keys %$seqs)
	  {
	    $seq = $seqs->{$frame};
	    $seq = $self->reverse_complement($seq) if $rc;
	    $seq = join ("\n", wrap("","",$seq)) if $col;
	    $fasta .= $head. " $type frame $frame\n".$seq."\n";
	  }
      }
    else
      {
	$seq = join ("\n", wrap("","",$seq)) if $col;
	$fasta = $head."\n".$seq."\n";
      }
    return $fasta;
  }

sub reverse_complement
  {
    my $self = shift;
    my $seq = shift;# || $self->genomic_sequence;
    my $rcseq = reverse($seq);
    $rcseq =~ tr/ATCGatcg/TAGCtagc/; 
    return $rcseq;
  }

#this sub determines the correct directory structure for storing the sequence files for a dataset group
#By Dan Hembry
#The idea is to build a dir structure that holds large amounts of files, and is easy to lookup based on dataset_group ID number.
	#The strucuture is three levels of directorys, and each dir holds 1000 files and/or directorys.
	#Thus:
	#./0/0/0/ will hold files 0-999
	#./0/0/1/ will hold files 1000-1999
	#./0/0/2/ will hold files 2000-2999
	#./0/1/0 will hold files 1000000-1000999 (I think).
	#./level0/level1/level2/

sub get_path
  {
    my $self = shift;
    my $dataset_group_id = $self->id;
    my $level0 = floor($dataset_group_id / 1000000000) % 1000;
    my $level1 = floor($dataset_group_id / 1000000) % 1000;
    my $level2 = floor($dataset_group_id / 1000) % 1000;
    my $path = catdir($level0,$level1,$level2, $dataset_group_id); #adding dataset_group_id for final directory.  blast's formatdb will be run on the faa file and this will help keep that stuff organized
    return $path;
  }


1;
