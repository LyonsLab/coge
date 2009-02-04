package CoGeX::DatasetGroup;

# Created by DBIx::Class::Schema::Loader v0.03009 @ 2006-12-01 18:13:38

use strict;
use warnings;
use base 'DBIx::Class';

__PACKAGE__->load_components("PK::Auto", "ResultSetManager", "Core");
__PACKAGE__->table("dataset_group");
__PACKAGE__->add_columns(
  "dataset_group_id",{ data_type => "INT", default_value => undef, is_nullable => 0, size => 10 },
  "name",{ data_type => "VARCHAR", default_value => "", is_nullable => 0, size => 50 },
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
  "organism_id",{ data_type => "INT", default_value => 0, is_nullable => 0, size => 10 },
  "genomic_sequence_type_id",{ data_type => "INT", default_value => 0, is_nullable => 0, size => 10 },
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


sub get_genomic_sequence_new {
  my $self = shift;
  my %opts = @_;
  my $start = $opts{start} || $opts{begin};
  my $stop = $opts{stop} || $opts{end};
  my $chr = $opts{chr} || $opts{chromosome};
  my $strand = $opts{strand};
  my $skip_length_check = $opts{skip_length_check} || 0;
  my $str = "";
  if (defined $start && defined $stop)
    {
      $chr = "1" unless defined $chr;
      $start = 1 if $start < 1;
      $stop = $start unless $stop;
      return undef unless ($start =~ /^\d+$/ and  $stop =~ /^\d+$/);
      ($start, $stop) = ($stop, $start) if $stop < $start;

      if (! $skip_length_check)
	{
	  my $last = $self->last_chromosome_position($chr);
	  $stop = $last if $stop > $last;
	}
    }
  my $SEQDIR = "/opt/apache/CoGe/data/genomic_sequence";
  my $chr_total = $self->genomic_sequences->count();
  my $chr_level=1;
  if ($chr_total > 1000)
    {
      my $chr_count=0;
      foreach my $testchr (sort $self->chromosomes)
	{
	  last if $chr eq $testchr;
	  $chr_count++;
	  $chr_level++ unless $chr_count%1000;
	}
    }
  my $file = join ("/", $SEQDIR,ceil($self->organism->id/1000), "org_".$self->organism->id, "seqtype_".$self->sequence_type->id, "dataset_".$self->id,$chr_level,$chr);
  return $file;
}

sub get_genomic_sequence {
  my $self = shift;
  my %opts = @_;
  my $start = $opts{start} || $opts{begin};
  my $stop = $opts{stop} || $opts{end};
  my $chr = $opts{chr} || $opts{chromosome};
  my $strand = $opts{strand};
  my $skip_length_check = $opts{skip_length_check} || 0;
  my $str = "";
  if (defined $start && defined $stop)
    {
      $chr = "1" unless defined $chr;
      $start = 1 if $start < 1;
      $stop = $start unless $stop;
      return undef unless ($start =~ /^\d+$/ and  $stop =~ /^\d+$/);
      ($start, $stop) = ($stop, $start) if $stop < $start;

      if (! $skip_length_check)
	{
	  my $last = $self->last_chromosome_position($chr);
	  $stop = $last if $stop > $last;
	}
      # make sure two numbers were sent in
      #my $fstart = $start%10000 ? $start - ($start % 10000) + 1 : ($start -1)- (($start-1) % 10000) +1;
      my $fstart = $start - (($start -1) % 10000);
#      print STDERR "start: $start, fstart: $fstart\n";
      my @starts;
#      push (@starts, $fstart) if $fstart == $stop;
      for(my $i=$fstart;$i<=$stop;$i+=10000){
        push(@starts,$i);
      }
      return unless @starts;
      my @seqs = $self->genomic_sequences(
					  {chromosome=>$chr,
					   -and=>[ start=> { 'in', \@starts } ],
					  },
					  {order_by=>"start asc"}
					 )->all;
      return unless @seqs;
#      print STDERR Dumper \@starts;
      $str = join ("", map{$_->sequence_data} @seqs);  
#      print STDERR "!!",$str,"\n";;
      $str = $self->trim_sequence( $str, $seqs[0]->start, $seqs[-1]->stop, $start, $stop );
#      print STDERR "!!",$str,"\n";;
    } 
  elsif ( $chr ) 
    {
    $str = join("",map { $_->sequence_data } $self->genomic_sequences( { 'chromosome' => $chr},
					      {order_by=>"start asc"}
					    ));
    } 
  else 
    {                 # entire sequence
      my $allseqs = $self->genomic_sequences({},{order_by=>"start asc"});
      while ( my $g = $allseqs->next ) {
	$str .= $g->sequence_data;
      }
    }
  $str = $self->reverse_complement($str) if $strand && $strand =~ /-/;
  return $str;
}

sub get_genome_sequence
  {
    return shift->get_genomic_sequence(@_);
  }
sub genomic_sequence
  {
    return shift->get_genomic_sequence(@_);
  }

sub trim_sequence {
  my $self = shift;
  my( $seq, $seqstart, $seqend, $newstart, $newend ) = @_;
  
  my $start = $newstart-$seqstart;
  my $stop = length($seq)-($seqend-$newend)-1;  
#  print STDERR join ("\t", $seqstart, $seqend, $newstart, $newend),"\n";
#  print STDERR join ("\t", length ($seq), $start, $stop, $stop-$start+1),"\n";
  $seq = substr($seq, $start, $stop-$start+1);
#  print STDERR "final seq lenght: ",length($seq),"\n";
  return($seq);
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


 sub last_chromosome_position_new
   {
     my $self = shift;
     my $chr = shift;
     my ($item) =  $self->genomic_sequences(
					  {
					   chromosome=>"$chr",
					  },
#					 )->get_column('stop')->max;
					 );
     my $stop = $item->stop();
     unless ($stop)
      {
        warn "No genomic sequence for ",$self->name," for chr $chr\n";
        return;
      }
     $stop;
   }

sub last_chromosome_position
   {
     my $self = shift;
     my $chr = shift;
     my $stop =  $self->genomic_sequences(
                                          {
                                           chromosome=>"$chr",
                                          },
                                         )->get_column('stop')->max;
     unless ($stop)
      {
        warn "No genomic sequence for ",$self->name," for chr $chr\n";
        return;
      }
     $stop;
   }



sub sequence_type
  {
    my $self = shift;
    my ($type) = $self->genomic_sequences->slice(0,0);
    return $type ? $type->genomic_sequence_type : undef;
  }
sub genomic_sequence_type
  {
    my $self = shift;
    return $self->sequence_type(@_);
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
#     unless (@data)
#       {
# 	foreach my $item ($self->features(
# 					    {
					     
# 					    },
# 					    {
# 					     select=>"locations.chromosome",
# 					     join=>"prefetch",
# 					     locations=>"locations",
# 					     distinct=>["locations.chromosome"],
# 					     prefetch=>["locations"],
# 					    },
# 					    ))
# 	  {
# 	    print $item->chromosome,"\n";
# 	  }
#       }
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
#    my $chr = $opts{chr};
    my $seq = $self->genomic_sequence(%opts);
    my $length = length $seq;
    return unless $length;
    my ($gc) = $seq =~ tr/GCgc/GCgc/;
    my ($at) = $seq =~ tr/ATat/ATat/;
    my ($n) = $seq =~ tr/nNxX/nNxX/;
    return ($gc,$at, $n) if $count;
    return sprintf("%.4f", $gc/$length),sprintf("%.4f", $at/$length),,sprintf("%.4f", $n/$length);
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
    my $stop = $opts{stop} || $self->last_chromosome_position($chr);
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

sub trans_type
  {
    my $self = shift;
    my $trans_type;
    foreach my $feat ($self->features)
      {
	next unless $feat->type->name =~ /cds/i;
	my ($code, $type) = $feat->genetic_code;
	($type) = $type =~/transl_table=(\d+)/ if $type =~ /transl_table/;
	return $type if $type;
      }
    return 1; #universal genetic code type;
  }

sub reverse_complement
  {
    my $self = shift;
    my $seq = shift;# || $self->genomic_sequence;
    my $rcseq = reverse($seq);
    $rcseq =~ tr/ATCGatcg/TAGCtagc/; 
    return $rcseq;
  }

1;
