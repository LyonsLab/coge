package CoGeX::Feature;

use strict;
use warnings;
use base 'DBIx::Class';

__PACKAGE__->load_components("PK::Auto", "ResultSetManager", "Core");
__PACKAGE__->table("feature");
__PACKAGE__->add_columns(
  "feature_id",
  { data_type => "INT", default_value => undef, is_nullable => 0, size => 10 },
  "feature_type_id",
  { data_type => "INT", default_value => 0, is_nullable => 0, size => 10 },
  "dataset_id",
  { data_type => "INT", default_value => 0, is_nullable => 0, size => 10 },
  "start",
  { data_type => "INT", default_value => 0, is_nullable => 1, size => 10 },
  "stop",
  { data_type => "INT", default_value => 0, is_nullable => 1, size => 10 },
  "chromosome",
  { data_type => "VARCHAR", default_value => 0, is_nullable => 1, size => 50 },
  "strand",
  { data_type => "VARCHAR", default_value => 0, is_nullable => 1, size => 2 },
);

__PACKAGE__->set_primary_key("feature_id");

# feature has many feature_names
__PACKAGE__->has_many( 'feature_names' => "CoGeX::FeatureName", 'feature_id');

# feature has many annotations
__PACKAGE__->has_many( 'annotations' => "CoGeX::Annotation", 'feature_id');

# feature has many locations
__PACKAGE__->has_many( 'locations' => "CoGeX::Location", 'feature_id');

# feature has many sequences - note this is non-genomic sequence (see
# genomic_sequence() and sequence
__PACKAGE__->has_many( 'sequences' => "CoGeX::Sequence", 'feature_id');

# feature_type has many features
__PACKAGE__->belongs_to("feature_type" => "CoGeX::FeatureType", 'feature_type_id');

# dataset has many features
__PACKAGE__->belongs_to("dataset" => "CoGeX::Dataset", 'dataset_id');


#__PACKAGE__->mk_group_accessors(['start', 'stop', 'chromosome', 'strand']);

sub esearch : ResultSet {
    my $self = shift;
    my $join = $_[1]{'join'};
    map { push(@$join, $_ ) } 
        ();


    my $prefetch = $_[1]{'prefetch'};
    map { push(@$prefetch, $_ ) } 
        ('feature_type','locations', 
            { 'dataset' => 'organism' }
        );

    $_[1]{'join'} = $join;
    $_[1]{'prefetch'} = $prefetch;
    my $rs = $self->search(@_);
    return $rs;

}
sub type
  {
    my $self = shift;
    return $self->feature_type();
  }

sub organism
  {
    my $self = shift;
    return $self->dataset->organism();
  }

sub org
  {
    my $self = shift;
    return $self->organism();
  }

sub names
  {
    my $self = shift;
    if ($self->{_names})
      {
	return wantarray ? @{$self->{_names}} : $self->{_names};
      }
    my @names =  $self->feature_names()->get_column('name')->all;
    $self->{_names}=\@names;
    return wantarray ? @names : \@names;
  }

sub locs
  {
    my $self = shift;
    return $self->locations();
  }

sub seqs
  {
    my $self = shift;
    return $self->sequences();
  }

sub eannotations
  {
    my $self = shift;
    return $self->annotations(undef,{prefetch=>['annotation_type']});
  }

sub annos
  {
    shift->eannotations(@_);
  }

################################################ subroutine header begin ##

=head2 annotation_pretty_print

 Usage     : my $pretty_annotation = $feat->annotation_pretty_print
 Purpose   : returns a string with information and annotations about a feature
             in a nice format with tabs and new-lines and the like.
 Returns   : returns a string
 Argument  : none
 Throws    : 
 Comments  : uses Coge::Genome::Accessory::Annotation to build the annotations,
           : specifying delimters, and printing to string.   Pretty cool object.

See Also   : CoGe::Genome::Accessory::Annotation

=cut

################################################## subroutine header end ##


sub annotation_pretty_print
  {
    my $self = shift;
    my $anno_obj = new CoGe::Genome::Accessory::Annotation(Type=>"anno");
    $anno_obj->Val_delimit("\n");
    $anno_obj->Val_delimit("\n");
    $anno_obj->Add_type(0);
    $anno_obj->String_end("\n");
    my $start = $self->start;
    my $stop = $self->stop;
    my $chr = $self->chr;
    my $strand = $self->strand;
    #look into changing this to set_id
    my $info_id = $self->dataset->id;
    my $location = "Chr ".$chr." ";
    $location .= join (", ", map {$_->start."-".$_->stop} $self->locs);
    $location .="(".$strand.")";
    #my $location = "Chr ".$chr. "".$start."-".$stop.""."(".$strand.")";
    $anno_obj->add_Annot(new CoGe::Genome::Accessory::Annotation(Type=>"Location", Values=>[$location], Type_delimit=>": ", Val_delimit=>" "));
    my $anno_type = new CoGe::Genome::Accessory::Annotation(Type=>"Name(s)");
    $anno_type->Type_delimit(": ");
    $anno_type->Val_delimit(", ");
    foreach my $name ($self->names)
      {
	$anno_type->add_Annot($name);
      }
    
    $anno_obj->add_Annot($anno_type);
    foreach my $anno (sort {$b->type->name cmp $a->type->name} $self->annos)
      {
	my $type = $anno->type();
	my $group = $type->group();
	my $anno_type = new CoGe::Genome::Accessory::Annotation(Type=>$type->name);
	$anno_type->Val_delimit("\n");

	$anno_type->add_Annot($anno->annotation);
	if (ref ($group) =~ /group/i)
	  {
	    my $anno_g = new CoGe::Genome::Accessory::Annotation(Type=>$group->name);
	    $anno_g->add_Annot($anno_type);
	    $anno_g->Type_delimit(": ");
	    $anno_g->Val_delimit(", ");
	    $anno_obj->add_Annot($anno_g);
	  }
	else
	  {
	    $anno_type->Type_delimit(": ");
	    $anno_obj->add_Annot($anno_type);
	  }
      }
    return $anno_obj->to_String;
  }


################################################ subroutine header begin ##

=head2 annotation_pretty_print_html

 Usage     : my $pretty_annotation_html = $feat->annotation_pretty_print_html
 Purpose   : returns a string with information and annotations about a feature
             in a nice html format with breaks and class tags (called "annotation")
 Returns   : returns a string
 Argument  : none
 Throws    : 
 Comments  : uses Coge::Genome::Accessory::Annotation to build the annotations,
           : specifying delimters, and printing to string.   Pretty cool object.

See Also   : CoGe::Genome::Accessory::Annotation

=cut

################################################## subroutine header end ##


sub annotation_pretty_print_html
  {
    my $self = shift;
    my %opts = @_;
    my $loc_link = $opts{loc_link};
    my $anno_obj = new CoGe::Genome::Accessory::Annotation(Type=>"anno");
    $anno_obj->Val_delimit("<BR/>");
    $anno_obj->Add_type(0);
    $anno_obj->String_end("<BR/>");
    my $start = $self->start;
    my $stop = $self->stop;
    my $chr = $self->chr;
    my $strand = $self->strand;
    my $dataset_id = $self->dataset->id;
    my $anno_type = new CoGe::Genome::Accessory::Annotation(Type=>"<span class=\"title4\">"."Name(s):"."</span>");
    $anno_type->Type_delimit("");
    $anno_type->Val_delimit(", ");
    foreach my $name ($self->names)
      {
	$anno_type->add_Annot("<a class=\"data\" href=\"FeatView.pl?accn=".$name."\">".$name."</a>");
      }
    
    $anno_obj->add_Annot($anno_type);
    foreach my $anno (sort {$b->type->name cmp $a->type->name} $self->annos)
      {
	my $type = $anno->type();
	my $group = $type->group();
	my $anno_name = $type->name;
	$anno_name = "<span class=\"title4\">". $anno_name."</span>" unless ref($group) =~ /group/i;
	
	my $anno_type = new CoGe::Genome::Accessory::Annotation(Type=>$anno_name);
	$anno_type->Val_delimit(", ");

	$anno_type->add_Annot("<span class=\"data\">".$anno->annotation."</span>");
	if (ref ($group) =~ /group/i)
	  {
	    my $anno_g = new CoGe::Genome::Accessory::Annotation(Type=>"<span class=\"title4\">".$group->name."</span>");
	    $anno_g->add_Annot($anno_type);
	    $anno_g->Type_delimit(": ");
	    $anno_g->Val_delimit(", ");
#	    $anno_g->Val_delimit(" ");
	    $anno_obj->add_Annot($anno_g);
	  }
	else
	  {
	    $anno_type->Type_delimit(": ");
	    $anno_obj->add_Annot($anno_type);
	  }
      }
    my $location = "Chr ".$chr." ";
    $location .= join (", ", map {$_->start."-".$_->stop} $self->locs);
    $location .="(".$strand.")";
    my $featid = $self->id;
    $location = qq{<a href="$loc_link?featid=$featid&start=$start&stop=$stop&chr=$chr&dsid=$dataset_id&strand=$strand" target=_new>}.$location."</a>" if $loc_link;
    $location = qq{<span class="data">$location</span>};
    $anno_obj->add_Annot(new CoGe::Genome::Accessory::Annotation(Type=>"<span class=\"title4\">Location</span>", Values=>[$location], Type_delimit=>": ", Val_delimit=>" "));
    return $anno_obj->to_String;
  }




################################################ subroutine header begin ##

=head2 genbank_location_string

 Usage     : my $genbank_loc = $feat->genbank_location_string
 Purpose   : generates a genbank location string for the feature in genomic coordinates or
           : based on a recalibration number that is user specified
           : e.g.: complement(join(10..100,200..400))
 Returns   : a string
 Argument  : hash:  recalibrate => number of positions to subtract from genomic location
 Throws    : none
 Comments  : 
           : 

See Also   : 

=cut

################################################## subroutine header end ##


sub genbank_location_string
  {
    my $self = shift;
    my %opts = @_;
    my $recal = $opts{recalibrate};
    my $string;
    my $count= 0;
    my $comp = 0;
    foreach my $loc (sort {$a->start <=> $b->start}  $self->locs())
      {
  #?
	# $comp = 1 if $loc->strand =~ "-";
	$comp = 1 if $loc->strand == "-1";
	$string .= "," if $count;
	$string .= $recal ? ($loc->start-$recal+1)."..".($loc->stop-$recal+1): $loc->start."..".$loc->stop;
	$count++;
      }
    $string = "join(".$string.")" if $count > 1;
    $string = "complement(".$string.")" if $comp;
    return $string;
  }


################################################ subroutine header begin ##

=head2 start

 Usage     : my $feat_start = $feat->start
 Purpose   : returns the start of the feature (does not take into account the strand on which
             the feature is located)
 Returns   : a string, number usually
 Argument  : none
 Throws    : 
 Comments  : this simply calles $feat->locs, sorts them based on their starting position, and
           : returns the smallest position

See Also   : 

=cut

################################################## subroutine header end ##


#sub start
#  {
#    my $self = shift;
#    return $self->{_start} if $self->{_start};
#    my @loc =  $self->locations({},
#				 {
#				  order_by=>'start asc',
#				 });
#    $self->{_start}=($loc[0]->start);
#    $self->{_stop}=($loc[-1]->stop);
#    $self->{_strand}=($loc[0]->strand);
#    $self->{_chromosome}=($loc[0]->chromosome);
#    return $self->{_start};
#  }

################################################ subroutine header begin ##

=head2 stop

 Usage     : my $feat_end = $feat->stop
 Purpose   : returns the end of the feature (does not take into account the strand on which
             the feature is located)
 Returns   : a string, number usually
 Argument  : none
 Throws    : 
 Comments  : this simply calles $feat->locs, sorts them based on their ending position, and
           : returns the largest position

See Also   : 

=cut

################################################## subroutine header end ##


#sub stop
#  {
#    my $self = shift;
#    return $self->{_stop} if $self->{_stop};
#    my @loc =  $self->locations({},
#				 {
#				  order_by=>'stop desc',
#				 });
#    $self->{_start}=($loc[-1]->start);
#    $self->{_stop}=($loc[0]->stop);
#    $self->{_strand}=($loc[0]->strand);
#    $self->{_chromosome}=($loc[0]->chromosome);
#    return $self->{_stop};
#  }
################################################ subroutine header begin ##

=head2 chromosome

 Usage     : my $chr = $feat->chromosome
 Purpose   : return the chromosome of the feature
 Returns   : a string
 Argument  : none
 Throws    : none
 Comments  : returns $self->locs->next->chr
           : 

See Also   : 

=cut

################################################## subroutine header end ##


#sub chromosome
#  {
#    my $self = shift;
#    return $self->{_chromosome} if $self->{_chromosome};
#    $self->start;
#    return $self->{_chromosome};
#  }
#
################################################# subroutine header begin ##

=head2 chr

 Usage     : my $chr = $feat->chr
 Purpose   : alias for $feat->chromosome

=cut

################################################## subroutine header end ##

sub chr
  {
    my $self = shift;
    return $self->chromosome;
  }

################################################ subroutine header begin ##

=head2 strand

 Usage     : my $strand = $feat->strand
 Purpose   : return the chromosome strand of the feature
 Returns   : a string (usally something like 1, -1, +, -, etc)
 Argument  : none
 Throws    : none
 Comments  : returns $self->locs->next->strand
           : 

See Also   : 

=cut

################################################## subroutine header end ##


#sub strand
#  {
#    my $self = shift;    
#    return $self->{_strand} if $self->{_strand};
#    $self->start;
#    return $self->{_strand};
#  }


################################################ subroutine header begin ##

=head2 version

 Usage     : my $version = $feat->version
 Purpose   : return the dataset version of the feature
 Returns   : an integer
 Argument  : none
 Throws    : none
 Comments  : returns $self->dataset->version
           : 

See Also   : 

=cut

################################################## subroutine header end ##


sub version
  {
    my $self = shift;
    return $self->dataset->version();
  }

################################################ subroutine header begin ##

=head2 genomic_sequence

 Usage     : my $genomic_seq = $feat->genomic_sequence
 Purpose   : gets the genomic seqence for a feature
 Returns   : a string
 Argument  : none
 Comments  : This method simply creates a CoGe::Genome object and calls:
             get_genomic_sequence_for_feature($self)
See Also   : CoGe::Genome

=cut

################################################## subroutine header end ##

sub genomic_sequence {
  my $self = shift;
  my $dataset = $self->dataset();
  my @sequences;
  my $lociter = $self->locations();
  while ( my $loc = $lociter->next() ) {
    my $fseq = $dataset->get_genome_sequence(
                                             chromosome=>$loc->chromosome(),
                                             skip_length_check=>1,
                                             start=>$loc->start,
                                             stop=>$loc->stop );
    if ( $loc->strand == -1 ) {
      push @sequences, $self->reverse_complement($fseq);
    } else {
      push @sequences, $fseq;
    }
  }
  return wantarray ? @sequences : join( "", @sequences );
}

sub genome_sequence
  {
   shift->genomic_sequence(@_);
  }

sub has_genomic_sequence
  {
    my $self = shift;
    return 1 if $self->dataset->has_genomic_sequence;
    return 0;
  }

################################################ subroutine header begin ##

=head2 blast_bit_score

 Usage     : my $bit_score = $feature->blast_bit_score();
 Purpose   : returns the blast bit score for the feature's self-self identical hit
 Returns   : an int -- the blast bit score
 Argument  : optional hash
             match    => the score for a nucleotide match. DEFAULT: 1
             mismatch => the score for a nucleotide mismatch.  DEFAULT: -3
 Throws    : 
 Comments  : 
           : 

See Also   : 

=cut

################################################## subroutine header end ##

sub blast_bit_score
  {
    my $self = shift;
    my %opts = @_;
    my $match = $opts{match} || 1;
    my $mismatch = $opts{mismatch} || -3;
    my $lambda = $self->_estimate_lambda(match=>$match, mismatch=>$mismatch);
    my $seq = $self->genomic_sequence();
    warn "No genomic sequence could be obtained for this feature object.  Can't calculate a blast bit score.\n" unless $seq;
    my $bs = sprintf("%.0f", $lambda*length($seq)*$match/log(2));
    return $bs;
  }

################################################ subroutine header begin ##

=head2 _estimate_lambda

 Usage     : my $lambda = $feature->_estimate_lambda
 Purpose   : estimates lambda for calculating blast bit scores.  Lambda is
             a matrix-specific constant for normalizing raw blast scores 
 Returns   : a number, lambda
 Argument  : optional hash
             match    => the score for a nucleotide match. DEFAULT: 1
             mismatch => the score for a nucleotide mismatch.  DEFAULT: -3
             precision=> the different between the high and low estimate 
                         of lambda before lambda is returned.  
                         DEFAULT: 0.001
 Throws    : a warning if there is a problem with the calcualted expected_score
             or the match score is less than 0;
 Comments  : Assumes an equal probability for each nucleotide.
           : this routine is based on example 4-1 from 
           : BLAST: An essential guide to the Basic Local Alignment Search Tool 
           : by Korf, Yandell, and Bedell published by O'Reilly press.

See Also   : 

=cut

################################################## subroutine header end ##


sub _estimate_lambda
  {
    #this routine is based on example 4-1 from BLAST: An essential guide to the Basic Local Alignment Search Tool by Korf, Yandell, and Bedell published by O'Reilly press.
    my $self = shift;
    my %opts = @_;
    my $match = $opts{match} || 1;
    my $mismatch = $opts{mismatch} || -3;
    my $precision = $opts{precision} || 0.001;
      
    use constant Pn => 0.25; #prob of any nucleotide
    my $expected_score = $match * 0.25 + $mismatch * 0.75; 
    if ($match <= 0 or $expected_score >= 0)
      {
	warn qq{
Problem with scores.  Match: $match (should be greater than 0).
             Expected score: $expected_score (should be less than 0).
};
	return 0;
      }
    # calculate lambda 
    my ($lambda, $high, $low) = (1, 2, 0); # initial estimates 
    while ($high - $low > $precision) 
      {         # precision 
	# calculate the sum of all normalized scores 
	my $sum = Pn * Pn * exp($lambda * $match) * 4 
	  + Pn * Pn * exp($lambda * $mismatch) * 12; 
	# refine guess at lambda 
	if ($sum > 1) 
	  { 
	    $high = $lambda;
	    $lambda = ($lambda + $low)/2; 
	  } 
	else 
	  {
	  $low = $lambda; 
	  $lambda = ($lambda + $high)/2; 
	}
      }
    # compute target frequency and H 
    my $targetID = Pn * Pn * exp($lambda * $match) * 4; 
    my $H = $lambda * $match    *     $targetID 
      + $lambda * $mismatch * (1 -$targetID); 
    # output 
#    print "expscore: $expected_score\n"; 
#    print "lambda:   $lambda nats (", $lambda/log(2), " bits)\n"; 
#    print "H:        $H nats (", $H/log(2), " bits)\n"; 
#    print "%ID:      ", $targetID * 100, "\n"; 

    return $lambda;
  }


sub reverse_complement
  {
    my $self = shift;
    my $seq = shift;# || $self->genomic_sequence;
    if (ref($seq) =~ /Feature/)
      {
	$seq = $self->genomic_sequence unless $seq; #self seq unless we have a seq
      }
    else #we were passed a sequence without invoking self
      {
	$seq = $self unless $seq;
      }
    my $rcseq = reverse($seq);
    $rcseq =~ tr/ATCGatcg/TAGCtagc/; 
    return $rcseq;
  }

sub reverse_comp {
  shift->reverse_complement(@_);
}


sub protein_sequence {
  my $self = shift;
  my $siter = $self->sequences();
  my @sequence_objects;
  while ( my $seq = $siter->next() ) {
    push @sequence_objects, $seq;
  }

  if (@sequence_objects == 1) {
    return $sequence_objects[0]->sequence_data();
  } elsif ( @sequence_objects > 1 )  {
    return \@sequence_objects;
  } else {
    return undef;
  }
}



1;


################################################ subroutine header begin ##

=head2 

 Usage     : 
 Purpose   : 
 Returns   : 
 Argument  : 
 Throws    : 
 Comments  : 
           : 

See Also   : 

=cut

################################################## subroutine header end ##

