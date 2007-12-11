package CoGe::Accessory::GenBank;

use strict;
use Data::Dumper;
use LWP::UserAgent;
use base qw(Class::Accessor);
use CoGe::Accessory::GenBank::Feature;
use CoGeX::Feature;


__PACKAGE__->mk_accessors qw(id locus accn seq_length moltype division date definition verison keywords data_source dataset organism sequence srcfile dir anntoation features start stop chromosome add_gene_models);

sub new
  {
    my $self = shift;
    $self = __PACKAGE__->SUPER::new(@_);
    $self->features([]) unless $self->features;
    return $self;
  }

sub get_genbank_from_ncbi
  {
    my $self = shift;
    my $id = shift;
    my $rev = shift;
    my $ua = new LWP::UserAgent;
    my $url = "http://www.ncbi.nlm.nih.gov/entrez/viewer.fcgi?db=nucleotide&qty=1&c_start=1&list_uids=$id&dopt=gb&send=Send&sendto=t&from=begin&to=end&extrafeatpresent=1&ef_MGC=16";
    my $search = $ua->get($url);
    unless ($search->is_success)
      {
	print STDERR "Trouble retrieving $id: ", $search->status_link,"\n";
      }

    return $self->parse_genbank($search->content, $rev);    

  }

sub parse_genbank
  {
    my $self = shift;
    my $gb = shift;
    my $rev = shift;
    my @gb = split /\n/, $gb;
    my %moltype_list = ( 'mRNA' => '1', 'RNA' => '1', 'tRNA' => '1',
			 'rRNA' => '1', 'DNA' => '1', 'scRNA' => '1',
			 'ds-mRNA' => '1', 'ss-DNA' => '1', 'snRNA' => '1'    );
    my($pos) = 0;
    my($locus,$length,$moltype,$division,$date) = "";
    my($accn,$definition,$version,$source,$organism) = "";
    my($keywords,$sequence) = "";
    my(@features) = ();

    unless ($gb[0] && $gb[0] =~ /^LOCUS/)
      {
	my $error =  "Doesn't appear to be a genbank file.  Returning. . .$gb\n";
	print STDERR $error;
	return (0,$error);
      }
    foreach my $line ( @gb ) 
      {
	if ( $line =~ /^LOCUS\s+(\S+)\s+(\S+)\s+bp\s+(.*)/ ) 
	  {
	    $locus = $1; $length = $2;
	    my @temp = split( /\s+/, $3 );
	    if ( @temp >= 3 ) 
	      {
		if ( exists($moltype_list{ $temp[0] } ) ) 
		  {
		    $moltype = $temp[0];
		  }
		else 
		  {
		    $moltype = "?";
		  }
		$division = $temp[2];
		$date =  $temp[3] ;
	      }
	    elsif ( @temp == 2 ) 
	      {
		if ( exists($moltype_list{ $temp[0] } ) ) 
		  {
		    $moltype = $temp[0];
		    $division = "?";
		  }
		else 
		  {
		    $moltype = "?";
		    $division = $temp[0];
		  }
		$date = $temp[2] ;
	      }
	    else 
	      {
		warn "  GBFile::Locus regexp - Problems with input entry\n";
		return 0;
	      }
	  }
	elsif ( $line =~ /^DEFINITION\s+(.*)/ ) {
	  $definition .= $definition ? " $1" : $1;
	  if ($gb[ $pos+1 ] !~ /^ACCESSION/ ) {
	    $gb[ $pos+1 ] =
	      "DEFINITION" . $gb[ $pos+1 ];
	  }
	}
	elsif ( $line =~ /^ACCESSION\s+(\S+)/ ) {
	  $accn = $1;
	}
	elsif ( $line =~ /^VERSION\s+(.*)/ ) {
	  my @temp = split(/\s+/,$1);
	  $version = join(" ", @temp);
	}
	elsif ( $line =~ /^KEYWORDS\s+(.*)\./ ) {
	  $keywords = $1;
	}
	elsif ( $line =~ /^SOURCE\s+(.*)/ ) {
	  $source = $1;
	  $source =~ s/\.$//; #get rid of trailing ".", if present
	}
	elsif ( $line =~ /^\s+ORGANISM\s+(.*)/ ) {
	  $organism = $1;
	}
	elsif ( $line =~ /^FEATURES/ ) 
	  {
	    my $x = $pos+1;
	    my $iso;
	    my %feature;
	    while ( $pos < $x ) 
	      {
		if  ($gb[ $x ] =~ /^ORIGIN/)
		  {
		    my %tmp = %feature;
		    push @features, \%tmp if keys %tmp;
		    last;
		  }
		if ($gb[$x] =~ /^\s+(\S+)\s\s+(.*)$/)
		  {
		    my %tmp = %feature;
		    push @features, \%tmp if keys %tmp;
		    %feature = ();
		    $feature{type} = $1;
		    $iso = 'location';
		    $feature{$iso} = $2;
		  }
		elsif($gb[$x] =~ /\/(.*?)=(.*)$/)
		  {
		    $iso = $1;
		    $feature{$iso} = $2;
		  }
	        else
		  {
		    my $line = $gb[$x];
		    $line =~ s/\s+//g;
		    $feature{$iso} .= $line if $iso;
		  }
	        $x++;
	      }
	    }
	elsif ( $line =~ /^ORIGIN/ ) {
	  # get the sequence
	  my $x = $pos+1;
	  $sequence = "";
	  while ( $pos < $x ) {
	    if ( $gb[ $x ] !~ /^\/\// ) {
	      foreach (split(/\s+/, $gb[ $x ] ))
		{
		  next unless $_;
		  next if /^\d+/;
		  $sequence .= $_;
		}
	      $x++;
	    }
	    else {
	      $x = 0;
	    }
	  }
	}
	elsif ( $line =~ /^\/\// ) {
	  # last line of buffer
	  $accn = $locus unless $accn;
	  $self->locus($locus);
	  $self->seq_length($length);
	  $self->moltype($moltype);
	  $self->division($division);
	  $self->verison($version);
	  $self->date($date);
	  $self->accn($accn);
	  $self->keywords($keywords);
	  $self->data_source($source);
	  $self->organism($organism);
	  $self->definition($definition);
	  $sequence = CoGeX::Feature->reverse_complement($sequence) if $rev;
	  $self->sequence( $sequence );
	  foreach my $feat (@features)
	    {
	      #  print STDERR "$locus\n",Dumper $feat;
	      my %quals;
	      my @names;
	      my $anno;
	      foreach (keys %$feat)
		{
		  my ($qual, $val) = ($_, $feat->{$_});
		  $val =~ s/"//g;
		  next if $qual eq "type" || $qual eq "location";
		  push @{$quals{$qual}},$val;
		  $anno .= "$qual: $val\n";
		  push @names, $val if $qual =~ /gene/;
		}
	      $quals{names} = \@names;
	      $feat->{location} = $self->reverse_genbank_location(loc=>$feat->{location}) if $rev;
	      my $strand = $feat->{location} =~ /complement/i ? -1 : 1;
	      $self->add_feature(
				 type=>$feat->{type},
				 location=>$feat->{location},
				 qualifiers=>\%quals,
				 annotation=>$anno,
				 strand=>$strand,
				);
	    }
	  $self->_check_for_gene_models if $self->add_gene_models;
	  return $self;
	}
	$pos++; # incr position counter
      }
    
  }

sub _check_for_gene_models
  {
    my $self = shift;
    my %data;
    foreach my $feat (sort {$a->start <=> $b->start} $self->get_features)
      {
	next unless $feat->type =~ /cds/i || $feat->type =~ /rna/i;
	push @{$data{$feat->strand}{lc($feat->type)}}, [$feat->start,$feat->stop];
      }
    foreach my $strand (keys %data)
      {
	foreach my $type (keys %{$data{$strand}})
	  {
	    next if $type =~ /gene/i;
	    foreach my $item (@{$data{$strand}{$type}})
	      {
		my ($start, $stop) = @$item;
		  {

		    my $location = $start."..".$stop;
		    $location  = "complement(".$location.")" if $strand =~ /-/;
		    $self->add_feature(
				       type=>"gene",
				       location=>$location,
				       qualifiers=>{names=>["autogenerated"]},
				       annotation=>"gene model automatically added because none was annotated",
				      );
		  }
	      }
	  }
      }
  }

## mask_exons: mask exons over entire sequence: 
## find all "features" and replace them with "N".  Uses
## "CDS" to determine mask, since many genomic GB entries lack mRNA
## entries.
sub mask_exons 
    {
      my $self = shift;
      my $seq = shift;
      foreach my $type (qw (CDS tRNA rRNA snRNA snoRNA))
	{
	  foreach my $block (@{$self->get_blocks(type=>$type)})
	    {
	      my $seglength = $block->[1] - $block->[0];
	      next if $block->[0] > length($seq);
	      my $maskstring = "X" x $seglength;
	      substr( $seq, $block->[0], $seglength ) = $maskstring;
	    }
	}
      return( $seq );
    }

#mask non coding sequences
sub mask_ncs 
  {
    my $self = shift;
    my $seq = shift;
    my $tmp_seq = "X"x length $seq;
    foreach my $type (qw (CDS tRNA rRNA snRNA snoRNA))
      {
	foreach my $block (@{$self->get_blocks(type=>$type)})
	  {
	    my $seglength = $block->[1] - $block->[0];
	    next if $block->[0] > length($seq);
	    my $coding_seq = substr($seq, $block->[0], $seglength);
	    substr( $tmp_seq, $block->[0], $seglength ) = $coding_seq;
	  }
      }
    return( $tmp_seq );
  }

# process_location: given a location line from a gb file, try to pull
# out the coordinates
sub process_location 
    {
      my $self = shift;
      my $locline = shift;
      my @locdata = ();
      $locline =~ s/\s+//g;                   #remove leading spaces
      $locline =~ s/complement\((.*)\)/$1/g;  #forget about "complement"
      $locline =~ s/join\((.*)\)/$1/g;        #forget about "join"
      $locline =~ s/[<>]//g;                  #forget about "<" or ">"
      if ( $locline =~ /,/ ) 
	{
	  my @pairs = split( /,/, $locline );
	  foreach my $pair ( @pairs ) 
	    {
	      my @temp = split( /\.\./, $pair );
	      if ( not defined $temp[1] ) 
		{
		  # need this for single nucleotide, which sometimes happens
		  $temp[1] = $temp[0];
		}
	      push(@locdata, [ $temp[0], $temp[1] ] );
	    }
	} 
      else 
	{
	  my @temp = split( /\.\./, $locline );
	  if ( not defined $temp[1] ) 
	    {
	      # need this for single nucleotide, which sometimes happens
	      $temp[1] = $temp[0];
	    }
	  push(@locdata, [ $temp[0], $temp[1] ] );
	}
      return(\@locdata);
    }

sub get_blocks 
  {
    my $self = shift;
    my %opts = @_;
    my $type = $opts{type};
    my $start = $opts{start};
    my $stop = $opts{stop};
    my @blocks;
    
    if ( $start && $stop ) 
      {
	foreach my $feature ( @{$self->features} ) 
	  {
	    if (!$type || $feature->type() eq $type ) 
	      {
		my $locdata = $self->process_location( $feature->location() );
		foreach my $subblock ( @{ $locdata } ) 
		  {
		    if (($subblock->[0] >= $start) and ($subblock->[1] <= $stop)) 
		      {
			push @blocks, [ 
				       $subblock->[0] - $start,
				       $subblock->[1]  - $start
				      ];
		      } 
		    elsif (($subblock->[1] >= $start) and ($subblock->[1] <= $stop)) 
		      {
			# begins up of region, but stops in region
			push @blocks, [
				       0, 
				       $subblock->[1] - $start
				      ];
		      } 
		    elsif (($subblock->[0] >= $start) and ($subblock->[0] <= $stop)) 
		      {
			# begins in the region, but stops outside region
			push @blocks, [
				       $subblock->[0] - $start, 
				       $stop - $start 
				      ];
		      }
		  }
	      }
	  }
      } 
    else 
      {
	# return all blocks
	foreach my $feature ( @{$self->features} ) 
	  {
	    if (!$type || $feature->type eq $type) 
	      {
		# this is a match
		my $locdata = $self->process_location( $feature->location() );
		foreach my $subblock ( @{ $locdata } ) 
		  {
		    push @blocks, [
				   $subblock->[0], $subblock->[1] 
				  ];
		  }
	      }
	  }
      }
    return( \@blocks );
  }

sub get_blocks_all
  {
    my $self = shift;
    return $self->get_blocks(@_);
  }

sub get_features
  {
    my $self = shift;
    my %opts = @_;
    my $start = $opts{start};
    my $stop = $opts{stop};
#    print "<pre>",Dumper ($self),"</pre>";
    my @features;
    if ( $start && $stop ) 
      {
	# only return blocks in this range
	foreach my $feature (@{$self->features})
	  {
	    my $locdata = $self->process_location( $feature->location() );
	    my @blocks;
	    foreach my $subblock ( @{ $locdata } ) 
	      {
		if (($subblock->[0] >= $start) and ($subblock->[1] <= $stop)) 
		  {
		    push @blocks, [ 
				   $subblock->[0] - $start,
				   $subblock->[1]  - $start
				  ];
		  } 
		elsif (($subblock->[1] > $start) and ($subblock->[1] < $stop)) 
		  {
		    # begins up of region, but ends in region
		    push @blocks,[
				  0, 
				  $subblock->[1] - $start
				 ];
		  } 
		elsif (($subblock->[0] >= $start) and ($subblock->[0] < $stop)) 
		  {
		    # begins in the region, but ends outside region
		    push @blocks, [
				   $subblock->[0] - $start, 
				   $stop - $start 
				  ];
		  }
		elsif ( ($subblock->[0] < $start) and ($subblock->[0] > $stop) )  #begins and ends outside of region
		  {
		    push @blocks, [
				   0, 
				   $stop-$start, 
				  ];
		  
		  }
	      }
	    $feature->blocks( \@blocks );
	    push @features, $feature if @blocks;
	  }
      } 
    else 
      {
	# return all blocks
	foreach my $feature (@{$self->features})
	  {
	    my $locdata = $self->process_location( $feature->location() );
	    my @blocks;
	    foreach my $subblock ( @{ $locdata } ) 
	      {
		push @blocks, [
			       $subblock->[0], 
			       $subblock->[1] 
			      ];
	      }
	    $feature->blocks(\@blocks);
	    push @features, $feature;
	  }
      }
    return( wantarray ? @features : \@features );
  }

sub add_feature 
    {
      my $self = shift;
      my %options = @_;
      my $strand = $options{location} =~ /complement/i ? "-1" : "1";
      my %info;
      $info{type}=$options{type} if $options{type};
      $info{location}=$options{location} if $options{location};
      $info{qualifiers}=$options{qualifiers} if $options{qualifiers};
      $info{annotation}=$options{annotation} if $options{annotation};
      $info{strand}=$options{strand} if $options{strand};
      my $feature = new CoGe::Accessory::GenBank::Feature(\%info);
      push @{$self->features}, $feature;
      return;
}


# write out a header for a fasta dump;  if any extra parameters are
# sent in they will get tacked on at the end.
sub get_headerfasta {
	my $self = shift;
	my $extra = "";
	if ( @_ ) {
		$extra = join(" ", @_);
	}
	my $header = ">";
	$header .= $self->accn() . " ";
	$header .= $self->definition() . " " if $self->definition;
	if ( $extra ne "" ) {
		$header .= $extra;
	}
	return( $header );
}

sub subsequence {
	my $self = shift;
	my $begin = shift;
	my $end = shift;

	# need to decrement by one, since substr works with arrays
	$begin--;
	$end--;

	my $seq = "";
	if ( @_ ) {
		$seq = shift;
	} else {
		$seq = $self->sequence();
	}
	my $newseq = substr( $seq, $begin, ( $end - $begin+1) );
	$self->sequence($newseq);
	$self->seq_length(length($newseq));
	return( $newseq	);
}

sub reverse_genbank_location
  {
    my $self = shift;
    my %opts = @_;
    my $loc = $opts{location} || $opts{loc};
    my $length = $opts{length} || length($self->sequence);
    my $complement = $loc=~/complement/ ? 0 : 1; #get the opposite
    $loc =~ s/complement//;
    $loc =~ s/join//;
    $loc =~ s/\(|\)|<|>//g;

    my $new_loc;
    my $count = 0;
    
    foreach my $set (split /,/, $loc)
      {
	my ($start, $stop) = split /\.\./, $set;
	$new_loc .= "," if $new_loc;
	$new_loc .= ($length-$stop)."..".($length-$start);
	$count++;
      }
    $new_loc = "join($new_loc)" if $count > 1;
    $new_loc = "complement($new_loc)" if $complement;
    return $new_loc;
  }

1;
