package CoGe::Accessory::GenBank;

use strict;
use Data::Dumper;
#use LWP::Simple;
use LWP::UserAgent;
use base qw(Class::Accessor);
use CoGe::Accessory::GenBank::Feature;
use CoGeX::Feature;
use Roman;

__PACKAGE__->mk_accessors qw(id locus accn seq_length moltype division date definition version gi keywords data_source dataset organism sequence srcfile dir anntoation features start stop chromosome add_gene_models _has_genes);



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
    my %opts = @_;
    my $id = $opts{accn};
    my $rev = $opts{rev};
    my $start = $opts{start};
    my $length = $opts{length};
    my $file = $opts{file};
    my $url = "http://www.ncbi.nlm.nih.gov/entrez/viewer.fcgi?db=nucleotide&qty=1&c_start=1&list_uids=$id&dopt=gb&send=Send&sendto=t&from=begin&to=end&extrafeatpresent=1&ef_MGC=16";
    my $ua = new LWP::UserAgent;
    if ($file)
      {
	unless (-r $file)
	  {
	    $ua->mirror($url, $file);
	  }
	return $self->parse_genbank_file(file=>$file, rev=>$rev, start=>$start, length=>$length);
      }

    my $search = $ua->get($url);
    unless ($search->is_success)
      {
	print STDERR "Trouble retrieving $id: ", $search->status_link,"\n";
      }

    return $self->parse_genbank(content=>$search->content, rev=>$rev, start=>$start);    

  }


sub parse_genbank_file
  {
    my $self = shift;
    my %opts = @_;
    my $file = $opts{file};
    my $rev = $opts{rev};
    my $start = $opts{start};
    my $length = $opts{length};
    $/="\n";

    open (IN, $file) || die "Can't open $file for reading";
    seek IN, 0, 0;
    my $working_line;
    my $feature_flag = 0;
    my $seq_flag = 0;
    while (my $line = <IN>)
      {
	chomp $line;
	next unless $line;
	$line =~ s/^\s+// if ($line =~ /ORGANISM/);
	if ($line=~/^FEATURES/)
	  {
	    $self->process_line(line=>$working_line, rev=>$rev, start=>$start, feature_flag=>$feature_flag);
	    $working_line = "";;
	    $feature_flag = 1;
	    next;
	  }
	if ($line=~/^ORIGIN/)
	  {
	    $self->process_line(line=>$working_line, rev=>$rev, start=>$start, feature_flag=>$feature_flag, length=>$length);
	    $working_line = "";
	    $seq_flag = 1;
	    $feature_flag = 0;
	    next;
	  }
	if ($feature_flag && $line =~ /^\s\s\s\s\s\w/)
	  {
	    $self->process_line(line=>$working_line, rev=>$rev, start=>$start, feature_flag=>$feature_flag, length=>$length);
	    $line =~ s/^\s+//;
	    $working_line = $line;
	    next;
	  }
	elsif ($feature_flag)
	  {
	    $line =~ s/^\s+//;
	    if ($line =~ /^\//)
	      {
		$working_line .= "\n".$line;
	      }
	    else
	      {
		$working_line .= " ".$line;
	      }
	    next;
	  }

	if ($seq_flag)
	  {
	    if ($line =~ /^\/\//)
	      {
		$working_line = CoGeX::Feature->reverse_complement($working_line) if $rev;
		$self->sequence($working_line);
		$seq_flag = 0;
		next;
	      }
	    $line =~ s/^\s+//;
	    my ($num, $seq) = split/\s+/,$line,2;
	    $seq =~ s/\s//g;
	    $working_line .= $seq;
	    next;
	  }
	if ($line =~ /^\S/)
	  {
	    $self->process_line(line=>$working_line, rev=>$rev, start=>$start, feature_flag=>$feature_flag);
	    $working_line = $line;
	  }
	else
	  {
	    $line =~ s/^\s+/ /;
	    $working_line .= $line;
	  }
      }
    close (IN);
    if ($seq_flag) #terminal // was missing from entry, process
      {
	$self->sequence($working_line);
      }
    unless ($self->chromosome && $self->chromosome eq "X") #dunno if it is a sex chromosome
      {
	$self->chromosome(arabic( $self->chromosome )) if $self->chromosome && isroman($self->chromosome);
      }
    $self->_check_for_gene_models() if $self->add_gene_models;
    return $self;
  }


sub process_line
  {
    my $self = shift;
    my %opts = @_;
    my $line = $opts{line};
    my $rev = $opts{rev};
    my $start = $opts{start};
    my $length = $opts{length};
    my $feature_flag = $opts{feature_flag};
    return unless $line;
    if ($line =~ /^LOCUS/)
      {
	my @tmp = split /\s+/, $line;
	$self->locus($tmp[1]);
	$self->seq_length($tmp[2]);
	$self->date($tmp[-1]);
	$tmp[4]="?" unless $tmp[4];
	$tmp[5]="?" unless $tmp[5];
	$self->moltype($tmp[4]);
	$self->division($tmp[5]);
      }
    elsif ( $line =~ /^DEFINITION\s+(.*)/ ) {
      $self->definition($1);
      if ($self->definition() =~ /(linkage group \w+),?/i || $self->definition =~ /chromosome (\w+),?/)
	{
	  $self->chromosome($1);
	}
    }
    elsif ( $line =~ /^ACCESSION\s+(\S+)/ ) {
      $self->accn($1);
    }
    elsif ( $line =~ /^VERSION\s+.*?\.(\d+)/ ) {
      my @temp = split(/\s+/,$1);
      $self->version(join(" ", @temp));
      if ($line =~ /GI:(\d+)/)
	{
	  $self->gi($1);
	}
    }
    elsif ( $line =~ /^KEYWORDS\s+(.*)\./ ) {
      $self->keywords($1);
    }
    elsif ( $line =~ /^SOURCE\s+(.*)/ ) {
      my $tmp = $1;
      $tmp =~ s/\.$//; #get rid of trailing ".", if present
      $self->data_source($tmp);
    }
    elsif ( $line =~ /^ORGANISM\s+(.*)/ ) {
      $self->organism($1)
    }
    elsif ( $feature_flag ) 
      {
	my $iso;
	my %feature;
	foreach my $item (split/\n/,$line)
	  {
	    if ($item !~ /^\// && $item =~ /^(\S+)\s\s+(.*)$/)
	      {
		$feature{type} = $1;
		$feature{location} = $2;
		$feature{location}=~ s/\s//g;
	      }
	    elsif($item =~ /\//)#(.*?)=?(.*)$/)
	      {
		$item =~ s/\///;
		my $val;
		($iso,$val) = split/=/,$item,2;
		if ($iso =~ /^pseudo/i)
		  {
		    $feature{type} = "pseudogene";
		    next;
		  }
		$val =1 unless $val;
		$val =~ s/"//g;
		$feature{$iso} = $val;
		$self->chromosome($feature{$iso}) if $iso eq "chromosome";
	      }
	    else
	      {
		$item =~ s/\s+/ /g;
		$feature{$iso} .= $line if $iso;
	      }
	  }
	my %quals;
	my %names;
	my $anno;
	foreach (keys %feature)
	  {
	    my ($qual, $val) = ($_, $feature{$_});
	    $val =~ s/"//g;
	    next if $qual eq "type" || $qual eq "location";
	    push @{$quals{$qual}},$val;
	    $anno .= "$qual: $val\n";
	    $names{$val}=1 if $qual =~ /gene/;
	    $names{$val}=1 if $qual =~ /synonym/;
	    $names{$val}=1 if $qual =~ /locus/;
	    $names{$val}=1 if $qual =~ /transcript_id/;
	    $names{$val}=1 if $qual =~ /protein_id/;
	  }
	$quals{names} = [sort keys %names];
	$feature{location} = $self->reverse_genbank_location(loc=>$feature{location}) if $rev;
	my $strand = $feature{location} =~ /complement/i ? -1 : 1;
	$start = $self->seq_length -$start+1 if $rev;
	$start -= $length if $rev && $length;
	$self->add_feature(
			   type=>$feature{type},
			   location=>$feature{location},
			   qualifiers=>\%quals,
			   annotation=>$anno,
			   strand=>$strand,
			   start=>$start,
			   length=>$length,
			  );
      }
  }



sub parse_genbank
  {
    my $self = shift;
    my %opts = @_;
    my $gb = $opts{content};
    my $rev = $opts{rev};
    my $start = $opts{start};
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
		unless ($gb[$x])
		  {
		    $x++;
		    next;
		  }
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
		    $line =~ s/\s+/ /g;
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
	  $self->version($version);
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
				 start=>$start,
				);
	    }
	  $self->_check_for_gene_models() if $self->add_gene_models;
	  return $self;
	}
	$pos++; # incr position counter
      }
    
  }

sub _check_for_gene_models
  {
    my $self = shift;
    return if $self->_has_genes;
    my %opts = @_;
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
      my %type = (CDS=>1, tRNA=>1, rRNA=>1, snRNA=>1, snoRNA=>1);
      foreach my $feat (@{$self->features})
	{
	  next unless $type{$feat->type};
	  foreach my $block (@{$feat->blocks})
	    {
	      next if $block->[0] > length($seq);
	      next if $block->[1] < 1;
	      my ($start, $stop) = @$block;
	      $start = 1 if $start < 1;
	      my $seglength = $stop - $start+1;
	      my $maskstring = "X" x $seglength;
	      substr( $seq, $start-1, $seglength ) = $maskstring;
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
    my %type = (CDS=>1, tRNA=>1, rRNA=>1, snRNA=>1, snoRNA=>1);
    foreach my $feat (@{$self->features})
      {
	next unless $type{$feat->type};
	foreach my $block (@{$feat->blocks})
	  {
	    next if $block->[0] > length($seq);
	    next if $block->[1] < 1;
	    my ($start, $stop) = @$block;
	    $start = 1 if $start < 1;
	    my $seglength = $stop - $start+1;
	    my $coding_seq = substr($seq, $start-1, $seglength);
	    substr( $tmp_seq, $start-1, $seglength ) = $coding_seq;
	  }
      }
    return( $tmp_seq );
  }

# process_location: given a location line from a gb file, try to pull
# out the coordinates
sub process_location 
    {
      my $self = shift;
      my %opts = @_;
      my $locline = $opts{loc};
      my $start = $opts{start};
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
	      $temp[0] = $temp[0]-$start+1 if $start;
	      $temp[1] = $temp[1]-$start+1 if $start;
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
	  $temp[0] = $temp[0]-$start+1 if $start;
	  $temp[1] = $temp[1]-$start+1 if $start;
	  push(@locdata, [ $temp[0], $temp[1] ] );
	}
      return(\@locdata);
    }

sub get_features
  {
    my $self = shift;
    return( wantarray ? @{$self->features} : $self->features );
  }

sub add_feature 
  {
    my $self = shift;
    my %options = @_;
    my $type=$options{type} if $options{type};
    my $location=$options{location} if $options{location};
    $self->_has_genes(1) if $type =~ /^gene$/i;
    my $qualifiers = $options{qualifiers};
    $qualifiers = {} unless $qualifiers;
    $qualifiers=$qualifiers;
    my$annotation=$options{annotation} if $options{annotation};
    $annotation .= "location: ".$location if $location;
    my$strand=$options{strand} if $options{strand};
    my $start = $options{start};
    my $length = $options{length};
    unless ($strand)
      {
	$strand = $location=~/complement/i ? -1 : 1;
      }
    my $feature = new CoGe::Accessory::GenBank::Feature({
							 type=>$type,
							 location=>$location,
							 qualifiers=>$qualifiers,
							 annotation=>$annotation,
							 strand=>$strand,
							});
    $feature->blocks($self->process_location( loc=>$feature->location, start=>$start)) if $feature->location;
    return unless $feature->stop && $feature->start;
    return if $feature->stop < 0;
    return if $start && $start > 1 && $length && $feature->start > $length;
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
    my $length = $opts{length} || $self->seq_length || length($self->sequence);
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
	$new_loc .= ($length-$stop+1)."..".($length-$start+1);
	$count++;
      }
    $new_loc = "join($new_loc)" if $count > 1;
    $new_loc = "complement($new_loc)" if $complement;
    return $new_loc;
  }


sub accession
  {
    shift->accn(@_);
  }

1;
