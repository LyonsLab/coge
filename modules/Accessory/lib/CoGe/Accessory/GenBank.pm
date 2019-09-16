package CoGe::Accessory::GenBank;

use strict;
use Data::Dumper;

#use LWP::Simple;
use LWP::UserAgent;
use base qw(Class::Accessor);
use CoGe::Accessory::GenBank::Feature;
use CoGeX::Result::Feature;
use Roman;

__PACKAGE__->mk_accessors(
	'id',              'locus',
	'accn',            'seq_length',
	'moltype',         'division',
	'date',            'definition',
	'version',         'gi',
	'keywords',        'data_source',
	'dataset',         'organism',
	'sequence',        'srcfile',
	'dir',             'anntoation',
	'features',        'start',
	'stop',            'chromosome',
	'add_gene_models', '_has_genes',
	'wgs',             'wgs_scafld',
	'wgs_data',        'strain',
	'substrain',       'genomic_sequence_type_id',
	'debug',           'no_wgs',
	'requested_id',    'other_stuff',
	'ncbi_link',       'max_entries',
	'entries_count',   'logfile'
);

sub new {
	my $self = shift;
	$self = __PACKAGE__->SUPER::new(@_);
	$self->features( [] ) unless $self->features;
	$self->wgs_data( [] ) unless $self->wgs_data;
	return $self;
}

sub get_genbank_from_ncbi {
	my $self = shift;
	my %opts = @_;
	my $id   = $opts{accn};
	my $rev    = $opts{rev};
	my $start  = $opts{start};
	my $length = $opts{length};
	my $file   = $opts{file};
	my $reload = $opts{reload};    #option to force reloading of genbank file
	my $retmax = $opts{retmax};    #number of sequences to return
	my $complexity = $opts{complexity};    #what to return from genbank
	$complexity = 0                   unless defined $complexity;

	print STDERR "Fetching $id\n" if $self->debug();
	print {$self->logfile} "log: Fetching $id\n" if $self->logfile;

	$id         = $self->requested_id unless $id;
	$file       = $self->srcfile      unless $file;
	$self->requested_id($id);
	$self->srcfile($file);

#    my $url = "http://www.ncbi.nlm.nih.gov/entrez/viewer.fcgi?db=nucleotide&qty=1&c_start=1&dopt=gbwithparts&send=Send&sendto=t&from=begin&to=end&extrafeatpresent=1&ef_CDD=8&ef_MGC=16&ef_HPRD=32&ef_STS=64&ef_tRNA=128&ef_microRNA=256&ef_Exon=512&list_uids=";
	my $url = "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&rettype=gbwithparts&retmode=text";
	$url .= "&complexity=$complexity";
	$url .= "&retmax=$retmax" if $retmax;
	$url .= "&id=";
	print STDERR $url . $id, "\n";
	my $ua = new LWP::UserAgent;
	$self->ncbi_link("http://www.ncbi.nlm.nih.gov/nuccore/$id");

	unless ($file) {
		print STDERR "must specify a file where the gb file exists or where it will be written using file=>$file!\n";
		return;
	}

	my @files = ( [ $file, $id ] );

#    if ($id =~ /^NZ_\w\w\w\w00000000/i)
#      {
#we have shotgun sequence.  There are inconsistencies with how these can be stored.  Often there accession NZ_XXXX00000000 and XXXX00000000 are both available.  Sometimes one of these accessions has a more useful sequence by either having annotations, or having better contig assembly.  We'll try to figure that out!
#	my $tmp = $file;
#	$tmp =~ s/NZ_//i;
#	my $tmpid = $id;
#	$tmpid =~ s/NZ_//i;
#	push @files, [$tmp, $tmpid];
#      }
#    print Dumper \@files;
	foreach my $tmp (@files) {
		my ( $fileloc, $accn ) = @$tmp;
		unless ( -r $fileloc && !$reload ) {
			my $response = $ua->get( $url . "$accn" );

			my $entry = $response->content;
			$entry = "" unless $response->is_success();
			my $tries = 10;
			my $count = 0;
		  try_loop: while ( !$entry
				|| $entry =~ /(temporarily unavailable)/i
				|| $entry =~ /error/i )
			{
				last try_loop if $entry =~ /^LOCUS/;
				last if $count >= $tries;
				print STDERR "Warning: $url".$accn." failed\nFetching from NCBI yielding '$entry'\nRetrying\n";
				sleep 10;
				$response = $ua->get( $url . "$accn" );
				$entry    = $response->content;
				$count++;
			}
			open( OUT, ">" . $fileloc )
			  || die "Can't open $fileloc for writing: $!";
			print OUT $entry;
			close OUT;
		}
	}
	if ( @files == 1 ) {
		return $self->parse_genbank_file(
			file   => $file,
			rev    => $rev,
			start  => $start,
			length => $length
		);
	}
	my @gbs;
	foreach my $tmp (@files) {
		my ( $fileloc, $accn ) = @$tmp;
		my $gb = $self->new();
		$gb->debug( $self->debug ) if $self->debug;
		$gb->srcfile($fileloc);
		$gb->parse_genbank_file(
			file   => $fileloc,
			rev    => $rev,
			start  => $start,
			length => $length
		);
		push @gbs, $gb;
	}

	#try to find the one with the most features
	my $best_gb   = $gbs[0];
	my $max_count = 0;
	foreach my $gb (@gbs) {
		my $feat_count = 0;
		foreach my $item ( @{ $gb->wgs_data } ) {
			$feat_count += scalar @{ $item->features } - 1;
			$feat_count
			  --; #subtract one for each dataset -- this gives an advantage to that dataset with fewer contigs
		}
		if ( $feat_count > $max_count ) {
			$best_gb   = $gb;
			$max_count = $feat_count;
		}
	}
	$self->parse_genbank_file(
		file   => $best_gb->srcfile,
		rev    => $rev,
		start  => $start,
		length => $length
	);            #reparse and go!
	return $self;
}

sub parse_genbank_file {
	my $self = shift;
	my %opts = @_;
	my $file = $opts{file};
	my $rev    = $opts{rev};
	my $start  = $opts{start};
	my $length = $opts{length};
	$/ = "\n";

	print STDERR "parsing $file\n" if $self->debug;

	open( IN, $file ) || die "Can't open $file for reading";
	seek IN, 0, 0;
	my $working_line;
	my $feature_flag = 0;
	my $seq_flag     = 0;
	while ( my $line = <IN> ) {
		chomp $line;
		next unless $line;
		if ( $line =~ /ORGANISM/ ) {
			$line =~ s/^\s+//;
			$line .= "::";
		}
		if ( $line =~ /^FEATURES/ ) {
			$self->process_line(
				line         => $working_line,
				rev          => $rev,
				start        => $start,
				feature_flag => $feature_flag
			);
			$working_line = "";
			$feature_flag = 1;
			next;
		}
		if ( $line =~ /^ORIGIN/ ) {
			$self->process_line(
				line         => $working_line,
				rev          => $rev,
				start        => $start,
				feature_flag => $feature_flag,
				length       => $length
			);
			$working_line = "";
			$seq_flag     = 1;
			$feature_flag = 0;
			next;
		}
		if ( $line =~ /^WGS/ ) {
			$self->process_line(
				line         => $working_line,
				rev          => $rev,
				start        => $start,
				feature_flag => $feature_flag,
				length       => $length
			);
			$feature_flag = 0;
		}
		if ( $feature_flag && $line =~ /^\s\s\s\s\s\S/ ) {
			$self->process_line(
				line         => $working_line,
				rev          => $rev,
				start        => $start,
				feature_flag => $feature_flag,
				length       => $length
			);
			$line =~ s/^\s+//;
			$working_line = $line;
			next;
		}
		elsif ($feature_flag) {
			$line =~ s/^\s+//;
			if ( $line =~ /^\// ) {
				$working_line .= "\n" . $line;
			}
			else {
				$working_line .= " " . $line;
			}
			next;
		}

		if ($seq_flag) {
			if ( $line =~ /^\/\// ) {
				$working_line =
				  CoGeX::Result::Feature->reverse_complement($working_line)
				  if $rev;
				$self->sequence($working_line);
				$seq_flag = 0;
				next;
			}
			$line =~ s/^\s+//;
			my ( $num, $seq ) = split( /\s+/, $line, 2 );
			$seq =~ s/\s//g;
			$working_line .= $seq;
			next;
		}
		if ( $line =~ /^\S/ ) {
			$self->process_line(
				line         => $working_line,
				rev          => $rev,
				start        => $start,
				feature_flag => $feature_flag
			);
			$working_line = $line;
		}
		else {
			$line =~ s/^\s+/ /;
			$working_line .= $line;
		}
	}
	close(IN);
	if ($seq_flag)    #terminal // was missing from entry, process
	{
		$self->sequence($working_line);
	}

#    unless ($self->chromosome && ($self->chromosome eq "X" || $self->chromosome=~/\s/)) #dunno if it is a sex chromosome
#      {
#	print "---\n";
#	print STDERR $self->chromosome,"\n";
#	$self->chromosome(arabic( $self->chromosome )) if $self->chromosome && isroman($self->chromosome);
#	print STDERR $self->chromosome,"!!!\n";
#      }
	$self->_check_for_gene_models() if $self->add_gene_models;
	$self->process_wgs              if ( $self->wgs );
	return $self;
}

sub process_line {
	my $self         = shift;
	my %opts         = @_;
	my $line         = $opts{line};
	my $rev          = $opts{rev};
	my $start        = $opts{start};
	my $length       = $opts{length};
	my $feature_flag = $opts{feature_flag};
	return unless $line;

	if ( $line =~ /^LOCUS/ ) {
		my @tmp = split( /\s+/, $line );
		$self->locus( $tmp[1] );
		$self->seq_length( $tmp[2] );
		$self->date( $tmp[-1] );
		$tmp[4] = "?" unless $tmp[4];
		$tmp[5] = "?" unless $tmp[5];
		$self->moltype( $tmp[4] );
		$self->division( $tmp[5] );
	}
	elsif ( $line =~ /^DEFINITION\s+(.*)/ ) {
		$self->definition($1);
		if (   $self->definition() =~ /(linkage group .+),?/i
			|| $self->definition =~ /chromosome ([^,]*),?/
			|| $self->definition =~ /(plasmid[^,]*),?/i
			|| $self->definition =~ /(extrachromosomal[^,]*),?/i
			|| $self->definition =~ /(scaffold[^,]*),?/i
			|| $self->definition =~ /(segment[^,]*),?/i
			|| $self->definition =~ /\s(part [^,]*),?/i
			|| $self->definition =~ /\s(clone [^,]*),?/i
			|| $self->definition =~ /(mitochondrion)/i
			|| $self->definition =~ /(chloroplast)/i
			|| $self->definition =~ /(apicoplast)/i
			|| $self->definition =~ /(plastid)/i )
		{
			my $tmp = $1;
			$tmp =~ s/clone/contig/;
			$tmp =~ s/\.$//;
			$tmp =~ s/complete//;
			$tmp =~ s/sequence//;
			$tmp =~ s/genomic scaffold//;
			$tmp =~ s/\s+$//;
			$tmp =~ s/^\s+//;
			$self->chromosome($tmp);
		}
	}
	elsif ( $line =~ /^ACCESSION\s+(\S+)/ ) {
		$self->accn($1);
	}
	elsif ( $line =~ /^VERSION\s+.*?\.(\d+)/ ) {
		my @temp = split( /\s+/, $1 );
		$self->version( join( " ", @temp ) );
		if ( $line =~ /GI:(\d+)/ ) {
			$self->gi($1);
		}
	}
	elsif ( $line =~ /^KEYWORDS\s+(.*)\./ ) {
		$self->keywords($1);
	}
	elsif ( $line =~ /^SOURCE\s+(.*)/ ) {
		my $tmp = $1;
		$tmp =~ s/\.$//;    #get rid of trailing ".", if present
		$self->data_source($tmp);
	}
	elsif ( $line =~ /^ORGANISM\s+(.*)/ ) {
		$self->organism($1);
	}
	elsif ( $line =~ /^WGS\s+(.*)/ ) {
		$self->wgs($1);
	}
	elsif ( $line =~ /^WGS_SCAFLD\s+(.*)/ ) {
		if ( $self->wgs_scafld && $self->wgs_scafld ne $1 ) {
			return;
			my $prev = $self->wgs_scafld;
			$prev .= "," . $1;
			$self->wgs_scafld($prev);
		}
		else {
			$self->wgs_scafld($1);
		}
	}
	elsif ($feature_flag) {
		my $iso;
		my %feature;
		foreach my $item ( split(/\n/, $line) ) {
			if ( $item !~ /^\// && $item =~ /^(\S+)\s\s+(.*)$/ ) {
				$feature{type}     = $1;
				$feature{location} = $2;
				$feature{location} =~ s/\s//g;
			}
			elsif ( $item =~ /\// )    #(.*?)=?(.*)$/)
			{
				$item =~ s/\///;
				my $val;
				( $iso, $val ) = split( /=/, $item, 2 );
				if ( $iso =~ /^pseudo/i ) {
					$feature{type} = "pseudogene";
					next;
				}
				$val = 1 unless $val;
				$val =~ s/"//g;
				$feature{$iso} = $val;

				#funky ways for getting chromosome info
				if (
					   $feature{type}
					&& $feature{type} eq "source"
					&& (   $iso eq "segment"
						|| $iso eq "chromosome"
						|| $iso eq "plasmid" )
				  )
				{
					my $chr = $val;
					$chr = $iso . " " . $val
					  unless $iso eq "chromosome"
					  || $chr =~ /\s/
					  || $chr =~ /$iso/i;
					$self->chromosome($chr) unless $self->chromosome;
					$self->chromosome($chr)
					  if length($chr) > length( $self->chromosome );
				}
				if (   $feature{type}
					&& $feature{type} eq "source"
					&& $iso           eq "strain" )
				{
					$self->strain($val);
				}
				if (   $feature{type}
					&& $feature{type} eq "source"
					&& $iso           eq "substrain" )
				{
					$self->substrain($val);
				}
				if (   $feature{type}
					&& $feature{type} eq "source"
					&& $iso           eq "sub_strain" )
				{
					$self->substrain($val);
				}
			}
			else {
				$item =~ s/\s+/ /g;
				$feature{$iso} .= $line if $iso;
			}
		}
		my %quals;
		my %names;
		my $anno;
		foreach ( keys %feature ) {
			my ( $qual, $val ) = ( $_, $feature{$_} );
			$val =~ s/"//g;
			next if $qual eq "type" || $qual eq "location";
			push @{ $quals{$qual} }, $val;
			$anno .= "$qual: $val\n";
			$names{$val} = 1 if $qual =~ /gene/;
			$names{$val} = 1 if $qual =~ /synonym/;
			$names{$val} = 2 if $qual =~ /locus/;
			$names{$val} = 1 if $qual =~ /transcript_id/;
			$names{$val} = 1 if $qual =~ /protein_id/;
		}
		$quals{names} =
		  [ sort { $names{$b} <=> $names{$a} || $a cmp $b } keys %names ];
		$feature{location} =
		  $self->reverse_genbank_location( loc => $feature{location} )
		  if $rev;
		my $strand =
		  ( $feature{location} && $feature{location} =~ /complement/i )
		  ? -1
		  : 1;
		$start = $self->seq_length - $start + 1 if $rev;
		$start -= $length if $rev && $length;
		$self->add_feature(
			type       => $feature{type},
			location   => $feature{location},
			qualifiers => \%quals,
			annotation => $anno,
			strand     => $strand,
			start      => $start,
			length     => $length,
		);
	}
}

sub parse_genbank {
	my $self         = shift;
	my %opts         = @_;
	my $gb           = $opts{content};
	my $rev          = $opts{rev};
	my $start        = $opts{start};
	my @gb           = split(/\n/, $gb);
	my %moltype_list = (
		'mRNA'    => '1',
		'RNA'     => '1',
		'tRNA'    => '1',
		'rRNA'    => '1',
		'DNA'     => '1',
		'scRNA'   => '1',
		'ds-mRNA' => '1',
		'ss-DNA'  => '1',
		'snRNA'   => '1'
	);
	my ($pos) = 0;
	my ( $locus, $length,     $moltype, $division, $date )     = "";
	my ( $accn,  $definition, $version, $source,   $organism ) = "";
	my ( $keywords, $sequence ) = "";
	my (@features) = ();

	unless ( $gb[0] && $gb[0] =~ /^LOCUS/ ) {
		my $error = "Doesn't appear to be a genbank file.  Returning. . .$gb\n";
		print STDERR $error;
		return ( 0, $error );
	}
	foreach my $line (@gb) {
		if ( $line =~ /^LOCUS\s+(\S+)\s+(\S+)\s+bp\s+(.*)/ ) {
			$locus  = $1;
			$length = $2;
			my @temp = split( /\s+/, $3 );
			if ( @temp >= 3 ) {
				if ( exists( $moltype_list{ $temp[0] } ) ) {
					$moltype = $temp[0];
				}
				else {
					$moltype = "?";
				}
				$division = $temp[2];
				$date     = $temp[3];
			}
			elsif ( @temp == 2 ) {
				if ( exists( $moltype_list{ $temp[0] } ) ) {
					$moltype  = $temp[0];
					$division = "?";
				}
				else {
					$moltype  = "?";
					$division = $temp[0];
				}
				$date = $temp[2];
			}
			else {
				if ( exists ($moltype_list { $temp[0] } ) ) {
					$moltype = $temp[0];
					$division = "?";
					$date = "?";
				}
				else {
					$moltype = "?";
					$division = "?";
					$date = "?";
				}
				#warn "  GBFile::Locus regexp - Problems with input entry\n";
				#return 0;
			}
		}
		elsif ( $line =~ /^DEFINITION\s+(.*)/ ) {
			$definition .= $definition ? " $1" : $1;
			if ( $gb[ $pos + 1 ] !~ /^ACCESSION/ ) {
				$gb[ $pos + 1 ] = "DEFINITION" . $gb[ $pos + 1 ];
			}
		}
		elsif ( $line =~ /^ACCESSION\s+(\S+)/ ) {
			$accn = $1;
		}
		elsif ( $line =~ /^VERSION\s+(.*)/ ) {
			my @temp = split( /\s+/, $1 );
			$version = join( " ", @temp );
		}
		elsif ( $line =~ /^KEYWORDS\s+(.*)\./ ) {
			$keywords = $1;
		}
		elsif ( $line =~ /^SOURCE\s+(.*)/ ) {
			$source = $1;
			$source =~ s/\.$//;    #get rid of trailing ".", if present
		}
		elsif ( $line =~ /^\s+ORGANISM\s+(.*)/ ) {
			$organism = $1;
		}
		elsif ( $line =~ /^FEATURES/ ) {
			my $x = $pos + 1;
			my $iso;
			my %feature;
			while ( $pos < $x ) {
				unless ( $gb[$x] ) {
					$x++;
					next;
				}
				if ( $gb[$x] =~ /^ORIGIN/ ) {
					my %tmp = %feature;
					push @features, \%tmp if keys %tmp;
					last;
				}
				if ( $gb[$x] =~ /^\s+(\S+)\s\s+(.*)$/ ) {
					my %tmp = %feature;
					push @features, \%tmp if keys %tmp;
					%feature       = ();
					$feature{type} = $1;
					$iso           = 'location';
					$feature{$iso} = $2;
				}
				elsif ( $gb[$x] =~ /\/(.*?)=(.*)$/ ) {
					$iso = $1;
					$feature{$iso} = $2;
				}
				else {
					my $line = $gb[$x];
					$line =~ s/\s+/ /g;
					$feature{$iso} .= $line if $iso;
				}
				$x++;
			}
		}
		elsif ( $line =~ /^ORIGIN/ ) {

			# get the sequence
			my $x = $pos + 1;
			$sequence = "";
			while ( $pos < $x ) {
				if ( $gb[$x] !~ /^\/\// ) {
					foreach ( split( /\s+/, $gb[$x] ) ) {
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
			$sequence = CoGeX::Result::Feature->reverse_complement($sequence)
			  if $rev;
			$self->sequence($sequence);

			foreach my $feat (@features) {
				my %quals;
				my @names;
				my $anno;
				foreach ( keys %$feat ) {
					my ( $qual, $val ) = ( $_, $feat->{$_} );
					$val =~ s/"//g;
					next if $qual eq "type" || $qual eq "location";
					push @{ $quals{$qual} }, $val;
					$anno .= "$qual: $val\n";
					push @names, $val if $qual =~ /gene/;
				}
				$quals{names} = \@names;
				$feat->{location} =
				  $self->reverse_genbank_location( loc => $feat->{location} )
				  if $rev;
				my $strand = $feat->{location} =~ /complement/i ? -1 : 1;
				$self->add_feature(
					type       => $feat->{type},
					location   => $feat->{location},
					qualifiers => \%quals,
					annotation => $anno,
					strand     => $strand,
					start      => $start,
				);
			}
			$self->_check_for_gene_models() if $self->add_gene_models;
			return $self;
		}
		$pos++;    # incr position counter
	}

}

sub _check_for_gene_models {
	my $self = shift;
	return if $self->_has_genes;
	my %opts = @_;
	my %data;
	foreach my $feat ( sort { $a->start <=> $b->start } $self->get_features ) {
		next unless $feat->type =~ /cds/i || $feat->type =~ /rna/i;
		push @{ $data{ $feat->strand }{ lc( $feat->type ) } },
		  [ $feat->start, $feat->stop ];
	}
	foreach my $strand ( keys %data ) {
		foreach my $type ( keys %{ $data{$strand} } ) {
			next if $type =~ /gene/i;
			foreach my $item ( @{ $data{$strand}{$type} } ) {
				my ( $start, $stop ) = @$item;
				{

					my $location = $start . ".." . $stop;
					$location = "complement(" . $location . ")"
					  if $strand =~ /-/;
					$self->add_feature(
						type       => "gene",
						location   => $location,
						qualifiers => { names => ["autogenerated"] },
						annotation => "gene model automatically added because none was annotated",
					);
				}
			}
		}
	}
}

sub mask_cds {
	my $self = shift;
	my $seq  = shift;
	my %type = ( CDS => 1 );
	foreach my $feat ( @{ $self->features } ) {
		next unless $type{ $feat->type };
		foreach my $block ( @{ $feat->blocks } ) {
			next if $block->[0] > length($seq);
			next if $block->[1] < 1;
			my ( $start, $stop ) = @$block;
			$start = 1 if $start < 1;
			my $seglength  = $stop - $start + 1;
			my $maskstring = "X" x $seglength;
			substr( $seq, $start - 1, $seglength ) = $maskstring;
		}
	}
	return ($seq);
}

sub mask_exons {
	my $self = shift;
	warn "please use sub mask_cds\n";
	$self->mask_cds(@_);
}

#mask rna sequences
sub mask_rna {
	my $self = shift;
	my $seq  = shift;
	foreach my $feat ( @{ $self->features } ) {
		next unless $feat->type =~ /rna/i;
		foreach my $block ( @{ $feat->blocks } ) {
			next if $block->[0] > length($seq);
			next if $block->[1] < 1;
			my ( $start, $stop ) = @$block;
			$start = 1 if $start < 1;
			my $seglength  = $stop - $start + 1;
			my $maskstring = "X" x $seglength;
			substr( $seq, $start - 1, $seglength ) = $maskstring;
		}
	}
	return ($seq);
}

#mask non coding sequences
sub mask_ncs {
	my $self    = shift;
	my $seq     = shift;
	my $tmp_seq = "X" x length $seq;
	my %type    = ( CDS => 1, tRNA => 1, rRNA => 1, snRNA => 1, snoRNA => 1 );
	foreach my $feat ( @{ $self->features } ) {
		next unless $type{ $feat->type };

		foreach my $block ( @{ $feat->blocks } ) {
			next if $block->[0] > length($seq);
			next if $block->[1] < 1;
			my ( $start, $stop ) = @$block;
			$start = 1 if $start < 1;
			my $seglength = $stop - $start + 1;
			my $coding_seq = substr( $seq, $start - 1, $seglength );
			substr( $tmp_seq, $start - 1, $seglength ) = $coding_seq;
		}
	}
	return ($tmp_seq);
}

#mask non genic sequence (everything but start to stop of genes
sub mask_ngene {
	my $self    = shift;
	my $seq     = shift;
	my $tmp_seq = "X" x length $seq;
	my %type    = ( gene => 1 );
	foreach my $feat ( @{ $self->features } ) {
		next unless $type{ lc( $feat->type ) };
		foreach my $block ( @{ $feat->blocks } ) {
			next if $block->[0] > length($seq);
			next if $block->[1] < 1;
			my ( $start, $stop ) = @$block;
			$start = 1 if $start < 1;
			my $seglength = $stop - $start + 1;
			my $coding_seq = substr( $seq, $start - 1, $seglength );
			substr( $tmp_seq, $start - 1, $seglength ) = $coding_seq;
		}
	}
	return ($tmp_seq);
}

# process_location: given a location line from a gb file, try to pull
# out the coordinates
sub process_location {
	my $self    = shift;
	my %opts    = @_;
	my $locline = $opts{loc};
	my $start   = $opts{start};
	my @locdata = ();
	$locline =~ s/\s+//g;                     #remove leading spaces
	$locline =~ s/complement\((.*)\)/$1/g;    #forget about "complement"
	$locline =~ s/join\((.*)\)/$1/g;          #forget about "join"
	$locline =~ s/order\((.*)\)/$1/g;         #forget about "order"
	$locline =~ s/[<>]//g;                    #forget about "<" or ">"

	if ( $locline =~ /,/ ) {
		my @pairs = split( /,/, $locline );
		foreach my $pair (@pairs) {
			my @temp = split( /\.\./, $pair );
			if ( not defined $temp[1] ) {

				# need this for single nucleotide, which sometimes happens
				$temp[1] = $temp[0];
			}
			$temp[0] = $temp[0] - $start + 1 if $start;
			$temp[1] = $temp[1] - $start + 1 if $start;
			push( @locdata, [ $temp[0], $temp[1] ] );
		}
	}
	else {
		my @temp;

		# mdb added 11/14/13 - add handling of variation location syntax
		if ($locline =~ /\d+\^\d+/) {
			@temp = split( /\^/, $locline )
		}
		else {
			@temp = split( /\.\./, $locline );
			if ( not defined $temp[1] ) {
				# need this for single nucleotide, which sometimes happens
				$temp[1] = $temp[0];
			}
		}
		$temp[0] = $temp[0] - $start + 1 if $start;
		$temp[1] = $temp[1] - $start + 1 if $start;
		push( @locdata, [ $temp[0], $temp[1] ] );
	}
	return ( \@locdata );
}

sub get_features {
	my $self = shift;
	return ( wantarray ? @{ $self->features } : $self->features );
}

sub add_feature {
	my $self     = shift;
	my %options  = @_;
	my $type     = $options{type} if $options{type};
	my $location = $options{location} if $options{location};
	$self->_has_genes(1) if $type && $type =~ /^gene$/i;
	my $qualifiers = $options{qualifiers};
	$qualifiers = {} unless $qualifiers;
	$qualifiers = $qualifiers;
	my $annotation;
	$annotation = $options{annotation} if $options{annotation};

	#    $annotation .= "location: ".$location if $location;
	my $strand = $options{strand} if $options{strand};
	my $start  = $options{start};
	my $length = $options{length};
	my $force  =
	  $options{force}; #force addition of feature even if it is outside of sequence range;
	unless ($strand) {
		$strand = $location =~ /complement/i ? -1 : 1;
	}
	my $feature = new CoGe::Accessory::GenBank::Feature(
		{
			type       => $type,
			location   => $location,
			qualifiers => $qualifiers,
			annotation => $annotation,
			strand     => $strand,
		}
	);
	$feature->blocks(
		$self->process_location( loc => $feature->location, start => $start ) )
	  if $feature->location;
	$feature->blocks( [] ) unless $feature->blocks;
	return unless $force || ( $feature->stop && $feature->start );
	return if !$force && $feature->stop < 0;
	return
	  if !$force
	  && $start
	  && $start > 1
	  && $length
	  && $feature->start > $length;
	push @{ $self->features }, $feature;
	return;
}

# write out a header for a fasta dump;  if any extra parameters are
# sent in they will get tacked on at the end.
sub get_headerfasta {
	my $self  = shift;
	my $extra = "";
	if (@_) {
		$extra = join( " ", @_ );
	}
	my $header = ">";
	$header .= $self->accn() . " ";
	$header .= $self->definition() . " " if $self->definition;
	if ( $extra ne "" ) {
		$header .= $extra;
	}
	return ($header);
}

sub subsequence {
	my $self  = shift;
	my $begin = shift;
	my $end   = shift;

	# need to decrement by one, since substr works with arrays
	$begin--;
	$end--;

	my $seq = "";
	if (@_) {
		$seq = shift;
	}
	else {
		$seq = $self->sequence();
	}
	my $newseq = substr( $seq, $begin, ( $end - $begin + 1 ) );
	$self->sequence($newseq);
	$self->seq_length( length($newseq) );
	return ($newseq);
}

sub reverse_genbank_location {
	my $self   = shift;
	my %opts   = @_;
	my $loc    = $opts{location} || $opts{loc};
	my $length = $opts{length}
	  || $self->seq_length
	  || length( $self->sequence );
	my $complement = $loc =~ /complement/ ? 0 : 1;    #get the opposite
	$loc =~ s/complement//;
	$loc =~ s/join//;
	$loc =~ s/\(|\)|<|>//g;

	my $new_loc;
	my $count = 0;

	foreach my $set ( split( /,/, $loc ) ) {
		my ( $start, $stop ) = split(/\.\./, $set);
		$new_loc .= "," if $new_loc;
		$new_loc .= ( $length - $stop + 1 ) . ".." . ( $length - $start + 1 );
		$count++;
	}
	$new_loc = "join($new_loc)"       if $count > 1;
	$new_loc = "complement($new_loc)" if $complement;
	return $new_loc;
}

sub accession {
	shift->accn(@_);
}

sub process_wgs {
	my $self = shift;
	return if $self->no_wgs;
	my $wgs = $self->wgs;
	$wgs = $self->wgs_scafld if $self->wgs_scafld;    #this is a better assembly
	return unless $wgs;
	my @ids;
	foreach my $item ( split(/,/, $wgs) ) {
		if ( $item =~ /-/ ) {
			my ( $id1, $id2 ) = split(/-/, $item);
			my ($head) = $id1 =~ /^(\D+)/;
			my ($num1) = $id1 =~ /(\d+)/;
			my ($num2) = $id2 =~ /(\d+)/;
			my $num_len = length($num1);
			for ( my $i = $num1 ; $i <= $num2 ; $i++ ) {
				my $num = $i;
				while ( length($num) < $num_len ) {
					$num = "0" . $num;
				}
				push @ids, $head . $num;
			}
		}
		else {
			push @ids, $item;
		}
	}
	$self->wgs_data( [] ) unless $self->wgs_data();
	my $file = $self->srcfile;
	$file =~ s/[^\/]*$//;
	print "Have WGS:  " . scalar(@ids) . " entries\n";
	$self->entries_count( scalar(@ids) );
	if ( $self->max_entries && scalar(@ids) > $self->max_entries ) {
		print "Exceeded maximum entries for recursive processing.\n";
		print "Max: " . $self->max_entries, "\n";
		print "Have " . scalar(@ids), "\n";
		return;
	}

	foreach my $id (@ids) {
		my $gb = $self->new;
		$gb->debug( $self->debug ) if $self->debug;
		$gb->logfile( $self->logfile ) if $self->logfile;
		$gb->no_wgs(1);    #don't recursively get wgs data
		$gb->get_genbank_from_ncbi( accn => $id, file => $file . $id . ".gbk" );
		push @{ $self->wgs_data }, $gb;
	}
}

sub gstid {
	shift->genomic_sequence_type_id(@_);
}

1;
