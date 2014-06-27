# hacked by BCT
# Based on GBlite, by Korf
# Significant changes...
#   remove "DataBrowser" dependency
#   switch feature handling to hash-based struct
#   allow parsing of gzip'd files
##################################
package CoGe::Accessory::GBlite;
use strict;
use Compress::Zlib;
use base ('Class::Accessor');
use Data::Dumper;
# Month - for converting GenBank date format to numeric format
  my %Month = ( JAN=>'01', FEB=>'02', MAR=>'03', APR=>'04',
		MAY=>'05', JUN=>'06', JUL=>'07', AUG=>'08', SEP=>'09',
		OCT=>'10', NOV=>'11', DEC=>'12');
#  __PACKAGE__->mk_accessors(
#"organism",
#"organism_long",
#"locus",
#"accn",
#"version",

#);

sub new {
	my ($class, $file) = @_;
	if ( not defined $file ) {
		die "Must instantiate GBlite object with a File Name, Handle, or GZIP'd File Name!\n";
	}
	my $self = bless {};
	if (ref $file =~ /GLOB/) {
		$self->{FILETYPE} = "GLOB";
	} else {
		if ( $file =~ /\.gz$/ ) {
			$self->{FILETYPE} = "GZIP";
		} else {
			$self->{FILETYPE} = "FLAT";
		}
	}
	# open the file and get the file handle...
	if ( $self->{FILETYPE} eq "FLAT" ) {
		open( IN, "< $file" ) or die "Couldn't open $file\n";
		$self->{FH} = *IN;
	} elsif ( $self->{FILETYPE} eq "GZIP" ) {
		$self->{FH} = gzopen( $file, "rb" );
	}

	# initialize these vars
	$self->{LASTLINE} = "";
	$self->{DONE} = 0;
	return $self;
}

sub DESTROY {
	my ($self) = @_;
	if ( $self->{FILETYPE} eq "GZIP" ) {
		$self->{FH}->gzclose();
	} else {
		close( $self->{FH} );
	}
}

sub nextEntry {
  my ($self) = @_;
  $self->_fastForward or return 0;
  my $FH = $self->{FH};

  # These are the fields that will be kept
  my ($locus, $mol_type, $division, $date, $definition, $accession, $version,
      $gi, $keywords, $organism, $organism_long, $features, $sequence, $source);

  # get LOCUS, MOL_TYPE, DIVISION, DATE from LOCUS line
  my $locus_line = $self->{LASTLINE};
  my @field = split(/\s+/, $locus_line);
  $date = $field[-1];
  my ($day, $month, $year) = split(/\-/, $date);
  $locus    = $field[1];
  $date     = "$year-$Month{$month}-$day";
  $mol_type = $field[4];
  $division = $field[-2];

  # get DEFINITION, which may span several lines
  my @def_line;
  if ( $self->{FILETYPE} eq "GZIP" ) {
    while ($FH->gzreadline($_) > 0) {
      if (/^ACCESSION/) {
	$self->{LASTLINE} = $_;
	last;
      } else {
	push @def_line, $_;
      }
    }
  } else {
    while (<$FH>) {
      if (/^ACCESSION/) {
	$self->{LASTLINE} = $_;
	last;
      } else {
	push @def_line, $_;
      }
    }
  }
  $definition = join("", @def_line);
  $definition =~ s/\s+/ /g;
  $definition = substr($definition, 11);

  # get ACCESSION, VERSION, and GI from the VERSION line
  if ( $self->{FILETYPE} eq "GZIP" ) {
    while ($FH->gzreadline($_) > 0) {
      last if ( /^VERSION/ );
    }
  } else {
    while (<$FH>) {
      last if (/^VERSION/);
    }
  }
  my $versionline = $_;
  ($accession, $version, $gi) =
    $versionline =~ /^\S+\s+(\w+)\.(\d+)\s+GI:(\d+)/;
  if (not defined $gi) {
    warn "no gi identified >>> $versionline";
  }
  if (not defined $accession) {
    warn "no accession identified acc >> $versionline";
    $accession = $self->{LASTLINE};
    $accession =~ s/ACCESSION\s+//;
    $accession =~ s/\s+$//;
  }
  if (not defined $version) {
    warn "no version identified >> $versionline";
  }

  # parse the KEYWORDS, which may span several lines
  my %keyword;
  if ( $self->{FILETYPE} eq "GZIP" ) {
    while ($FH->gzreadline($_) > 0) {
      if (/^SOURCE/ or /^SEGMENT/ ) {
	$self->{LASTLINE} = $_;
	last;
      } else {
	$_ =~ s/[\.;]//g;	# remove punctuation
	my @words = split;
	foreach my $word (@words) {
	  $keyword{$word}++;
	}
      }
    }
  } else {
    while (<$FH>) {
      if (/^SOURCE/ or /^SEGMENT/ ) {
	$self->{LASTLINE} = $_;
	last;
      } else {
	$_ =~ s/[\.;]//g;	# remove punctuation
	my @words = split;
	foreach my $word (@words) {
	  $keyword{$word}++;
	}
      }
    }
  }
  delete $keyword{KEYWORDS};
  $keywords = [keys %keyword];

  # parse the SOURCE
  if ( $self->{LASTLINE} =~ /^SOURCE\s+(.+)/ ) {
    $source = $1;
  } else {
    if ( $self->{FILETYPE} eq "GZIP" ) {
      while ($FH->gzreadline($_) > 0) {
	if ( /^SOURCE/ ) {
	  $self->{LASTLINE} = $_;
	  last;
	}
      }
    } else {
      while (<$FH>) {
	if ( /^SOURCE/ ) {
	  $self->{LASTLINE} = $_;
	  last;
	}
      }
    }
    my $source_line = $_;
    ($source) = $source_line =~ /SOURCE\s+(.+)/;
  }
  $source =~ s/[\.;]//g;	# remove punctuation

  # parse the ORGANISM
  my @lines = ();
  if ( $self->{FILETYPE} eq "GZIP" ) {
    while ($FH->gzreadline($_) > 0) {
      if ( /^\s+ORGANISM/ )
	{
	  push( @lines, $_ );
	  last;
	}
    }
  } else {
    while (<$FH>) {
      if (/^\s+ORGANISM/)
	{
	  push( @lines, $_ );
	  last;
	}
    }
  }
  if ( $self->{FILETYPE} eq "GZIP" ) {
    while ($FH->gzreadline($_) > 0) {
      chomp;
      last if ( /^\w/);
      push( @lines, $_ );
    }
  } else {
    while (<$FH>) {
      chomp;
      last if ( /^\w/);
      push( @lines, $_ );
    }
  }
  my $orgline = shift(@lines);
  ($organism) = $orgline =~ /ORGANISM\s+(.+)/;
  $organism_long = join(" ", @lines);
  $organism_long =~ s/\s+/ /g;
  # parse the FEATURES
  if ( $self->{FILETYPE} eq "GZIP" ) {
    while ($FH->gzreadline($_) > 0) {
      last if ( /^FEATURES/ );
    }
  } else {
    while (<$FH>) {
      last if ( /^FEATURES/);
    }
  }
  @lines = ();
  $features = [];		# array of Feature objects
  if ( $self->{FILETYPE} eq "GZIP" ) {
    while ($FH->gzreadline($_) > 0) {
      chomp;
      last if ( /^BASE COUNT/ or /^ORIGIN/);
      if ( substr($_, 5, 1 ) ne ' ' and @lines ) {
	push @$features, CoGe::Accessory::GBlite::Feature::new(\@lines);
	@lines = ( $_ );
      } else {
	push @lines, $_;
      }
    }
  } else {
    while (<$FH>) {
      chomp;
      last if ( /^BASE COUNT/ or /^ORIGIN/);
      if (substr($_, 5, 1) ne ' ' and @lines) {
	push @$features, CoGe::Accessory::GBlite::Feature::new(\@lines);
	@lines = ($_);
      } else {
	push @lines, $_;
      }
    }
  }
  if ( @lines > 0 ) {
    push @$features, CoGe::Accessory::GBlite::Feature::new(\@lines);
  }
  if (@$features == 0) {
    die "unexpected fatal parsing error\n";
  }

  # parse the SEQUENCE
  my @seq;
  if ( $self->{FILETYPE} eq "GZIP" ) {
    while ($FH->gzreadline($_) > 0) {
      last if ( /^\/\// );
      next if /^ORIGIN/;
      $_ =~ s/\d+//g;
      $_ =~ s/\s+//g;
      push( @seq, $_ );
    }
  } else {
    while (<$FH>) {
      last if /^\/\//;
      next if /^ORIGIN/;
      $_ =~ s/\d+//g;
      $_ =~ s/\s+//g;
      push @seq, $_;
    }
  }
  $sequence = join("", @seq);
  $sequence = uc $sequence;

  $self->{LASTLINE} = $_;
#  print join (",\t","l: $locus", "mt: $mol_type", "D: $division", "DA: $date" ,"DEF: $definition", "ACC: $accession", "V: $version", "GI: $gi", "KW: $keywords", "S: $source", "O: $organism", "OL: $organism_long", "F: $features"),"\n";
  my $entry = CoGe::Accessory::GBlite::Entry::new($locus, $mol_type, $division, $date,$definition, $accession, $version, $gi, $keywords, $source, $organism, $organism_long, $features, $sequence);

  return $entry;
}

sub _fastForward {
  my ($self) = @_;
  return 0 if $self->{DONE};
  return 1 if $self->{LASTLINE} =~ /^LOCUS/;
  my $FH = $self->{FH};
  if ( $self->{FILETYPE} eq "GZIP" ) {
    while ( $FH->gzreadline($_) > 0 ) {
      if ($_ =~ /^LOCUS/) {
	$self->{LASTLINE} = $_;
	return 1;
      }
    }
  } else {
    while (<$FH>) {
      if ($_ =~ /^LOCUS/) {
	$self->{LASTLINE} = $_;
	return 1;
      }
    }
  }
  #	print STDERR "INSIDE _fastForward ($_)\n";
  return 0 if not defined $_;
  #warn "Possible parse error in _fastForward in GBlite.pm\n", $_;
}

################################################################################
# GBlite::Entry
################################################################################
package CoGe::Accessory::GBlite::Entry;
use strict;

# Field - these are the fields that will be parsed for every GenBank entry

my @Field = qw(
	LOCUS
	MOL_TYPE
	DIVISION
	DATE
	DEFINITION
	ACCESSION
	VERSION
	GI
	KEYWORDS
	SOURCE
	ORGANISM
	ORGANISM_LONG
	FEATURES
	SEQUENCE
);

sub new {
	my $entry = bless {};

	($entry->{LOCUS}, $entry->{MOL_TYPE}, $entry->{DIVISION}, $entry->{DATE},
		$entry->{DEFINITION}, $entry->{ACCESSION}, $entry->{VERSION},
		$entry->{GI}, $entry->{KEYWORDS}, $entry->{SOURCE},
    $entry->{ORGANISM}, $entry->{ORGANISM_LONG}, $entry->{FEATURES},
    $entry->{SEQUENCE}) = @_;

	my $CONSTRUCTOR_ERROR = 0;
	foreach my $name (@Field) {
    # this can only get filled in the Feature parsing stage of a script, so just skip it
	  next if $name eq "CHROMOSOME";
		if (not defined $entry->{$name}) {
			$CONSTRUCTOR_ERROR++;
			warn "CoGe::Accessory::GBlite::Entry constructor error, $name undefined\n";
		}
	}
	if ($CONSTRUCTOR_ERROR) {
		warn "Problems with entry object construction";
	}

	$entry->{DIR} = set_dir( $entry );

	return $entry;
}

sub locus           {shift->{LOCUS}}
sub mol_type        {shift->{MOL_TYPE}}
sub division        {shift->{DIVISION}}
sub date            {shift->{DATE}}
sub definition      {shift->{DEFINITION}}
sub accession       {shift->{ACCESSION}}
sub version         {shift->{VERSION}}
sub gi              {shift->{GI}}
sub keywords        {shift->{KEYWORDS}}
sub source          {shift->{SOURCE}}
sub organism        {shift->{ORGANISM}}
sub organism_long   {shift->{ORGANISM_LONG}}
sub features        {shift->{FEATURES}}
sub sequence        {shift->{SEQUENCE}}
sub length          {length(shift->{SEQUENCE})}
sub dir             {shift->{DIR}}
sub chromosome {
  my $self = shift;
  if ( @_ ) {
    $self->{CHROMOSOME} = shift;
  }
  return $self->{CHROMOSOME};
}
sub set_dir {
	my $self = shift;
	# reads definition line and tries to determine
	# the orientation of self sequence.
	my $defn = $self->{DEFINITION};
	my $division = $self->{DIVISION};
	my $dir = "";

	if ( $division eq "EST" ) {
		if ( $defn =~ /5'/ ) {
			$dir = "5'";
		} elsif ( $defn =~ /3'/ ) {
			$dir="3'";
		} else {
			$dir = "?'";
		}
	} else {
		# these seqs must have come from a PRI or ROD, etc.
		# Thus they are "knowns", and if moltype=mRNA, then they are 5'
		$dir = "5'";
	}

	return $dir;
}

################################################################################
# GBlite::Feature
################################################################################
package CoGe::Accessory::GBlite::Feature;
use strict;

sub new {
  my $feature = bless {};
  my ($lines) = @_;
  if ( not @$lines > 0 ) {
    return 0;
  }
  my $string = join("", @$lines);  # join all lines
  my @part = split(/\s{10,30}\/(?=\S)/x, $string); # split qualifiers from key/location
  $string =~ s/\s+/ /g;            # trim multiple spaces to one space
  my $key_loc = shift @part;
  my ($key, $location) = $key_loc =~ /^\s+(\S+)\s+(.+)/;
  $location =~ s/\s+//g;
  my $qualifiers = {};
  foreach my $qual (@part) {
    if ($qual =~ /\S=\S/) {
      my ($key, $value) = $qual =~ /^(\S+?)=(.+)/;
      $value =~ s/"//g;  # remove "
      $value =~ s/\s+$//g; # remove trailing whitespace
      $value =~ s/\s+/ /g; # substitute multiple spaces for one
      if (not defined $key) {
	print "$key --> $qual\n";
	print "@$lines\n";
	die "CoGe::Accessory::GBlite::Feature constructor error\n";
      }
      if ($key eq 'translation') {$value =~ s/\s+//g}
      $qualifiers->{$key} .= "$value ";
    } else {
      if ( $qual eq "pseudo" ) {
	$qualifiers->{"pseudo"} = 1;
      }
      next if ( $key eq "gene" );
      $qualifiers->{$key} = "";
    }
  }

  if ( keys %$qualifiers > 0 ) {
    $feature->{QUALIFIERS} = $qualifiers;
  } else {
    $feature->{QUALIFIERS} = {};
  }

  $feature->{KEY}        = $key;
  $feature->{LOCATION}   = $location;

  $feature->{STRAND} = $location =~ /complement/ ? "-1" : 1;
  $location =~ s/complement//;
  $location =~ s/join//;
  $location =~ s/\(|\)//g;
  my @temp = split( /\.\./, $location );
  $feature->{START} = $temp[0];
  $feature->{STOP} = $temp[1];

  return $feature;
}

sub key        {shift->{KEY}}
sub location   {shift->{LOCATION}}
sub start      {shift->{START}}
sub stop       {shift->{STOP}}
sub strand     {shift->{STRAND}}
sub qualifiers {shift->{QUALIFIERS}}

1;

__END__

=head1 NAME

GBlite.pm

=head1 SYNOPSIS

 use GBlite;
 my $genbank = new GBlite(\*STDIN);
 while(my $entry = $genbank->nextEntry) {
   $entry->locus;
   $entry->mol_type;
   $entry->division;
   $entry->date;           # yyyy-mm-dd
   $entry->definition;
   $entry->accession;
   $entry->version;
   $entry->gi;
   $entry->keywords;       # reference to ARRAY
   $entry->organism;
   $entry->organism_long;
   $entry->features;       # reference to ARRAY
   $entry->sequence;
   $entry->chromosome;

   foreach my $feature (@{$entry->features}) {
     $feature->key;
     $feature->location;
     $feature->qualifiers; # reference to HASH
   }
 }

=head1 DESCRIPTION

GBlite is a package for parsing concatenated GenBank flat files. The GenBank
format is a common format for bioinformatics. Its specification is complicated,
and anyone using this module should at least skim the GenBank release.notes and
the DDJB/EMBL/GenBank feature table specification. These documents are
available from the NCBI.

=head1 AUTHOR

Ian Korf (ikorf@sapiens.wustl.edu, http://sapiens.wustl.edu/~ikorf)

=head1 ACKNOWLEDGEMENTS

This software was developed at Washington Univeristy, St. Louis, MO.

=head1 COPYRIGHT

Copyright (C) 2000 Ian Korf. All Rights Reserved.

=head1 DISCLAIMER

This software is provided "as is" without warranty of any kind.

=cut
