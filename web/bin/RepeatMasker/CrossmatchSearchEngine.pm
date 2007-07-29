#!/usr/bin/perl -w
##---------------------------------------------------------------------------##
##  File:
##      @(#) CrossmatchSearchEngine.pm
##  Author:
##      Robert M. Hubley   rhubley@systemsbiology.org
##  Description:
##      An implementation of SearchEngineI for the
##      the cross_match search engine.
##
#******************************************************************************
#* Copyright (C) Institute for Systems Biology 2003-2004 Developed by
#* Arian Smit and Robert Hubley.
#*
#* This work is licensed under the Open Source License v2.1.  To view a copy
#* of this license, visit http://www.opensource.org/licenses/osl-2.1.php or
#* see the license.txt file contained in this distribution.
#*
#******************************************************************************
# Implementation Details:
#
# bless(
#      'CrossmatchSearchEngine' );
#
###############################################################################
# ChangeLog
#
#     $Log: CrossmatchSearchEngine.pm,v $
#     Revision 1.56  2007/05/17 21:01:53  rhubley
#     Cleanup before a distribution
#
#
###############################################################################
# To Do:
#
#

=head1 NAME

CrossmatchSearchEngine

=head1 SYNOPSIS

use CrossmatchSearchEngine

  my $cEngine = CrossmatchSearchEngine->new( 
                       pathToEngine=>"/usr/local/bin/crossmatch" );

  $cEngine->setMatrix( "/users/bob/simple.matrix" );
  $cEngine->setQuery( "/users/bob/query.fasta" );
  $cEngine->setSubject( "/users/bob/subject.fasta" );
  my $searchResults = $cEngine->search();

Usage: 

=head1 DESCRIPTION

An encapsulation class for the CrossMatch search engine.  It is capable
for running searches and returning the results as a SearchResultCollection.

The CrossMatch parser built into this object captures several types of i
data from the voluminous cross_match output stream.  The first is 
the high-scoring pair line or "hit" line:

Foward Strand:

 SW    perc perc perc qry   qry   qry  qry   subj           subj  subj subj 
 score div. del. ins. seq   begin end (left) seq            begin end (left) OV 
 ------------------------------------------------------------------------------
 2334  8.44 0.00 3.25 Human 127   737 (8222) AluSx#SINE/Alu  1    298 (14)   *  

Reverse Strand:

 SW    perc perc perc qry   qry   qry  qry     subj       subj  subj subj 
 score div. del. ins. seq   begin end (left) C seq       (left) end  begin OV 
 ------------------------------------------------------------------------------
 2334  8.44 0.00 3.25 Human 127   737 (8222) C AluSx#SINE (14)  298  1       


  SW score    = smith-waterman score of the match (complexity-adjusted, 
                by default).
  perc div.   = %substitutions in matching region.
  perc del.   = %deletions (in query seq rel to subject) in matching region.
  perc ins.   = %insertions (in query seq rel to subject) in matching region.
  qry seq     = id of query sequence.
  qry begin   = starting position of match in query sequence.
  qry end     = ending position of match in query sequence.
  qry (left)  = no. of bases in query sequence past the ending position of 
                match (so 0 means that the match extended all the way to 
                the end of the query sequence).
  C           = match is with the Complement of subject sequence.
  subj seq    = id of the subject sequence.
  subj (left) = The remaining bases in (complement of) subject sequence 
                prior to beginning of the match.
  subj end    = starting position of match in subject sequence (using 
                top-strand numbering).
  subj begin  = ending position of match in subject sequence.
  OV          = A "*" in this field indicates that there is a higher-scoring 
                match whose domain partly includes the domain of this match.

RepeatMasker has added two new fields to this format for it's standard
annotation output:

 SW    perc perc perc qry qry   qry  qry   subj subj  subj subj     Lineage
 score div. del. ins. seq begin end (left) seq  begin end (left) Id   Id    OV 
 ------------------------------------------------------------------------------

  Id         = A unique identifier for the repeat which corresponds
               to the repeat in the .align files.
  LineageId  = ?

The second type of data this object collects (optionally) is the
alignment data.  If the cross_match output contains alignment
data; typically in the form:

 NT_004321_1       1047 CACCCACATGCACACACACACGCGCGCACACACGCACACGCACACACATG 1096
                           v    ii           i i i       i     i        ii
 (CA)n#Simple_re      1 CACACACACACACACACACACACACACACACACACACACACACACACACA 50

 NT_004321_1       1097 CACACACGCGCAC--ACACGCACACATATGCACACACAAACGCACA 1142
                               i i         i      i ii        v  i
 (CA)n#Simple_re     51 CACACACACACACACACACACACACACACACACACACACACACACA 96

This object will parse these lines and include the
query and subject alignment sequences in the object.  NOTE: This can
greatly increase memory usage as the sequences are not compressed and
will include the gap characters ("-").

If the full cross_match output is available the matrix filename will be 
read from the cross_match invocation line.  If the input file is from 
RepeatMasker and contains the:

  Assumed background GC level in scoring matrices is 49 %

line, the GC attribute will be set.

Lastly, if cross_match was run with alignments turned on the 
transitions and transversion information will be read from
the following line:

  Transitions / transversions = 2.00 (42 / 21)


=head1 SEE ALSO

=over 4

SearchEngineI, SearchResultCollection

=back

=head1 COPYRIGHT

Copyright 2004 Institute for Systems Biology

=head1 AUTHOR

Robert Hubley <rhubley@systemsbiology.org>

=head1 INSTANCE METHODS

=cut 

package CrossmatchSearchEngine;
use strict;
use POSIX qw(:sys_wait_h);
use SearchEngineI;
use SearchResultCollection;
use Data::Dumper;
use Carp;
use FileHandle;
use IPC::Open3;
use vars qw(@ISA @EXPORT @EXPORT_OK %EXPORT_TAGS $VERSION);

require Exporter;

@ISA = qw(Exporter SearchEngineI);

@EXPORT = qw();

@EXPORT_OK = qw();

%EXPORT_TAGS = ( all => [ @EXPORT_OK ] );

my $CLASS = "CrossmatchSearchEngine";
my $DEBUG = 0;

##-------------------------------------------------------------------------##
## Constructor
##-------------------------------------------------------------------------##
sub new {
  my $class          = shift;
  my %nameValuePairs = @_;

  croak $CLASS
      . "::new: Missing path to search engine!\n\n"
      . "use \$searchEngine = $CLASS->new( pathToEngine=>\"/usr/local/"
      . "bin/cross_match\")\n"
      if ( not defined $nameValuePairs{'pathToEngine'} );

  # Create ourself as a hash
  my $this = {};

  # Bless this hash in the name of the father, the son...
  bless $this, $class;

  $this->setPathToEngine( $nameValuePairs{'pathToEngine'} );

  return $this;
}

##-------------------------------------------------------------------------##
## Get and Set Methods
##-------------------------------------------------------------------------##

##-------------------------------------------------------------------------##

=head2 get_setPathToEngine()

  Use: my $value    = getPathToEngine( );
  Use: my $oldValue = setPathToEngine( $value );

  Get/Set the fully qualified path to the search engine
  binary file.

=cut

##-------------------------------------------------------------------------##
sub getPathToEngine {
  my $this = shift;

  return $this->{'pathToEngine'};
}

sub setPathToEngine {
  my $this  = shift;
  my $value = shift;

  croak $CLASS. "::setPathToEngine( $value ): Program does not exist!"
      if ( not -x $value || `which $value` );

  my $result = `$value 2>&1`;
  if ( $result =~ /.*cross_match version ([\.\d]+).*/ ) {
    $this->{'version'} = $1;
  }

  my $oldValue = $this->{'pathToEngine'};
  $this->{'pathToEngine'} = $value;

  return $oldValue;
}

##-------------------------------------------------------------------------##
## General Object Methods
##-------------------------------------------------------------------------##

##-------------------------------------------------------------------------##

=head2 search()

  Use: my ( $resultCode, $SearchResultCollectionI ) = search( );
 or
  Use: my ( $resultCode, $SearchResultCollectionI )
                          = search( matrix=>"7p16g.matrix",
                                    ...
                                  );

  Run the search and return a SearchResultsCollectionI.

=cut

##-------------------------------------------------------------------------##
sub search {
  my $this           = shift;
  my %nameValuePairs = @_;

  if ( %nameValuePairs ) {
    while ( my ( $name, $value ) = each( %nameValuePairs ) ) {
      my $method = "set" . _ucFirst( $name );
      unless ( $this->can( $method ) ) {
        croak( $CLASS . "::search: Instance variable $name doesn't exist." );
      }
      $this->$method( $value );
    }
  }

  # Test if engine is available
  my $engine = $this->getPathToEngine();
  if ( !defined $engine || !-f "$engine" ) {
    croak $CLASS
        . "::search: The path to the search engine is undefined or\n"
        . "is set incorrectly: $engine\n";
  }

  # Generate parameter line
  my $parameters;
  my $value;
  if ( ( $value = $this->getScoreMode() ) ) {
    if ( $value == SearchEngineI::basicScoreMode ) {
      $parameters .= " -raw";
    }
  }
  if ( ( $value = $this->getWordRaw() ) ) {
    $parameters .= " -word_raw" if ( $value > 0 );
  }
  if ( ( $value = $this->getGenerateAlignments() ) ) {
    $parameters .= " -alignments" if ( $value > 0 );
  }
  if ( defined( $value = $this->getGapInit() )
       && $value =~ /\d+/ )
  {
    $parameters .= " -gap_init $value";
  }
  if ( defined( $value = $this->getInsGapExt() )
       && $value =~ /\d+/ )
  {
    $parameters .= " -ins_gap_ext $value";
  }
  if ( defined( $value = $this->getDelGapExt() )
       && $value =~ /\d+/ )
  {
    $parameters .= " -del_gap_ext $value";
  }
  if ( ( $value = $this->getMinMatch() ) ) {
    $parameters .= " -minmatch $value" if ( $value > 0 );
  }
  if ( ( $value = $this->getMinScore() ) ) {
    $parameters .= " -minscore $value" if ( $value > 0 );
  }
  if ( ( $value = $this->getBandwidth() ) ) {
    $parameters .= " -bandwidth $value" if ( $value > 0 );
  }
  if ( ( $value = $this->getMaskLevel() ) ) {
    $parameters .= " -masklevel $value" if ( $value > 0 );
  }
  if ( ( $value = $this->getMatrix() ) ) {

    # test if matrix exists
    if ( -f $value ) {
      $parameters .= " -matrix $value";
    }
    else {
      croak $CLASS. "::search: Error...matrix ($value) does not exist!\n";
    }
  }
  if ( ( $value = $this->getQuery() ) ) {
    if ( -f $value ) {
      $parameters .= " $value";
    }
    else {
      croak $CLASS. "::search: Error...query ($value) does not exist!\n";
    }
  }
  else {
    croak $CLASS. "::search: Error query undefined!\n";
  }
  if ( ( $value = $this->getSubject() ) ) {
    if ( -f $value ) {
      $parameters .= " $value";
    }
    else {
      croak $CLASS. "::search: Error...subject ($value) does not exist!\n";
    }
  }
  else {
    croak $CLASS. "::search: Error subject undefined!\n";
  }

  # Invoke engine and handle errors
  #print "Did run parameters=$parameters\n";
  my $PINPUT  = new FileHandle;
  my $POUTPUT = new FileHandle;
  my $PERROR  = new FileHandle;
  my $pid;

  #eval { $pid = open3( $PINPUT, $POUTPUT, $PERROR, "$engine $parameters" ) };
  #$pid = open3( $PINPUT, $POUTPUT, $PERROR, "$engine $parameters" );
  # NOTE: Here I wanted to use open3.  The problem with open3 is that
  #       programs writing output to both stderr and stdout will block
  #       if either buffer becomes full while the other is being read from.
  #       So....IO::Select needs to be used to alternate between reading
  #       from each one.

  print $CLASS
      . "::search(): Running crossmatch as:\n    "
      . "$engine $parameters 2>/dev/null |\n"
      if ( $DEBUG );

  $pid = open( $POUTPUT, "$engine $parameters 2>/dev/null |" );

  my $resultCode = 0;
  my $searchResultsCollection;

  # Create SearchResultCollection object from
  # the engine results.
  if ( $DEBUG ) {
    ## Create a debug file
    my $outFile;
    do {
      $outFile = "cmResults-" . time() . ".out";
    } while ( -f $outFile );
    open OUT, ">$outFile";
    while ( <$POUTPUT> ) {
      print OUT $_;
    }
    close OUT;
    close $POUTPUT;
    $resultCode = $?;
    $searchResultsCollection = parseOutput( searchOutput => $outFile );
  }
  else {
    $searchResultsCollection = parseOutput( searchOutput => $POUTPUT );
    close $POUTPUT;
    $resultCode = $?;
  }

  return ( $resultCode, $searchResultsCollection );
}

##---------------------------------------------------------------------##

=head1 CLASS METHODS

=cut

##---------------------------------------------------------------------##

##---------------------------------------------------------------------##

=head2 parseOutput()

  Use: my $SearchResultCollection = parseOutput(
                                     searchOutput => $filename|$FH,
                                     [excludeAlignments => 1],
                                     [scoreHighThresh => #],
                                     [scoreLowThresh => #],
                                     [subjPattern => ""],
                                     [queryPattern => ""] 
                                               );

  Parse the result of a search and return a SearchResultCollection.

=cut

##---------------------------------------------------------------------##
sub parseOutput {
  my %nameValueParams = @_;

  croak $CLASS. "::parseOutput() missing searchOutput parameter!\n"
      if ( !exists $nameValueParams{'searchOutput'} );

  my $CMFILE;
  if ( ref( $nameValueParams{'searchOutput'} ) !~ /GLOB|FileHandle/ ) {
    print $CLASS
        . "::parseOutput() Opening file "
        . $nameValueParams{'searchOutput'} . "\n"
        if ( $DEBUG );
    open $CMFILE, $nameValueParams{'searchOutput'}
        or die $CLASS
        . "::parseOutput: Unable to open "
        . "results file: $nameValueParams{'searchOutput'} : $!";
  }
  else {
    $CMFILE = $nameValueParams{'searchOutput'};
  }

  #
  #
  #
  my @outfileFwdStrandKeys = qw( score pctDiverge pctDelete pctInsert queryName
      queryStart queryEnd queryRemaining orientation subjName subjType
      subjStart subjEnd subjRemaining id lineageId
      overlap );
  my @outfileRevStrandKeys = qw( score pctDiverge pctDelete pctInsert
      queryName queryStart queryEnd queryRemaining orientation
      subjName subjType subjRemaining subjEnd
      subjStart id lineageId overlap );
  my @crossmatchFwdStrandKeys = qw( score pctDiverge pctDelete pctInsert
      queryName queryStart queryEnd queryRemaining subjName
      subjStart subjEnd subjRemaining id lineageId
      overlap );
  my @crossmatchRevStrandKeys = qw( score pctDiverge pctDelete pctInsert
      queryName queryStart queryEnd queryRemaining orientation
      subjName subjRemaining subjEnd
      subjStart id lineageId overlap );

  my $seqName;
  my $querySeq;
  my $transI;
  my $transV;
  my $subjSeq;
  my $result;
  my $alignPos          = 0;
  my $sourceFileIndex   = 0;    # The index of the alignment
                                #  within the source file.
  my $queryBackgroundGC = 0;
  my $subjComplemented  = 0;
  my $matrix            = "";
  my @hdrLineArray      = ();
  my $alignmts;

  my $resultColl = SearchResultCollection->new();

  while ( <$CMFILE> ) {

    #
    # RepeatMasker adds this to the crossmatch alignment output
    #  TODO: REMOVE: I don't believe this is the case anymore
    $queryBackgroundGC = $1 if ( /^Assumed background.*is (\d+) %/ );

    #
    # If we have complete output we can also grab the matrix
    # used.
    #
    if ( /Score matrix\s+(\S+)/ ) {
      my @path = split( /[\\\/]/ );
      $matrix = $path[ $#path ];
      $matrix =~ s/[\n\r]//g;
    }

    #
    # Look for transition/transversion ratio
    #
    if ( /Transitions.*\((\d+)\s*\/\s*(\d+)\)/ ) {
      $transI = $1;
      $transV = $2;
    }

    # Look for a header line (e.g.):
    if ( /^\s*\d+\s+\d+(\.\d+)?/ ) {
      @hdrLineArray = split;    # Create an array out of this line
           # Check to see if we should filter out low scoring hits
           #unless (  ( exists $this->{scoreHighThresh} &&
           #            $hdrLineArray[0] >= $this->{scoreHighThresh}) ||
           #          ( exists $this->{scoreLowThresh} &&
           #            $hdrLineArray[0] <= $this->{scoreLowThresh}) ||
           #          ( exists $this->{subjPattern} &&
           #           (($hdrLineArray[8] =~ /$this->{subjPattern}/o)<1 &&
           #           ($hdrLineArray[9] =~ /$this->{subjPattern}/o)<1)) ||
           #          ( exists $this->{queryPattern} &&
           #            ($hdrLineArray[4] =~ /$this->{queryPattern}/o)<1 ) ||
           #          ( exists $this->{gcLowThresh} &&
           #           $queryBackgroundGC < $this->{gcLowThresh}) ||
           #          ( exists $this->{gcHighThresh} &&
           #           $queryBackgroundGC > $this->{gcHighThresh})) {
      if ( $#hdrLineArray > 10 ) {

        # break up the header line into a hash using @fwdStrandKeys and
        # @revStrandKeys depending the existance of a "C"
        my @dataArray      = @hdrLineArray;
        my @fieldKeys      = ();
        my %nameValuePairs = ();

        # Is this a crossmatch line or an outfile line?
        #   Forward strand results have a "+" orientation field
        #   which is not the case in crossmatch files.
        #   Out files also have a non-numeric type field following
        #   the subject name.
        if ( $hdrLineArray[ 8 ] eq "+"
             || !( $hdrLineArray[ 10 ] =~ /^[\(\)\d]+$/ ) )
        {

          # Definately an out file
          if ( $hdrLineArray[ 8 ] eq "+" ) {
            $nameValuePairs{'orientation'} = "";
            @fieldKeys = @outfileFwdStrandKeys;
          }
          else {
            $nameValuePairs{'orientation'} = "C";
            @fieldKeys = @outfileRevStrandKeys;
          }
        }
        else {

          # Probably a crossmatch file
          if ( $hdrLineArray[ 8 ] eq "C" ) {
            $nameValuePairs{'orientation'} = "C";
            @fieldKeys = @crossmatchRevStrandKeys;
          }
          else {
            $nameValuePairs{'orientation'} = "";
            @fieldKeys = @crossmatchFwdStrandKeys;
          }
        }

        map {
          my $datum = shift( @dataArray );
          my $field = $_;
          if ( $field ne "orientation" ) {
            $datum =~ s/[\(\)]//g
                if (    $field eq "subjRemaining"
                     || $field eq "queryRemaining" );
            if ( $#dataArray == 0 && $datum eq "*" ) {
              $field = 'id';
            }
            $nameValuePairs{$field} = $datum;
          }
        } @fieldKeys;

        $nameValuePairs{'matrixName'} = $matrix;

        $result = SearchResult->new( %nameValuePairs );
        $resultColl->add( $result );
        $querySeq = "";
        $subjSeq  = "";
      }

      @hdrLineArray = ();

    }

    #
    # Concatenate the alignment sequences
    #
    if ( /^(C?)\s+(\S+)\s+\d+\s+(\S+)\s+\d+\s*$/
         && !exists $nameValueParams{'excludeAlignments'} )
    {
      $subjComplemented = 1 if ( $1 eq "C" );
      if ( $alignPos == 0 ) {

        # This is the query sequence
        $querySeq .= $3;
      }
      else {

        # This is the subj sequence
        $subjSeq .= $3;
      }
      $alignPos ^= 1;
    }

    #
    # Look for a signal for the end of an alignment
    #
    if ( /Gap_init rate/ && defined $result ) {

      # Store this alignment in our data structure.
      if ( $subjComplemented ) {

        # Crossmatch does not complement the subject
        # sequence in it's alignments ( only the query
        # sequence).  However our SearchResult object
        # requires that the query sequence be in
        # the forward direction.  So...do a simple
        # reversal of the sequence so that we can
        # store it in the object.
        $querySeq = reverse $querySeq;
        $querySeq =~ tr/ACGTYRMKHBVD/TGCARYKMDVBH/;    # complement
        $subjSeq = reverse $subjSeq;
        $subjSeq =~ tr/ACGTYRMKHBVD/TGCARYKMDVBH/;     # complement
      }
      $result->setQueryString( $querySeq );
      $result->setSubjString( $subjSeq );

      $result           = undef;
      $querySeq         = "";
      $subjSeq          = "";
      $transV           = 0;
      $transI           = 0;
      $alignPos         = 0;
      $subjComplemented = 0;
      @hdrLineArray     = ();
    }
  }
  close $CMFILE;

  return $resultColl;
}

##-------------------------------------------------------------------------##
## Private Methods
##-------------------------------------------------------------------------##

##-------------------------------------------------------------------------##
## Use: my _ucFirst( $string );
##
##   Uppercases the first character in a string and returns it.
##
##-------------------------------------------------------------------------##
sub _ucFirst {
  my $string = shift;

  if ( defined $string && $string ne "" ) {
    substr( $string, 0, 1 ) = uc( substr( $string, 0, 1 ) );
  }
  return $string;
}

##-------------------------------------------------------------------------##
## Use:  my _systemint( $cmd );
##
##     Interruptible system call routine.
##
##  Returns
##
##  Globals Used: None
##-------------------------------------------------------------------------##?
sub _systemint {
  my ( $cmd ) = @_;
  my ( $pid );
  my ( $flag ) = 0;

  local $SIG{INT}  = sub { &_handler( @_ ) if ( $flag ) };    #^C
  local $SIG{QUIT} = sub { &_handler( @_ ) if ( $flag ) };    #^\
  local $SIG{TERM} =
      sub { &_handler( @_ ) if ( $flag ) };    #kill command or system crash
  local $SIG{HUP} = sub { &_handler( @_ ) if ( $flag ) };

FORK: {
    if ( $pid = fork ) {
      $flag = 1;
      waitpid( $pid, 0 );                      #Waits for child to finish...
      my ( $status ) = $?;
      if ( WIFSTOPPED( $status ) ) {
        my ( $signal ) = WSTOPSIG( $status );
        print "\nforksys:  Program terminated by a signal $signal.\n";
        print "The executing command was:  $cmd\n";
        return 1;
      }
      if ( WIFEXITED( $status ) ) {
        my ( $temp ) = WEXITSTATUS( $status );
        return $temp;
      }
      if ( WIFSIGNALED( $status ) ) {
        my ( $signal ) = WTERMSIG( $status );
        return $signal;
      }
      if ( WIFSIGNALED( $status ) ) {
        my ( $signal ) = WTERMSIG( $status );
        print "\nforksys:  Program terminated by a signal $signal.\n";
        print "The executing command was:  $cmd\n";
        return 1;
      }
    }
    elsif ( defined $pid ) {
      exec( "$cmd" ) or die "Exec $cmd failed\n";
    }
    elsif ( $! =~ /No more process/o ) {
      print "$!\n";
      sleep 5;
      redo FORK;
    }
    else {
      die "Can't fork...\n";
    }
  }
}
##-------------------------------------------------------------------------##
## Use: my _handler( $sig );
##
##  Interrupt handler used by systemint() ###
##
##  Returns
##
##  Globals Used: None
##-------------------------------------------------------------------------##?
sub _handler {
  my ( $sig ) = @_;

  print $CLASS. "_handler(): Aborting with a SIG$sig\n";
  exit( -1 );
}

##-------------------------------------------------------------------------##
## Serialization & Debug Routines
##-------------------------------------------------------------------------##

##-------------------------------------------------------------------------##
## Use: my $string = toString([$this]);
##
##      $this         : Normally passed implicitly
##
##  Returns
##
##      Uses the Data::Dumper to create a printable reprentation
##      of a data structure.  In this case the object data itself.
##
##-------------------------------------------------------------------------##
sub toString {
  my $this = shift;
  my $data_dumper = new Data::Dumper( [ $this ] );
  $data_dumper->Purity( 1 )->Terse( 1 )->Deepcopy( 1 );
  return $data_dumper->Dump();
}

##-------------------------------------------------------------------------##
## Use: my serializeOUT( $filename );
##
##	  $filename	: A filename to be created
##
##  Returns
##
##	Uses the Data::Dumper module to save out the data
##	structure as a text file.  This text file can be
##	read back into an object of this type.
##
##-------------------------------------------------------------------------##
sub serializeOUT {
  my $this     = shift;
  my $fileName = shift;

  my $data_dumper = new Data::Dumper( [ $this ] );
  $data_dumper->Purity( 1 )->Terse( 1 )->Deepcopy( 1 );
  open OUT, ">$fileName";
  print OUT $data_dumper->Dump();
  close OUT;
}

##-------------------------------------------------------------------------##
## Use: my serializeIN( $filename );
##
##	$filename	: A filename containing a serialized object
##
##  Returns
##
##	Uses the Data::Dumper module to read in data
##	from a serialized PERL object or data structure.
##
##-------------------------------------------------------------------------##
sub serializeIN {
  my $this         = shift;
  my $fileName     = shift;
  my $fileContents = "";
  my $oldSep       = $/;
  undef $/;
  my $in;
  open $in, "$fileName";
  $fileContents = <$in>;
  $/            = $oldSep;
  close $in;
  return eval( $fileContents );
}

1;
