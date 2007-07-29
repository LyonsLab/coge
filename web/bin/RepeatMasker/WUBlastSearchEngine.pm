
#!/usr/bin/perl -w
##---------------------------------------------------------------------------##
##  File:
##      @(#) WUBlastSearchEngine.pm
##  Author:
##      Robert M. Hubley   rhubley@systemsbiology.org
##  Description:
##      An implementation of SearchEngineI for the
##      the WU-Blast search engine.
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
#      'WUBlastSearchEngine' );
#
###############################################################################
# ChangeLog
#
#     $Log: WUBlastSearchEngine.pm,v $
#     Revision 1.84  2007/06/14 17:33:41  rhubley
#      - Lou Shaoke reported a problem with the latest version of wublast
#        ( 04-May-2006 ).  The new version of wublast recognizes NCBI "gi"
#        identifiers and omits them in the output now by default. To get
#        the sequence identifier back in the output you must use "-gi"
#        option by default.  This change appears to be backwards compatable
#        with previous version of wublast.
#
#     Revision 1.83  2007/05/17 21:01:54  rhubley
#     Cleanup before a distribution
#
#
###############################################################################
# To Do:
#
#

=head1 NAME

WUBlastSearchEngine

=head1 SYNOPSIS

use WUBlastSearchEngine

Usage: 

  use SearchEngineI;
  use WUBlastSearchEngine;
  use SearchResultCollection;

  my $wuEngine = WUBlastSearchEngine->new( 
                    pathToEngine=>"/usr/local/wublast/blastp" );

  $wuEngine->setMatrix( "/users/bob/simple.matrix" );
  $wuEngine->setQuery( "/users/bob/query.fasta" );
  $wuEngine->setSubject( "/users/bob/subject.fasta" );
  my $searchResults = $wuEngine->search();

=head1 DESCRIPTION

  A concrete implementation of the abstract class / interface SearchEngineI
  which use the WUBlast sequence search engine.

=head1 INSTANCE METHODS

=cut 

package WUBlastSearchEngine;
use strict;
use SearchEngineI;
use SearchResultCollection;
use Data::Dumper;
use FileHandle;
use File::Basename;
use Carp;
use vars qw(@ISA @EXPORT @EXPORT_OK %EXPORT_TAGS $VERSION);

require Exporter;

@ISA = qw(Exporter SearchEngineI);

@EXPORT = qw();

@EXPORT_OK = qw();

%EXPORT_TAGS = ( all => [ @EXPORT_OK ] );

#
# Version
#
my $VERSION = 0.1;
my $CLASS   = "WUBlastSearchEngine";
my $DEBUG   = 0;

##-------------------------------------------------------------------------##
## Constructor
##-------------------------------------------------------------------------##
sub new {
  my $class          = shift;
  my %nameValuePairs = @_;

  croak $CLASS
      . "::new: Missing path to search engine!\n\n"
      . "use \$searchEngine = $CLASS->new( pathToEngine=>\"/usr/local/"
      . "bin/blastp\")\n"
      if ( not defined $nameValuePairs{'pathToEngine'} );

  # Create ourself as a hash
  my $this = {};

  # Bless this hash in the name of the father, the son...
  bless $this, $class;

  $this->setPathToEngine( $nameValuePairs{'pathToEngine'} );

  # TODO: Figure out a better design
  $this->setUseDustSeg( 1 );

  return $this;
}

##-------------------------------------------------------------------------##
## Get and Set Methods
##-------------------------------------------------------------------------##

##-------------------------------------------------------------------------##

=head2 get_setOverrideParameters()

  Use: my $value    = getOverrideParameters( );
  Use: my $oldValue = setOverrideParameters( $value );

  Get/Set the the override paramters.  These are used instead
  of all the SearchEngineI default parameters if set.

=cut

##-------------------------------------------------------------------------##
sub getOverrideParameters {
  my $this = shift;

  return $this->{'overrideParameters'};
}

sub setOverrideParameters {
  my $this  = shift;
  my $value = shift;

  my $oldValue = $this->{'overrideParameters'};
  $this->{'overrideParameters'} = $value;

  return $oldValue;
}

##-------------------------------------------------------------------------##

=head2 MaskLevelSequence()

  Use: my $value    = getMaskLevelSequence( );
  Use: my $oldValue = setMaskLevelSequence( $value );

  Get/Set the MaskLevelSequence paramter.  This is the
  sequence ( SearchResult::QueryStrand or SearchResult::SubjectStrand )
  which will be considered when applying the mask level.  The
  default is SearchResult::QueryStrand.

=cut

##-------------------------------------------------------------------------##
sub getMaskLevelSequence {
  my $this = shift;

  return $this->{'maskLevelSequence'};
}

sub setMaskLevelSequence {
  my $this  = shift;
  my $value = shift;

  croak $CLASS
      . "::setMaskLevelSequence: Invalid value ( $value ). "
      . "Should be either SearchResult::QueryStrand or "
      . "SearchResult::SubjectStrand\n"
      if (    $value != SearchResult::QueryStrand
           && $value != SearchResult::SubjectStrand );

  my $oldValue = $this->{'maskLevelSequence'};
  $this->{'maskLevelSequence'} = $value;

  return $oldValue;
}

##-------------------------------------------------------------------------##

=head2 UseDustSeg()

  Use: my $value    = getUseDustSeg( );
  Use: my $oldValue = setUseDustSeg( $value );

  Turn on / off the Dust/Seg screening of words.

=cut

##-------------------------------------------------------------------------##
sub getUseDustSeg {
  my $this = shift;

  return $this->{'useDustSeg'};
}

sub setUseDustSeg {
  my $this  = shift;
  my $value = shift;

  my $oldValue = $this->{'useDustSeg'};
  $this->{'useDustSeg'} = $value;

  return $oldValue;
}

##-------------------------------------------------------------------------##

=head2 TempDir()

  Use: my $value    = getTempDir( );
  Use: my $oldValue = setTempDir( $value );

  Set the directory to use as a temp directory for a search.  
  The default is to use the directory which contains the query sequence.

=cut

##-------------------------------------------------------------------------##
sub getTempDir {
  my $this = shift;

  return $this->{'tempDir'};
}

sub setTempDir {
  my $this  = shift;
  my $value = shift;

  my $oldValue = $this->{'tempDir'};
  $this->{'tempDir'} = $value;

  return $oldValue;
}

##-------------------------------------------------------------------------##

=head2 get_setAdditionalParameters()

  Use: my $value    = getAdditionalParameters( );
  Use: my $oldValue = setAdditionalParameters( $value );

  Get/Set the additional paramters.  These are used in addition
  to the existing parameter set. 

=cut

##-------------------------------------------------------------------------##
sub getAdditionalParameters {
  my $this = shift;

  return $this->{'additionalParameters'};
}

sub setAdditionalParameters {
  my $this  = shift;
  my $value = shift;

  my $oldValue = $this->{'additionalParameters'};
  $this->{'additionalParameters'} = $value;

  return $oldValue;
}

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
  if ( $result =~ /^(BLAST[PN]) (\S+) \[.*/ ) {
    $this->{'version'}    = $2;
    $this->{'engineName'} = $1;
  }
  else {
    croak $CLASS
        . "::setPathToEngine( $value ): Cannot determine "
        . "engine variant and version!\n";
  }

  my $oldValue = $this->{'pathToEngine'};
  $this->{'pathToEngine'} = $value;

  return $oldValue;
}

##-------------------------------------------------------------------------##
## Instance Methods
##-------------------------------------------------------------------------##

##-------------------------------------------------------------------------##

=head2 search()

  Use: my ( $resultCode, $SearchResultCollectionI ) = search( );
 or                                                                            
  Use: my ( $resultCode, $SearchResultCollectionI )
                          = search( matrix=>"7p16g.matrix",
                                    ...
                                  );

  Run the search and return a SearchResultCollection.

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
  my $parameters    = "";
  my $spanParameter = "";
  my $value;
  if ( ( $value = $this->getSubject() ) ) {

    # The -span1 and -span2 parameters ( from the wublast
    # manual for -span1 ):
    #  "This option relaxes the criteria for judging whether an HSP
    #   spans another, prior to discarding one of them if spanning is
    #   detected. With this option, it is merely a matter of either the query
    #   seg-ment or the database segment (or both) spans the corresponding
    #   segment(s) in the other HSP, whereas the -span2 option requires
    #   that both segments be spanned. The -span1 option may be useful in
    #   suppressing reports of HSPs when the query or a database sequence
    #   contains internal repeats."
    # MaskerAid used -span1 for dealing with simple repeats.  The
    # only way to catch that is to check for the name of the lib.
    # This is a hack and I am not very happy about it.
    # Decide when to use the span1 parameter:
    if ( $value =~ /simple.lib/ || $value =~ /at.lib/ ) {
      $spanParameter = " -span1";
    }

    # Make sure we have the compressed form of the database handy
    if ( -f "$value.ahd" || -f "$value.xns" || -f "$value.xps" ) {
      $parameters .= " $value";
    }
    else {
      croak $CLASS
          . "::search: Error...compressed subject "
          . "database ($value) does not exist!\n";
    }
  }
  else {
    croak $CLASS. "::search: Error subject undefined!\n";
  }

  my $outputDirName;
  if ( ( $value = $this->getQuery() ) ) {
    if ( -f $value ) {
      $parameters .= " $value";

      # TODO: Figure out a better design
      if ( $value =~ /simple.lib/ || $value =~ /at.lib/ ) {
        $spanParameter = " -span1";
      }
    }
    else {
      croak $CLASS. "::search: Error...query ($value) does not exist!\n";
    }
    if ( defined $this->getTempDir() && -d $this->getTempDir() ) {
      $outputDirName = $this->getTempDir();
    }
    else {
      $outputDirName = dirname( $value );
    }
  }
  else {
    croak $CLASS. "::search: Error query undefined!\n";
  }

  # Constant parameters
  if ( $this->{'engineName'} eq "BLASTP" ) {
    $parameters .= " -warnings -T=1000000 -p=1 -hspmax=2000000000 -gi ";
    $parameters .= " V=0 B=100000000";
  }
  elsif ( $this->{'engineName'} eq "BLASTN" ) {
    $parameters .= " -warnings -kap";
    if ( $this->{'useDustSeg'} >= 1 ) {
      $parameters .= " wordmask=dust wordmask=seg";
    }
    $parameters .= " maskextra=10 -hspmax=0 V=0 B=100000000";
  }
  else {
    croak $CLASS
        . "::search: Could not determine default parameters "
        . " for engineName ( "
        . $this->{'engineName'} . " )!\n";
  }
  $parameters .= $spanParameter;

  if ( defined( $value = $this->getGapInit() )
       && $value =~ /\d+/ )
  {
    $parameters .= " Q=" . abs( $value );
  }
  else {
    $parameters .= " Q=12";
  }

  if ( defined( $value = $this->getInsGapExt() )
       && $value =~ /\d+/ )
  {
    $parameters .= " R=" . abs( $value );
  }
  else {
    $parameters .= " R=2";
  }

  if ( defined( $value = $this->getMinMatch() )
       && $value =~ /\d+/ )
  {
    $parameters .= " W=$value";
  }
  else {
    $parameters .= " W=14";
  }

  #     S:  Overall score threshold...not sure how this relates to the others
  #    S2:  Score threshold for ungapped HSPs
  # gapS2:  Score threshold for gapped HSPs
  #     X:  Dropoff score for ungapped HSPs
  #  gapX:  Dropoff score for gapped HSPs
  if ( defined( $value = $this->getMinScore() )
       && $value =~ /\d+/ )
  {
    $parameters .= " S=$value";
    $parameters .= " gapS2=$value";
    $parameters .= " S2=" . int( $value / 2 );
    $parameters .= " X=$value";
    $parameters .= " gapX=" . ( $value * 2 );
  }
  else {
    $parameters .= " S=30";
    $parameters .= " gapS2=30";
    $parameters .= " S2=15";
    $parameters .= " X=30";
    $parameters .= " gapX=60";
  }

  if ( defined( $value = $this->getBandwidth() )
       && $value =~ /\d+/ )
  {
    $parameters .= " gapW=" . ( ( $value * 3 ) + 1 );
  }

  my ( $matrixRef, $matrixAlphabet, $matrixFreqRef ) = ();

  if ( defined( $value = $this->getMatrix() ) ) {

    # Test if matrix exists
    if ( -f $value ) {

      # Read in the matrix if we need to
      # adjust scores
      if ( $this->getScoreMode() == SearchEngineI::complexityAdjustedScoreMode )
      {
        ( $matrixRef, $matrixAlphabet, $matrixFreqRef ) = _readMatrix( $value );
      }

      # WUBLAST requires that the matrix filename parameter
      # be relative to a directory path specified in
      # environment variables.
      my @path   = split( /[\\\/]/, $value );
      my $matrix = pop @path;

      # WUBLAST expects the path to be above "aa" or "nt" sub
      # directories
      if ( $path[ $#path ] eq "aa" || $path[ $#path ] eq "nt" ) {
        pop @path;
      }

      # Set the environment
      $ENV{BLASTMAT}   = join( "/", @path );
      $ENV{WUBLASTMAT} = $ENV{BLASTMAT};

      $parameters .= " -matrix=$matrix";

    }
    else {
      croak $CLASS. "::search: Error...matrix ($value) does not exist!\n";
    }
  }

  # Invoke engine and handle errors
  my $POUTPUT = new FileHandle;
  my $errFile;
  do {
    $errFile = $outputDirName . "/wuResults-" . time() . ".err";
  } while ( -f $errFile );
  my $pid;
  my $runParameters;
  if ( defined $this->{'overrideParameters'}
       && $this->{'overrideParameters'} ne "" )
  {
    $runParameters = $this->{'overrideParameters'};
  }
  else {
    $runParameters = $parameters . " ";
  }

  if ( defined $this->{'additionalParameters'}
       && $this->{'additionalParameters'} ne "" )
  {
    $runParameters .= " " . $this->{'additionalParameters'};
  }

  print $CLASS
      . "::search() Invoking search engine as: $engine "
      . "$runParameters 2>$errFile |\n"
      if ( $DEBUG );

  $pid = open( $POUTPUT, "$engine $runParameters 2>$errFile |" );

  # Create SearchResultCollection object from
  # the engine results.
  my $searchResultsCollection;
  if ( $DEBUG ) {
    ## Create a debug file
    my $outFile;
    do {
      $outFile = $outputDirName . "/wuResults-" . time() . ".out";
    } while ( -f $outFile );
    open OUT, ">$outFile";
    while ( <$POUTPUT> ) {
      print OUT $_;
    }
    close OUT;
    close $POUTPUT;
    $searchResultsCollection = parseOutput( searchOutput => $outFile );
  }
  else {

    # Just pipe to the parser
    $searchResultsCollection = parseOutput( searchOutput => $POUTPUT );
  }
  close $POUTPUT;

  my $resultCode = ( $? >> 8 );
  print "WUBlast returned a the following result code >$resultCode<\n"
      if ( $DEBUG );

  # This result code can occur if wublast encounters a sequence
  # (at the end of a query collection) which is short enough
  # that it would never score higher than the threshold. An
  # exit code of 16 | 17 occurs with a printed error of
  # "FATAL:   There are no valid contexts in the requested search."
  # Look for this and ignore.  NOTE: In WUBLast 2.0 ( non free version )
  # you can use the -novalidctxok option to warn but not return a bad exit
  # code.
  if ( $resultCode != 0 ) {
    open ERR, "<$errFile";
    my $errOk = 1;
    while ( <ERR> ) {
      if ( /^FATAL:/ ) {
        if ( !/context|cpus|P option|uses P|shorter/i ) {
          $errOk = 0;
          print STDERR $CLASS . "::search: $_\n";
        }
      }
    }
    close ERR;
    $resultCode = 0 if ( $errOk == 1 );
  }
  unlink $errFile unless ( $DEBUG );

  #
  # Postprocess the results
  #
  if ( defined $searchResultsCollection
       && $searchResultsCollection->size() > 0 )
  {

    #
    # Calculate the complexity adjusted score if necessary
    #
    my $minScore = $this->getMinScore();
    if ( defined $this->getScoreMode()
        && $this->getScoreMode() == SearchEngineI::complexityAdjustedScoreMode )
    {

      # We need a defined matrix in order to calculate the complexity adjusted
      # score.  If the user overrides the parameters then we cannot be sure
      # which matrix was used.
      if ( !defined $matrixRef ) {
        croak $CLASS
            . "::search: Evidently you overrode the parameters to "
            . "this SearchEngine but tried to use complexityAdjustedScoreMode "
            . "bad bad bad!\n";
      }
      print $CLASS
          . "::search: "
          . $searchResultsCollection->size()
          . " hits "
          . "before complexity adjustment and minScore filtering\n"
          if ( $DEBUG );
      my $lambda = _calculateLambda( $matrixRef, $matrixFreqRef );

      for ( my $i = $searchResultsCollection->size() - 1 ; $i >= 0 ; $i-- ) {
        my $adjScore = _complexityAdjust(
                          $searchResultsCollection->get( $i )->getScore(),
                          $searchResultsCollection->get( $i )->getQueryString(),
                          $searchResultsCollection->get( $i )->getSubjString(),
                          $lambda,
                          $matrixAlphabet,
                          $matrixRef,
                          $matrixFreqRef
        );

        if ( defined $minScore && $adjScore < $minScore ) {
          $searchResultsCollection->remove( $i );
        }
        else {
          $searchResultsCollection->get( $i )->setScore( $adjScore );
        }
      }
      print $CLASS
          . "::search: "
          . $searchResultsCollection->size()
          . " hits "
          . "after complexity adjustment and minScore filtering\n"
          if ( $DEBUG );
    }

    #
    # Mask level filtering
    #
    my $maskLevel;
    my $strand;
    if ( defined( $maskLevel = $this->getMaskLevel() ) ) {
      if ( defined( $strand = $this->getMaskLevelSequence() ) ) {
        $searchResultsCollection->maskLevelFilter( value  => $maskLevel,
                                                   strand => $strand );
      }
      else {
        $searchResultsCollection->maskLevelFilter( value => $maskLevel );
      }
      print $CLASS
          . "::search: "
          . $searchResultsCollection->size()
          . " hits "
          . "after masklevel filtering\n"
          if ( $DEBUG );
    }

    # The final result collection should be sorted by
    #  queryname and secondarily by query start position.
    $searchResultsCollection->sort(
      sub ($$) {
        $_[ 0 ]->getQueryName cmp $_[ 1 ]->getQueryName()
            || $_[ 0 ]->getQueryStart() <=> $_[ 1 ]->getQueryStart();
      }
    );
  }

  return ( $resultCode, $searchResultsCollection );
}

##-------------------------------------------------------------------------##

=head1 Class Methods

=cut

##-------------------------------------------------------------------------##

##-------------------------------------------------------------------------##

=head2 parseOutput()

  Use: my $SearchResultCollection = WUBLASTSearchEngine::parseOutput(
                                     searchOutput => $filename|$FH,
                                     [excludeAlignments => 1],
                                     [scoreHighThresh => #],
                                     [scoreLowThresh => #],
                                     [subjPattern => ""],
                                     [queryPattern => ""] );

  Parse the result of a search and return a SearchResultCollection.

=cut

##-------------------------------------------------------------------------##
sub parseOutput {
  my %nameValueParams = @_;

  croak $CLASS. "::parseOutput() missing searchOutput parameter!\n"
      if ( !exists $nameValueParams{'searchOutput'} );

  my $WUFILE;
  if ( ref( $nameValueParams{'searchOutput'} ) !~ /GLOB|FileHandle/ ) {
    print $CLASS
        . "::parseOutput() Opening file "
        . $nameValueParams{'searchOutput'} . "\n"
        if ( $DEBUG );
    open $WUFILE, $nameValueParams{'searchOutput'}
        or die $CLASS
        . "::parseOutput: Unable to open "
        . "results file: $nameValueParams{'searchOutput'} : $!";
  }
  else {
    $WUFILE = $nameValueParams{'searchOutput'};
  }

  my $inAlignState = 0;
  my $sbjID        = "";
  my $qryID        = "";
  my $absIndex     = 0;
  my $score        = 0;
  my $adjScore     = 0;
  my $sbjSeq       = "";
  my $qrySeq       = "";
  my $sbjOrient    = "";
  my $qryOrient    = "";    # TODO: Check to see if this ever happens
  my $sbjStart     = 0;
  my $sbjEnd       = 0;
  my $qryStart     = 0;
  my $qryEnd       = 0;
  my $qryLength    = 0;
  my $sbjLength    = 0;
  my $matrix;

  my $resultColl = SearchResultCollection->new();

  while ( <$WUFILE> ) {

    #
    # Conditions for the end of a hit record:
    #   o Must have seen a score
    #   o Must be in the alignment state (ie. have seen Query: and
    #     Subj: recently)
    #   o Must see something which isn't either "Query:" "Subj:" or " "
    #     or must see the end of the file
    #
    if ( $inAlignState ) {
      if ( !/^(Query:|Sbjct:|\s{8}|\n|\r)/ || eof ) {

        #
        # Reorient if this is a reverse strand
        # hit.
        #
        my $orientation = "";
        if ( $qryOrient eq "C" ) {

          # Fix the sequence orientation so that it
          # matches the SearchResult.pm convention of
          # the query being in the forward direction.
          $qrySeq = reverse $qrySeq;
          $qrySeq =~ tr/ACGTYRMKHBVD/TGCARYKMDVBH/;    # complement
          $sbjSeq = reverse $sbjSeq;
          $sbjSeq =~ tr/ACGTYRMKHBVD/TGCARYKMDVBH/;    # complement
          $orientation = "C";
        }
        elsif ( $sbjOrient eq "C" ) {

          # Fix the subject coordinates since they
          # are based on a pre-reversed sequence.
          my $tmpEnd = $sbjEnd;
          $sbjEnd      = $sbjLength - $sbjStart + 1;
          $sbjStart    = $sbjLength - $tmpEnd + 1;
          $orientation = "C";
        }

        #
        # Calculate percent divergence
        #           percent insertions
        #           percent deletions
        #
        my %baseFreq = ();
        my $mismatch = 0;
        for ( my $i = 0 ; $i < length( $qrySeq ) ; $i++ ) {
          my $qryBase = substr( $qrySeq, $i, 1 );
          my $sbjBase = substr( $sbjSeq, $i, 1 );
          next if ( $qryBase =~ /-|x/i || $sbjBase =~ /-|x/i );
          $baseFreq{$qryBase}++;
          $mismatch++ if ( $qryBase ne $sbjBase );
        }
        my $percDiv =
            sprintf( "%4.2f", $mismatch * 100 / ( $qryEnd + 1 - $qryStart ) );
        my $qgap = $qrySeq =~ tr/-/-/;
        my $sgap = $sbjSeq =~ tr/-/-/;
        my $percIns =
            sprintf( "%4.2f", $sgap * 100 / ( ( $sbjEnd + 1 ) - $sbjStart ) );
        my $percDel =
            sprintf( "%4.2f", $qgap * 100 / ( ( $qryEnd + 1 ) - $qryStart ) );

        my $result = SearchResult->new(
                                     queryName      => $qryID,
                                     queryStart     => $qryStart,
                                     queryEnd       => $qryEnd,
                                     queryRemaining => ( $qryLength - $qryEnd ),
                                     queryString    => $qrySeq,
                                     subjString     => $sbjSeq,
                                     orientation    => $orientation,
                                     subjName       => $sbjID,
                                     subjStart      => $sbjStart,
                                     subjEnd        => $sbjEnd,
                                     subjRemaining  => ( $sbjLength - $sbjEnd ),
                                     queryString    => $qrySeq,
                                     pctDiverge     => $percDiv,
                                     pctInsert      => $percIns,
                                     pctDelete      => $percDel,
                                     matrixName     => $matrix,
                                     score          => $score
        );

        $resultColl->add( $result );

        $score        = "";
        $sbjSeq       = "";
        $qrySeq       = "";
        $qryOrient    = "";
        $sbjStart     = 0;
        $sbjEnd       = 0;
        $qryStart     = 0;
        $qryEnd       = 0;
        $inAlignState = 0;
      }
    }

    #
    # Query alignment
    #
    if ( /Query:\s+(\d+)\s+(\S+)\s+(\d+)/ ) {
      $qrySeq .= $2;
      if ( $qryStart < 1 ) {
        $qryStart = _min( $1, $3 );
        $qryEnd   = _max( $1, $3 );
      }
      else {
        $qryStart = _min( _min( $1, $3 ), $qryStart );
        $qryEnd   = _max( _max( $1, $3 ), $qryEnd );
      }
      $qryOrient = "C" if ( $1 > $3 );
      $inAlignState = 1;
    }

    #
    # Subject alignment
    #
    if ( /Sbjct:\s+(\d+)\s+(\S+)\s+(\d+)/ ) {
      $sbjSeq .= $2;
      my $leftNum  = $1;
      my $rightNum = $3;
      if ( $leftNum > $rightNum ) {
        if ( ( $leftNum - $rightNum ) == 1 ) {

          # There is currently a bug in blast.  It occasionally
          # mis-indexes the alignment data so that the last
          # line contains an index which is descending.
          $leftNum  = $3;
          $rightNum = $1;
        }
        else {
          croak $CLASS
              . "::parseOutput: Strange blast alignment format.  "
              . "The subject sequence should not contain descending "
              . "indexes: ( $_ )!\n";
        }
      }
      if ( $sbjStart < 1 ) {
        $sbjStart = _min( $leftNum, $rightNum );
        $sbjEnd   = _max( $leftNum, $rightNum );
      }
      else {
        $sbjStart = _min( _min( $leftNum, $rightNum ), $sbjStart );
        $sbjEnd   = _max( _max( $leftNum, $rightNum ), $sbjEnd );
      }
      $inAlignState = 1;
    }

    #
    # Query source name
    #
    if ( /^Query=\s+(\S+)/ ) {
      $qryID = $1;
    }

    #
    # Query length
    #
    if ( /^\s+\(([,\d]+) letters/ ) {
      $qryLength = $1;
      $qryLength =~ s/,//g;
    }

    #
    # TODO: Database source name
    #
    #$this->{databaseName} = $1 if ( /Database:\s+(\S+)/ );

    #
    # Length of database ( can have comas in it! )
    #
    if ( /^\s+Length = ([\d,]+)\s*$/ ) {
      $sbjLength = $1;
      $sbjLength =~ s/,//g;
    }

    #
    # Score
    #
    if ( /Score = (\d+)/ ) {
      $score = $1;
    }

    #
    # Hit description line
    #
    if ( /^>(.*)/ ) {
      $sbjID = $1;
      ## Take note!  This is a very special case here.
      ## In order to search DNA with an expanded alphabet
      ## we have used blastp to search the DNA sequence
      ## as if it were amino acids.  Since blastp does
      ## not search both strands we seed the database with
      ## both the forward and reverse strands labeling the
      ## reverse copies with the string "(anti)".  If this
      ## is located in the results treat this hit as
      ## occuring on the revsere complement of the database
      ## entry.
      if ( $sbjID =~ /anti/ ) {
        $sbjOrient = "C";
      }
      else { $sbjOrient = ""; }
      ( $sbjID ) = $sbjID =~ /(\S+)/;
    }

    #
    # Grab the matrix from the parameter list
    #
    if ( /matrix=\s*(\S+)/ ) {
      my @path = split( /[\\\/]/ );
      $matrix = $path[ $#path ];
      $matrix =~ s/[\n\r]//;

      # TODO: Go back and set matrix for all previous hits
    }

  }
  close $WUFILE;

  return $resultColl;
}

##-------------------------------------------------------------------------##
## Private Methods
##-------------------------------------------------------------------------##

##-------------------------------------------------------------------------##
## Use: my _complexityAdjust( $swnScore, $qrySeq, $sbjSeq, $matLambda,
##                            $matAlphabet, $matScoresRef, $matFreqsRef );
##
##      $swnScore    : Smith-Waterman Raw Alignment Score
##      $qrySeq      : Query string from the alignment
##      $sbjSeq      : Subject string from the alignment
##      $matLambda   : Matrix Lambda parameter
##      $matAlphabet : Matrix alphabet string (in column order)
##      $matScoresRef: Score Matrix
##      $matFreqsRef : Matrix alphabet freq vector
##
## Returns
##      $adj_score : The complexity adjusted score according to
##                   to Phil Green's swat/cross_match program.
##-------------------------------------------------------------------------##
sub _complexityAdjust {
  my ( $swnScore, $qrySeq, $sbjSeq, $matLambda, $matAlphabet, $matScoresRef,
       $matFreqsRef )
      = @_;
  my $mCol      = 0;
  my $t_factor  = 0;
  my $t_sum     = 0;
  my $t_counts  = 0;
  my $n_letters = 0;
  my $baseIndex = 0;
  my $qBase     = "";
  my $sBase     = "";
  my $pos_score = 0;
  my @matCounts = ();
  my $adj_score = 0;

  #
  # Recalculates the raw score ( $pos_score ) based on the matrix we
  # have been handed.
  #
  # Creates a vector of base counts from the query sequence.  It
  # ignores base instances when they are part of a deletion or
  # are insertion characters "-".
  #
  for ( $baseIndex = 0 ; $baseIndex < length( $qrySeq ) ; $baseIndex++ ) {
    $qBase = substr( $qrySeq, $baseIndex, 1 );
    $sBase = substr( $sbjSeq, $baseIndex, 1 );
    if ( $qBase ne "-" && $sBase ne "-" ) {    ## && $qBase eq $sBase) {
      $pos_score +=
          $$matScoresRef[ index( $matAlphabet, $qBase ) ]
          [ index( $matAlphabet, $sBase ) ];
      $matCounts[ index( $matAlphabet, $qBase ) ]++
          if ( ( index $matAlphabet, $qBase ) >= 0 );
      ## TODO: TEST: DO NOT LEAVE IN PRODUCTION CODE!
      ##$matCounts[ index( $matAlphabet, $sBase ) ]++
      ##    if ( ( index $matAlphabet, $sBase ) >= 0 );
      ## TODO: TEST
    }
  }

  #
  #
  #
  for ( $mCol = 0 ; $mCol < length( $matAlphabet ) ; $mCol++ ) {
    if ( defined $matCounts[ $mCol ] && $matCounts[ $mCol ] > 0 ) {
      if ( $$matFreqsRef[ $mCol ] > 0 && log( $$matFreqsRef[ $mCol ] ) != 0 ) {
        my $count = $matCounts[ $mCol ];
        ## TODO: TEST: DO NOT LEAVE IN PRODUCTION CODE!
        ##my $count = $matCounts[ $mCol ] / 2;
        ## TODO: TEST
        $t_factor += $count * log( $count );
        $t_sum    += $count * log( $$matFreqsRef[ $mCol ] );
        $t_counts += $count;
        $n_letters++;
      }
    }
  }

  #
  #
  #
  $t_factor -= $t_counts * log( $t_counts );
  $t_sum    -= $t_factor;

  #
  # Looks like Phil changed his mind here
  #
  #my $complexity_factor = 0.25;
  #if ( $n_letters > ( 1 / $complexity_factor )) {
  #  $t_factor /= $t_counts * log(1 / $n_letters);
  #}else {
  #  $t_factor /= $t_counts * log($complexity_factor);
  #}
  #
  #$old_adj_score = $swnScore + ($t_factor * $pos_score) - $pos_score + .5;

  #
  #
  #
  $adj_score = sprintf( "%0.0d", $swnScore + $t_sum / $matLambda + .999 );

  $adj_score = 0 if ( !( $adj_score =~ /\d+/ ) || $adj_score < 0 );
  return ( $adj_score );
}

##-------------------------------------------------------------------------##
## Use: my _calculateLambda( $matScoresRef, $matFreqsRef );
##
##      $matScoresRef: Score Matrix
##      $matFreqsRef : Matrix alphabet freq vector
##
## Returns
##      $lambda : The lambda parameter derived from the matrix
##                and the matrix alphabet frequencies.  This
##                is derived from Phil Green's swat/cross_match
##                programs.
##-------------------------------------------------------------------------##
sub _calculateLambda {
  my ( $matScoresRef, $matFreqsRef ) = @_;
  my $lambda_upper = 0;
  my $lambda_lower = 0;
  my $lambda       = 0.5;
  my $sum          = 0;

  do {
    $sum = _getSum( $lambda, $matScoresRef, $matFreqsRef );
    if ( $sum < 1.0 ) {
      $lambda_lower = $lambda;
      $lambda *= 2.0;
    }
  } while ( $sum < 1.0 );

  $lambda_upper = $lambda;

  while ( $lambda_upper - $lambda_lower > .00001 ) {
    $lambda = ( $lambda_lower + $lambda_upper ) / 2.0;
    $sum = _getSum( $lambda, $matScoresRef, $matFreqsRef );
    if ( $sum >= 1.0 ) {
      $lambda_upper = $lambda;
    }
    else {
      $lambda_lower = $lambda;
    }
  }
  return ( $lambda );
}

##-------------------------------------------------------------------------##
## Use: my _getSum( $lambda, $matScoresRef, $matFreqsRef );
##
##      $matScoresRef: Score Matrix
##      $matFreqsRef : Matrix alphabet frequency vector
##
## Returns
##      $sum : Good question....???  This is used by the
##             lambda estimation procedure.  _getSum is
##             derived from Phil Green's swat/cross_match
##             programs.
##-------------------------------------------------------------------------##
sub _getSum {
  my ( $lambda, $matScoresRef, $matFreqsRef ) = @_;
  my $check = 0;
  my $total = 0;
  my $i     = 0;
  my $j     = 0;

  for ( $i = 0 ; $i <= $#$matFreqsRef ; $i++ ) {
    for ( $j = 0 ; $j <= $#$matFreqsRef ; $j++ ) {
      if ( $$matFreqsRef[ $i ] && $$matFreqsRef[ $j ] ) {
        $total += $$matFreqsRef[ $i ] * $$matFreqsRef[ $j ] *
            exp( $lambda * $$matScoresRef[ $i ][ $j ] );
        $check += $$matFreqsRef[ $i ] * $$matFreqsRef[ $j ];
      }
    }
  }
  die "error in _getSum!! check=$check\n"
      if (    $check > 1.001
           || $check < .999 );

  return ( $total );
}

## TODO: Use SequenceSimilarityMatrix.pm here
##-------------------------------------------------------------------------##
## Use: my ( \@matrix, $alphabet, \@freqArray ) =
##                                        _readMatrix( $matrixFileName );
##
##      $matFileName : Name of the file containing the WUBlast Matrix
##
## Returns
##      This is a very specific matrix reader for WUBlast matrices
##      which have been created for RepeatMasker.  This routine
##      assumes that the standard crossmatch matrix frequency line
##      is included in the comment section. Ie.:
##
##        # FREQS A 0.325 C 0.175 G 0.175 T 0.325
##
##      The routine returns the matrix stored as a 2-d array.  It also
##      returns the alphabet as a compressed string and the frequencies
##      as an array ( sizeof alphabet ).
##
##-------------------------------------------------------------------------##
sub _readMatrix {
  my ( $matrixFileName ) = @_;
  my %freqHash           = ();
  my @matrix             = ();
  my @rowValues          = ();
  my $row                = 0;
  my $alphabet           = "";
  my $i                  = 0;
  my @freqArray          = ();

  open MATRIX, "<$matrixFileName" || die "Can't open $matrixFileName!\n";

  while ( <MATRIX> ) {
    $_ = uc( $_ );
    chomp;
    if ( /FREQS\s+(.*)/ ) {
      %freqHash = split " ", $1;
    }
    if ( /^\s*[A-Z]\s+[A-Z]\s+[A-Z]\s+[A-Z]\s+/ ) {
      s/ //g;
      $alphabet = $_;
    }
    elsif ( $alphabet
            && /^\s*\S\s+[\d-]+\s+[\d-]+\s+[\d-]+\s+[\d-]+\s+/ )
    {
      @rowValues = split;
      shift @rowValues;
      $matrix[ $row++ ] = [ @rowValues ];
    }
  }
  close MATRIX;

  if ( scalar( keys %freqHash ) == 0 ) {
    $freqHash{"A"} = 0.25;
    $freqHash{"C"} = 0.25;
    $freqHash{"G"} = 0.25;
    $freqHash{"T"} = 0.25;
  }

  for ( $i = 0 ; $i < length( $alphabet ) ; $i++ ) {
    $freqArray[ $i ] = $freqHash{ substr( $alphabet, $i, 1 ) } || 0;
  }

  return ( \@matrix, $alphabet, \@freqArray );
}

##-------------------------------------------------------------------------##
## Use: my _min( $num1, $num2 );
##
##              $num1   :       A number to be compared
##              $num2   :       A number to be comprared
##
##      Returns:                The minimum of the two numbers
##
##-------------------------------------------------------------------------##
sub _min {
  my ( $num1, $num2 ) = @_;
  if ( $num1 < $num2 ) {
    return ( $num1 );
  }
  else {
    return ( $num2 );
  }
}

##-------------------------------------------------------------------------##
## Use: my _max( $num1, $num2 );
##
##              $num1   :       A number to be compared
##              $num2   :       A number to be comprared
##
##      Returns:                The maximum of the two numbers
##
##-------------------------------------------------------------------------##
sub _max {
  my ( $num1, $num2 ) = @_;
  if ( $num1 < $num2 ) {
    return ( $num2 );
  }
  else {
    return ( $num1 );
  }
}

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
