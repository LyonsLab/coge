#!/usr/local/bin/perl -w
##---------------------------------------------------------------------------##
##  File:
##      @(#) SearchResult.pm
##  Author:
##      Robert M. Hubley   rhubley@systemsbiology.org
##  Description:
##      An for holding a generic biological sequence similarity
##      search result.
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
# bless( {
#          'qryGCBackground' => '43',
#          'querySeq' => 'AGCAA..TGTAAA',
#          'sbjEnd' => '2760',
#          'qryEnd' => '319032',
#          'qryBegin' => '318751',
#          'percDel' => '4.26',
#          'sbjBegin' => '2484',
#          'qryName' => 'ctg12382',
#          'score' => '1005',
#          'sbjName' => 'Charlie1',
#          'percIns' => '6.03',
#          'sbjOrient' => 'C',
#          'qryLeft' => '4582070',
#          'sbjLeft' => '1',
#          'id' => '5',
#          'uniqId' => '1383820',
#          'percDiv' => '16.31',
#          'subjSeq' => 'AGCGGT...AA',
#          'matrix' => '35p40g.matrix',
#          'transV' => '3',
#          'transI' => '1',
#          'overlap' => '*',
#        }
#     }, 'SearchResult' );
#
###############################################################################
# ChangeLog
#
#     $Log: SearchResult.pm,v $
#     Revision 1.47  2007/05/17 21:01:54  rhubley
#     Cleanup before a distribution
#
#
###############################################################################
# To Do:
#    - Expose overlap, transitions and transversions attributes!
#
#

=head1 NAME

SearchResult

=head1 SYNOPSIS

use SearchResult

Usage: 

    $SearchResultsCollection = SearchResult->new();

  or 

    $SearchResultsCollection = SearchResult->new( 
                                 queryName=>value, subjName=>value,
                                 pctInsert=>value, pctDelete=>value, 
                                 queryStart=>value, queryEnd=>value, 
                                 score=>value, pctDiverge=>value, 
                                 subjRemaining=>value, subjType=>value,
                                 queryRemaining=>value, id=>value,
                                 orientation=>value, queryString=>value,
                                 subjString=>value, matrix=>value,
                                 id=>value, lineageId=>value );

=head1 DESCRIPTION

A class for storing a result from the crossmatch search engine.

=head1 SEE ALSO

=over 4

SearchResultCollection

=back

=head1 COPYRIGHT

Copyright 2004 Institute for Systems Biology

=head1 AUTHOR

Robert Hubley <rhubley@systemsbiology.org>

=head1 INSTANCE METHODS

=cut 

package SearchResult;
use strict;
use Data::Dumper;
use Carp;
use vars qw(@ISA @EXPORT @EXPORT_OK %EXPORT_TAGS $VERSION);

use constant NoAlign => 1;

# Indicates that alignments should always be reported with
# the query sequences all in the forward direction
use constant AlignWithQuerySeq => 2;

# Indicates that alignments should always be reported with
# the subj sequences all in the forward direction
use constant AlignWithSubjSeq   => 3;
use constant OutFileFormat      => 4;
use constant CompressedAlignCSV => 5;
use constant PSL                => 6;

use constant QueryStrand   => 1;
use constant SubjectStrand => 2;

require Exporter;

@ISA = qw(Exporter);

@EXPORT = qw();

@EXPORT_OK = qw();

%EXPORT_TAGS = ( all => [ @EXPORT_OK ] );

my $CLASS = "SearchResult";

##-------------------------------------------------------------------------##
## Constructor:
##-------------------------------------------------------------------------##
sub new {
  my $class          = shift;
  my %nameValuePairs = @_;

  # Create ourself as a hash
  my $this = {};

  # Bless this hash in the name of the father, the son...
  bless $this, $class;

  # Allow import of values
  if ( %nameValuePairs ) {
    while ( my ( $name, $value ) = each( %nameValuePairs ) ) {
      my $method = "set" . _ucFirst( $name );
      unless ( $this->can( $method ) ) {
        croak(
             "SearchResult::add: Instance variable $name doesn't exist." . "" );
      }
      $this->$method( $value );
    }
  }

  return $this;
}

##-------------------------------------------------------------------------##

=head2 clone()

  Use: my $newObj = $obj->clone();

  Clone a SearchResult *duplicating* all the values of the old
  object in the new one.

=cut

##-------------------------------------------------------------------------##
sub clone {
  my $this = shift;

  my %newHash = %{$this};
  my $newObj  = \%newHash;

  bless $newObj, ref( $this );

  return $newObj;
}

##-------------------------------------------------------------------------##
## Get and Set Methods
##-------------------------------------------------------------------------##

##-------------------------------------------------------------------------##

=head2 get_setMatrixName()

  Use: my $value    = getMatrixName( );
  Use: my $oldValue = setMatrixName( $value );

  Get/Set the name of the matrix.

=cut

##-------------------------------------------------------------------------##
sub getMatrixName {
  my $obj = shift;

  my $value = $obj->{'matrixName'};

  return $value;
}

sub setMatrixName {
  my $obj      = shift;
  my $value    = shift;
  my $oldValue = undef;

  $oldValue = $obj->{'matrixName'};
  $obj->{'matrixName'} = $value;

  return $oldValue;
}

##-------------------------------------------------------------------------##

=head2 get_setLineageId()

  Use: my $value    = getLineageId( );
  Use: my $oldValue = setLineageId( $value );

  Get/Set the lineage id.

=cut

##-------------------------------------------------------------------------##
sub getLineageId {
  my $obj = shift;

  my $value = $obj->{'lineageId'};

  return $value;
}

sub setLineageId {
  my $obj      = shift;
  my $value    = shift;
  my $oldValue = undef;

  $oldValue = $obj->{'lineageId'};
  $obj->{'lineageId'} = $value;

  return $oldValue;
}

##-------------------------------------------------------------------------##

=head2 get_setQueryName()

  Use: my $value    = getQueryName();
  Use: my $oldValue = setQueryName( $value );

  Get/Set the name of the query sequence.

=cut

##-------------------------------------------------------------------------##
sub getQueryName {
  my $obj = shift;

  my $value = $obj->{'qryName'};

  return $value;
}

sub setQueryName {
  my $obj      = shift;
  my $value    = shift;
  my $oldValue = undef;

  $oldValue = $obj->{'qryName'};
  $obj->{'qryName'} = $value;

  return $oldValue;
}

##-------------------------------------------------------------------------##

=head2 get_setSubjName()

  Use: my $value    = getSubjName();
  Use: my $oldValue = setSubjName( $value );

  Get/Set the name of the subject sequence.

=cut

##-------------------------------------------------------------------------##
sub getSubjName {
  my $obj = shift;

  my $value = $obj->{'sbjName'};

  return $value;
}

sub setSubjName {
  my $obj      = shift;
  my $value    = shift;
  my $oldValue = undef;

  $oldValue = $obj->{'sbjName'};
  $obj->{'sbjName'} = $value;

  return $oldValue;
}

##-------------------------------------------------------------------------##

=head2 get_setSubjType()

  Use: my $value    = getSubjType();
  Use: my $oldValue = setSubjType( $value );

  Get/Set the type of the subject sequence ( for out files only ).

=cut

##-------------------------------------------------------------------------##
sub getSubjType {
  my $obj = shift;

  my $value = $obj->{'sbjType'};

  return $value;
}

sub setSubjType {
  my $obj      = shift;
  my $value    = shift;
  my $oldValue = undef;

  $oldValue = $obj->{'sbjType'};
  $obj->{'sbjType'} = $value;

  return $oldValue;
}

##-------------------------------------------------------------------------##

=head2 get_setPctInsert()

  Use: my $value    = getPctInsert();
  Use: my $oldValue = getPctInsert( $value );

  Get/Set the percent insert value.

=cut

##-------------------------------------------------------------------------##
sub getPctInsert {
  my $obj = shift;

  my $value = $obj->{'percIns'};

  return $value;
}

sub setPctInsert {
  my $obj      = shift;
  my $value    = shift;
  my $oldValue = undef;

  $oldValue = $obj->{'percIns'};
  $obj->{'percIns'} = $value;

  return $oldValue;
}

##-------------------------------------------------------------------------##

=head2 get_setPctDelete()

  Use: my $value    = getPctDelete();
  Use: my $oldValue = setPctDelete( $value );

  Get/Set the percent deletion value.

=cut

##-------------------------------------------------------------------------##
sub getPctDelete {
  my $obj = shift;

  my $value = $obj->{'percDel'};

  return $value;
}

sub setPctDelete {
  my $obj      = shift;
  my $value    = shift;
  my $oldValue = undef;

  $oldValue = $obj->{'percDel'};
  $obj->{'percDel'} = $value;

  return $oldValue;
}

##-------------------------------------------------------------------------##

=head2 get_setSubjStart()

  Use: my $value    = getSubjStart();
  Use: my $oldValue = setSubjStart( $value );

  Get/Set the subject start position (1 based).

=cut

##-------------------------------------------------------------------------##
sub getSubjStart {
  my $obj = shift;

  my $value = $obj->{'sbjBegin'};

  return $value;
}

sub setSubjStart {
  my $obj      = shift;
  my $value    = shift;
  my $oldValue = undef;

  $oldValue = $obj->{'sbjBegin'};
  $obj->{'sbjBegin'} = $value;

  return $oldValue;
}

##-------------------------------------------------------------------------##

=head2 get_setQueryStart()

  Use: my $value    = getQueryStart();
  Use: my $oldValue = setQueryStart( $value );

  Get/Set the query start position (1 based).

=cut

##-------------------------------------------------------------------------##
sub getQueryStart {
  my $obj = shift;

  my $value = $obj->{'qryBegin'};

  return $value;
}

sub setQueryStart {
  my $obj      = shift;
  my $value    = shift;
  my $oldValue = undef;

  $oldValue = $obj->{'qryBegin'};
  $obj->{'qryBegin'} = $value;

  return $oldValue;
}

##-------------------------------------------------------------------------##

=head2 get_setQueryEnd()

  Use: my $value    = getQueryEnd();
  Use: my $oldValue = setQueryEnd( $value );

  Get/Set the query end position (1 based).

=cut

##-------------------------------------------------------------------------##
sub getQueryEnd {
  my $obj = shift;

  my $value = $obj->{'qryEnd'};

  return $value;
}

sub setQueryEnd {
  my $obj      = shift;
  my $value    = shift;
  my $oldValue = undef;

  $oldValue = $obj->{'qryEnd'};
  $obj->{'qryEnd'} = $value;

  return $oldValue;
}

##-------------------------------------------------------------------------##

=head2 get_setSubjEnd()

  Use: my $value    = getSubjEnd();
  Use: my $oldValue = setSubjEnd( $value );

  Get/Set the subject end position (1 based).

=cut

##-------------------------------------------------------------------------##
sub getSubjEnd {
  my $obj = shift;

  my $value = $obj->{'sbjEnd'};

  return $value;
}

sub setSubjEnd {
  my $obj      = shift;
  my $value    = shift;
  my $oldValue = undef;

  $oldValue = $obj->{'sbjEnd'};
  $obj->{'sbjEnd'} = $value;

  return $oldValue;
}

##-------------------------------------------------------------------------##

=head2 get_setScore()

  Use: my $value    = getScore();
  Use: my $oldValue = setScore( $value );

  Get/Set the score value.

=cut

##-------------------------------------------------------------------------##
sub getScore {
  my $obj = shift;

  my $value = $obj->{'score'};

  return $value;
}

sub setScore {
  my $obj      = shift;
  my $value    = shift;
  my $oldValue = undef;

  $oldValue = $obj->{'score'};
  $obj->{'score'} = $value;

  return $oldValue;
}

##-------------------------------------------------------------------------##

=head2 get_setPValue()

  Use: my $value    = getPValue();
  Use: my $oldValue = setPValue( $value );

  Get/Set the PValue.

=cut

##-------------------------------------------------------------------------##
sub getPValue {
  my $obj = shift;

  my $value = $obj->{'PValue'};

  return $value;
}

sub setPValue {
  my $obj      = shift;
  my $value    = shift;
  my $oldValue = undef;

  $oldValue = $obj->{'PValue'};
  $obj->{'PValue'} = $value;

  return $oldValue;
}

##-------------------------------------------------------------------------##

=head2 get_setPctDiverge()

  Use: my $value    = getPctDiverge();
  Use: my $oldValue = setPctDiverge( $value );

  Get/Set the percent divergence value.

=cut

##-------------------------------------------------------------------------##
sub getPctDiverge {
  my $obj = shift;

  my $value = $obj->{'percDiv'};

  return $value;
}

sub setPctDiverge {
  my $obj      = shift;
  my $value    = shift;
  my $oldValue = undef;

  $oldValue = $obj->{'percDiv'};
  $obj->{'percDiv'} = $value;

  return $oldValue;
}

##-------------------------------------------------------------------------##

=head2 get_setSubjRemaining()

  Use: my $value    = getSubjRemaining();
  Use: my $oldValue = setSubjRemaining( $value );

  Get/Set the subject remaining length value.

=cut

##-------------------------------------------------------------------------##
sub getSubjRemaining {
  my $obj = shift;

  my $value = $obj->{'sbjLeft'};

  return $value;
}

sub setSubjRemaining {
  my $obj      = shift;
  my $value    = shift;
  my $oldValue = undef;

  $oldValue = $obj->{'sbjLeft'};
  $obj->{'sbjLeft'} = $value;

  return $oldValue;
}

##-------------------------------------------------------------------------##

=head2 get_setQueryRemaining()

  Use: my $value    = getQueryRemaining();
  Use: my $oldValue = setQueryRemaining( $value );

  Get/Set the query remaining length value.

=cut

##-------------------------------------------------------------------------##
sub getQueryRemaining {
  my $obj = shift;

  my $value = $obj->{'qryLeft'};

  return $value;
}

sub setQueryRemaining {
  my $obj      = shift;
  my $value    = shift;
  my $oldValue = undef;

  $oldValue = $obj->{'qryLeft'};
  $obj->{'qryLeft'} = $value;

  return $oldValue;
}

##-------------------------------------------------------------------------##

=head2 get_setOverlap()

  Use: my $value    = getOverlap();
  Use: my $oldValue = setOverlap( $value );

  Get/Set the the overlap value.

=cut

##-------------------------------------------------------------------------##
sub getOverlap {
  my $obj = shift;

  my $value = $obj->{'overlap'};

  return $value;
}

sub setOverlap {
  my $obj      = shift;
  my $value    = shift;
  my $oldValue = undef;

  $oldValue = $obj->{'overlap'};
  $obj->{'overlap'} = $value;

  return $oldValue;
}

##-------------------------------------------------------------------------##

=head2 get_setId()

  Use: my $value    = getId();
  Use: my $oldValue = setId( $value );

  Get/Set the the identifier value.

=cut

##-------------------------------------------------------------------------##
sub getId {
  my $obj = shift;

  my $value = $obj->{'id'};

  return $value;
}

sub setId {
  my $obj      = shift;
  my $value    = shift;
  my $oldValue = undef;

  $oldValue = $obj->{'id'};
  $obj->{'id'} = $value;

  return $oldValue;
}

##-------------------------------------------------------------------------##

=head2 get_setOrientation()

  Use: my $value    = getOrientation();
  Use: my $oldValue = setOrientation( $value );

  Get/Set the orientation of the sequence.  The orientation is 
  interpreted as the orientation of the subject sequence.
  The query is always assumed to be in the forward direction.

=cut

##-------------------------------------------------------------------------##
sub getOrientation {
  my $obj = shift;

  my $value = $obj->{'sbjOrient'};

  return $value;
}

sub setOrientation {
  my $obj      = shift;
  my $value    = shift;
  my $oldValue = undef;

  $oldValue = $obj->{'sbjOrient'};
  $obj->{'sbjOrient'} = $value;

  return $oldValue;
}

##-------------------------------------------------------------------------##

=head2 get_setQueryString()

  Use: my $value    = getQueryString();
  Use: my $oldValue = setQueryString( $value );

  Get/Set the query portion of the alignment string.

=cut

##-------------------------------------------------------------------------##
sub getQueryString {
  my $obj = shift;

  my $value = $obj->{'querySeq'};

  return $value;
}

sub setQueryString {
  my $obj      = shift;
  my $value    = shift;
  my $oldValue = undef;

  $oldValue = $obj->{'querySeq'};
  $obj->{'querySeq'} = $value;

  return $oldValue;
}

##-------------------------------------------------------------------------##

=head2 get_setSubjString()

  Use: my $value    = getSubjString();
  Use: my $oldValue = setSubjString( $value );

  Get/Set the subject portion of the alignment string.

=cut

##-------------------------------------------------------------------------##
sub getSubjString {
  my $obj = shift;

  my $value = $obj->{'subjSeq'};

  return $value;
}

sub setSubjString {
  my $obj      = shift;
  my $value    = shift;
  my $oldValue = undef;

  $oldValue = $obj->{'subjSeq'};
  $obj->{'subjSeq'} = $value;

  return $oldValue;
}

##-------------------------------------------------------------------------##
## General Object Methods
##-------------------------------------------------------------------------##

##-------------------------------------------------------------------------##

=head2

  Use: toStringFormatted( $format );

    $format        : SearchResult::NoAlign
                     SearchResult::AlignWithQuerySeq
                     SearchResult::AlignWithSubjSeq
                     SearchResult::OutFileFormat

  Create an string representation of a single result.

=cut

##-------------------------------------------------------------------------##
sub toStringFormatted {
  my $obj    = shift;
  my $format = shift;

  $format = SearchResult::NoAlign
      if ( !defined $format );

  if (    $format == SearchResult::NoAlign
       || $format == SearchResult::AlignWithQuerySeq
       || $format == SearchResult::AlignWithSubjSeq )
  {
    return $obj->_toCrossMatchFormat( $format );
  }
  elsif ( $format == SearchResult::OutFileFormat ) {
    return $obj->_toOUTFileFormat();
  }
  elsif ( $format == SearchResult::CompressedAlignCSV ) {
    return $obj->_toCSVFormat();
  }
  else {
    croak $CLASS . "::toStringFormatted: Unknown format " . "( $format )\n";
  }
}

##
sub _toCrossMatchFormat {
  my $obj           = shift;
  my $alignmentMode = shift;

  $alignmentMode = SearchResult::NoAlign
      if ( !defined $alignmentMode );
  croak $CLASS
      . "::toStringFormatted: Unknown alignment mode "
      . "( $alignmentMode )\n"
      if (    $alignmentMode != SearchResult::NoAlign
           && $alignmentMode != SearchResult::AlignWithQuerySeq
           && $alignmentMode != SearchResult::AlignWithSubjSeq );

  my $retStr = "";
  my $sbjDir;
  my $sbjDirStr;
  my $alignIndex = 0;
  my $alignCol;

  my ( $qryName ) = ( $obj->{qryName} =~ /(\S+).*/ );
  $retStr .=
        "$obj->{score} $obj->{percDiv} $obj->{percDel} "
      . "$obj->{percIns} $qryName $obj->{qryBegin} "
      . "$obj->{qryEnd} ($obj->{qryLeft}) ";
  my ( $sbjName ) = ( $obj->{sbjName} =~ /(\S+).*/ );
  if ( $obj->{sbjOrient} eq "C" ) {
    $retStr .=
        "C $sbjName ($obj->{sbjLeft}) " . "$obj->{sbjEnd} $obj->{sbjBegin}";
  }
  else {
    $retStr .=
        "$sbjName $obj->{sbjBegin} $obj->{sbjEnd} " . "($obj->{sbjLeft})";
  }
  if ( defined $obj->{id} ) {
    $retStr .= " $obj->{id}";
  }
  if ( defined $obj->{lineageId} ) {
    $retStr .= " $obj->{lineageId}";
  }
  if ( defined $obj->{overlap} ) {
    $retStr .= " $obj->{overlap}";
  }
  $retStr .= "\n";

  if ( $alignmentMode != SearchResult::NoAlign
       && defined $obj->{querySeq} )
  {
    my $qMasked;
    my $sMasked;
    my $insertions = 0;
    my $deletions  = 0;

    my $query   = $obj->{querySeq};
    my $subject = $obj->{subjSeq};

    $retStr .= "\n";

    if (    $obj->{sbjOrient} eq "C"
         && $alignmentMode == SearchResult::AlignWithSubjSeq )
    {
      $query = reverse $query;
      $query =~ tr/ACGTYRMKHBVD/TGCARYKMDVBH/;    # complement
      $subject = reverse $subject;
      $subject =~ tr/ACGTYRMKHBVD/TGCARYKMDVBH/;    # complement
    }

    my $qStart = $obj->{qryBegin};
    if (    $obj->{sbjOrient} eq "C"
         && $alignmentMode == SearchResult::AlignWithSubjSeq )
    {
      $qStart = $obj->{qryEnd};
    }
    my $qEnd   = 0;
    my $sStart = $obj->{sbjBegin};
    if (    $obj->{sbjOrient} eq "C"
         && $alignmentMode == SearchResult::AlignWithQuerySeq )
    {
      $sStart = $obj->{sbjEnd};
    }
    my $sEnd;
    while ( $query ) {
      $query =~ s/^(.{1,50})//;
      my $qSeq = $1;
      $subject =~ s/^(.{1,50})//;
      my $sSeq = $1;

      $insertions = ( $qSeq =~ tr/-/-/ );
      $deletions  = ( $sSeq =~ tr/-/-/ );

      if (    $obj->{'sbjOrient'} eq "C"
           && $alignmentMode == SearchResult::AlignWithSubjSeq )
      {
        $qEnd = $qStart - length( $qSeq ) + 1 + $insertions;
        $retStr .= "C ";
      }
      else {
        $qEnd = $qStart + length( $qSeq ) - 1 - $insertions;
        $retStr .= "  ";
      }
      $qEnd = $qStart if ( length( $qSeq ) == $insertions );

      $retStr .= substr( $obj->{'qryName'}, 0, 13 )
          . " " x (
           13 - (
             length( $obj->{'qryName'} ) < 13 ? length( $obj->{'qryName'} ) : 13
           )
          );
      $retStr .= " " x ( 9 - length( $qStart ) ) . $qStart . " $qSeq $qEnd\n";
      $retStr .= " " x 25;
      for ( my $j = 0 ; $j < length( $qSeq ) ; $j++ ) {
        my $qChar = substr( $qSeq, $j, 1 );
        my $sChar = substr( $sSeq, $j, 1 );
        if ( $qChar eq $sChar ) {
          $retStr .= " ";
        }
        elsif ( $qChar eq "-" || $sChar eq "-" ) {
          $retStr .= "-";
        }
        elsif (    $qChar . $sChar eq "CT"
                || $qChar . $sChar eq "TC"
                || $qChar . $sChar eq "AG"
                || $qChar . $sChar eq "GA" )
        {
          $retStr .= "i";
        }
        elsif (    $qChar . $sChar eq "GT"
                || $qChar . $sChar eq "TG"
                || $qChar . $sChar eq "GC"
                || $qChar . $sChar eq "CG"
                || $qChar . $sChar eq "CA"
                || $qChar . $sChar eq "AC"
                || $qChar . $sChar eq "AT"
                || $qChar . $sChar eq "TA" )
        {
          $retStr .= "v";
        }
        elsif (    ( $qChar =~ /[BDHVRYKMSWNX]/ )
                || ( $sChar =~ /[BDHVRYKMSWNX]/ ) )
        {
          $retStr .= "?";
        }
        else {
          $retStr .= " ";
        }
      }
      $retStr .= "\n";
      if (    $obj->{'sbjOrient'} eq "C"
           && $alignmentMode == SearchResult::AlignWithQuerySeq )
      {
        $retStr .= "C ";
        $sEnd = $sStart - length( $sSeq ) + 1 + $deletions;
      }
      else {
        $retStr .= "  ";
        $sEnd = $sStart + length( $sSeq ) - 1 - $deletions;
      }

      # Handle 3 ---- 3 case
      $sEnd = $sStart if ( length( $sSeq ) == $deletions );
      $retStr .= substr( $obj->{'sbjName'}, 0, 13 )
          . " " x (
           13 - (
             length( $obj->{'sbjName'} ) < 13 ? length( $obj->{'sbjName'} ) : 13
           )
          );
      $retStr .= " " x ( 9 - length( $sStart ) ) . $sStart . " $sSeq $sEnd\n";
      $retStr .= "\n";
      $sStart = $sEnd + 1;
      if ( $obj->{sbjOrient} eq "C" ) {
        if ( $alignmentMode == SearchResult::AlignWithSubjSeq ) {
          $qStart = $qEnd - 1;
          $sStart = $sEnd + 1;
        }
        else {
          $qStart = $qEnd + 1;
          $sStart = $sEnd - 1;
        }
      }
      else {
        $qStart = $qEnd + 1;
        $sStart = $sEnd + 1;
      }
    }
    ## TESTING
    $retStr .= "Matrix = " . $obj->getMatrixName() . "\n";
    ##
    $retStr .= "Transitions / transversions = ";
    my ( $mismatches, $transitions, $transversions, $numGaps, $totGapLen ) =
        $obj->_getAlignmentStats();
    if ( defined $transitions ) {
      if ( $transversions > 0 ) {
        $retStr .= sprintf( "%0.2f", ( $transitions / $transversions ) );
      }
      else {
        $retStr .= "0.0";
      }
      $retStr .= " ($transitions / $transversions)\n";
      if ( ( $obj->getQueryEnd() - $obj->getQueryStart() ) > 0 ) {
        $retStr .= "Gap_init rate = "
            . sprintf( "%0.2f",
                       $numGaps /
                           ( $obj->getQueryEnd() - $obj->getQueryStart() ) )
            . " ($numGaps / "
            . ( $obj->getQueryEnd() - $obj->getQueryStart() ) . ")";
      }
      else {
        $retStr .= "Gap_init rate = 0.0 ( $numGaps / 0 )";
      }
      if ( $numGaps > 0 ) {
        $retStr .=
              ", avg. gap size = "
            . sprintf( "%0.2f", ( $totGapLen / $numGaps ) )
            . " ($totGapLen / $numGaps)\n\n";
      }
      else {
        $retStr .= ", avg. gap size = 0.0 (0 / 0)\n\n";
      }

    }
    else {
      $retStr .= "Transitions / transversions = Unknown\n";
      $retStr .= "Gap_init rate = Unknown, avg. gap size = Unknown\n\n";
    }
  }

  return $retStr;
}

##
## TODO: Document
##
sub _toOUTFileFormat {
  my $obj = shift;

  my $outStr    = "";
  my $orient    = "+";
  my $sbjCoord1 = $obj->{sbjBegin};
  my $sbjCoord2 = $obj->{sbjEnd};
  my $sbjCoord3 = "(" . $obj->{sbjLeft} . ")";
  if ( $obj->{sbjOrient} eq "C" ) {
    $orient    = "C";
    $sbjCoord1 = "(" . $obj->{sbjLeft} . ")";
    $sbjCoord2 = $obj->{sbjEnd};
    $sbjCoord3 = $obj->{sbjBegin};
  }
  $outStr = sprintf(
                     "%6d %4.1f %4.1f %4.1f %-17s %8d %8d %8s %1s %-15s %-15s "
                         . "%7s %7s %7s %-5s %3s %3s\n",
                     $obj->{score},   $obj->{percDiv},
                     $obj->{percDel}, $obj->{percIns},
                     $obj->{qryName}, $obj->{qryBegin},
                     $obj->{qryEnd},  "(" . $obj->{qryLeft} . ")",
                     $orient,         $obj->{sbjName},
                     $obj->{sbjType}, $sbjCoord1,
                     $sbjCoord2,      $sbjCoord3,
                     $obj->{id},      $obj->{lineageId},
                     $obj->{overlap}
  );

  return $outStr;

}

##
## TODO: Document
##
sub parseFromCSVFormat {
  my $record = shift;

  my @flds   = split( /\,/, $record );
  my $rawSeq = $flds[ 16 ];

  my $inSub  = 0;
  my $inDel  = 0;
  my $inIns  = 0;
  my $qrySeq = "";
  my $sbjSeq = "";
  while ( $rawSeq =~ s/^(\S)// ) {
    my $char = $1;
    if ( $char eq "/" ) {
      $inSub = 1;
    }
    elsif ( $char eq "-" ) {
      $inDel ^= 1;
    }
    elsif ( $char eq "+" ) {
      $inIns ^= 1;
    }
    else {
      if ( $inSub == 1 ) {
        substr( $sbjSeq, length( $sbjSeq ) - 1, 1 ) = $char;
        $inSub = 0;
      }
      elsif ( $inDel == 1 ) {
        $qrySeq .= $char;
        $sbjSeq .= "-";
      }
      elsif ( $inIns == 1 ) {
        $qrySeq .= "-";
        $sbjSeq .= $char;
      }
      else {
        $qrySeq .= $char;
        $sbjSeq .= $char;
      }
    }
  }

  my $orient = "";
  $orient = "C" if ( $flds[ 13 ] == 1 );
  my $retVal = SearchResult->new(
                                  score          => $flds[ 0 ],
                                  pctDiverge     => $flds[ 1 ],
                                  pctDelete      => $flds[ 2 ],
                                  pctInsert      => $flds[ 3 ],
                                  queryName      => $flds[ 4 ],
                                  queryStart     => $flds[ 5 ],
                                  queryEnd       => $flds[ 6 ],
                                  queryRemaining => $flds[ 7 ],
                                  subjName       => $flds[ 8 ],
                                  subjType       => $flds[ 9 ],
                                  subjStart      => $flds[ 10 ],
                                  subjEnd        => $flds[ 11 ],
                                  subjRemaining  => $flds[ 12 ],
                                  overlap        => $flds[ 14 ],
                                  id             => $flds[ 15 ],
                                  orientation    => $orient,
                                  queryString    => $qrySeq,
                                  subjString     => $sbjSeq
  );
  return $retVal;
}

##
## TODO: Document
##
sub _toCSVFormat {
  my $obj = shift;

  # AAGAA
  #   |
  # AACAA
  #
  # Encode as:  AAG/CAA
  my $qryPos    = 0;
  my $sbjPos    = 0;
  my $insStart  = -1;
  my $delStart  = -1;
  my $indCount  = 0;
  my $indRec    = "";
  my $outSeq    = "";
  my $totIndels = 0;

  my $qrySeq = $obj->getQueryString();
  my $sbjSeq = $obj->getSubjString();
  while ( $qrySeq ne "" ) {
    $qrySeq =~ s/^(\S)//;
    my $qryChar = $1;
    $sbjSeq =~ s/^(\S)//;
    my $sbjChar = $1;
    if ( $qryChar eq "-" ) {
      if ( $insStart == -1 ) {
        $outSeq .= "+";
        $insStart = 1;
      }
      $outSeq .= $sbjChar;
      $sbjPos++;
    }
    elsif ( $sbjChar eq "-" ) {
      if ( $delStart == -1 ) {
        $outSeq .= "-";
        $delStart = 1;
      }
      $outSeq .= $qryChar;
      $qryPos++;
    }
    else {
      if ( $delStart != -1 ) {
        $outSeq .= "-";
        $delStart = -1;
      }
      elsif ( $insStart != -1 ) {
        $outSeq .= "+";
        $insStart = -1;
      }
      if ( $qryChar eq $sbjChar ) {
        $outSeq .= $qryChar;
      }
      else {
        $outSeq .= $qryChar . "/" . $sbjChar;
      }
      $qryPos++;
      $sbjPos++;
    }
  }

  my $cRec =
        $obj->getScore() . ","
      . $obj->getPctDiverge() . ","
      . $obj->getPctDelete() . ","
      . $obj->getPctInsert() . ","
      . $obj->getQueryName() . ","
      . $obj->getQueryStart() . ","
      . $obj->getQueryEnd() . ","
      . $obj->getQueryRemaining() . ","
      . $obj->getSubjName() . ","
      . $obj->getSubjType() . ","
      . $obj->getSubjStart() . ","
      . $obj->getSubjEnd() . ","
      . $obj->getSubjRemaining() . ",";

  if ( $obj->getOrientation =~ /C|c/ ) {
    $cRec .= "1,";
  }
  else {
    $cRec .= "0,";
  }

  $cRec .= $obj->getOverlap() . "," . $obj->getId() . "," . $outSeq;

  return $cRec;
}

##-------------------------------------------------------------------------##
## Private Methods
##-------------------------------------------------------------------------##

##-------------------------------------------------------------------------##
## Use: my ($mismatches, $transitions, $transversion, $numGaps, $totGapLen) =
##                                                 _getAlignmentStats();
##
##   Given an alignment calculate the number of mismatches, transitions
##   transversions, gap count, average gap size etc.
##
##-------------------------------------------------------------------------##
sub _getAlignmentStats {
  my $obj = shift;

  my $querySeq = $obj->getQueryString();
  my $subjSeq  = $obj->getSubjString();

  return if ( !defined $querySeq || $querySeq eq "" );

  my $transitions   = 0;
  my $transversions = 0;
  my $gaps          = 0;
  my $totGapLen     = 0;
  my $mismatches    = 0;
  my $inGap         = 0;

  while ( $querySeq ne "" ) {
    $querySeq =~ s/^(\S)//;
    my $qryChar = $1;
    $subjSeq =~ s/^(\S)//;
    my $sbjChar = $1;
    if ( $qryChar eq $sbjChar ) {
      $inGap = 0;
      next;
    }
    if ( $qryChar eq "-" || $sbjChar eq "-" ) {
      if ( $inGap == 0 ) {
        $gaps++;
        $inGap = 1;
      }
      $totGapLen++;
    }
    else {
      $inGap = 0;
      my $basePair = uc( $qryChar . $sbjChar );
      if ( $basePair =~ /CT|TC|GA|AG/ ) {
        $transitions++;
      }
      elsif ( $basePair =~ /GT|TG|TA|AT|CA|AC|CG|GC/ ) {
        $transversions++;
      }
      $mismatches++;
    }
  }

  return ( $mismatches, $transitions, $transversions, $gaps, $totGapLen );
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
