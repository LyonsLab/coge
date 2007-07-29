#!/usr/bin/perl -w
##---------------------------------------------------------------------------##
##  File:
##      @(#) SearchResultCollection.pm
##  Author:
##      Robert M. Hubley   rhubley@systemsbiology.org
##  Description:
##      An datastructure for holding search results
##      from various search engines.
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
# bless( { 'results' => [ $SearchResultRef1,
#                         $SearchResultRef2,
#                         ..
#                       ]
#        }
#     }, 'SearchResultCollection' );
#
###############################################################################
# ChangeLog
#
#     $Log: SearchResultCollection.pm,v $
#     Revision 1.39  2007/05/17 21:01:54  rhubley
#     Cleanup before a distribution
#
#
###############################################################################
# To Do:
#
#

=head1 NAME

SearchResultCollection

=head1 SYNOPSIS

use SearchResultCollection

Usage: 

    $SearchResultsCollection = SearchResultCollection->new();

=head1 DESCRIPTION

A class for storing the results from a sequence search engine.
NOTE: This is basically an ArrayList with a specialized write
method.  See ArrayList.pm for accessors methods.

=head1 INSTANCE METHODS

=cut 

package SearchResultCollection;
use strict;
use SearchResult;
use Data::Dumper;
use Carp;
use vars qw(@ISA @EXPORT @EXPORT_OK %EXPORT_TAGS $VERSION);

require ArrayList;
require Exporter;

@ISA = qw(ArrayList Exporter);

@EXPORT = qw();

@EXPORT_OK = qw();

%EXPORT_TAGS = ( all => [ @EXPORT_OK ] );

my $CLASS = "SearchResultCollection";
my $DEBUG = 0;

##-------------------------------------------------------------------------##
## Constructor
##-------------------------------------------------------------------------##
sub new {
  my $class = shift;

  # TODO: Make this handle filenames!!!!
  #       ie. move parser back!
  my $this = $class->SUPER::new( @_ );

  return $this;
}

##-------------------------------------------------------------------------##
## Get and Set Methods
##-------------------------------------------------------------------------##

##-------------------------------------------------------------------------##
## General Object Methods
##-------------------------------------------------------------------------##

##-------------------------------------------------------------------------##

=head2 write()

  Use: $obj->write( $file, $alignmentMode );

    $file : String
    $alignmentMode : SearchResult::NoAlign
                     SearchResult::AlignWithQuerySeq
                     SearchResult::AlignWithSubjSeq

  Write object to a file.  If $printAlignments
  then include the alignment data.

=cut

##-------------------------------------------------------------------------##
## TODO: Make this able to write to stdout as well as to a file
##       using the filehandle check below:
##  my $OUT;
##  if ( ref( $destination ) !~ /GLOB|FileHandle/ ) {
##    print $CLASS
##        . "::writeHTMLMultipleAlignments() Opening file "
##        . $destination . "\n"
##        if ( $DEBUG );
##    open $OUT, $destination
##        or die $CLASS
##        . "::writeHTMLMultipleAlignments(): Unable to open "
##        . "results file: $desitnation : $!";
##  }
##  else {
##    $OUT = $desitination;
##  }
##
sub write {
  my $this          = shift;
  my $filename      = shift;
  my $alignmentMode = shift;

  open OUT, ">$filename"
      or croak "SearchResultCollection::write - "
      . "Error...could not open file $filename "
      . "for writing!\n";

  for ( my $i = 0 ; $i < $this->size() ; $i++ ) {
    print OUT $this->get( $i )->toStringFormatted( $alignmentMode );
  }
  close OUT;
}

##-------------------------------------------------------------------------##

=head2 maskLevelFilter()

  Use: $numRemoved = $obj->maskLevelFilter(
                                -value => 80,
                                [-strand => SearchResult::QueryStrand
                                            SearchResult::SubjectStrand] );

    value : Integer - masklevel
    strand: Integer - SearchResult::QueryStrand (default) or 
                      SearchResult::SubjectStrand.  

  The masklevel controls the reporting of matches based on the 
  overlap of aligned bases.  Typically a match is reported only
  if at least (100 - masklevel)% of the bases in its "domain"
  (the part of the query that is aligned) are not contained within
  the domain of any higher-scoring match.

  Returns the count of entries filtered out. Does *not* alter
  the sort order of the remaining SearchResultCollection.

=cut

##-------------------------------------------------------------------------##
sub maskLevelFilter {
  my $this       = shift;
  my %parameters = @_;

  croak $CLASS
      . "::maskLevelFilter: Missing or invalid mask level value "
      . "parameter ( value = "
      . $parameters{'value'} . " )!\n"
      if ( !( $parameters{'value'} =~ /\d+/ ) );
  my $maskLevel = $parameters{'value'};

  croak $CLASS
      . "::maskLevelFilter: Invalid strand parameter "
      . "parameter ( strand = "
      . $parameters{'strand'}
      . ")!\n"
      if (
           defined $parameters{'strand'}
           && (    $parameters{'strand'} != SearchResult::QueryStrand
                || $parameters{'strand'} != SearchResult::SubjectStrand )
      );

  #
  # Filter set using "masklevel" concept
  #
  return ( 0 ) if ( $maskLevel >= 101 );

  my $strandName  = "getQueryName";
  my $strandStart = "getQueryStart";
  my $strandEnd   = "getQueryEnd";
  if ( defined $parameters{'strand'}
       && $parameters{'strand'} == SearchResult::SubjectStrand )
  {
    $strandName  = "getSubjName";
    $strandStart = "getSubjStart";
    $strandEnd   = "getSubjEnd";

  }

  # Sort by strand and then scores high to low
  my @indices = sort {
           $this->get( $a )->$strandName() cmp $this->get( $b )->$strandName()
        || $this->get( $b )->getScore <=> $this->get( $a )->getScore
  } ( 0 .. $this->size() - 1 );

  my %deleteHash = ();    # Hash to hold indexes of filtered entries
                          # Coverage Linked List
                          #    [ start, end, next ]
                          #
  my $coverageRanges;
  my $prevID = "";
  for ( my $i = 0 ; $i <= $#indices ; $i++ ) {
    my $refResult = $this->get( $indices[ $i ] );
    my $refID     = $refResult->$strandName();
    my $refStart  = $refResult->$strandStart();
    my $refEnd    = $refResult->$strandEnd();
    my $refScore  = $refResult->getScore();

    if ( $refID ne $prevID ) {
      $coverageRanges = [
                          $this->get( $indices[ $i ] )->$strandStart(),
                          $this->get( $indices[ $i ] )->$strandEnd(),
                          0
      ];
      $prevID = $refID;
      next;
    }

    my $overlapBP = 0;

    if ( $DEBUG ) {
      print "Considering reference hit:\n";
      print "    Subject: " . $refResult->getSubjName() . "\n";
      print "    Range: $refStart - $refEnd\n";
      print "    Score: $refScore\n";
    }

    my $rangePtr = $coverageRanges;
    print "Checking for overlap with any coverage ranges:\n" if ( $DEBUG );
    while ( $rangePtr != 0 ) {
      my $rangeStart = $rangePtr->[ 0 ];
      my $rangeEnd   = $rangePtr->[ 1 ];
      print "  Range: $rangeStart - $rangeEnd, next=$rangePtr->[2]\n"
          if ( $DEBUG );
      if ( $rangeEnd < $refStart ) {
        $rangePtr = $rangePtr->[ 2 ];
        next;
      }
      last if ( $rangeStart > $refEnd );
      my $overlapStart = $refStart;
      $overlapStart = $rangeStart if ( $overlapStart < $rangeStart );
      my $overlapEnd = $refEnd;
      $overlapEnd = $rangeEnd if ( $overlapEnd > $rangeEnd );
      $overlapBP += $overlapEnd - $overlapStart + 1;
      print "  Overlap:  $overlapStart - $overlapEnd = $overlapBP\n"
          if ( $DEBUG );
      $rangePtr = $rangePtr->[ 2 ];
    }

    my $percOverlap = sprintf(
                               "%0.0f",
                               int(
                                  ( $overlapBP / ( $refEnd - $refStart + 1 ) ) *
                                      100
                               )
    );
    print "Done....overlap perc = $percOverlap\n" if ( $DEBUG );

    if ( $percOverlap < $maskLevel ) {

      # keep it!
      # Find it's place in the list
      print "Keeping it...finding its place in the coverage list\n"
          if ( $DEBUG );
      my $rangePtrA     = $coverageRanges;
      my $prevRangePtrA = 0;
      while ( $rangePtrA != 0 ) {
        print "Considering segment: "
            . $rangePtrA->[ 0 ] . "-"
            . $rangePtrA->[ 1 ]
            . ", next="
            . $rangePtrA->[ 2 ] . "\n"
            if ( $DEBUG );

        if ( $rangePtrA->[ 2 ] == 0 && $refStart > $rangePtrA->[ 1 ] ) {
          print "  Not in list...adding after\n" if ( $DEBUG );
          $rangePtrA->[ 2 ] = [ $refStart, $refEnd, $rangePtrA->[ 2 ] ];
        }
        elsif ( $refEnd < $rangePtrA->[ 0 ] ) {
          print "  Not in list...adding before\n" if ( $DEBUG );
          if ( $prevRangePtrA != 0 ) {
            $prevRangePtrA->[ 2 ] = [ $refStart, $refEnd, $rangePtrA ];
          }
          else {
            $coverageRanges = [ $refStart, $refEnd, $rangePtrA ];
          }
          last;
        }
        else {
          print
"  In list...Extending existing: rangePtrA: $rangePtrA->[0] - $rangePtrA->[1]\n"
              if ( $DEBUG );
          if (   $refStart < $rangePtrA->[ 0 ]
              || $refEnd > $rangePtrA->[ 1 ] && $refStart <= $rangePtrA->[ 1 ] )
          {
            if ( $refStart < $rangePtrA->[ 0 ] ) {
              $rangePtrA->[ 0 ] = $refStart;
              print "      rangePtrA-start after = $rangePtrA->[0]\n"
                  if ( $DEBUG );
            }

            if ( $refEnd > $rangePtrA->[ 1 ] && $refStart <= $rangePtrA->[ 1 ] )
            {
              $rangePtrA->[ 1 ] = $refEnd;
              my $rangePtrB = $rangePtrA->[ 2 ];
              while ( $rangePtrB != 0 && $rangePtrB->[ 2 ] != 0 ) {
                last if ( $rangePtrB->[ 0 ] > $refEnd );

                # Combine ranges
                $rangePtrA->[ 1 ] = $rangePtrB->[ 1 ]
                    if ( $refEnd < $rangePtrB->[ 1 ] );
                $rangePtrA->[ 2 ] = $rangePtrB->[ 2 ];
                $rangePtrB = $rangePtrB->[ 2 ];
              }
            }

            last;
          }
        }
        print
"  In list...Extending after: rangePtrA: $rangePtrA->[0] - $rangePtrA->[1]\n"
            if ( $DEBUG );
        last if ( $refEnd < $rangePtrA->[ 0 ] );
        $prevRangePtrA = $rangePtrA;
        $rangePtrA     = $rangePtrA->[ 2 ];
      }

    }
    else {

      # Mark for deletion
      print "Tossing it ( marking )\n" if ( $DEBUG );
      $deleteHash{ $indices[ $i ] } = 1;
    }
  }

  # Remove all hits which were filtered above
  my $numRemoved = scalar( keys( %deleteHash ) );
  if ( keys( %deleteHash ) ) {
    foreach my $index ( sort { $b <=> $a } keys( %deleteHash ) ) {
      $this->remove( $index );
    }
  }
  undef %deleteHash;

  return ( $numRemoved );

}    # sub maskLevelFilter

##-------------------------------------------------------------------------##
## Private Methods
##-------------------------------------------------------------------------##
##-------------------------------------------------------------------------##

=head2 old...maskLevelFilter()

  Using Arians suggestion

  Use: $numRemoved = $obj->maskLevelFilter(
                                -value => 80,
                                [-strand => SearchResult::QueryStrand
                                            SearchResult::SubjectStrand] );

    value : Integer - masklevel
    strand: Integer - SearchResult::QueryStrand (default) or 
                      SearchResult::SubjectStrand.  

  The masklevel controls the reporting of matches based on the 
  overlap of aligned bases.  Typically a match is reported only
  if at least (100 - masklevel)% of the bases in its "domain"
  (the part of the query that is aligned) are not contained within
  the domain of any higher-scoring match.

  Returns the count of entries filtered out. Does *not* alter
  the sort order of the remaining SearchResultCollection.

=cut

##-------------------------------------------------------------------------##
sub maskLevelFilterArian {
  my $this       = shift;
  my %parameters = @_;

  croak $CLASS
      . "::maskLevelFilter: Missing or invalid mask level value "
      . "parameter ( value = "
      . $parameters{'value'} . " )!\n"
      if ( !( $parameters{'value'} =~ /\d+/ ) );
  my $maskLevel = $parameters{'value'};

  croak $CLASS
      . "::maskLevelFilter: Invalid strand parameter "
      . "parameter ( strand = "
      . $parameters{'strand'}
      . ")!\n"
      if (
           defined $parameters{'strand'}
           && (    $parameters{'strand'} != SearchResult::QueryStrand
                || $parameters{'strand'} != SearchResult::SubjectStrand )
      );

  #
  # Filter set using "masklevel" concept
  #
  return ( 0 ) if ( $maskLevel >= 101 );

  my $strandName  = "getQueryName";
  my $strandBegin = "getQueryStart";
  my $strandEnd   = "getQueryEnd";
  if ( defined $parameters{'strand'}
       && $parameters{'strand'} == SearchResult::SubjectStrand )
  {
    $strandName  = "getSubjName";
    $strandBegin = "getSubjStart";
    $strandEnd   = "getSubjEnd";

  }

  # endIndices = ( collIndex#1, collIndex#2, collIndex#3 .. )
  my @endIndices = sort {
           $this->get( $a )->$strandName() cmp $this->get( $b )->$strandName()
        || $this->get( $a )->$strandEnd <=> $this->get( $b )->$strandEnd
  } ( 0 .. $this->size() - 1 );

# startIndices = ( [ collIndex#1, endIndex#1 ], [ collIndex#2, endIndex#2 ],... )
  my $index = 0;
  my @startIndices = map { [ $_, $index++ ] } @endIndices;
  @startIndices = sort {
    $this->get( $a->[ 0 ] )->$strandName()
        cmp $this->get( $b->[ 0 ] )->$strandName()
        || $this->get( $a->[ 0 ] )->$strandBegin <=> $this->get( $b->[ 0 ] )
        ->$strandBegin
  } @startIndices;

  my %deleteHash = ();    # Hash to hold indexes of filtered entries
  for ( my $i = 0 ; $i <= $#startIndices ; $i++ ) {
    next if ( $deleteHash{ $startIndices[ $i ]->[ 0 ] } == 1 );
    my $refResult = $this->get( $startIndices[ $i ]->[ 0 ] );
    my $refID     = $refResult->$strandName();
    my $refBegin  = $refResult->$strandBegin();
    my $refEnd    = $refResult->$strandEnd();
    my $refScore  = $refResult->getScore();

    if ( $DEBUG ) {
      print "Considering reference hit:\n";
      print "    Subject: " . $refResult->getSubjName() . "\n";
      print "    Range: $refBegin - $refEnd\n";
      print "    Score: $refScore\n";
    }

    my @overlapRanges = ();

    #
    # Search Left
    #
    my $endIndex = $startIndices[ $i ]->[ 1 ];
    for ( my $j = $endIndex - 1 ; $j >= 0 ; $j-- ) {

      #
      #next if ( $deleteHash{ $endIndices[ $j ] } == 1 );
      my $neighborResult = $this->get( $endIndices[ $j ] );
      my $neighborID     = $neighborResult->$strandName();
      last if ( $neighborID ne $refID );
      my $neighborBegin = $neighborResult->$strandBegin();
      my $neighborEnd   = $neighborResult->$strandEnd();
      my $neighborScore = $neighborResult->getScore();

      if ( $DEBUG ) {
        print "      Considering Left Range: $neighborBegin - "
            . "$neighborEnd, $neighborScore"
            . $neighborResult->getSubjName() . "\n";
      }

      last if ( $neighborEnd < $refBegin );

      if ( $neighborScore < $refScore ) {
        if ( $maskLevel < 101 && $neighborBegin >= $refBegin ) {

          # This is a lower score *and* it's completely
          # contained inside the ref....mark this for deletion!
          #print "Deleting! " .  $neighborResult->getSubjName() . "\n";
          $deleteHash{ $endIndices[ $j ] } = 1;
        }
        next;
      }
      elsif ( $neighborScore == $refScore ) {

        #print "Deleting2! " .  $neighborResult->getSubjName() . "\n";
        $deleteHash{ $endIndices[ $j ] } = 1;
        next;
      }

      # Calc overlap range
      my $overlapBegin = $refBegin;
      my $overlapEnd   = $refEnd;
      $overlapBegin = $neighborBegin if ( $neighborBegin > $overlapBegin );
      $overlapEnd   = $neighborEnd   if ( $neighborEnd < $overlapEnd );

      print "           Overlaps!!!\n" if ( $DEBUG );

      # Store range
      push @overlapRanges, ( [ $overlapBegin, $overlapEnd ] );
    }

    #
    # Search Right
    #
    for ( my $j = $i + 1 ; $j <= $#startIndices ; $j++ ) {

      #
      #next if ( $deleteHash{ $startIndices[ $j ]->[0] } == 1 );
      my $neighborResult = $this->get( $startIndices[ $j ]->[ 0 ] );
      my $neighborID     = $neighborResult->$strandName();
      last if ( $neighborID ne $refID );
      my $neighborBegin = $neighborResult->$strandBegin();
      my $neighborEnd   = $neighborResult->$strandEnd();
      my $neighborScore = $neighborResult->getScore();

      if ( $DEBUG ) {
        print "      Considering Right Range: $neighborBegin - "
            . "$neighborEnd, $neighborScore"
            . $neighborResult->getSubjName() . "\n";
      }

      last if ( $neighborBegin > $refEnd );
      if ( $neighborScore < $refScore ) {
        if ( $maskLevel < 101 && $neighborEnd <= $refEnd ) {

          # This is a lower score *and* it's completely
          # contained inside the ref....mark this for deletion!
          #print "Deleting! " .  $neighborResult->getSubjName() . "\n";
          $deleteHash{ $startIndices[ $j ]->[ 0 ] } = 1;
        }
        next;
      }
      elsif ( $neighborScore == $refScore ) {

        #print "Deleting2! " .  $neighborResult->getSubjName() . "\n";
        $deleteHash{ $startIndices[ $j ]->[ 0 ] } = 1;
        next;
      }
      next if ( $neighborScore < $refScore );

      # Calc overlap range
      my $overlapBegin = $refBegin;
      my $overlapEnd   = $refEnd;
      $overlapBegin = $neighborBegin if ( $neighborBegin > $overlapBegin );
      $overlapEnd   = $neighborEnd   if ( $neighborEnd < $overlapEnd );

      print "           Overlaps!!!\n" if ( $DEBUG );

      # Store range
      push @overlapRanges, ( [ $overlapBegin, $overlapEnd ] );

      #
    }

    @overlapRanges = sort { $a->[ 0 ] <=> $b->[ 0 ] } @overlapRanges;

    my $bpCovered = 0;
    my $lastEnd   = 0;
    print "Coverage analysis\n" if ( $DEBUG );
    foreach my $range ( @overlapRanges ) {
      print "Considering range: $range->[0] - $range->[1]\n" if ( $DEBUG );
      if ( $bpCovered == 0 && $lastEnd == 0 ) {
        print "First one\n" if ( $DEBUG );
        $bpCovered = $range->[ 1 ] - $range->[ 0 ] + 1;
        $lastEnd   = $range->[ 1 ];
      }
      elsif ( $range->[ 1 ] > $lastEnd ) {
        if ( $range->[ 0 ] > $lastEnd ) {

          # outside of last hit
          print "Adding it's bp\n" if ( $DEBUG );
          $bpCovered += $range->[ 1 ] - $range->[ 0 ] + 1;
        }
        else {

          # overlaps last hit
          print "Adding it's non-overlapping bp\n" if ( $DEBUG );
          $bpCovered += $range->[ 1 ] - $lastEnd;
        }
        $lastEnd = $range->[ 1 ];
      }
      else {

        # Contained in previous hit
      }
    }

    print "  - BP Covered = $bpCovered\n" if ( $DEBUG );
    print "  - % of ref = "
        . sprintf( "%0.0f", ( $bpCovered / ( $refEnd - $refBegin + 1 ) ) * 100 )
        . "\n"
        if ( $DEBUG );
    my $percCovered = ( $bpCovered / ( $refEnd - $refBegin + 1 ) ) * 100;

    if ( $percCovered > $maskLevel ) {
      print "Marking for deletion!\n" if ( $DEBUG );
      $deleteHash{ $startIndices[ $i ]->[ 0 ] } = 1;
    }

  }

  # Remove all hits which were filtered above
  my $numRemoved = scalar( keys( %deleteHash ) );
  if ( keys( %deleteHash ) ) {
    foreach my $index ( sort { $b <=> $a } keys( %deleteHash ) ) {
      $this->remove( $index );
    }
  }
  undef %deleteHash;

  return ( $numRemoved );

}    # sub maskLevelFilter

##-------------------------------------------------------------------------##

=head2 old....maskLevelFilter()

  Using MaskerAid method

  Use: $numRemoved = $obj->maskLevelFilter(
                                -value => 80,
                                [-strand => SearchResult::QueryStrand
                                            SearchResult::SubjectStrand] );

    value : Integer - masklevel
    strand: Integer - SearchResult::QueryStrand (default) or 
                      SearchResult::SubjectStrand.  

  The masklevel controls the reporting of matches based on the 
  overlap of aligned bases.  Typically a match is reported only
  if at least (100 - masklevel)% of the bases in its "domain"
  (the part of the query that is aligned) are not contained within
  the domain of any higher-scoring match.

  Returns the count of entries filtered out. Does *not* alter
  the sort order of the remaining SearchResultCollection.

=cut

##-------------------------------------------------------------------------##
sub maskLevelFilterMaskerAid {
  my $this       = shift;
  my %parameters = @_;

  croak $CLASS
      . "::maskLevelFilter: Missing or invalid mask level value "
      . "parameter ( value = "
      . $parameters{'value'} . " )!\n"
      if ( !( $parameters{'value'} =~ /\d+/ ) );
  my $maskLevel = $parameters{'value'};

  croak $CLASS
      . "::maskLevelFilter: Invalid strand parameter "
      . "parameter ( strand = "
      . $parameters{'strand'}
      . ")!\n"
      if (
           defined $parameters{'strand'}
           && (    $parameters{'strand'} != SearchResult::QueryStrand
                || $parameters{'strand'} != SearchResult::SubjectStrand )
      );

  #
  # Filter set using "masklevel" concept
  #
  return ( 0 ) if ( $maskLevel >= 101 );

  my $strandName  = "getQueryName";
  my $strandStart = "getQueryStart";
  my $strandEnd   = "getQueryEnd";
  if ( defined $parameters{'strand'}
       && $parameters{'strand'} == SearchResult::SubjectStrand )
  {
    $strandName  = "getSubjName";
    $strandStart = "getSubjStart";
    $strandEnd   = "getSubjEnd";

  }

  # Sort by strand and then scores high to low
  my @indices = sort {
           $this->get( $a )->$strandName() cmp $this->get( $b )->$strandName()
        || $this->get( $b )->getScore <=> $this->get( $a )->getScore
  } ( 0 .. $this->size() - 1 );

  my %deleteHash = ();    # Hash to hold indexes of filtered entries
  for ( my $i = 0 ; $i <= $#indices ; $i++ ) {
    next if ( $deleteHash{ $indices[ $i ] } == 2 );

    for ( my $j = $i + 1 ; $j <= $#indices ; $j++ ) {
      next if ( defined $deleteHash{ $indices[ $j ] } );

      my $result1 = $this->get( $indices[ $i ] );
      my $result2 = $this->get( $indices[ $j ] );

      last if ( $result1->$strandName() ne $result2->$strandName() );

      # Get members
      my $begin1 = $result1->$strandStart();
      my $begin2 = $result2->$strandStart();
      my $end1   = $result1->$strandEnd();
      my $end2   = $result2->$strandEnd();

      # Check if they overlap
      next if ( $begin2 > $end1 || $end2 < $begin1 );
      next if ( $begin2 > $end1 || $end2 < $begin1 );

      # Calc overlap extension
      my $extension = 0;
      $extension = $begin1 - $begin2 if ( $begin2 < $begin1 );
      $extension += $end2 - $end1 if ( $end2 > $end1 );

      # % of the low scoring HSP outside the domain of the
      # high scoring HSP
      my $perc = ( $extension / ( $end2 - $begin2 + 1 ) ) * 100;
      if ( $perc < ( 100 - $maskLevel ) ) {
        if ( $perc == 0 ) {
          $deleteHash{ $indices[ $j ] } = 2;
        }
        else {
          $deleteHash{ $indices[ $j ] } = 1;
        }
      }
    }    # for ( my $j = $i + 1...
  }    # for ( my $i = 0...
  undef @indices;

  # Remove all hits which were filtered above
  my $numRemoved = scalar( keys( %deleteHash ) );
  if ( keys( %deleteHash ) ) {
    foreach my $index ( sort { $b <=> $a } keys( %deleteHash ) ) {
      $this->remove( $index );
    }
  }
  undef %deleteHash;

  return ( $numRemoved );

}    # sub maskLevelFilter

1;
