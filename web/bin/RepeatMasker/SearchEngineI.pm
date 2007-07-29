#!/usr/local/bin/perl -w
##---------------------------------------------------------------------------##
##  File:
##      @(#) SearchEngineI.pm
##  Author:
##      Robert M. Hubley   rhubley@systemsbiology.org
##  Description:
##      An interface which defines how RepeatMasker interacts
##      with a search engine. As with most perl interfaces
##      this also includes a simple implementation of a few
##      methods.  In particular this implements the getters/setters
##      for convenience.
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
#      'SearchEngineI' );
#
###############################################################################
# ChangeLog
#
#     $Log: SearchEngineI.pm,v $
#     Revision 1.40  2007/05/17 21:01:54  rhubley
#     Cleanup before a distribution
#
#
###############################################################################
#
# To Do:
#
#

=head1 NAME

SearchEngineI

=head1 SYNOPSIS

use SearchEngineI

Usage: 

  use SearchEngineI;
  use CrossmatchSearchEngine;

  my $cEngine = CrossmatchSearchEngine->new( 
                         pathToEngine=>"/usr/local/bin/crossmatch" );
  $cEngine->setScoreMode( SearchEngineI::complexityAdjustedScoreMode );


=head1 DESCRIPTION

  An interface for sequence search engines such as CrossMatch and WUBlast. 
  This interface is really a cross between an interface and abstract class.  
  It is not intended for use on it's own but in conjunction with a concrete
  implementation of all of it's methods such as CrossmatchSearchEngine.pm.

=head1 INSTANCE METHODS

=cut 

package SearchEngineI;
use strict;
use SearchEngineI;
use SearchResultCollection;
use Data::Dumper;
use Carp;

#
# Constants
#
use constant basicScoreMode              => 1;
use constant complexityAdjustedScoreMode => 2;

use vars qw(@ISA @EXPORT @EXPORT_OK %EXPORT_TAGS $VERSION);

require Exporter;

@ISA = qw(Exporter);

@EXPORT = qw();

@EXPORT_OK = qw();

%EXPORT_TAGS = ( all => [ @EXPORT_OK ] );

#
# Version
#
my $VERSION = 0.1;
my $CLASS   = "SearchEngineI";

##-------------------------------------------------------------------------##
## Constructor
##-------------------------------------------------------------------------##
sub new {
  my $class = shift;

  # Create ourself as a hash
  my $this = {};

  # Bless this hash in the name of the father, the son...
  bless $this, $class;

  return $this;
}

##-------------------------------------------------------------------------##
## Get and Set Methods
##-------------------------------------------------------------------------##

##-------------------------------------------------------------------------##

=over 4

=item Use: my $value = getPathToEngine( );

=item Use: my $oldValue = setPathToEngine( $value );

Get/Set the fully qualified path to the search engine
binary file.

=back 

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

  my $oldValue = $this->{'pathToEngine'};
  $this->{'pathToEngine'} = $value;

  return $oldValue;
}

##-------------------------------------------------------------------------##

=over 4

=item Use: my $value = getScoreMode( );

=item Use: my $oldValue = setScoreMode( $value );

Get/Set the score mode.  Typically a search engine has it's own
scoring method.  This is an attempt to define the normalization
process so that one scheme can be used for many search engines.
Unfortunately both these schemes require some knowledge of the
alignment of the hit.  If this is not available via the primary
search engine the implmementor may need to perform a local alignment
themselves once the region of the hit has been determined.

The parameter passed is a string constant:

    $value :  SearchEngineI::basicScoreMode
                 This mode uses the supplied matrix, gap penalties
                 and the alignment itself to compute a basic score.

              SearchEngineI::complexityAdjustedScoreMode
                This mode uses Phil Green's sequence complexity
                adjustment to compute the score.

=back 

=cut

##-------------------------------------------------------------------------##
sub getScoreMode {
  my $this = shift;

  return $this->{'scoreMode'};
}

sub setScoreMode {
  my $this  = shift;
  my $value = shift;

  if (    $value != SearchEngineI::basicScoreMode
       && $value != SearchEngineI::complexityAdjustedScoreMode )
  {
    croak $CLASS
        . "::setScoreMode( $value ): Value does not match\n"
        . "either SearchEngineI::basicScoreMode or "
        . "SearchEngineI::complexityAdjustedScoreMode\n";
  }

  my $oldValue = $this->{'scoreMode'};
  $this->{'scoreMode'} = $value;

  return $oldValue;
}

##-------------------------------------------------------------------------##

=over 4

=item Use: my $value = getGenerateAlignments( );

=item Use: my $oldValue = setGenerateAlignments( $value );

Get/Set the generate alignments attribute.  This....

=back

=cut

##-------------------------------------------------------------------------##
sub getGenerateAlignments {
  my $this = shift;

  return $this->{'generateAlignments'};
}

sub setGenerateAlignments {
  my $this  = shift;
  my $value = shift;

  my $oldValue = $this->{'generateAlignments'};
  $this->{'generateAlignments'} = $value;

  return $oldValue;
}

##-------------------------------------------------------------------------##

=over 4

=item Use: my $value = getWordRaw( );

=item Use: my $oldValue = setWordRaw( $value );

Get/Set the word_raw attribute.  This....

=back 

=cut

##-------------------------------------------------------------------------##
sub getWordRaw {
  my $this = shift;

  return $this->{'wordRaw'};
}

sub setWordRaw {
  my $this  = shift;
  my $value = shift;

  my $oldValue = $this->{'wordRaw'};
  $this->{'wordRaw'} = $value;

  return $oldValue;
}

##-------------------------------------------------------------------------##

=over 4

=item Use: my $value = getGapInit( );

=item Use: my $oldValue = setGapInit( $value );

Get/Set the gap_init attribute.  This is the penalty
imposed (subtracted from the score of the hit) on the candidate 
match when a gap is initiated in the query or subject sequence. 

  $value :  Integer >= 0 

=back 

=cut

##-------------------------------------------------------------------------##
sub getGapInit {
  my $this = shift;

  return $this->{'gapInit'};
}

sub setGapInit {
  my $this  = shift;
  my $value = shift;

  my $oldValue = $this->{'gapInit'};
  $this->{'gapInit'} = $value;

  return $oldValue;
}

##-------------------------------------------------------------------------##

=over 4

=item Use: my $value = getInsGapExt( );

=item Use: my $oldValue = setInsGapExt( $value );

Get/Set the insertion extension type gap penalty attribute.  
This is the penalty imposed (subtracted from the score of the hit) 
on the candidate match when an already initiated insertion gap is
extended by one base. 

  $value :  Integer >= 0 

=back 

=cut

##-------------------------------------------------------------------------##
sub getInsGapExt {
  my $this = shift;

  return $this->{'insGapExt'};
}

sub setInsGapExt {
  my $this  = shift;
  my $value = shift;

  my $oldValue = $this->{'insGapExt'};
  $this->{'insGapExt'} = $value;

  return $oldValue;
}

##-------------------------------------------------------------------------##

=over 4

=item Use: my $value = getDelGapExt( );

=item Use: my $oldValue = setDelGapExt( $value );

Get/Set the deletion extension type gap penalty attribute.  
This is the penalty imposed (subtracted from the score of the hit) 
on the candidate match when an already initiated deletion gap is
extended by one base. 

  $value :  Integer >= 0 

=back 

=cut

##-------------------------------------------------------------------------##
sub getDelGapExt {
  my $this = shift;

  return $this->{'delGapExt'};
}

sub setDelGapExt {
  my $this  = shift;
  my $value = shift;

  my $oldValue = $this->{'delGapExt'};
  $this->{'delGapExt'} = $value;

  return $oldValue;
}

##-------------------------------------------------------------------------##

=over 4

=item Use: my $value = getMinMatch( );

=item Use: my $oldValue = setMinMatch( $value );

Get/Set the minimum match length.  This is the minimum
length of a hit that should be returned by search
engine (in bp).

  $value :  Integer >= 0 

=back 

=cut

##-------------------------------------------------------------------------##
sub getMinMatch {
  my $this = shift;

  return $this->{'minMatch'};
}

sub setMinMatch {
  my $this  = shift;
  my $value = shift;

  my $oldValue = $this->{'minMatch'};
  $this->{'minMatch'} = $value;

  return $oldValue;
}

##-------------------------------------------------------------------------##

=over 4

=item Use: my $value = getMinScore( );

=item Use: my $oldValue = setMinScore( $value );

Get/Set the minimum score value.  This is the minimum
score of a hit returned by search engine (in bp).

  $value :  Integer >= 0 

=back 

=cut

##-------------------------------------------------------------------------##
sub getMinScore {
  my $this = shift;

  return $this->{'minScore'};
}

sub setMinScore {
  my $this  = shift;
  my $value = shift;

  my $oldValue = $this->{'minScore'};
  $this->{'minScore'} = $value;

  return $oldValue;
}

##-------------------------------------------------------------------------##

=over 4

=item Use: my $value = getBandwidth( );

=item Use: my $oldValue = setBandwidth( $value );

Get/Set the bandwidth value.  This is the maximum width
of sequence surrounding a short exact match which is often
used to search for the best local alignment to entire sequence.
Lowering the bandwidth increases speed at the cost of
sensitivity.

  $value :  Integer > 0 

=back 

=cut

##-------------------------------------------------------------------------##
sub getBandwidth {
  my $this = shift;

  return $this->{'bandwidth'};
}

sub setBandwidth {
  my $this  = shift;
  my $value = shift;

  my $oldValue = $this->{'bandwidth'};
  $this->{'bandwidth'} = $value;

  return $oldValue;
}

##-------------------------------------------------------------------------##

=over 4

=item Use: my $value = getMaskLevel( );

=item Use: my $oldValue = setMaskLevel( $value );

Get/Set the masklevel value.  The masklevel controls the reporting
of matches based on the overlap of aligned bases.  Typically a match is 
reported only if at least (100 - masklevel)% of the bases in its "domain" 
(the part of the query that is aligned) are not contained within the domain 
of any higher-scoring match.

 Special cases:
    masklevel 0     report only the single highest scoring match for each query                                                                                
    masklevel 100   report any match whose domain is not completely contained
                        within a higher scoring match
    masklevel 101   report all matches


  $value :  0 <= Integer <= 101

=back 

=cut

##-------------------------------------------------------------------------##
sub getMaskLevel {
  my $this = shift;

  return $this->{'maskLevel'};
}

sub setMaskLevel {
  my $this  = shift;
  my $value = shift;

  my $oldValue = $this->{'maskLevel'};
  $this->{'maskLevel'} = $value;

  return $oldValue;
}

##-------------------------------------------------------------------------##

=over 4

=item Use: my $value = getMatrix( );

=item Use: my $oldValue = setMatrix( $value );

Get/Set the matrix value.  

=back 

=cut

##-------------------------------------------------------------------------##
sub getMatrix {
  my $this = shift;

  return $this->{'matrixName'};
}

sub setMatrix {
  my $this  = shift;
  my $value = shift;

  my $oldValue = $this->{'matrixName'};
  $this->{'matrixName'} = $value;

  return $oldValue;
}

##-------------------------------------------------------------------------##

=over 4

=item Use: my $value = getSubject( );

=item Use: my $oldValue = setSubject( $value );

Get/Set the subject database value.  

=back 

=cut

##-------------------------------------------------------------------------##
sub getSubject {
  my $this = shift;

  return $this->{'subject'};
}

sub setSubject {
  my $this  = shift;
  my $value = shift;

  my $oldValue = $this->{'subject'};
  $this->{'subject'} = $value;

  return $oldValue;
}

##-------------------------------------------------------------------------##

=over 4

=item Use: my $value = getQuery( );

=item Use: my $oldValue = setQuery( $value );

Get/Set the query database value.  

=back 

=cut

##-------------------------------------------------------------------------##
sub getQuery {
  my $this = shift;

  return $this->{'query'};
}

sub setQuery {
  my $this  = shift;
  my $value = shift;

  my $oldValue = $this->{'query'};
  $this->{'query'} = $value;

  return $oldValue;
}

##-------------------------------------------------------------------------##
## General Object Methods
##-------------------------------------------------------------------------##

##-------------------------------------------------------------------------##

=over 4

=item Use: my $value = getVersion( );

Get the search engine binary version.  

=back 

=cut

##-------------------------------------------------------------------------##
sub getVersion {
  my $this = shift;

  return $this->{'version'};
}

##-------------------------------------------------------------------------##

=over 4

=item Use: my ( $resultCode, $SearchResultCollectionI ) = search( );

=item Use: my ( $resultCode, $SearchResultCollectionI )
                          = search( matrix=>"7p16g.matrix",
                                    ...
                                  );

Run the search and return a SearchResultsCollectionI.

=back 

=cut

##-------------------------------------------------------------------------##
sub search {
  croak $CLASS. "::search() not implemented!\n";
}

##-------------------------------------------------------------------------##

=over 4

=item Use: my $SearchResultCollection = parseOutput( 
                                     searchOutput => $filename|$FH,
                                     [excludeAlignments => 1],
                                     [scoreHighThresh => #],
                                     [scoreLowThresh => #],
                                     [subjPattern => ""],
                                     [queryPattern => ""] );

Parse the result of a search and return a SearchResultCollection.

=back 

=cut

##-------------------------------------------------------------------------##
sub parseOutput {
  croak $CLASS. "::parseOutput() not implemented!\n";
}

##-------------------------------------------------------------------------##
## Private Methods
##-------------------------------------------------------------------------##

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
