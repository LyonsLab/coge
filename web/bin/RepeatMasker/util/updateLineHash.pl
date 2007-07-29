#!/usr/bin/perl -w
##---------------------------------------------------------------------------##
##  File:
##      @(#) updateLineHash.pl
##  Author:
##      Robert M. Hubley   rhubley@systemsbiology.org
##  Description:
##      Update the RepeatAnnotationData.pm file.
##
#******************************************************************************
#*  This software is provided ``AS IS'' and any express or implied            *
#*  warranties, including, but not limited to, the implied warranties of      *
#*  merchantability and fitness for a particular purpose, are disclaimed.     *
#*  In no event shall the authors or the Institute for Systems Biology        *
#*  liable for any direct, indirect, incidental, special, exemplary, or       *
#*  consequential damages (including, but not limited to, procurement of      *
#*  substitute goods or services; loss of use, data, or profits; or           *
#*  business interruption) however caused and on any theory of liability,     *
#*  whether in contract, strict liability, or tort (including negligence      *
#*  or otherwise) arising in any way out of the use of this software, even    *
#*  if advised of the possibility of such damage.                             *
#*                                                                            *
#******************************************************************************
#
# ChangeLog
#
#     $Log: updateLineHash.pl,v $
#     Revision 1.3  2007/05/17 21:01:54  rhubley
#     Cleanup before a distribution
#
#
###############################################################################
#
# To Do:
#

=head1 NAME

updateLineHash.pl - Update the RepeatAnnnotationData.pm

=head1 SYNOPSIS

  updateLineHash.pl [-version]

=head1 DESCRIPTION

The options are:

=over 4

=item -version

Displays the version of the program

=back

=head1 SEE ALSO

=head1 COPYRIGHT

Copyright 2006 Robert Hubley, Institute for Systems Biology

=head1 AUTHOR

Robert Hubley <rhubley@systemsbiology.org>

=cut

#
# Module Dependence
#
use strict;
use FindBin;
use lib $FindBin::RealBin;
use lib "$FindBin::RealBin/../";
use Getopt::Long;
use RepeatAnnotationData;
use CrossmatchSearchEngine;
use Data::Dumper;

#
# Magic numbers/constants here
#  ie. my $PI = 3.14159;
#
my $DEBUG   = 0;
my $Version = "0.1";

#
# Option processing
#  e.g.
#   -t: Single letter binary option
#   -t=s: String parameters
#   -t=i: Number paramters
#
my @getopt_args = (
                    '-version',    # print out the version and exit
                    '-help',
                    '-q=s',
                    '-equate',
                    '-replace',
                    '-equiv'
);

my %options = ();
Getopt::Long::config( "noignorecase", "bundling_override" );
unless ( GetOptions( \%options, @getopt_args ) ) {
  usage();
}

sub usage {
  print "$0 - $Version\n";
  exec "pod2text $0";
  exit;
}

if ( $options{'version'} ) {
  print "$Version\n";
  exit;
}

&usage() if ( $options{'help'} );

my $lineHash = \%RepeatAnnotationData::lineHash;

if ( $options{'equiv'} ) {
  my $modifyElementName = shift @ARGV;

  # If doesn't exist add one
  if ( !defined $lineHash->{$modifyElementName} ) {
    $lineHash->{$modifyElementName} = {};
  }
  elsif ( $options{'replace'} ) {
    foreach my $equiv ( keys( %{ $lineHash->{$modifyElementName} } ) ) {
      delete $lineHash->{$equiv}->{$modifyElementName};
    }
    $lineHash->{$modifyElementName} = {};
  }

  foreach my $equiv ( @ARGV ) {
    print "Making $modifyElementName equiv to $equiv\n";
    $lineHash->{$modifyElementName}->{$equiv} = { 'overThresh' => 0 };

    # Now the reciprocal
    $lineHash->{$equiv}->{$modifyElementName} = { 'overThresh' => 0 };
  }
  print $modifyElementName . ":\n";
  foreach my $equiv ( keys( %{ $lineHash->{$modifyElementName} } ) ) {
    print "    $equiv";
    if ( $lineHash->{$modifyElementName}->{$equiv}->{'overThresh'} > 0 ) {
      print "   ** if overlap > "
          . $lineHash->{$modifyElementName}->{$equiv}->{'overThresh'};
    }
    print "\n";
  }

##
## Output results
##
  open OUT, ">RepeatAnnotationDataModified.pm";
  print OUT "package RepeatAnnotationData;\n";
  print OUT "require Exporter;\n";
  print OUT "\@EXPORT_OK = qw( \%repeatDB \%lineHash \%preProcData );\n";
  print OUT "\%EXPORT_TAGS = ( all => [ \@EXPORT_OK ] );\n";
  print OUT "\@ISA         = qw(Exporter);\n";
  print OUT "\n";
  print OUT "BEGIN {\n";
  print OUT "\n";
  print OUT "  \%repeatDB = (\n";
  my $output = Dumper( \%RepeatAnnotationData::repeatDB );
  $output =~ s/^\$VAR1 = \{//;
  $output =~ s/\};//;
  print OUT $output;
  print OUT "\n";
  print OUT "  );\n";
  print OUT "  \%lineHash = (\n";
  $output = Dumper( \%RepeatAnnotationData::lineHash );
  $output =~ s/^\$VAR1 = \{//;
  $output =~ s/\};//;
  print OUT $output;
  print OUT "\n";
  print OUT "  );\n";
  print OUT "  \%preProcData = (\n";
  $output = Dumper( \%RepeatAnnotationData::preProcData );
  $output =~ s/^\$VAR1 = \{//;
  $output =~ s/\};//;
  print OUT $output;
  print OUT "\n";
  print OUT "  );\n";
  print OUT "}\n";
  close OUT;

}
##
##  Add/Edit element and make it equivalent to another
##  element.
##
elsif ( $options{'equate'} ) {
  my $modifyElementName   = $ARGV[ 0 ];
  my $templateElementName = $ARGV[ 1 ];

  # If doesn't exist add one
  if ( !defined $lineHash->{$modifyElementName} ) {
    $lineHash->{$modifyElementName} = {};
  }
  elsif ( $options{'replace'} ) {
    foreach my $equiv ( keys( %{ $lineHash->{$modifyElementName} } ) ) {
      delete $lineHash->{$equiv}->{$modifyElementName};
    }
    $lineHash->{$modifyElementName} = {};
  }

  if ( defined $lineHash->{$templateElementName} ) {
    foreach my $equiv ( keys( %{ $lineHash->{$templateElementName} } ) ) {
      if ( defined $lineHash->{$modifyElementName}->{$equiv} ) {
        if ( $lineHash->{$modifyElementName}->{$equiv}->{'overThresh'} !=
             $lineHash->{$templateElementName}->{$equiv}->{'overThresh'} )
        {
          warn "Relationship $modifyElementName->$equiv already defined but "
              . "with a different overThresh parameter.";
        }
        else {
          warn "Relationship $modifyElementName->$equiv already defined";
        }
      }
      $lineHash->{$modifyElementName}->{$equiv} =
          { 'overThresh' =>
            $lineHash->{$templateElementName}->{$equiv}->{'overThresh'} };

      # Now the reciprocal
      $lineHash->{$equiv}->{$modifyElementName} =
          { 'overThresh' =>
            $lineHash->{$templateElementName}->{$equiv}->{'overThresh'} };
    }
  }
  else {
    warn "Template element $templateElementName does not exist "
        . "in the lineHash!\n";
  }

  print $modifyElementName . ":\n";
  foreach my $equiv ( keys( %{ $lineHash->{$modifyElementName} } ) ) {
    print "    $equiv";
    if ( $lineHash->{$modifyElementName}->{$equiv}->{'overThresh'} > 0 ) {
      print "   ** if overlap > "
          . $lineHash->{$modifyElementName}->{$equiv}->{'overThresh'};
    }
    print "\n";
  }

##
## Output results
##
  open OUT, ">RepeatAnnotationDataModified.pm";
  print OUT "package RepeatAnnotationData;\n";
  print OUT "require Exporter;\n";
  print OUT "\@EXPORT_OK = qw( \%repeatDB \%lineHash \%preProcData );\n";
  print OUT "\%EXPORT_TAGS = ( all => [ \@EXPORT_OK ] );\n";
  print OUT "\@ISA         = qw(Exporter);\n";
  print OUT "\n";
  print OUT "BEGIN {\n";
  print OUT "\n";
  print OUT "  \%repeatDB = (\n";
  my $output = Dumper( \%RepeatAnnotationData::repeatDB );
  $output =~ s/^\$VAR1 = \{//;
  $output =~ s/\};//;
  print OUT $output;
  print OUT "\n";
  print OUT "  );\n";
  print OUT "  \%lineHash = (\n";
  $output = Dumper( \%RepeatAnnotationData::lineHash );
  $output =~ s/^\$VAR1 = \{//;
  $output =~ s/\};//;
  print OUT $output;
  print OUT "\n";
  print OUT "  );\n";
  print OUT "  \%preProcData = (\n";
  $output = Dumper( \%RepeatAnnotationData::preProcData );
  $output =~ s/^\$VAR1 = \{//;
  $output =~ s/\};//;
  print OUT $output;
  print OUT "\n";
  print OUT "  );\n";
  print OUT "}\n";
  close OUT;
}
elsif ( $options{'q'} ) {

  #print "" . Dumper( $lineHash->{ $options{'q'} } ) . "\n";
  print $options{'q'} . ":\n";
  foreach my $equiv ( keys( %{ $lineHash->{ $options{'q'} } } ) ) {
    print "    $equiv";
    if ( $lineHash->{ $options{'q'} }->{$equiv}->{'overThresh'} > 0 ) {
      print "   ** if overlap > "
          . $lineHash->{ $options{'q'} }->{$equiv}->{'overThresh'};
    }
    print "\n";
  }
}

## Done
exit;

######################## S U B R O U T I N E S ############################

##-------------------------------------------------------------------------##
## Use: my getOverlapSize( $range1Begin, $range1End,
##                         $range2Begin, $range2End );
##
##      $range1Begin       : Range1 Start Position
##      $range1End         : Range1 End Position
##      $range2Begin       : Range2 Start Position
##      $range2End         : Range2 End Position
##
##  Returns
##      The number of positions that the ranges overlap or
##      0 if they don't.
##
##-------------------------------------------------------------------------##
sub getOverlapSize {
  my $range1Begin = shift;
  my $range1End   = shift;
  my $range2Begin = shift;
  my $range2End   = shift;

  my $overlap = 0;
  if (    $range1Begin >= $range2Begin
       && $range1Begin <= $range2End )
  {

    #      -------
    #   ------
    # or
    #     -----
    #   --------
    if ( $range1End <= $range2End ) {

      #     -----
      #   --------
      $overlap = $range1End - $range1Begin + 1;
    }
    else {

      #      -------
      #   ------
      $overlap = $range2End - $range1Begin + 1;
    }
  }
  elsif (    $range1End >= $range2Begin
          && $range1End <= $range2End )
  {

    #   -------
    #      ------
    # or
    #   --------
    #    -----
    if ( $range1End <= $range2End ) {

      #   --------
      #    -----
      $overlap = $range2End - $range2Begin + 1;
    }
    else {

      #   -------
      #      ------
      $overlap = $range1End - $range2Begin + 1;
    }
  }
  return $overlap;
}

1;
