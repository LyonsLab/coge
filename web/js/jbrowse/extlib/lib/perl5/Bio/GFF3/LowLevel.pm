package Bio::GFF3::LowLevel;
BEGIN {
  $Bio::GFF3::LowLevel::AUTHORITY = 'cpan:RBUELS';
}
{
  $Bio::GFF3::LowLevel::VERSION = '1.5';
}
# ABSTRACT: fast, low-level functions for parsing and formatting GFF3

use strict;

use Scalar::Util ();
use URI::Escape ();


require Exporter;
our @ISA = qw(Exporter);
our @EXPORT_OK = qw(
  gff3_parse_feature
  gff3_parse_attributes
  gff3_parse_directive
  gff3_format_feature
  gff3_format_attributes
  gff3_escape
  gff3_unescape
);

my @gff3_field_names = qw(
    seq_id
    source
    type
    start
    end
    score
    strand
    phase
    attributes
);


sub gff3_parse_feature {
  my ( $line ) = @_;
  no warnings 'uninitialized';

  my @f = split /\t/, $line;
  for( @f ) {
      if( $_ eq '.' ) {
          $_ = undef;
      }
  }

  # unescape only the ref and source columns
  $f[0] =~ s/%([0-9A-Fa-f]{2})/chr(hex($1))/eg;
  $f[1] =~ s/%([0-9A-Fa-f]{2})/chr(hex($1))/eg;

  $f[8] = gff3_parse_attributes( $f[8] );
  my %parsed;
  @parsed{@gff3_field_names} = @f;
  return \%parsed;
}


sub gff3_parse_attributes {
    my ( $attr_string ) = @_;

    return {} if !defined $attr_string || $attr_string eq '.';

    $attr_string =~ s/\r?\n$//;

    my %attrs;
    for my $a ( split ';', $attr_string ) {
        no warnings 'uninitialized';
        my ( $name, $values ) = split '=', $a, 2;
        next unless defined $values;
        push @{$attrs{$name}}, map { s/%([0-9A-Fa-f]{2})/chr(hex($1))/eg; $_ } split ',', $values;
    }

    return \%attrs;
}


sub gff3_parse_directive {
    my ( $line ) = @_;

    my ( $name, $contents ) = $line =~ /^ \s* \#\# \s* (\S+) \s* (.*) $/x
        or return;

    my $parsed = { directive => $name };
    if( length $contents ) {
        $contents =~ s/\r?\n$//;
        $parsed->{value} = $contents;
    }

    # do a little additional parsing for sequence-region and genome-build directives
    if( $name eq 'sequence-region' ) {
        my ( $seqid, $start, $end ) = split /\s+/, $contents, 3;
        s/\D//g for $start, $end;
        @{$parsed}{qw( seq_id start end )} = ( $seqid, $start, $end );
    }
    elsif( $name eq 'genome-build' ) {
        my ( $source, $buildname ) = split /\s+/, $contents, 2;
        @{$parsed}{qw(source buildname)} = ( $source, $buildname );
    }

    return $parsed;
}


sub gff3_format_feature {
    my ( $f ) = @_;

    my $attr_string = $f->{attributes};
    $attr_string = '.' unless defined $attr_string;

    $attr_string = gff3_format_attributes( $attr_string )
        if    ref( $attr_string ) eq 'HASH'
           && ! Scalar::Util::blessed( $attr_string );

    return join( "\t",
                 ( map { defined $_ ? gff3_escape($_) : '.' }
                   @{$f}{@gff3_field_names[0..7]}
                 ),
                 $attr_string
               )."\n";
}


my %force_attr_first = (
    ID     => 1,
    Name   => 2,
    Alias  => 3,
    Parent => 4,
);
sub _cmp_attr_names {
    no warnings 'uninitialized';
    my ( $fa, $fb ) = @force_attr_first{ $a, $b };
    return $fa <=> $fb if $fa && $fb;

    return -1 if  $fa && !$fb;
    return  1 if !$fa &&  $fb;

    return $a cmp $b;
}

sub gff3_format_attributes {
  my ( $attr ) = @_;

  return '.' unless defined $attr;

  my $astring = join ';' => (
    map {
      my $key = $_;
      my $val = $attr->{$key};
      no warnings 'uninitialized';
      $val = join( ',', map gff3_escape($_), ref $val eq 'ARRAY' ? @$val : $val );
      if( length $val ) {
          "$key=$val"
      } else {
          ()
      }
    }
    sort _cmp_attr_names
    keys %$attr
 );

  return length $astring ? $astring : '.';
}


sub gff3_escape {
    URI::Escape::uri_escape( $_[0], '\n\r\t;=%&,\x00-\x1f\x7f-\xff' )
}


*gff3_unescape = \&URI::Escape::uri_unescape;


__END__
=pod

=encoding utf-8

=head1 NAME

Bio::GFF3::LowLevel - fast, low-level functions for parsing and formatting GFF3

=head1 SYNOPSIS

  use Bio::GFF3::LowLevel qw/ gff3_parse_feature /;

  open my $gff3_fh, 'myfile.gff3' or die;
  while( <$gff3_fh> ) {
    next if /^#/;
    my $feat = gff3_parse_feature( $_ );
  }

=head1 DESCRIPTION

These are low-level, fast functions for parsing GFF version 3 files.
All they do is convert back and forth between low-level Perl data
structures and GFF3 text.

Sometimes this is what you need when you are just doing simple
transformations on GFF3.  I found myself writing these functions over
and over again, until I finally got fed up enough to just package them
up properly.

These functions do no validation, do not reconstruct feature
hierarchies, or anything like that.  If you want that, use
L<Bio::FeatureIO>.

All of the functions in this module are EXPORT_OK, meaning that you
can add their name after using this module to make them available in
your namespace.

=head1 FUNCTIONS

=head2 gff3_parse_feature( $line )

Given a string containing a GFF3 feature line (i.e. not a comment),
parses it and returns a hashref of its information, of the form:

    {
        seq_id => 'chr02',
        source => 'AUGUSTUS',
        type   => 'transcript',
        start  => '23486',
        end    => '48209',
        score  => '0.02',
        strand => '+',
        phase  => undef,
        attributes => {
            ID => [
                'chr02.g3.t1'
              ],
            Parent => [
                'chr02.g3'
              ],
          },
    }

Note that all values are simple scalars, except for C<attributes>,
which is a hashref as returned by L</gff3_parse_attributes> below.

Unescaping is performed according to the GFF3 specification.

=head2 gff3_parse_attributes( $attr_string )

Given a GFF3 attribute string, parse it and return a hashref of its
data, of the form:

    {
      'attribute_name' => [ value, value, ... ],
      ...
    }

Always returns a hashref.  If the passed attribute string is
undefined, or ".", the hashref returned will be empty.  Attribute
values are always arrayrefs, even if they have only one value.

=head2 gff3_parse_directive( $line )

Parse a GFF3 directive/metadata line.  Returns a hashref as:

  {  directive => 'directive-name',
     value     => 'the contents of the directive'
  }

Or nothing if the line could not be parsed as a GFF3 directive.

In addition, C<sequence-region> and C<genome-build> directives are
parsed further.  C<sequence-region> hashrefs have additional
C<seq_id>, C<start>, and C<end> keys, and C<genome-build> hashrefs
have additional C<source> and C<buildname> keys

=head2 gff3_format_feature( \%fields )

Given a hashref of feature information in the same format returned by
L</gff3_parse_feature> above, constructs a correctly-escaped line of
GFF3 encoding that information.

The line ends with a single newline character, a UNIX-style line
ending, regardless of the local operating system.

=head2 gff3_format_attributes( \%attrs )

Given a hashref of GFF3 attributes in the same format returned by
L</gff3_parse_attributes> above, returns a correctly formatted and
escaped GFF3 attribute string (the 9th column of a GFF3 feature line)
encoding those attributes.

For convenience, single-valued attributes can have simple scalars as
values in the passed hashref.  For example, if a feature has only one
C<ID> attribute (as it should), you can pass C<{ ID =E<gt> 'foo' }>
instead of C<{ ID =E<gt> ['foo'] }}>.

=head2 gff3_escape( $string )

Given a string, escapes special characters in that string according to
the GFF3 specification.

=head2 gff3_unescape( $string )

Unescapes a GFF3-escaped string.

=head1 AUTHOR

Robert Buels <rmb32@cornell.edu>

=head1 COPYRIGHT AND LICENSE

This software is copyright (c) 2012 by Robert Buels.

This is free software; you can redistribute it and/or modify it under
the same terms as the Perl 5 programming language system itself.

=cut

