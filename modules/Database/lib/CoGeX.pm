package CoGeX;

use strict;
use warnings;

use vars qw( $VERSION );

$VERSION = 0.01;

use base 'DBIx::Class::Schema';

__PACKAGE__->load_classes;

=head1 NAME

CoGeX - CoGeX

=head1 SYNOPSIS

  use CoGeX;
  blah blah blah


=head1 DESCRIPTION

Primary object for interacting with CoGe database system.

=head1 USAGE

  use CoGeX;

  my $connstr = 'dbi:mysql:genomes:biocon:3306';
  my $s = CoGeX->connect($connstr, 'cnssys', 'CnS' ); # Biocon's ro user

  my $rs = $s->resultset('Feature')->search(
                {
                  'organism.name' => "Arabidopsis thaliana"
                },
                { join => [ 'dataset', 'organism' ] }
  );

=head1 BUGS


=head1 SUPPORT


=head1 AUTHOR

=head1 COPYRIGHT

This program is free software; you can redistribute
it and/or modify it under the same terms as Perl itself.

The full text of the license can be found in the
LICENSE file included with this module.


=head1 SEE ALSO

perl(1).

=cut

1;
