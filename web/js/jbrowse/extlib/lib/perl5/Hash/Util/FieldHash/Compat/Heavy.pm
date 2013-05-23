#!/usr/bin/perl

package Hash::Util::FieldHash::Compat;

use strict;
use warnings;

use Tie::RefHash::Weak;

use vars qw(@ISA %EXPORT_TAGS @EXPORT_OK);

@ISA = qw(Exporter);

%EXPORT_TAGS = (
    'all' => [ qw(
        fieldhash
        fieldhashes
        idhash
        idhashes
        id
        id_2obj
        register
    )],
);

@EXPORT_OK = ( @{ $EXPORT_TAGS{'all'} } );

sub fieldhash (\%) {
	my $hash = shift;

	tie %$hash, 'Hash::Util::FieldHash::Compat::Tie::FieldHash', %$hash;

	return $hash;
}

sub fieldhashes { map { &fieldhash($_) } @_ }

sub idhash (\%) {
	tie %{$_[0]}, 'Hash::Util::FieldHash::Compat::Tie::IdHash', %{$_[0]};
	$_[0];
}

sub idhashes { map { &idhash($_) } @_ }

sub id ($) {
	my $obj = shift;

	if ( defined ( my $refaddr = Tie::RefHash::refaddr($obj) ) ) {
		return $refaddr;
	} else {
		return $obj;
	}
}

tie my %registry, 'Tie::RefHash::Weak';

sub id_2obj {
	my $id = shift;

	my $registry_by_id = tied(%registry)->[0];

	if ( my $record = $registry_by_id->{$id} ) {
		return $record->[0]; # first slot is the key
	}

	return;
}

sub register {
	my ( $obj, @args ) = @_;
	( $registry{$obj} ||= Hash::Util::FieldHash::Compat::Destroyer->new($obj) )->register(@args);
}

{
	package Hash::Util::FieldHash::Compat::Tie::IdHash;
	use Tie::Hash ();
	use vars qw(@ISA);
	@ISA = qw(Tie::StdHash);

	# this class always stringifies using id().
	
	sub TIEHASH {
		my ( $class, @args ) = @_;
		my $self = bless {}, $class;

		while ( @args ) {
			my ( $key, $value ) = splice @args, 0, 2;
			$self->STORE($key, $value);
		}

		$self;
	}

	BEGIN {
		foreach my $method ( qw(STORE FETCH DELETE EXISTS) ) {
			eval 'sub '.$method.' {
				my ( $self, $key, @args ) = @_;
				$self->SUPER::'.$method.'( Hash::Util::FieldHash::Compat::id($key), @args );
			}';
		}
	}

	package Hash::Util::FieldHash::Compat::Tie::FieldHash;
	use vars qw(@ISA);
	@ISA = qw(Tie::RefHash::Weak);

	# this subclass retains weakrefs to the objects in the keys, but pretends
	# the keys are actually strings

	BEGIN {
		# always return strings from keys

		foreach my $method ( qw(FIRSTKEY NEXTKEY) ) {
			eval 'sub '.$method.' {
				my ( $self, @args ) = @_;
				Hash::Util::FieldHash::Compat::id($self->SUPER::'.$method.'(@args));
			}';
		}

		sub EXISTS {
			my ( $self, $key ) = @_;
			my $str_key = Hash::Util::FieldHash::Compat::id($key);
			exists $_->{$str_key} and return 1 for @{ $self }[0, 1];
			return;
		}

		sub FETCH {
			my($self, $key) = @_;

			my $str_key = Hash::Util::FieldHash::Compat::id($key);

			if ( exists $self->[0]{$str_key} ) {
				return $self->[0]{$str_key}[1];
			} else {
				$self->[1]{$str_key};
			}
		}

		sub STORE {
			my ( $self, $key, $value ) = @_;

			my $str_key = Hash::Util::FieldHash::Compat::id($key);

			delete $self->[1]{$str_key};

			$self->SUPER::STORE( $key, $value );
		}

		sub DELETE {
			my ( $self, $key ) = @_;

			foreach my $key ( $key, Hash::Util::FieldHash::Compat::id($key) ) {
				if ( defined ( my $ret = $self->SUPER::DELETE($key) ) ) {
					return $ret;
				}
			}
		}
	}

	package Hash::Util::FieldHash::Compat::Destroyer;
	use Scalar::Util qw(weaken);

	sub new {
		my ( $class, $obj ) = @_;

		tie my %hashes, 'Tie::RefHash::Weak';

		my $self = bless {
			object => $obj,
			hashes => \%hashes,
		}, $class;

		weaken($self->{object});

		$self;
	}

	sub register {
		my ( $self, @hashes ) = @_;
		$self->{hashes}{$_}++ for @hashes;
	}

	sub DESTROY {
		my $self = shift;
		my $object = $self->{object};
		delete $_->{Hash::Util::FieldHash::Compat::id($object)} for keys %{ $self->{hashes} };
	}
}

__PACKAGE__

__END__

=pod

=head1 NAME

Hash::Util::FieldHash::Compat::Heavy - Emulate Hash::Util::FieldHash using
L<Tie::RefHash> etc.

=head1 SYNOPSIS

	# this module will be used automatically by L<Hash::Util::FieldHash::Compat> if necessary

=head1 DESCRIPTION

See L<Hash::Util::FieldHash::Compat> for the documentation

=cut


