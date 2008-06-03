package CoGeX::User;

# Created by DBIx::Class::Schema::Loader v0.03009 @ 2006-12-01 18:13:38

use strict;
use warnings;

use base 'DBIx::Class';

__PACKAGE__->load_components("PK::Auto", "Core");
__PACKAGE__->table("user");
__PACKAGE__->add_columns(
  "user_id",
  { data_type => "INT", default_value => undef, is_nullable => 0, size => 10 },
  "user_name",
  { data_type => "VARCHAR", default_value => "", is_nullable => 0, size => 10 },
  "first_name",
  { data_type => "VARCHAR", default_value => "", is_nullable => 0, size => 10 },
  "last_name",
  { data_type => "VARCHAR", default_value => "", is_nullable => 0, size => 10 },
  "email",
  {
    data_type => "VARCHAR",
    default_value => undef,
    is_nullable => 1,
    size => 50,
  },
  "description",
  {
    data_type => "VARCHAR",
    default_value => undef,
    is_nullable => 1,
    size => 255,
  },
  "passwd",
  { data_type => "VARCHAR", default_value => "", is_nullable => 0, size => 255 },
);
__PACKAGE__->set_primary_key("user_id");
__PACKAGE__->has_many('sessions'=>"CoGeX::UserSession",'user_id');


sub generate_passwd
  {
    my $self = shift;
    my %opts = @_;
    my $pwd = $opts{passwd} || $opts{pwd};
    my $crypt_pwd = crypt( $pwd, "12" );
  }

#sub to check if entered password matches user passwd

sub check_passwd
  {
    my $self = shift;
    my %opts = @_;
    my $pwd = $opts{passwd} || $opts{pwd};
    return crypt($pwd, $self->passwd) eq $self->passwd;
  }
1;

