package CoGeX::Result::Log;

use strict;
use warnings;

use base 'DBIx::Class::Core';


=head1 NAME

My::Schema::Result::Log

=cut

__PACKAGE__->table("log");

=head1 ACCESSORS

=cut

__PACKAGE__->add_columns(
	"time",
	{ data_type => "TIMESTAMP", default_value => undef, is_nullable => 0 },
	"user_id",
	{ data_type => "INT", default_value => 0, is_nullable => 0, size => 11 },
#   "target", # table name
#	{ data_type => "VARCHAR", default_value => undef, is_nullable => 1, size => 30 },
	"page",
	{ data_type => "VARCHAR", default_value => undef, is_nullable => 0, size => 255 },
	"description",
	{ data_type => "VARCHAR", default_value => undef, is_nullable => 1, size => 255 },
	"link",
	{ data_type => "VARCHAR", default_value => undef, is_nullable => 1, size => 255 },
);
__PACKAGE__->belongs_to('user' => "CoGeX::Result::User", 'user_id');



1;
