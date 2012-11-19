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
	"log_id",
	{ data_type => "INT", default_value => undef, is_nullable => 0, size => 11 },	
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
	"status",
	{ data_type => "INT", default_value => 0, is_nullable => 0, size => 1 },	
	"comment",
	{ data_type => "VARCHAR", default_value => undef, is_nullable => 1, size => 255 },	
);
__PACKAGE__->set_primary_key("log_id");

__PACKAGE__->belongs_to('user' => "CoGeX::Result::User", 'user_id');


################################################ subroutine header begin ##

=head2 info

 Usage     : $self->info
 Purpose   : generate a string of information about the log entry
 Returns   : a string
 Argument  : None
 Throws    : None
 Comments  : 

=cut

################################################## subroutine header end ##

sub info {
	my $self = shift;
	my $user_name = ($self->user ? $self->user->user_name : '');
	return $self->time.' '.$user_name.' '.$self->page.' '.$self->description.' '.$self->link.' '.$self->comment; 
}


1;
