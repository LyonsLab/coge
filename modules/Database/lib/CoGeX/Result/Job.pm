package CoGeX::Result::Job;

use strict;
use warnings;

use base 'DBIx::Class::Core';


=head1 NAME

My::Schema::Result::Log

=cut

__PACKAGE__->table("job");

=head1 ACCESSORS

=cut

__PACKAGE__->add_columns(
	"job_id",
	{ data_type => "INT", default_value => undef, is_nullable => 0, size => 11 },	
	"user_id",
	{ data_type => "INT", default_value => 0, is_nullable => 0, size => 11 },
	"page",
	{ data_type => "VARCHAR", default_value => undef, is_nullable => 0, size => 255 },
	"start_time",
	{ data_type => "TIMESTAMP", default_value => undef, is_nullable => 0 },
	"link",
	{ data_type => "VARCHAR", default_value => undef, is_nullable => 1, size => 255 },
	"status",
	{ data_type => "TINYINT", default_value => 0, is_nullable => 0, size => 5 },	
	"type",
	{ data_type => "TINYINT", default_value => 0, is_nullable => 0, size => 5 },	
	"process_id",
	{ data_type => "INT", default_value => 0, is_nullable => 0, size => 5 },	
	"log_id",
	{ data_type => "INT", default_value => 0, is_nullable => 0, size => 11 },	
);
__PACKAGE__->set_primary_key("job_id");

__PACKAGE__->belongs_to('user' => "CoGeX::Result::User", 'user_id');
__PACKAGE__->belongs_to('log' => "CoGeX::Result::Log", 'log_id');


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

1;
