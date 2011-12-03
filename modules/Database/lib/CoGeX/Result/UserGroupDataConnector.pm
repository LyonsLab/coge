package CoGeX::Result::UserGroupDataConnector;

# Created by DBIx::Class::Schema::Loader
# DO NOT MODIFY THE FIRST PART OF THIS FILE

use strict;
use warnings;

use base 'DBIx::Class::Core';


=head1 NAME

My::Schema::Result::GroupDataConnector

=cut

__PACKAGE__->table("user_group_data_connector");

=head1 ACCESSORS

=head2 group_dataset_connector_id

  data_type: 'bigint'
  is_auto_increment: 1
  is_nullable: 0

=head2 dataset_group_id

  data_type: 'bigint'
  is_nullable: 0

=head2 dataset_id

  data_type: 'bigint'
  is_nullable: 0

=head2 group_id

  data_type: 'bigint'
  is_nullable: 0

=cut

__PACKAGE__->add_columns(
  "user_group_dataset_connector_id",  { data_type => "INT", default_value => undef, is_nullable => 0, size => 11},
  "dataset_group_id",  {data_type => "INT", default_value => undef, is_nullable => 1, size => 11},
  "dataset_id",  { data_type => "INT", default_value => undef, is_nullable => 1, size => 11 },
  "user_group_id",  { data_type => "INT", default_value => undef, is_nullable => 0, size => 11 },
);
__PACKAGE__->set_primary_key("user_group_dataset_connector_id");
__PACKAGE__->belongs_to("genome" => "CoGeX::Result::DatasetGroup", "dataset_group_id");
__PACKAGE__->belongs_to("dataset" => "CoGeX::Result::Dataset", "dataset_id");
__PACKAGE__->belongs_to("user_group" => "CoGeX::Result::UserGroup","user_group_id");


# Created by DBIx::Class::Schema::Loader v0.07002 @ 2011-08-29 09:28:12
# DO NOT MODIFY THIS OR ANYTHING ABOVE! md5sum:lcDzeTaaZSsmL3e0hILR8w


# You can replace this text with custom content, and it will be preserved on regeneration
1;
