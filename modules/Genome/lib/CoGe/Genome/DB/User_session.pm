package CoGe::Genome::DB::User_session;
use strict;
use base 'CoGe::Genome::DB';

BEGIN {
    use Exporter ();
    use vars qw ($VERSION @ISA @EXPORT @EXPORT_OK %EXPORT_TAGS);
    $VERSION     = 0.1;
    @ISA         = (@ISA, qw (Exporter));
    #Give a hoot don't pollute, do not export more than needed by default
    @EXPORT      = qw ();
    @EXPORT_OK   = qw ();
    %EXPORT_TAGS = ();
    __PACKAGE__->set_up_table('user_session');
    __PACKAGE__->set_sql('add_user_session'=> qq{
INSERT INTO user_session (user_id, date)
 VALUES (?, NOW());
});
 }


########################################### main pod documentation begin ##
# Below is the stub of documentation for your module. You better edit it!


=head1 NAME

CoGe::Genome::DB::User_session

=head1 SYNOPSIS

=head1 DESCRIPTION

This object accesses the user_session table.  The feature_list contains a list of 
user_ids and session_ids that can be used for tracking users.


=head1 USAGE



=head1 BUGS



=head1 SUPPORT



=head1 AUTHOR

	Eric Lyons
	elyons@nature.berkeley.edu

=head1 COPYRIGHT

This program is free software licensed under the...

	The Artistic License

The full text of the license can be found in the
LICENSE file included with this module.


=head1 SEE ALSO

CoGe::Genome
CoGe::Genome::DB
Class::DBI


perl(1).

=cut

############################################# main pod documentation end ##


=head2 Accessor Functions

 user_id          =>  database id for user
 u_id
 uid

 date             =>  timestamp of when the entry was created

 user_session_id  =>  database entry id
 id               =>  alias for feature_list_id

 

 new              =>  creates a new object (inherited from Class::Accessor)

=cut


sub desc
  {
    my $self = shift;
    return $self->description(@_);
  }

sub id
  {
    my $self = shift;
    return $self->user_session_id();
  }

sub log_user
  {
    my $self = shift;
    my %opts = @_;
    my $user = $opts{user};
    my $uid = ref($user) =~ /User/ ? $user->id : $user;
    unless ($uid =~ /^\d+$/)
      {
	warn "Error adding user_id to User_session.  Not a valid uid: $uid\n";
	return;
      }
    #FIRST REMOVE ALL ENTRIES FOR THIS USER
    foreach my $item ($self->retrieve(user_id=>$uid))
      {
	next unless $item;
	$item->delete;
      }
    #ADD NEW ENTRY
    my $sth = $self->sql_add_user_session;
    $sth->execute($uid);
    $sth->finish;
    my ($item) = $self->retrieve(user_id=>$uid);
    return $item->id;
  }

################################################ subroutine header begin ##

=head2 

 Usage     : 
 Purpose   : 
 Returns   : 
 Argument  : 
 Throws    : 
 Comments  : 

See Also   : 

=cut

################################################## subroutine header end ##

1; #this line is important and will help the module return a true value

