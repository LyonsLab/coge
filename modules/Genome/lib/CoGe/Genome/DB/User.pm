package CoGe::Genome::DB::User;
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
    __PACKAGE__->set_up_table('user');
 }


########################################### main pod documentation begin ##
# Below is the stub of documentation for your module. You better edit it!


=head1 NAME

CoGe::Genome::DB::User

=head1 SYNOPSIS

=head1 DESCRIPTION

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

 user_name        =>  name of user
 u_name
 uname
 name

 first_name
 f_name
 fname

 last_name
 l_name
 last_name

 description      =>  description
 desc             =>  alias for description

 passwd
 password
 pwd
 pw

 user_id          =>  database entry id
 id               =>  alias for sequence_type_id

 

 new              =>  creates a new object (inherited from Class::Accessor)

=cut

sub u_name
  {
    my $self = shift;
    return $self->user_name(@_);
  }

sub uname
  {
    my $self = shift;
    return $self->user_name(@_);
  }

sub name
  {
    my $self = shift;
    return $self->user_name(@_);
  }

sub f_name
  {
    my $self = shift;
    return $self->first_name(@_);
  }

sub fname
  {
    my $self = shift;
    return $self->first_name(@_);
  }

sub l_name
  {
    my $self = shift;
    return $self->last_name(@_);
  }

sub lname
  {
    my $self = shift;
    return $self->last_name(@_);
  }

sub password
  {
    my $self = shift;
    return $self->passwd(@_);
  }

sub pwd
  {
    my $self = shift;
    return $self->passwd(@_);
  }

sub pw
  {
    my $self = shift;
    return $self->passwd(@_);
  }


sub desc
  {
    my $self = shift;
    return $self->description(@_);
  }

sub id
  {
    my $self = shift;
    return $self->user_id();
  }

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
    return crypt($pwd, $self->pwd) eq $self->pwd;
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

