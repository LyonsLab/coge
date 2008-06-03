package CoGe::Accessory::LogUser;

use strict;
use CGI::Cookie;
use Data::Dumper;
use CoGeX;
use CoGeX::User;

use vars qw($cookie_name $coge);

$cookie_name = "CoGe";
$coge = CoGeX->dbconnect();
sub get_user
  {
    my $self = shift;
    my %cookies = fetch CGI::Cookie;
    my ($user, $uid, $session);# = "Public";
    if (ref $cookies{$cookie_name})
      {
	my %session = $cookies{$cookie_name}->value;
	$session = $session{session};
	my ($user_session) = $coge->resultset("UserSession")->find({session=>$session});
	$user = $user_session->user if $user_session;
      }
    unless ($user)
      {
	$user = new CoGeX::User;
	$user->user_name("public");
      }
    return ($user);
  }

sub gen_cookie
  {
    my $self = shift;
    my %opts = @_;
    my $session = $opts{session} || 0;
    my $exp = $opts{exp} || "+12M";
    my %params = (-name=>$cookie_name,
		 -path=>"/CoGe");

    $params{-expires} = $exp if $exp;
    $params{-values}={session=>$session} if $session;
    my $c = new CGI::Cookie(%params);
    return $c;
  }



1;
