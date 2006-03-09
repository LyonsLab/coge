package CoGe::Accessory::LogUser;

use strict;
use CGI::Cookie;

sub get_user
  {
    my %cookies = fetch CGI::Cookie;
    my $email = "not logged in";
    if (ref $cookies{"Session"})
      {
	my %session = $cookies{"Session"}->value;
	$email = $session{email} if $session{"email"};
      } 
    else
      {
	$email = $ENV{'REMOTE_USER'} if $ENV{'REMOTE_USER'};
      }
    return $email;
  }

1;
