package CoGe::Accessory::LogUser;

use strict;
use CGI::Cookie;
use Data::Dumper;

use vars qw($cookie_name);

$cookie_name = "CoGe";

sub get_user
  {
    my $self = shift;
    my %cookies = fetch CGI::Cookie;
    my ($user, $uid, $session);# = "Public";
    if (ref $cookies{$cookie_name})
      {
	my %session = $cookies{$cookie_name}->value;
	$user = $session{user_name} if $session{user_name};
	$uid = $session{uid} if $session{uid};
	$session = $session{session} if $session{session};
      }
    return ($user, $uid, $session);
  }

sub gen_cookie
  {
    my $self = shift;
    my %opts = @_;
    my $user_name = $opts{user_name};
    my $uid = $opts{uid};
    my $expires = $opts{exp} || "+24h";
    my $session = $opts{session};
    $expires = "+".$expires unless $expires =~ /^\+/;
    
    my $c = new CGI::Cookie(-name=>$cookie_name,
			    -expires=>$expires,
			    -values=>{
				      user_name=>$user_name,
				      uid => $uid,
				      session=>$session,
				     },
			    );
    return $c;
  }



1;
