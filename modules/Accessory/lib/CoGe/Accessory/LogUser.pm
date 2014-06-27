package CoGe::Accessory::LogUser;

use strict;
use CGI::Cookie;
use Data::Dumper;
use CoGeX::Result::User;

sub get_user {
	my $self        = shift;
	my %opts        = @_;
	my $cookie_name = $opts{cookie_name};
	my $coge        = $opts{coge};
	my %cookies     = fetch CGI::Cookie;
	my ( $user, $uid, $session );    # = "Public";

	#print STDERR "LogUser::get_user cookie=$cookie_name " . (defined $cookies{$cookie_name} ? 'exists' : '!exists') . "\n";
	if ( $cookie_name && ref $cookies{$cookie_name} ) {
		my %session = $cookies{$cookie_name}->value;
		$session = $session{session};
		my ($user_session) = $coge->resultset("UserSession")->find( { session => $session } );
		$user = $user_session->user if $user_session;
	}
	unless ($user) {
		$user = new CoGeX::Result::User;
		$user->user_name("public");
		$user->admin(0);
	}
	return ($user);
}

sub gen_cookie {
	my $self        = shift;
	my %opts        = @_;
	my $session     = $opts{session} || 0;
	my $exp         = $opts{exp} || '+7d'; # issue 48, this field must have lowercase-only
	my $cookie_name = $opts{cookie_name};
	my %params      = ( -name => $cookie_name, -path => "/" );
#	print STDERR "LogUser::gen_cookie session=$session cookie_name=$cookie_name\n";

	$params{-expires} = $exp if $exp;
	$params{ -values } = { session => $session } if $session;
	my $c = new CGI::Cookie(%params);
	return $c;
}

1;

=head1 NAME

LogUser

=head1 SYNOPSIS

use LogUser

=head1 DESCRIPTION

=head1 USAGE

=head1 BUGS

=head1 SUPPORT

=head1 AUTHOR

Eric Lyons

=head1 COPYRIGHT

This program is free software; you can redistribute
it and/or modify it under the same terms as Perl itself.

The full text of the license can be found in the
LICENSE file included with this module.

=head1 SEE ALSO

=cut
