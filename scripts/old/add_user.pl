#!/usr/bin/perl -w

use strict;
use Data::Dumper;
use CoGeX;

use Getopt::Long;

use vars qw($user_name $first_name $last_name $description $passwd $help $force $delete $email);

GetOptions (
	    "user_name|un=s"=>\$user_name,
	    "first_name|fn=s"=>\$first_name,
	    "last_name|ln=s"=>\$last_name,
	    "description|desc|d=s"=>\$description,
	    "password|passwd|ps|p=s"=>\$passwd,
	    "email|em=s"=>\$email,
	    "help|h"=>\$help,
            "force|f"=>\$force,
	    "delete"=>\$delete,
	   );
$help = 1 unless ($user_name);

help() if $help;

my $connstr = 'dbi:mysql:dbname=DB;host=HOST;port=PORT';
my $coge = CoGeX->connect($connstr, 'USER', 'PASSWORD' );
my ($user) = $coge->resultset('User')->find({user_name=>$user_name});
if (  ref($user)=~/user/i )
  {
    print "User account '$user_name' exists.\n";
    unless ($force || $delete)
      {
       print "Use the option '-force' to overwrite this account.\n";
       exit;
      }
    if ($delete)
      {
	print "Deleting $user_name.\n";
	$user->delete;
	exit;
      }
  }
unless ($passwd)
  {
    print "Need password to generate user account $user_name.\n";
    help();
  }
$first_name = $user->first_name if $user && !$first_name;
unless ($first_name)
  {
    print "Need first name to generate user account $user_name.\n";
    help();
  }
$last_name = $user->last_name if $user && !$last_name;
unless ($last_name)
  {
    print "Need last name to generate user account $user_name.\n";
    help();
  }
$description = $user->description if $user && !$description;
$email = $user->email if $user && !$email;
$user->delete if $user;
my $u = $coge->resultset('User')->create({
					  user_name=>$user_name,
					  first_name=>$first_name,
					  last_name=>$last_name,
					  description=>$description,
					  #					  passwd=>$cpwd,
					  email=>$email,
					 });
my $cpwd = $u->generate_passwd(pwd=>$passwd);
$u->passwd($cpwd);
$u->update;

print "User account ".$u->user_name." created\n";

sub help
  {
    print qq{
Welcome to $0.
  This program generates user accounts for coge.

Options

  -help | -h                     => print this message

  -user_name | -un               => user name in coge (required)

  -first_name | -fn              => user's first name (required)

  -last_name | -ln               => user's last name (required)

  -description | -desc | -d      => user's description

  -password | -passwd | -ps | -p => user account's password (required)

  -email | -em                   => email of user

  -force | -f                    => force update of user account (this will delete the old account)

};
    exit;
  }
