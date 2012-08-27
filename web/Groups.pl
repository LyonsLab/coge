#! /usr/bin/perl -w

use strict;
use CGI;

#use CGI::Ajax;
use JSON::XS;
use CoGe::Accessory::LogUser;
use CoGe::Accessory::Web;
use CoGeX;
use HTML::Template;
use Digest::MD5 qw(md5_base64);
use URI::Escape;
use Data::Dumper;
use File::Path;

no warnings 'redefine';

use vars
  qw($P $DBNAME $DBHOST $DBPORT $DBUSER $DBPASS $connstr $PAGE_NAME $TEMPDIR $USER $DATE $BASEFILE $coge $cogeweb %FUNCTION $COOKIE_NAME $FORM $URL $COGEDIR $TEMPDIR $TEMPURL);
$P = CoGe::Accessory::Web::get_defaults( $ENV{HOME} . 'coge.conf' );

$DATE = sprintf(
	"%04d-%02d-%02d %02d:%02d:%02d",
	sub { ( $_[5] + 1900, $_[4] + 1, $_[3] ), $_[2], $_[1], $_[0] }
	  ->(localtime)
);

$FORM = new CGI;

$DBNAME  = $P->{DBNAME};
$DBHOST  = $P->{DBHOST};
$DBPORT  = $P->{DBPORT};
$DBUSER  = $P->{DBUSER};
$DBPASS  = $P->{DBPASS};
$connstr =
  "dbi:mysql:dbname=" . $DBNAME . ";host=" . $DBHOST . ";port=" . $DBPORT;
$coge = CoGeX->connect( $connstr, $DBUSER, $DBPASS );

$COOKIE_NAME = $P->{COOKIE_NAME};
$URL         = $P->{URL};
$COGEDIR     = $P->{COGEDIR};
$TEMPDIR     = $P->{TEMPDIR} . "Groups/";
mkpath( $TEMPDIR, 0, 0777 ) unless -d $TEMPDIR;
$TEMPURL = $P->{TEMPURL} . "Groups/";

my ($cas_ticket) = $FORM->param('ticket');
$USER = undef;
($USER) = CoGe::Accessory::Web->login_cas(
	cookie_name => $COOKIE_NAME,
	ticket      => $cas_ticket,
	coge        => $coge,
	this_url    => $FORM->url()
  )
  if ($cas_ticket);
($USER) = CoGe::Accessory::LogUser->get_user(
	cookie_name => $COOKIE_NAME,
	coge        => $coge
  )
  unless $USER;

%FUNCTION = (
	gen_html                 => \&gen_html,
	get_groups_for_user      => \&get_groups_for_user,
	create_group             => \&create_group,
	delete_group             => \&delete_group,
	update_group             => \&update_group,
	modify_users             => \&modify_users,
	edit_group_info          => \&edit_group_info,
	add_user_to_group        => \&add_user_to_group,
	remove_user_from_group   => \&remove_user_from_group,
	add_list_to_group        => \&add_list_to_group,
	remove_list_from_group   => \&remove_list_from_group,
	add_genome_to_group      => \&add_genome_to_group,
	remove_genome_from_group => \&remove_genome_from_group,
);

dispatch();

sub dispatch {
	my %args  = $FORM->Vars;
	my $fname = $args{'fname'};
	if ($fname) {

		#my %args = $FORM->Vars;
		#print STDERR Dumper \%args;
		if ( $args{args} ) {
			my @args_list = split( /,/, $args{args} );
			print $FORM->header, $FUNCTION{$fname}->(@args_list);
		}
		else {
			print $FORM->header, $FUNCTION{$fname}->(%args);
		}
	}
	else {
		print $FORM->header, gen_html();
	}
}

sub gen_html {
	my $html;
	my $template =
	  HTML::Template->new( filename => $P->{TMPLDIR} . 'generic_page.tmpl' );
	$template->param( HELP => '/wiki/index.php?title=Groups' );
	my $name = $USER->user_name;
	$name = $USER->first_name if $USER->first_name;
	$name .= " " . $USER->last_name if $USER->first_name && $USER->last_name;
	$template->param( USER       => $name );
	$template->param( TITLE      => qq{Manage User Groups} );
	$template->param( PAGE_TITLE => qq{Manage Groups} );
	$template->param( LOGO_PNG   => "Groups-logo.png" );
	$template->param( LOGON      => 1 ) unless $USER->user_name eq "public";
	$template->param( DATE       => $DATE );
	$template->param( BODY       => gen_body() );
	$template->param( ADJUST_BOX => 1 );
	$name .= $name =~ /s$/ ? "'" : "'s";
	$template->param( BOX_NAME => $name . " groups" );
	$html .= $template->output;
}

sub gen_body {
	my $template =
	  HTML::Template->new( filename => $P->{TMPLDIR} . 'Groups.tmpl' );
	$template->param( PAGE_NAME => $FORM->url );
	$template->param( MAIN      => 1 );
	my $groups = get_groups_for_user();
	$template->param( MAIN_TABLE => $groups );
	my $roles = get_roles();
	$template->param( ROLE_LOOP => $roles );
	$template->param( ADMIN_AREA => 1 ) if $USER->is_admin;
	my $ugid;
	$ugid = $FORM->param('ugid') if defined $FORM->param('ugid');
	my $box_open = $ugid ? 'true' : 'false';
	$template->param( EDIT_BOX_OPEN    => $box_open );
	$template->param( ADDITIONAL_STUFF => qq{edit_group({ugid: $ugid});} )
	  if $ugid;
	return $template->output;
}

sub get_roles {
	my @roles;
	foreach my $role ( $coge->resultset('Role')->all() ) {
		next if $role->name =~ /admin/i && !$USER->is_admin;
		my $name = $role->name;
		$name .= ": " . $role->description if $role->description;
		push @roles, { RID => $role->id, NAME => $name };
	}
	return \@roles;
}

sub add_user_to_group {
	my %opts = @_;
	my $ugid = $opts{ugid};
	my $uid  = $opts{uid};

	#    return 1 if $uid == $USER->id;
	return "UGID and/or UID not specified" unless $ugid && $uid;
	my $group = $coge->resultset('UserGroup')->find($ugid);
	if ( $group->locked && !$USER->is_admin ) {
		return "This is a locked group.  Admin permission is needed to modify.";
	}
	my @res =
	  $coge->resultset('UserGroupConnector')
	  ->search( { user_id => $uid, user_group_id => $ugid } );
	return 1 if @res;    #previously added;
	my ($ugc) =
	  $coge->resultset('UserGroupConnector')
	  ->create( { user_id => $uid, user_group_id => $ugid } );
	return 1;
}

sub remove_user_from_group {
	my %opts = @_;
	my $ugid = $opts{ugid};
	my $uid  = $opts{uid};
	if ( $uid == $USER->id && !$USER->is_admin )    #only allow this for admins
	{
		return "Can't remove yourself from this group!";
	}
	my $group = $coge->resultset('UserGroup')->find($ugid);
	if ( $group->locked && !$USER->is_admin ) {
		return "This is a locked group.  Admin permission is needed to modify.";
	}

	foreach my $ugc ( $coge->resultset('UserGroupConnector')
		->search( { user_id => $uid, user_group_id => $ugid } ) )
	{
		$ugc->delete;
	}
	return 1;
}

sub add_genome_to_group {
	my %opts  = @_;
	my $ugid  = $opts{ugid};
	my $dsgid = $opts{dsgid};
	return "DSGID and/or UGID not specified" unless $dsgid && $ugid;
	my $is_owner = $USER->is_owner( dsg => $dsgid );
	$is_owner = 1 if $USER->is_admin;
	return "You are not the owner of this genome" unless $is_owner;
	my ($ugdc) =
	  $coge->resultset('UserGroupDataConnector')
	  ->create( { genome_id => $dsgid, user_group_id => $ugid } );
	return 1;
}

sub remove_genome_from_group {
	my %opts  = @_;
	my $ugid  = $opts{ugid};
	my $dsgid = $opts{dsgid};
	return "DSGID and/or UGID not specified" unless $dsgid && $ugid;
	my ($ugdc) =
	  $coge->resultset('UserGroupDataConnector')
	  ->search( { genome_id => $dsgid, user_group_id => $ugid } );
	$ugdc->delete;
	return 1;
}

sub modify_users {
	my %opts = @_;
	my $ugid = $opts{ugid};
	return 0 unless $ugid;
	my %data;
	my $group = $coge->resultset('UserGroup')->find($ugid);
	$data{title} = $group->name . "(id" . $group->id . ")";
	$data{title} .= ": " . $group->description if $group->description;
	my %users;
	my @users;

	foreach my $user (
		sort {
			     $a->last_name cmp $b->last_name
			  || $a->user_name cmp $b->user_name
		} $group->users
	  )
	{

		#          next if $user->user_name eq $USER->user_name; #skip self;
		my $name = $user->user_name;
		$name .= " : " . $user->first_name if $user->first_name;
		$name .= " " . $user->last_name
		  if $user->first_name && $user->last_name;
		push @users, { uid_name => $name, uid => $user->id };
		$users{ $user->id } = 1;
	}
	$data{users} = \@users;
	my @all_users;
	foreach my $user (
		sort {
			     $a->last_name cmp $b->last_name
			  || $a->user_name cmp $b->user_name
		} $coge->resultset('User')->all
	  )
	{
		next if $users{ $user->id };    #skip users we already have
		my $name = $user->user_name;
		$name .= " : " . $user->first_name if $user->first_name;
		$name .= " " . $user->last_name
		  if $user->first_name && $user->last_name;
		push @all_users, { uid_name => $name, uid => $user->id };
	}
	$data{all_users} = \@all_users;

#   my @genomes;
#   my %genomes;
#   foreach my $genome ($group->private_genomes)
#     {
#	my $name = $genome->organism->name;
#	$name .= " ".join (", ", map{$_->name} $genome->source) .": ";
#       $name .= $genome->name.", " if $genome->name;
#        $name .=  "v".$genome->version." ".$genome->type->name;#." ".commify($genome->length)."nt";
#	push @genomes, {dsgid_name=>$name, dsgid=>$genome->id};
#	$genomes{$genome->id}=1;
#     }
#   $data{genomes}=\@genomes;

#   my @genome_list;
#time to deal with admin privileges
#   if ($USER->is_admin)
#     {
#	@genome_list = $coge->resultset('Genome')->search({restricted=>1});
#     }
#   else
#     {
#	@genome_list = $USER->private_genomes;
#      }
#    my @all_genomes;
#    foreach my $genome (sort {$a->organism->name cmp $b->organism->name} @genome_list)
#      {
#	next if $genomes{$genome->id}; #skip users we already have
#	my $name = $genome->organism->name;
#	$name .= " ".join (", ", map{$_->name} $genome->source) .": ";
#	$name .= $genome->name.", " if $genome->name;
#	$name .=  "v".$genome->version." ".$genome->type->name;#." ".commify($genome->length)."nt";
#	push @all_genomes, {dsgid_name=>$name, dsgid=>$genome->id};
#      }
#    $data{all_genomes}=\@all_genomes;
	my $template =
	  HTML::Template->new( filename => $P->{TMPLDIR} . 'Groups.tmpl' );
	$template->param( MODIFY_USERS => 1 );
	$template->param( UGID         => $ugid );

	#    $template->param(NAME=>$data{name});
	#    $template->param(DESC=>$data{desc});
	$template->param( UGID_LOOP     => $data{users} );
	$template->param( ALL_UGID_LOOP => $data{all_users} );

	#    $template->param(DSGID_LOOP=>$data{genomes});
	#    $template->param(ALL_DSGID_LOOP=>$data{all_genomes});
	#add message if user is an admin
	#    if ($USER->is_admin)
	#      {
	#	$template->param(ADMIN=>"<span class=alert>Admin Mojo Invoked!</span>");
	#      }
	$data{output} = $template->output;
	return encode_json( \%data );
}

sub edit_group_info {
	my %opts = @_;
	my $ugid = $opts{ugid};
	return 0 unless $ugid;
	my %data;
	my $group = $coge->resultset('UserGroup')->find($ugid);
	$data{title} = $group->name . "(id" . $group->id . ")";
	$data{title} .= ": " . $group->description if $group->description;
	$data{name} = $group->name;
	$data{desc} = "";
	$data{desc} = $group->description if $group->description;
	my $template =
	  HTML::Template->new( filename => $P->{TMPLDIR} . 'Groups.tmpl' );
	$template->param( EDIT_GROUP_INFO => 1 );
	$template->param( UGID            => $ugid );
	$template->param( NAME            => $data{name} );
	$template->param( DESC            => $data{desc} );
	$data{output} = $template->output;
	return encode_json( \%data );
}

sub create_group {
	my %opts = @_;
	return "You need to be a registered user to create a user group!"
	  unless $USER->id;
	my $name = $opts{name};
	return "No specified name!" unless $name;
	my $desc = $opts{desc};
	my $rid  = $opts{rid};
	my ($role) = $coge->resultset('Role')->find( { name => "Reader" } );
	$rid = $role->id unless $rid;
	my $group =
	  $coge->resultset('UserGroup')
	  ->create( { name => $name, description => $desc, role_id => $rid } );
	my $ugc =
	  $coge->resultset('UserGroupConnector')
	  ->create( { user_id => $USER->id, user_group_id => $group->id } );
	return 1;
}

sub update_group {
	my %opts = @_;
	my $ugid = $opts{ugid};
	return 0 unless $ugid;
	my $name = $opts{name};
	return 0 unless $name;
	my $desc  = $opts{desc};
	my $group = $coge->resultset('UserGroup')->find($ugid);
	$group->name($name);
	$group->description($desc) if $desc;
	$group->update;
	return 1;
}

sub delete_group {
	my %opts  = @_;
	my $ugid  = $opts{ugid};
	my $group = $coge->resultset('UserGroup')->find($ugid);
	return unless $group;
	if ( $group->locked && !$USER->is_admin ) {
		return "This is a locked group.  Admin permission is needed to delete.";
	}
	$group->delete();
}

sub get_groups_for_user
{
	my %opts = @_;
	
	my @group_list;
	if ( $USER->is_admin ) {
		@group_list = $coge->resultset('UserGroup')->all();
	}
	else {
		@group_list = $USER->groups;
	}
	
	my @groups;
	foreach my $group (@group_list) {
		my %groups;
		$groups{NAME} = $group->name . " (id" . $group->id . ")";
		$groups{DESC} = $group->description if $group->description;
		my $role = $group->role->name;

#		$role .= ": ".$group->role->description if $group->role->description;
		$groups{ROLE} = $role;

#		my $perm = join (", ", map {$_->name} $group->role->permissions);
#		$groups{PERM}=$perm;
		my @users;
		foreach my $user ( sort { $a->last_name cmp $b->last_name } $group->users ) {
#	    	next if $user->user_name eq $USER->user_name; #skip self;
			my $name = $user->user_name;
			$name = $user->first_name if $user->first_name;
			$name .= " " . $user->last_name if $user->first_name && $user->last_name;
			push @users, $name;
		}

#		push @users, "Self only" unless @users;
		$groups{MEM} = join( ",<br>", @users );
		my %lists;
		foreach my $list (sort { $a->type->name cmp $b->type->name || $a->name cmp $b->name } $group->lists ) {
			my $name;
			$name .= "(P)" if !$list->restricted;    #list is open to the public (P)
			$name .= $list->name;
			$name = qq{<span class="link" onclick="window.open('ListView.pl?lid=} . $list->id . qq{');">} . $name . "</span>";
			push @{ $lists{ $list->type->name } }, $name;
		}
		my $lists = "<table class='small'>";
		foreach my $type ( sort keys %lists ) {
			$lists .= qq{<tr><td>$type<td>} . join( "<br>", @{ $lists{$type} } );
		}
		$lists .= "</table>";
		$groups{LISTS} = join( ",<br>", $lists );

###	my @genome;
#	push @genome, "Apotheosis" if $group->role->name =~ /admin/i;
#	foreach my $genome (sort {$a->organism->name cmp $b->organism->name || $a->version cmp $b->version || $a->genomic_sequence_type_id  <=> $b->genomic_sequence_type_id} $group->genomes)
#	  {
#	    my $name = $genome->organism->name;
#	    $name .= " (v".$genome->version.", ".$genome->type->name.")";
#	    $name = qq{<span class="link" onclick="window.open('OrganismView.pl?dsgid=}.$genome->id.qq{');">}.$name."</span>";
#	    push @genome, $name;
#	  }
#	$groups{GENOMES}=join (",<br>", @genome);

		$groups{UGID} = $group->id;
		push @groups, \%groups;
	}
	my $template = HTML::Template->new( filename => $P->{TMPLDIR} . 'Groups.tmpl' );
	$template->param( GROUP_TABLE => 1 );
	$template->param( GROUPS_LOOP => \@groups );
	return $template->output;
}

