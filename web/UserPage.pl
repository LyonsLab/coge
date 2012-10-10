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
no warnings 'redefine';

use vars qw($P $DBNAME $DBHOST $DBPORT $DBUSER $DBPASS $connstr $PAGE_NAME $TEMPDIR $USER $DATE $BASEFILE $coge $cogeweb $cgi %FUNCTION $COOKIE_NAME $FORM);
$P = CoGe::Accessory::Web::get_defaults($ENV{HOME}.'coge.conf');

$DATE = sprintf( "%04d-%02d-%02d %02d:%02d:%02d",
		sub { ($_[5]+1900, $_[4]+1, $_[3]),$_[2],$_[1],$_[0] }->(localtime));
$FORM = new CGI;

$DBNAME = $P->{DBNAME};
$DBHOST = $P->{DBHOST};
$DBPORT = $P->{DBPORT};
$DBUSER = $P->{DBUSER};
$DBPASS = $P->{DBPASS};
$connstr = "dbi:mysql:dbname=".$DBNAME.";host=".$DBHOST.";port=".$DBPORT;
$coge = CoGeX->connect($connstr, $DBUSER, $DBPASS );

$COOKIE_NAME = $P->{COOKIE_NAME};

my ($cas_ticket) =$FORM->param('ticket');
$USER = undef;
($USER) = CoGe::Accessory::Web->login_cas(cookie_name=>$COOKIE_NAME, ticket=>$cas_ticket, coge=>$coge, this_url=>$FORM->url()) if($cas_ticket);
($USER) = CoGe::Accessory::LogUser->get_user(cookie_name=>$COOKIE_NAME,coge=>$coge) unless $USER;

#my $pj = new CGI::Ajax(
%FUNCTION = (
    gen_html=>\&gen_html,
    profile=>\&profile,
    generate_worklist_header=>\&generate_worklist_header,
    view_all_work=>\&view_all_work,
    view_work_entry=>\&view_work_entry,
    remove_work_entry=>\&remove_work_entry,
    archive_work_entry=>\&archive_work_entry,
    generate_workflow=>\&generate_workflow,
    create_workflow=>\&create_workflow,
    remove_workflow=>\&remove_workflow,
    update_work_entry=>\&update_work_entry,
    show_workflow=>\&show_workflow,
    add_work_entry_to_workflow=>\&add_work_entry_to_workflow,
    remove_work_entry_from_workflow=>\&remove_work_entry_from_workflow,
    );
#print $pj->build_html($FORM, \&gen_html);
#print $FORM->header, gen_html();
#print STDERR Dumper \%ENV;
dispatch();

sub gen_html
  {
    my $html;    
    my $template = HTML::Template->new(filename=>$P->{TMPLDIR}.'generic_page.tmpl');
    $template->param(HELP=>'/wiki/index.php?title=UserProfile');
    # print STDERR "user is: ",$USER,"\n";
    my $name = $USER->user_name;
    $name = $USER->first_name if $USER->first_name;
    $name .= " ".$USER->last_name if $USER->first_name && $USER->last_name;
    $template->param(USER=>$name);
    $name .= $name =~ /s$/ ? "'" : "'s";
    $template->param(TITLE=>$name.qq{ User Settings});
    $template->param(PAGE_TITLE=>qq{User Settings});
    # $template->param(LOGO_PNG=>"FeatList-logo.png");
    $template->param(LOGON=>1) unless $USER->user_name eq "public";
    $template->param(DATE=>$DATE);
    $template->param(BODY=>gen_body());
    $template->param(ADJUST_BOX=>1);
    $html .= $template->output;
  }

sub dispatch
{
    my %args = $FORM->Vars;
    my $fname = $args{'fname'};
    if($fname)
    {
	#my %args = $FORM->Vars;
	#print STDERR Dumper \%args;
	if($args{args}){
	    my @args_list = split( /,/, $args{args} );
	    print $FORM->header, $FUNCTION{$fname}->(@args_list);
       	}
	else{
	    print $FORM->header, $FUNCTION{$fname}->(%args);
	}
    }
    else{
	print $FORM->header, gen_html();
    }
}

sub gen_body
  {
      my $template = HTML::Template->new(filename=>$P->{TMPLDIR}.'UserPage.tmpl');
      $template->param(STARTUP=>1);
      $template->param(MAIN=>1);
      return $template->output;
  }

sub profile
{
    my $name;
    $name = $USER->first_name if $USER->first_name;
    $name .= " ".$USER->last_name if $USER->first_name && $USER->last_name;
    $name = "No Name Specified" unless $name;
    my $template = HTML::Template->new(filename=>$P->{TMPLDIR}.'UserPage.tmpl');
    $template->param(PROFILE=>1);
    $template->param(USER_NAME=>$USER->user_name);
    $template->param(EMAIL=>$USER->email);
    $template->param(REAL_NAME=>$name);
    $template->param(DESC=>$USER->description);
    my @groups;
    foreach my $group ($USER->groups)
      {
	my %groups;
	$groups{GROUP_NAME}="<span class=link onclick=windown.open('GroupView.pl?ugid=".$group->id."');>".$group->name."</span>";
	$groups{GROUP_DESC}=$group->description if $group->description;
	my $role = $group->role->name;
	$role .= ": ".$group->role->description if $group->role->description;
	$groups{ROLE}=$role;
	my $perm = join (", ", map {$_->name} $group->role->permissions);
	$groups{GROUP_PERMISSION}=$perm;
	push @groups, \%groups;
      }
    $template->param(GROUPS_LOOP=>\@groups);
    #print STDERR $template->output;
    return $template->output;
}

sub generate_worklist_header
{
    my $template = HTML::Template->new(filename=>$P->{TMPLDIR}.'UserPage.tmpl');
    $template->param(WORK_TYPE_CHOICE=>1);
    my $content = view_all_work();
    my $html = $template->output; 
    $html .= '<div id=worklist>'.$content.'</div>';
    return $html;
}

sub view_all_work
{
    my %opts = @_;
   # print STDERR Dumper \%opts;
    my $type = $opts{type};
    my $message = $opts{message};
    my @works;
    my $count=0;
    my $template = HTML::Template->new(filename=>$P->{TMPLDIR}.'UserPage.tmpl');
    $template->param(WORK_LIST=>1);
    if($message)
    {
	$template->param(MESSAGE=>1);
	if($message =~ /warning/i) {$template->param(MESSAGE_CLASS=>"ui-state-alert ui-corner-all");}
	else{ $template->param(MESSAGE_CLASS=>"ui-state-highlight ui-corner-all");}
	$template->param(MESSAGE_CONTENT=>$message);
    }
    foreach my $work ($USER->works)
    {
	if ($type && $type =~ /tmp/) {next unless $work->archive != 1;}
	if ($type && $type =~ /archive/) {next unless $work->archive == 1;}
	my ($date) = $work->date =~ /^(\d{4}-\d{2}-\d{2}).*/;
	my $archive = $work->archive ? "Yes" : "No";
	my ($tool) = $work->page =~ /(\w+)\.pl/;
	my $name_html = "<a href=javascript:void(0) onClick=get_work_entry(".$work->id.")>".$work->name."</a>";
	push @works, ({WORK_NAME=>$name_html,
		       WORK_DESC=>$work->description,
		       WORK_DATE=>$date,
		     # WORK_ARCH=>$archive,#Disabled until Archive issue resolved
		       WORK_TOOL=>$tool,
		       WORK_ID=>$work->id,
		       WORK_VIEW=>qq{<a href='javascript:void(0)' title='View Work' onClick="get_work_entry('}.$work->id.qq{');" class='ui-button ui-button-icon-solo ui-corner-all'><span class='ui-icon ui-icon-info'></span></a>},
		       WORK_EDIT=>qq{<a href='javascript:void(0)' title='Edit Work' onClick="edit_work_entry(}.$work->id.qq{);" class='ui-button ui-button-icon-solo ui-corner-all'><span class='ui-icon ui-icon-pencil'></span></a>},
		       WORK_DELETE=>qq{<a href='javascript:void(0)' title='Delete Work' onClick="delete_work_entry(}.$work->id.qq{)" class='ui-button ui-button-icon-solo ui-corner-all'><span class='ui-icon ui-icon-trash'></span></a>},
		      });
	$count++;
    }
    $template->param(WORKS_LOOP=>\@works);
    my $html = $count ? $template->output : "No works exist";
    return $html;
}

sub view_work_entry
{
    my %opts = @_;
    my $work_id = $opts{'work_id'};
    my $edit = $opts{'edit'};
    my $template = HTML::Template->new(filename=>$P->{TMPLDIR}.'UserPage.tmpl');
    if($edit) {$template->param(WORK_EDIT=>1); }
    else {$template->param(WORK_VIEW=>1);}
    my ($work) = $coge->resultset("Work")->search({work_id=>$work_id,
						   user_id=>$USER->id,
						  });
   # $template->param(WORK_ID=>$work->id);
    my ($app) = $work->page =~ /(\w+)\.pl/;
    $template->param(WORK_NAME=>$work->name);
    $template->param(WORK_ID=>$work->id);
    $template->param(WORK_DESC=>$work->description);
    $template->param(WORK_PARAM=>$work->parameter.";ln=".$work->name);
    $template->param(WORK_APP=>$app);

    my $note = $work->note;

    $note =~ s/<br\/>/\n/g if $edit; #replace new lines with html breaks so text is formatted properly if editing
    
    $template->param(WORK_NOTE=>$note);
    my $workflow_list;
    my %workflows;
    foreach my $work_order ($work->work_orders)
    {
	my $wf_id = $work_order->workflow_id;
	my ($wf) = $coge->resultset("Workflow")->search({workflow_id=>$wf_id,
						     user_id=>$USER->id,
						    });
	next if $workflows{$wf->name};
	$workflow_list .= qq{<a onclick="\$('#work_dialog').dialog('close');get_workflow_header(}.$wf->id.qq{)" href="javascript:void(0)">}.$wf->name."</a><br/>";
	$workflows{$wf->name}++; #checks to make sure you don't list the same workflow more than once
    }
    $template->param(WORKFLOW_LIST=>$workflow_list || "None");
    my ($date) = $work->date =~ /^(\d{4}-\d{2}-\d{2}).*/;
    $template->param(WORK_DATE=>$date);
    
    unless ($edit) {
	my @workflows;
	foreach my $workflow ($USER->workflows)
	{
	    my $name = $workflow->name;
	    my $id = $workflow->id;
	    push @workflows, {WF_ID=>$id,WF_NAME=>$name};
	}
	if(@workflows)
	{
	    my $wf_names;
	    $template->param(WF_EXIST=>1);
	    $template->param(WF_SELECT_LOOP=>\@workflows);
	   # map { $wf_names .= $_->{WF_NAME}."<br/>" } @workflows;
	   # $template->param(WORKFLOW_LIST=>$wf_names);
	}
	else{
	  #  $template->param(WORKFLOW_LIST=>{WF_NAME=>"None"});
	}
    }
    
    my $json_ref = {html=>$template->output,work_name=>$work->name,work_id=>$work->id};
    my $json = encode_json $json_ref;
    return $json;
}

sub archive_work_entry
{
    my %opts = @_;
    my $work_ids = $opts{work_ids};
    chomp $work_ids;
    my $archive = $opts{archive};
    my $message;
    foreach my $work_id (split(/-/,$work_ids))
    {
	next unless $work_id;
	my ($work) = $coge->resultset("Work")->search({work_id=>$work_id,
						   user_id=>$USER->id,
						  });
	next unless $work;
	my $name = $work->name;
	if($archive) { $work->archive(1); }
	else { $work->archive(0); }
	$work->update;
	$message .= $archive ? "$name is now archived<br/>" : "$name is no longer archived<br/>";
    }

    return $message;
}


sub remove_work_entry
{
    my $work_ids = shift;
    chomp $work_ids;
   # print STDERR $work_ids,"\n";
    my $message;
    foreach my $work_id ( split(/-/, $work_ids) )
    {
#	print STDERR $work_id,"\n";
	next unless $work_id =~ /\d/;
#	print STDERR $work_id,"\n";
	my ($work) = $coge->resultset('Work')->search({user_id=>$USER->id,
						       work_id=>$work_id,
						      });
	if($work)
	{
	    my $name = $work->name;
	    $work->delete;
	    $message .= "Deletion of work $name successful!<br/>";
	}
	else{
	    print STDERR "Suspicious activity on account id $USER->id when trying to delete work entry id $work_id","\n";
	    $message .= "Warning: Deletion of work id $work_id unsuccessful.</br>";
	}
    }
    return $message;
}

# NOTE!!!! BIG SECURITY RISK!!! No check for malicious code here...someone
# can plant malicious code in the note, the entry is interpetted as HTML!
# If someone puts javascript in there, will probably be interpreted and
# run.
sub update_work_entry
{
    my %opts = @_;
   # print STDERR Dumper \%opts;
    my $name = $opts{name};
    my $desc = $opts{desc};
    my $notes = $opts{note};
    
    # SECURITY CHECK FOR JAVASCRIPT
    $notes =~ s/<\s*script.*?script\s*>//ig;

    chomp $notes;
    $notes =~ s/\n/<br\/>/g; #Keeps note formatting so new lines do not get lost

    my $work_id = $opts{work_id};
    my $wf_id = $opts{wf_id};
    my ($work) = $coge->resultset('Work')->search({user_id=>$USER->id,
						   work_id=>$work_id,
						  });
    $work->name($name);
    $work->description($desc);
    $work->note($notes);
    $work->update;
    
    my $message = $work ? "Update successful!" : "Update was unsuccesful, edits were not saved.";
    
    my $json_ref = {html=>$message,wf_id=>$wf_id};
    my $json = encode_json $json_ref;
    
    return $json;
}

sub add_work_entry_to_workflow
{
    my %opts = @_;
    my $work_id = $opts{work_id};
    my $wf_id = $opts{wf_id};
    my ($workflow) = $coge->resultset("Workflow")->search({workflow_id=>$wf_id,
							   user_id=>$USER->id,
							  });
   # print STDERR $workflow->id," is wf_id \n\n";
    my ($output) = $workflow->add_to_work_orders({ 'work_id'=>$work_id,
						   "work_order"=>23,
						 });
    my $message = $output ? "<p class='ui-state-highlight ui-corner-all'>Work $work_id successfully saved to workflow ".$workflow->name."</p>" : "<p class='ui-state-alert ui-corner-all'>Save unsuccessful :( </p>";
    return $message;
}

sub remove_work_entry_from_workflow
{
    my %opts = @_;
    print STDERR Dumper \%opts;
    my $work_id = $opts{work_id};
    my $wf_id = $opts{wf_id};
    my ($work_order) = $coge->resultset("WorkOrder")->search({workflow_id=>$wf_id,
							      work_id=>$work_id,
							     });
    # print STDERR $workflow->id," is wf_id \n\n";
    $work_order->delete if $work_order;
    #my $message = $output ? "<p class='ui-state-highlight ui-corner-all'>Work $work_id successfully saved to workflow</p>".$workflow->name : "<p class='ui-state-alert ui-corner-all'>Save unsuccessful :( </p>";
   # return $message;
    return;
}

sub generate_workflow
{
    my $wf_id=shift;
    my $count = 0;
    my $template = HTML::Template->new(filename=>$P->{TMPLDIR}.'UserPage.tmpl');
    my @workflows;
    foreach my $workflow ($USER->workflows)
    {
	my $name = $workflow->name;
	my $id = $workflow->id;
	push @workflows, ({WF_ID=>$id,WF_NAME=>$name});
    }
    if(@workflows)
    {
	my $wf_table = create_workflows_table();
	$template->param(WF_EXIST=>1);
	$template->param(WF_TABLE=>$wf_table);
	$template->param(WF_SELECT_LOOP=>\@workflows);
	#my $html = show_workflow($wf_id) if $wf_id;
	#$template->param(PRELOAD_WORKFLOW=>$html) if $html;
    }
    $template->param(WF_HEADER=>1);
    my $json_ref = {html=>$template->output,wf_id=>$wf_id};
    my $json = encode_json $json_ref;
    return $json;
}

sub show_workflow
{
    my $wf_id = shift;
    #print STDERR $wf_id,"\n";
    my $html;
    my @works;
    my $template = HTML::Template->new(filename=>$P->{TMPLDIR}.'UserPage.tmpl');
    $template->param(WF_CONTENT=>1);
    my ($workflow) = $coge->resultset("Workflow")->search({workflow_id=>$wf_id,
							   user_id=>$USER->id,
							  });
    $template->param(WF_NAME=>$workflow->name);
    $template->param(WF_DESC=>$workflow->description);
    my ($date) = $workflow->date =~ /^(\d{4}-\d{2}-\d{2}).*/;
    $template->param(WF_DATE=>$date);
    $template->param(WF_ID=>$workflow->id);
    foreach my $work ($workflow->works)
    {
	my ($tool) = $work->page =~ /(\w+)\.pl/;
	push @works, ({WORK_ID=>$work->id,WORK_NAME=>$tool." - ".$work->name,WORK_DESC=>$work->description,WORK_PARAM=>$work->parameter,WORK_NOTE=>$work->note,WORK_DATE=>$work->date,WF_ID=>$workflow->id});
    }
    if (@works)
    {
	$template->param(HAS_WORKS=>1);
	$template->param(WORKS_LOOP=>\@works);
    }
    return $template->output;
}

#sub refresh_work_in_workflow
#{
#    my $work_id = shift;
#    my $template = HTML::Template->new(filename=>$P->{TMPLDIR}.'UserPage.tmpl');
#    my ($tool) = $work->page =~ /(\w+)\.pl/;
#    my ($work) = $coge->resultset('Work')->search({user_id=>$USER->id,
#						   work_id=>$work_id,
#						  });
#    my $work_html = {WORK_ID=>$work->id,WORK_NAME=>$tool." - ".$work->name,WORK_DESC=>$work->description,WORK_PARAM=>$work->parameter,WORK_NOTE=>$work->note,WORK_DATE=>$work->date,WF_ID=>$workflow->id};
#    if ($work_html)
#    {
#	$template->param(HAS_WORKS=>1);
#	$template->param(WORKS_LOOP=>$work_html);
#    }
#    my $html = $template->output =~ /(<table>.+<\/table>)?/;
#    return $html;
#}

sub create_workflows_table
{
    my $html;
    my $i = 1;
    foreach my $workflow ($USER->workflows)
    {
	my $wf_name = $workflow->name;
	my $wf_id = $workflow->id;
	$html .= qq{<tr id="node-$i"><td><span class="folder" onClick='display_workflow($wf_id)' style="cursor:pointer">$wf_name</span></td></tr>};
	my $j = 1;
	foreach my $work ($workflow->works)
	{
	    my $workname = $work->name;
	    my $wid = $work->id;
	    $html .= qq{<tr id="node-$i-$j" class="child-of-node-$i"><td><span class="file" onClick='display_workflow($wf_id,$wid)' style="cursor:pointer">$workname</span></td></tr>};
	    $j++;
	}
	$i++;
    }
    return $html;
}

sub create_workflow
{
    my %opts = @_;
    #print STDERR Dumper \%opts;
    my $html;
    my $template = HTML::Template->new(filename=>$P->{TMPLDIR}.'UserPage.tmpl');
    $template->param(WF_CREATE=>1);
    if ($opts{name} && $opts{desc})
    {
	my $name=$opts{name};
	my $desc = $opts{desc};
	my ($workflow) = $USER->add_to_workflows({ 'name'=>$name,
						   'description'=>$desc,
			       });
	### Should perform check here to ensure $name eq $workflow->name ###

	$template->param(MESSAGE=>1);
	$template->param(SUCCESS=>1);
	$template->param(WF_NAME=>$name);
	$template->param(WF_ID=>$workflow->id);
	$html = $template->output;
    }
    else{
	if($opts{entry})
	{
	    $template->param(MESSAGE=>1);
	    $template->param(WF_NAME=>"<li>Workflow Name</li>") unless $opts{name};
	    $template->param(WF_DESC=>"<li>Workflow Description</li>") unless $opts{desc};
	}
	$template->param(WF_FORM=>1);
    }
    return $template->output;
}

sub remove_workflow
{
    my $wf_id = shift;
    my ($workflow) = $coge->resultset('Workflow')->search({user_id=>$USER->id,
							   workflow_id=>$wf_id,
						  });
    if($workflow)
    {
	my $name = $workflow->name;
	$workflow->delete;
	return "Deletion of workflow $name successful!";
    }
    else{
	print STDERR "Suspicious activity on account id $USER->id when trying to delete workflow entry id $workflow->id","\n";
	return "Warning: Delete unsuccessful";
    }
}

sub edit_workflow
{
    my %opts = @_;
    my $wf_name = $opts{wf_name};
    my $wf_desc = $opts{wf_desc};
    my $wf_id = $opts{wf_id};
    my ($workflow) = $coge->resultset('Workflow')->search({user_id=>$USER->id,
							   workflow_id=>$wf_id,
							  });
    if($workflow)
    {
	$workflow->name($wf_name);
	$workflow->description($wf_desc);
	$workflow->update;
	return "Update of workflow information successful!";
    }
    else{
	print STDERR "Suspicious activity on account id $USER->id when trying to update workflow entry id $workflow->id","\n";
	return "Warning: Update unsuccessful";
    }
    
    
}
