#! /usr/bin/perl -w

use strict;
use CGI;
use CoGeX;
use DBI;

use Data::Dumper;
use CoGe::Accessory::LogUser;
use CoGe::Accessory::Web;
use HTML::Template;
use JSON::XS;
use URI::Escape;
use Spreadsheet::WriteExcel;
use Digest::MD5 qw(md5_base64);
use DBIxProfiler;
use File::Path;
use Sort::Versions;
no warnings 'redefine';

use vars qw(
	$P $DBNAME $DBHOST $DBPORT $DBUSER $DBPASS $connstr $PAGE_TITLE 
	$TEMPDIR $USER $DATE $COGEDIR $coge $FORM $URL $TEMPURL $COOKIE_NAME 
	%FUNCTION
);

$P         = CoGe::Accessory::Web::get_defaults( $ENV{HOME} . 'coge.conf' );
$ENV{PATH} = $P->{COGEDIR};
$COGEDIR   = $P->{COGEDIR};
$URL       = $P->{URL};
$DATE = sprintf(
	"%04d-%02d-%02d %02d:%02d:%02d",
	sub { ( $_[5] + 1900, $_[4] + 1, $_[3] ), $_[2], $_[1], $_[0] }->(localtime)
);

$PAGE_TITLE = 'Genomes';

$TEMPDIR = $P->{TEMPDIR} . $PAGE_TITLE . '/';
mkpath( $TEMPDIR, 0, 0777 ) unless -d $TEMPDIR;
$TEMPURL = $P->{TEMPURL} . $PAGE_TITLE . '/';

$FORM = new CGI;

$DBNAME      = $P->{DBNAME};
$DBHOST      = $P->{DBHOST};
$DBPORT      = $P->{DBPORT};
$DBUSER      = $P->{DBUSER};
$DBPASS      = $P->{DBPASS};
$connstr     = "dbi:mysql:dbname=" . $DBNAME . ";host=" . $DBHOST . ";port=" . $DBPORT;
$coge        = CoGeX->connect( $connstr, $DBUSER, $DBPASS );
$COOKIE_NAME = $P->{COOKIE_NAME};

my ($cas_ticket) = $FORM->param('ticket');
$USER = undef;
($USER) = CoGe::Accessory::Web->login_cas( ticket => $cas_ticket, coge => $coge, this_url => $FORM->url() ) if ($cas_ticket);
($USER) = CoGe::Accessory::LogUser->get_user( cookie_name => $COOKIE_NAME, coge => $coge ) unless $USER;

#$SIG{'__WARN__'} = sub { };    # silence warnings

%FUNCTION = (
	generate_html			=> \&generate_html,
	create_genome			=> \&create_genome,
	# delete_genome			=> \&delete_genome,
	get_genomes_for_user	=> \&get_genomes_for_user,
);

if ( $FORM->param('jquery_ajax') ) {
	my %args  = $FORM->Vars;
	my $fname = $args{'fname'};
	if ($fname) {
		if ( $args{args} ) {
			my @args_list = split( /,/, $args{args} );
			print $FORM->header, $FUNCTION{$fname}->(@args_list);
		}
		else {
			print $FORM->header, $FUNCTION{$fname}->(%args);
		}
	}
}
else {
	print $FORM->header, "\n", generate_html();
}

# sub delete_genome {
# 	my %opts = @_;
# 	my $gid = $opts{gid};
# 	return "Must have valid genome id\n" unless ($gid);
	
# 	# Security check
# 	if (not $USER->is_owner(dsg => $gid)) {
# 		return 0;
# 	}

# 	# Delete the genome and associated connectors & datasets
# 	#FIXME add some error checking/logging here
# 	my $genome = $coge->resultset('Genome')->find($gid);
# 	return 0 unless ($genome);
# 	$genome->delete; #FIXME doesn't delete list connectors
	
# 	CoGe::Accessory::Web::log_history( db => $coge, user_id => $USER->id, page => "$PAGE_TITLE.pl", description => 'delete genome id' . $genome->id );

# 	return 1;
# }

sub get_genomes_for_user {
	#my %opts = @_;

	my @genomes;
#	if ( $USER->is_admin ) {
#		@genomes = $coge->resultset('Genome')->all();
#	}
#	else {
		push @genomes, $USER->owner_list->genomes;
#	}

	my @genome_info;
	foreach my $g (@genomes) {#( sort genomecmp @genomes ) {
		push @genome_info, 
		  { NAME  => qq{<span class="link" onclick='window.open("GenomeInfo.pl?gid=} . $g->id . qq{")'>} . $g->info . "</span>",
			VERSION  => $g->version,
			DATE =>  $g->date,
			EDIT_BUTTON => "<span class='link ui-icon ui-icon-gear' onclick=\"window.open('GenomeInfo.pl?gid=" . $g->id . "')\"></span>",
			DELETE_BUTTON => "<span class='link ui-icon ui-icon-trash' onclick=\"dialog_delete_genome({gid: '" . $g->id . "'});\"></span>"
		  };
	}

	my $template = HTML::Template->new( filename => $P->{TMPLDIR} . $PAGE_TITLE . '.tmpl' );
	$template->param( DO_GENOME_TABLE => 1 );
	$template->param( GENOME_LOOP  => \@genome_info );
#	$template->param( TYPE_LOOP => get_list_types() );

	return $template->output;
}

sub generate_html {
	my $html;
	my $template = HTML::Template->new( filename => $P->{TMPLDIR} . 'generic_page.tmpl' );
	$template->param( PAGE_TITLE => $PAGE_TITLE );
	$template->param( HELP       => '/wiki/index.php?title=' . $PAGE_TITLE . '.pl' );
	my $name = $USER->user_name;
	$name = $USER->first_name if $USER->first_name;
	$name .= ' ' . $USER->last_name if ( $USER->first_name && $USER->last_name );
	$template->param( USER     => $name );
	$template->param( LOGO_PNG => $PAGE_TITLE . "-logo.png" );
	$template->param( LOGON    => 1 ) unless $USER->user_name eq "public";
	$template->param( DATE     => $DATE );
	my $link = "http://" . $ENV{SERVER_NAME} . $ENV{REQUEST_URI};
	$link = CoGe::Accessory::Web::get_tiny_link( url => $link );

	$template->param( BODY => generate_body() );
	$template->param( ADJUST_BOX => 1 );

	$html .= $template->output;
}

sub generate_body {
	my $template = HTML::Template->new( filename => $P->{TMPLDIR} . $PAGE_TITLE . '.tmpl' );
	$template->param( MAIN => 1 );
	$template->param( PAGE_NAME => $PAGE_TITLE );
	$template->param( GENOME_TABLE  => get_genomes_for_user() );
	
	return $template->output;
}

# FIXME these comparison routines are duplicated elsewhere
sub genomecmp {
	no warnings 'uninitialized'; # disable warnings for undef values in sort
	$a->organism->name cmp $b->organism->name || versioncmp($b->version, $a->version) || $a->type->id <=> $b->type->id || $a->name cmp $b->name || $b->id cmp $a->id
}