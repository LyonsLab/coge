#! /usr/bin/perl -w

use strict;
use CGI;
use CoGeX;
use DBI;

use Data::Dumper;
use CoGe::Accessory::LogUser;
use CoGe::Accessory::Web;
use HTML::Template;
use URI::Escape;
use Spreadsheet::WriteExcel;
use Digest::MD5 qw(md5_base64);
use DBIxProfiler;
use File::Path;
no warnings 'redefine';

use vars
	qw( $P $DBNAME $DBHOST $DBPORT $DBUSER $DBPASS $connstr $PAGE_NAME 
		$TEMPDIR $USER $DATE $BASEFILE $COGEDIR $coge $cogeweb $FORM $URL 
		$TEMPURL $COOKIE_NAME );
		
$P = CoGe::Accessory::Web::get_defaults( $ENV{HOME} . 'coge.conf' );
print STDERR $ENV{HOME}, "\n";
$ENV{PATH} = $P->{COGEDIR};
$COGEDIR   = $P->{COGEDIR};
$URL       = $P->{URL};
$DATE      = sprintf(
	"%04d-%02d-%02d %02d:%02d:%02d",
	sub { ( $_[5] + 1900, $_[4] + 1, $_[3] ), $_[2], $_[1], $_[0] }->(localtime)
);
$PAGE_NAME = "ExperimentList.pl";

$TEMPDIR = $P->{TEMPDIR} . "ExperimentList/";
mkpath( $TEMPDIR, 0, 0777 ) unless -d $TEMPDIR;
$TEMPURL = $P->{TEMPURL} . "ExperimentList/";
$FORM    = new CGI;
$DBNAME  = $P->{DBNAME};
$DBHOST  = $P->{DBHOST};
$DBPORT  = $P->{DBPORT};
$DBUSER  = $P->{DBUSER};
$DBPASS  = $P->{DBPASS};
$connstr =
  "dbi:mysql:dbname=" . $DBNAME . ";host=" . $DBHOST . ";port=" . $DBPORT;
$coge        = CoGeX->connect( $connstr, $DBUSER, $DBPASS );
$COOKIE_NAME = $P->{COOKIE_NAME};

my ($cas_ticket) = $FORM->param('ticket');
$USER = undef;
($USER) = CoGe::Accessory::Web->login_cas(
	ticket   => $cas_ticket,
	coge     => $coge,
	this_url => $FORM->url()
  )
  if ($cas_ticket);
($USER) = CoGe::Accessory::LogUser->get_user(
	cookie_name => $COOKIE_NAME,
	coge        => $coge
  )
  unless $USER;

$SIG{'__WARN__'} = sub { };    #silence warnings

my %FUNCTION = (
	gen_html                => \&gen_html,
	get_feature_counts      => \&get_feature_counts,
	generate_excel_file     => \&generate_excel_file,
	generate_csv_file       => \&generate_csv_file,
	gen_data                => \&gen_data,
	send_to_ExperimentList  => \&send_to_ExperimentList,
	save_FeatList_settings  => \&save_FeatList_settings,
	add_to_user_history     => \&add_to_user_history,
);

if ( $FORM->param('jquery_ajax') ) {
	dispatch();
}
else {
	print $FORM->header, "\n", gen_html();
}

sub dispatch {
	my %args  = $FORM->Vars;
	my $fname = $args{'fname'};
	if ($fname) {

		#my %args = $cgi->Vars;
		#print STDERR Dumper \%args;
		if ( $args{args} ) {
			my @args_list = split( /,/, $args{args} );
			print $FORM->header, $FUNCTION{$fname}->(@args_list);
		}
		else {
			print $FORM->header, $FUNCTION{$fname}->(%args);
		}
	}
}

sub gen_html {
	my $html;
	my ($body) = gen_body();
	my $template =
	  HTML::Template->new( filename => $P->{TMPLDIR} . 'generic_page.tmpl' );
	$template->param( PAGE_TITLE => 'ExperimentList' );
	$template->param( HELP       => '/wiki/index.php?title=ExperimentList' );
	my $name = $USER->user_name;
	$name = $USER->first_name if $USER->first_name;
	$name .= " " . $USER->last_name if $USER->first_name && $USER->last_name;
	$template->param( USER     => $name );
	$template->param( LOGO_PNG => "ExperimentList-logo.png" );
	$template->param( LOGON    => 1 ) unless $USER->user_name eq "public";
	$template->param( DATE     => $DATE );
	my $link = "http://" . $ENV{SERVER_NAME} . $ENV{REQUEST_URI};
	$link = CoGe::Accessory::Web::get_tiny_link( url => $link );
	my $box_name = "Experiment List: ";
	my $list_name = $FORM->param('list_name') || $FORM->param('ln');
	$box_name .= " $list_name" if $list_name;
	$box_name .= "<span class=link onclick=window.open('$link');>$link</span>";
	$template->param( BOX_NAME   => $box_name );
	$template->param( BODY       => $body );
	$template->param( ADJUST_BOX => 1 );
	$html .= $template->output;
}

sub gen_body {
	my $template =
	  HTML::Template->new( filename => $P->{TMPLDIR} . 'ExperimentList.tmpl' );
	my $form = $FORM;
	my $no_values;
	$BASEFILE = $form->param('basename');
	my $sort_by_type     = $form->param('sort_type');
	my $sort_by_location = $form->param('sort_loc');
	my $prefs            = CoGe::Accessory::Web::load_settings(
		user => $USER,
		page => $PAGE_NAME,
		coge => $coge
	);
	$prefs = {} unless $prefs;
	$prefs->{display} = {
		NameDescD 	=> 1,
		TypeD 		=> 1,
		AnnotationD => 1,
		SourceD		=> 1,
		VerD		=> 1
	}
	unless $prefs->{display};

	foreach my $key ( keys %{ $prefs->{display} } ) {
		$template->param( $key => "checked" );
	}

	$template->param( SAVE_DISPLAY => 1 ) unless $USER->user_name eq "public";

	#Don't show button to save list data unless valid user
	$template->param( SAVE_DATA => 1 ) unless $USER->user_name eq "public";
	my %expids;
#	$dsgids = read_file() if $BASEFILE;    #: $opts{feature_list};
	foreach ( $form->param('dsgid') ) {
		foreach my $dsgid ( split /(::)|(,)/ ) {
			#next if $dsgid =~ /^\d+_?\d*$/;
			my $dsg  = $coge->resultset('Genome')->find($dsgid);
			next unless ($dsg);
			foreach my $exp ( $dsg->experiments ) {
				$expids{$exp->id}++;
			}
		}
	}
	
#	$expids = read_file() if $BASEFILE;    #: $opts{feature_list};
	foreach ( $form->param('expid') ) {
		foreach my $expid ( split /(::)|(,)/ ) {
			$expids{$expid}++;# if $expid =~ /^\d+_?\d*$/;
		}
	}

	my ( $table, $count ) =
	  generate_table( expids => [sort {$a<=>$b} keys %expids] );

	$template->param( 'EXPERIMENT_COUNT' => $count );

	if ($table) {
		$template->param( INFO => $table );
		return $template->output;
	}
	else {
		return "No dataset_group (genome) or experiment ids were specified.";
	}
}

sub add_to_user_history {
	my %opts = @_;

	#    my $url_params = $ENV{'QUERY_STRING'};
	my $url = $opts{url};
	if ( $opts{archive} ) {
		$USER->add_to_works(
			{
				'name'        => $opts{work_name},
				'archive'     => $opts{archive},
				'page'        => $PAGE_NAME,
				'parameter'   => $url,
				'description' => $opts{description},
				'note'        => $opts{note},
			}
		);
	}
	else {
		$USER->add_to_works(
			{
				'name'        => 'FeatList-' . $DATE,
				'archive'     => 0,
				'page'        => $PAGE_NAME,
				'parameter'   => $url,
				'description' => 'Feature List created on ' . $DATE,
			}
		);
	}
	return $opts{work_name};
}

sub generate_table {
	my %opts   = @_;
	my $expids = $opts{expids};
	return unless @$expids;
	my @table;
	my $count = 1;
	
	foreach my $expid (@$expids) {
		my $exp  = $coge->resultset('Experiment')->find($expid);
		next unless $exp;
		
		# Build sub-table of types
		my %types;
		foreach ($exp->types) {
			my $key = $_->name . ($_->desc ? ': ' . $_->desc : '');
			$types{$key}++;
		}
		my $type_tbl = '<table class="small"><tbody>';
		foreach (sort keys %types) {
			$type_tbl .= "<tr><td>$_</td></tr>";	
		}
		$type_tbl .= '</tbody></table>';
		
		# Build sub-table of annotations
		%types = ();
		foreach ($exp->annotations) {
			my $type = $_->annotation_type;
			my $group = $type->group;
			my $groupname = ($group ? $group->name : '');
			my $typename = $type->name;
			push @{$types{$groupname}{$typename}}, $_->annotation;
		}
		my $annot_tbl .= '<table class="small"><tbody>';
		foreach my $groupname (sort keys %types) {
			my $first1 = 1;
			foreach my $typename (sort keys %{$types{$groupname}}) {
				my $first2 = 1;
				foreach my $annot (sort @{$types{$groupname}{$typename}}) {
					$annot_tbl .= '<tr>';
					$annot_tbl .= '<td>' . ($first1 ? $groupname : '') . '</td>' if ($groupname);
					$annot_tbl .= '<td>' . ($first2 ? $typename  : '') . '</td>';
					$annot_tbl .= '<td></td>' if (not $groupname);
					$annot_tbl .= '<td>' . $annot . '</td></tr>';
					$first1 = $first2 = 0;
				}
			}
		}
		$annot_tbl .= '</tbody></table>';
		
		# Build source link
		my $src = $exp->source;
		my $link = $src->link;
		$link = 'http://' . $link if (not $link =~ /^http:\/\//);
		my $source = "<span class=link onclick=window.open('$link');>" . $src->name . ': ' . $src->desc . "</span>";
		
		push @table, 
		{
			COUNT		=> $count,
			EXPID		=> $expid,
			NAME		=> $exp->name,
			DESC		=> $exp->desc,
			TYPE		=> $type_tbl,
			ANNOTATION	=> $annot_tbl,
			SOURCE		=> $source,
			VER			=> $exp->version
		};
		$count++;
	}
	$count--;

	return \@table, $count;
}

sub get_feature_counts {
	my %opts  = @_;
	my $dsgid = $opts{dsgid};
	return "No specified dataset group id" unless $dsgid;
	my $dsg   = $coge->resultset('Genome')->find($dsgid);
	my $query = qq{
SELECT count(distinct(feature_id)), ft.name, ft.feature_type_id
  FROM feature
  JOIN feature_type ft using (feature_type_id)
  JOIN dataset_connector dc using (dataset_id)
 WHERE genome_id = $dsgid
  GROUP BY ft.name

};
	my $dbh = DBI->connect( $connstr, $DBUSER, $DBPASS );
	my $sth = $dbh->prepare($query);
	$sth->execute;
	my $feats = {};

	while ( my $row = $sth->fetchrow_arrayref ) {
		my $name = $row->[1];
		$name =~ s/\s+/&nbsp/g;
		$feats->{$name} = {
			count => $row->[0],
			id    => $row->[2],
			name  => $name,
		};
	}
	my $feat_string .= qq{<table class="small ui-widget-content ui-corner-all">
<thead></thead><tbody>};
	$feat_string .= "<tr valign=top>" . join(
		"\n<tr valign=top>",
		map {
			    "<td valign=top><div id=$_>"
			  . $feats->{$_}{name}
			  . "</div>"
			  . "<td valign=top align=right>"
			  . $feats->{$_}{count}
		  } sort { $a cmp $b } keys %$feats
	);
	$feat_string .= "</tbody></table>";
	$feat_string .= "None" unless keys %$feats;
	return $feat_string;
}

sub gen_data {
	my $message = shift;
	return qq{<font class="loading">$message. . .</font>};
}

sub read_file {
	my $file = "$TEMPDIR/$BASEFILE.featlist";
	my @featlist;
	unless ( -r $file ) {
		warn "unable to read file $file for feature ids\n";
		return \@featlist;
	}
	open( IN, $file ) || die "can't open $file for reading: $!";
	while (<IN>) {
		chomp;
		push @featlist, $_;
	}
	close IN;
	return \@featlist;
}

sub send_to_ExperimentList    #send to ExperimentList
{
	my %opts      = @_;
	my $accn_list = $opts{accn};
	$accn_list =~ s/^,//;
	$accn_list =~ s/,$//;
	my $url = $URL . "ExperimentList.pl?dsgid=$accn_list";
	return $url;
}

sub generate_excel_file { # FIXME mdb
	my %args      = @_;
	my $accn_list = $args{accn};
	$accn_list =~ s/^,//;
	$accn_list =~ s/,$//;
	$cogeweb = CoGe::Accessory::Web::initialize_basefile( tempdir => $TEMPDIR );
	my $basename = $cogeweb->basefilename;
	my $file     = "$TEMPDIR/Excel_$basename.xls";
	my $workbook = Spreadsheet::WriteExcel->new($file);
	$workbook->set_tempdir("$TEMPDIR");
	my $worksheet = $workbook->add_worksheet();
	my $i         = 1;

	$worksheet->write( 0, 0,  "Name" );
	$worksheet->write( 0, 1,  "Description" );
	$worksheet->write( 0, 2,  "Source" );
	$worksheet->write( 0, 3,  "Provenance" );
	$worksheet->write( 0, 4,  "Sequence Type" );
	$worksheet->write( 0, 5,  "Chr Count" );
	$worksheet->write( 0, 6,  "Length (bp)" );
	$worksheet->write( 0, 7,  "Percent GC" );
	$worksheet->write( 0, 8,  "Percent AT" );
	$worksheet->write( 0, 9,  "Percent N|X" );
	$worksheet->write( 0, 10, "OrganismView Link" );

	foreach my $dsgid ( split /,/, $accn_list ) {

		my ($dsg) = $coge->resultset("Genome")->find($dsgid);

		next unless $dsg;
		my $name = $dsg->name ? $dsg->name : $dsg->organism->name;
		my $desc =
		  $dsg->description ? $dsg->description : $dsg->organism->description;
		my ($ds_source) = $dsg->source;
		my $source = $ds_source->name;
		my $provenance = join( " ", map { $_->name } $dsg->datasets );
		my $length     = $dsg->length;
		my $chr_count  = $dsg->chromosome_count;
		my $type       = $dsg->type->name;
		my ( $gc, $at, $n ) = $dsg->percent_gc();
		$at *= 100;
		$gc *= 100;

		#my ($wgc, $wat) = $dsg->wobble_content;
		#$wat*=100;
		#$wgc*=100;

		$worksheet->write( $i, 0, $name );
		$worksheet->write( $i, 1, $desc );
		$worksheet->write( $i, 2, $source );
		$worksheet->write( $i, 3, $provenance );
		$worksheet->write( $i, 4, $type );
		$worksheet->write( $i, 5, $chr_count );
		$worksheet->write( $i, 6, $length );
		$worksheet->write( $i, 7, $gc . '%' );
		$worksheet->write( $i, 8, $at . '%' );
		$worksheet->write( $i, 9, $n . '%' );
		$worksheet->write( $i, 10,
			$P->{SERVER} . 'OrganismView.pl?dsgid=' . $dsgid );

		$i++;
	}
	$workbook->close() or die "Error closing file: $!";
	$file =~ s/$TEMPDIR/$TEMPURL/;
	return $file;
}

sub generate_csv_file { # FIXME mdb 
	my %args      = @_;
	my $accn_list = $args{accn};
	$accn_list =~ s/^,//;
	$accn_list =~ s/,$//;
	$cogeweb = CoGe::Accessory::Web::initialize_basefile( tempdir => $TEMPDIR );
	my $basename = $cogeweb->basefilename;
	my $file     = "$TEMPDIR/$basename.csv";
	open( OUT, ">$file" );
	print OUT join( "\t",
		"CoGe Genome ID",
		"Name",
		"Description",
		"Source",
		"Provenance",
		"Sequence Type",
		"Chr Count",
		"Length (bp)",
		"Percent GC",
		"Percent AT",
		"Percent N|X",
		"OrganismView Link" ),
	  "\n";

	foreach my $dsgid ( split /,/, $accn_list ) {
		next unless $dsgid;
		my ($dsg) = $coge->resultset("Genome")->find($dsgid);
		next unless $dsg;
		my $name = $dsg->name ? $dsg->name : $dsg->organism->name;
		my $desc =
		  $dsg->description ? $dsg->description : $dsg->organism->description;
		my ($ds_source) = $dsg->source;
		my $source = $ds_source->name;
		my $provenance = join( "||", map { $_->name } $dsg->datasets );
		my $chr_count  = $dsg->chromosome_count;
		my $length     = $dsg->length;
		my $type       = $dsg->type->name;
		my ( $gc, $at, $n ) = $dsg->percent_gc();
		$at *= 100;
		$gc *= 100;
		print OUT join( "\t",
			$dsgid,      $name,
			$desc,       $source,
			$provenance, $type,
			$chr_count,  $length,
			$gc,         $at,
			$n,          $P->{SERVER} . 'OrganismView.pl?dsgid=' . $dsgid ),
		  "\n";
	}
	close OUT;
	$file =~ s/$TEMPDIR/$TEMPURL/;
	return $file;
}

sub save_FeatList_settings {
	my %opts    = @_;
	my $display = $opts{display};
	my %save;
	if ($display) {
		my %settings = (
			1  => 'FeatNameD',
			2  => 'TypeD',
			3  => 'ChrD',
			4  => 'StartD',
			5  => 'StopD',
			6  => 'StrandD',
			7  => 'LengthD',
			8  => 'GCD',
			9  => 'ATD',
			10 => 'WGCD',
			11 => 'WATD',
			12 => 'OrgD',
			13 => 'AnnoD',
			14 => 'OtherD',
		);
		foreach my $index ( split /,/, $display ) {
			$save{display}{ $settings{$index} } = 1;
		}
	}
	CoGe::Accessory::Web::save_settings(
		opts => \%save,
		user => $USER,
		page => $PAGE_NAME,
		coge => $coge
	);
}

sub commify {
	my $text = reverse $_[0];
	$text =~ s/(\d\d\d)(?=\d)(?!\d*\.)/$1,/g;
	return scalar reverse $text;
}
