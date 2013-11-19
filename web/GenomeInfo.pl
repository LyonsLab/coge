#! /usr/bin/perl -w

use strict;
use CGI;
use CoGeX;
use CoGe::Accessory::Web;
use CoGe::Accessory::Utils;
use CoGe::Accessory::IRODS;
use HTML::Template;
use JSON::XS;
use Sort::Versions;
use File::Path qw(mkpath);
no warnings 'redefine';

use vars qw(
  $P $PAGE_TITLE $TEMPDIR $LOAD_ID $USER $CONFIGFILE $coge $FORM %FUNCTION $MAX_SEARCH_RESULTS $LINK $node_types
);

$PAGE_TITLE = 'GenomeInfo';

#EL: 10/31/13:  change to a global var
#my $node_types = CoGeX::node_types();
$node_types = CoGeX::node_types();


$FORM = new CGI;
( $coge, $USER, $P, $LINK ) = CoGe::Accessory::Web->init(
    cgi => $FORM,
    page_title => $PAGE_TITLE
);

$LOAD_ID = ( $FORM->Vars->{'load_id'} ? $FORM->Vars->{'load_id'} : get_unique_id() );
$TEMPDIR = $P->{SECTEMPDIR} . $PAGE_TITLE . '/' . $USER->name . '/' . $LOAD_ID . '/';
$CONFIGFILE = $ENV{COGE_HOME} . '/coge.conf';

$MAX_SEARCH_RESULTS = 100;

%FUNCTION = (
    get_genome_info    => \&get_genome_info,
    get_genome_data    => \&get_genome_data,
    edit_genome_info   => \&edit_genome_info,
    update_genome_info => \&update_genome_info,
    update_owner       => \&update_owner,
    search_organisms   => \&search_organisms,
    search_users       => \&search_users,
    delete_genome      => \&delete_genome,
    check_login        => \&check_login,
    copy_genome        => \&copy_genome,
	get_log            => \&get_log,
	export_fasta_irods => \&export_fasta_irods,
);

CoGe::Accessory::Web->dispatch( $FORM, \%FUNCTION, \&generate_html );

sub get_genome_info {
    my %opts   = @_;
    my $gid    = $opts{gid};
    my $genome = $opts{genome};
    return unless ( $gid or $genome );

    unless ($genome) {
        $genome = $coge->resultset('Genome')->find($gid);
        return unless ($genome);
    }

    my $template =
      HTML::Template->new( filename => $P->{TMPLDIR} . $PAGE_TITLE . '.tmpl' );
    $template->param(
        DO_GENOME_INFO => 1,
        ORGANISM       => $genome->organism->name,
        VERSION        => $genome->version,
        TYPE           => $genome->type->info,
        SOURCE         => join( ',', map { $_->name } $genome->source ),
        LINK           => $genome->link,
        RESTRICTED     => ( $genome->restricted ? 'Yes' : 'No' ),
        USERS_WITH_ACCESS => ( $genome->restricted ? join(', ', map { $_->display_name } $USER->users_with_access($genome))
                                                   : 'Everyone' ),
        NAME           => $genome->name,
        DESCRIPTION    => $genome->description,
        DELETED        => $genome->deleted
    );

    $template->param( GID => $genome->id );

    return $template->output;
}

sub edit_genome_info {
    my %opts = @_;
    my $gid  = $opts{gid};
    return unless ($gid);

    my $genome = $coge->resultset('Genome')->find($gid);
    return unless ($genome);

    my $template =
      HTML::Template->new( filename => $P->{TMPLDIR} . $PAGE_TITLE . '.tmpl' );
    $template->param(
        EDIT_GENOME_INFO => 1,
        ORGANISM         => $genome->organism->name,
        VERSION          => $genome->version,
        TYPE             => $genome->type->name,
        SOURCE           => join( ',', map { $_->name } $genome->source ),
        LINK             => $genome->link,
        RESTRICTED       => $genome->restricted,
        NAME             => $genome->name,
        DESCRIPTION      => $genome->description
    );

    $template->param(
        TYPES   => get_sequence_types( $genome->type->id ),
        SOURCES => get_sources()
    );

    return $template->output;
}

sub update_genome_info {
    my %opts        = @_;
    my $gid         = $opts{gid};
    my $name        = $opts{name};
    my $description = $opts{description};
    my $version     = $opts{version};
    my $type_id     = $opts{type_id};
    my $restricted  = $opts{restricted};
    my $org_name    = $opts{org_name};
    my $source_name = $opts{source_name};
    my $link        = $opts{link};
    my $timestamp   = $opts{timestamp};

# print STDERR "gid=$gid organism=$org_name version=$version source=$source_name\n";
    return "Error: missing params."
      unless ( $gid and $org_name and $version and $source_name );

    my $genome = $coge->resultset('Genome')->find($gid);
    return "Error: can't find genome." unless ($genome);

    my $organism = $coge->resultset('Organism')->find( { name => $org_name } );
    return "Error: can't find organism." unless ($organism);

    my $source =
      $coge->resultset('DataSource')->find( { name => $source_name } );
    return "Error: can't find source." unless ($source);

    $genome->organism_id( $organism->id );
    $genome->name($name);
    $genome->description($description);
    $genome->version($version);
    $genome->link($link);
    $genome->genomic_sequence_type_id($type_id);
    $genome->restricted( $restricted eq 'true' );

    foreach my $ds ( $genome->datasets ) {
        $ds->data_source_id( $source->id );
        $ds->update;
    }

    $genome->update;

    return;
}

sub search_organisms {
    my %opts        = @_;
    my $search_term = $opts{search_term};
    my $timestamp   = $opts{timestamp};

    #	print STDERR "$search_term $timestamp\n";
    return unless $search_term;

    # Perform search
    $search_term = '%' . $search_term . '%';
    my @organisms = $coge->resultset("Organism")->search(
        \[
            'name LIKE ? OR description LIKE ?',
            [ 'name',        $search_term ],
            [ 'description', $search_term ]
        ]
    );

    # Limit number of results displayed
    if ( @organisms > $MAX_SEARCH_RESULTS ) {
        return encode_json( { timestamp => $timestamp, items => undef } );
    }

    my %unique = map { $_->name => 1 } @organisms;
    return encode_json(
        { timestamp => $timestamp, items => [ sort keys %unique ] } );
}

sub search_users {
    my %opts        = @_;
    my $search_term = $opts{search_term};
    my $timestamp   = $opts{timestamp};

    #print STDERR "$search_term $timestamp\n";
    return unless $search_term;

    # Perform search
    $search_term = '%' . $search_term . '%';
    my @users = $coge->resultset("User")->search(
        \[
            'user_name LIKE ? OR first_name LIKE ? OR last_name LIKE ?',
            [ 'user_name',  $search_term ],
            [ 'first_name', $search_term ],
            [ 'last_name',  $search_term ]
        ]
    );

    # Limit number of results displayed
    # if (@users > $MAX_SEARCH_RESULTS) {
    # 	return encode_json({timestamp => $timestamp, items => undef});
    # }

    return encode_json(
        {
            timestamp => $timestamp,
            items     => [ sort map { $_->user_name } @users ]
        }
    );
}

sub update_owner {
    my %opts      = @_;
    my $gid       = $opts{gid};
    my $user_name = $opts{user_name};
    return unless $gid and $user_name;

    # Admin-only function
    return unless $USER->is_admin;

    # Make new user owner of genome
    my $user = $coge->resultset('User')->find( { user_name => $user_name } );
    unless ($user) {
        return "error finding user '$user_name'\n";
    }

    my $conn = $coge->resultset('UserConnector')->find_or_create(
        {
            parent_id   => $user->id,
            parent_type => $node_types->{user},
            child_id    => $gid,
            child_type  => $node_types->{genome},
            role_id     => 2                        # FIXME hardcoded
        }
    );
    unless ($conn) {
        return "error creating user connector\n";
    }

    # Remove admin user as owner
    $conn = $coge->resultset('UserConnector')->find(
        {
            parent_id   => $USER->id,
            parent_type => $node_types->{user},
            child_id    => $gid,
            child_type  => $node_types->{genome},
            role_id     => 2                        # FIXME hardcoded
        }
    );
    if ($conn) {
        $conn->delete;
    }

    return;
}

sub get_genome_data {
    my %opts   = @_;
    my $gid    = $opts{gid};
    my $genome = $opts{genome};
    return unless ( $gid or $genome );

    $gid = $genome->id if $genome;
    unless ($genome) {
        $genome = $coge->resultset('Genome')->find($gid);
        return unless ($genome);
    }

    my $template =
      HTML::Template->new( filename => $P->{TMPLDIR} . $PAGE_TITLE . '.tmpl' );
    $template->param(
        DO_GENOME_DATA   => 1,
        CHROMOSOME_COUNT => commify( $genome->chromosome_count() ),
        LENGTH           => commify( $genome->length ),
        GID              => $genome->id,
        LOGON            => ( $USER->user_name ne "public" )
    );

    return $template->output;
}

sub get_genome_download_links {
    my $genome = shift;

}

sub get_sequence_types {
    my $type_id = shift;

    my $html;
    foreach my $type ( sort { $a->info cmp $b->info }
        $coge->resultset('GenomicSequenceType')->all() )
    {
        $html .=
            '<option value="'
          . $type->id . '"'
          . ( defined $type_id && $type_id == $type->id ? ' selected' : '' )
          . '>'
          . $type->info
          . '</option>';
    }

    return $html;
}

sub get_sources {

    #my %opts = @_;

    my %unique;
    foreach ( $coge->resultset('DataSource')->all() ) {
        $unique{ $_->name }++;
    }

    return encode_json( [ sort keys %unique ] );
}

sub get_experiments {
    my %opts   = @_;
    my $gid    = $opts{gid};
    my $genome = $opts{genome};
    return unless ( $gid or $genome );

    unless ($genome) {
        $genome = $coge->resultset('Genome')->find($gid);
        return unless ($genome);
    }

    my @experiments = $genome->experiments;
    return "" unless @experiments;

    my @rows;
    foreach my $exp ( sort experimentcmp @experiments ) {
        next if ( $exp->deleted );

#next if ($exp->restricted && !$USER->is_admin && !$is_user && !$USER->has_access(experiment => $exp));

        my $id = $exp->id;
        my %row;
        $row{EXPERIMENT_INFO} =
qq{<span class="link" onclick='window.open("ExperimentView.pl?eid=$id")'>}
          . $exp->info
          . "</span>";

        push @rows, \%row;
    }

    my $template =
      HTML::Template->new( filename => $P->{TMPLDIR} . "$PAGE_TITLE.tmpl" );
    $template->param(
        DO_EXPERIMENTS  => 1,
        EXPERIMENT_LOOP => \@rows
    );
    return $template->output;
}

sub filter_dataset {    # detect chromosome-only datasets
    my $ds = shift;

    my @types = $ds->distinct_feature_type_ids;
    return ( @types <= 1 and shift(@types) == 4 );    #FIXME hardcoded type
}

sub get_datasets {
    my %opts        = @_;
    my $gid         = $opts{gid};
    my $genome      = $opts{genome};
    my $exclude_seq = $opts{exclude_seq};
    return unless ( $gid or $genome );

    unless ($genome) {
        $genome = $coge->resultset('Genome')->find($gid);
        return unless ($genome);
    }

    my @rows;
    foreach my $ds ( sort { $a->id <=> $b->id } $genome->datasets ) {

#next if ($exclude_seq && filter_dataset($ds)); #FIXME add dataset "type" field instead?
        push @rows, { DATASET_INFO => '<span>' . $ds->info . '</span>' };
    }
    return '' unless @rows;

    my $template =
      HTML::Template->new( filename => $P->{TMPLDIR} . "$PAGE_TITLE.tmpl" );
    $template->param(
        DO_DATASETS  => 1,
        DATASET_LOOP => \@rows
    );
    return $template->output;
}

sub delete_genome {
    my %opts = @_;
    my $gid  = $opts{gid};
    print STDERR "delete_genome $gid\n";
    return 0 unless $gid;

    my $genome = $coge->resultset('Genome')->find($gid);
    return 0 unless $genome;
    return 0 unless ( $USER->is_admin or $USER->is_owner( dsg => $gid ) );
    my $delete_or_undelete = ($genome->deleted ? 'undelete' : 'delete');
    #print STDERR "delete_genome " . $genome->deleted . "\n";
    $genome->deleted( !$genome->deleted ); # do undelete if already deleted
    $genome->update;
    
    # Record in log
    CoGe::Accessory::Web::log_history(
	    db          => $coge,
	    user_id     => $USER->id,
	    page        => $PAGE_TITLE,
	    description => "$delete_or_undelete genome id$gid"
	);

    return 1;
}

sub check_login {
	#print STDERR $USER->user_name . ' ' . int($USER->is_public) . "\n";
	return ($USER && !$USER->is_public);
}

sub copy_genome {
    my %opts = @_;
    my $gid  = $opts{gid};
    my $mask = $opts{mask};
    $mask = 0 unless $mask;

    print STDERR "copy_and_mask_genome: gid=$gid mask=$mask\n";

    if ($USER->is_public) {
    	return 'Not logged in';
    }

    # Setup staging area and log file
    my $stagepath = $TEMPDIR . '/staging/';
    mkpath $stagepath;

    my $logfile = $stagepath . '/log.txt';
    open( my $log, ">$logfile" ) or die "Error creating log file";
    print $log "Calling copy_load_mask_genome.pl ...\n";
    my $cmd =
        $P->{SCRIPTDIR} . "/copy_genome/copy_load_mask_genome.pl "
    	. "-gid $gid "
    	. "-uid " . $USER->id . " "
    	. "-mask $mask " 
    	. "-staging_dir $stagepath "
    	. "-conf_file $CONFIGFILE";
    print STDERR "$cmd\n";
    print $log "$cmd\n";
    close($log);

    if ( !defined( my $child_pid = fork() ) ) {
        return "Cannot fork: $!";
    }
    elsif ( $child_pid == 0 ) {
        print STDERR "child running: $cmd\n";
        `$cmd`;
        exit;
    }

    return;
}

sub get_log {
    #my %opts    = @_;
    #print STDERR "get_log $LOAD_ID\n";

    my $logfile = $TEMPDIR . "staging/log.txt";
    open( my $fh, $logfile ) or
      return encode_json( { status => -1, log => ["Error opening log file"] } );

    my @lines = ();
    my $gid;
    my $status = 0;
    while (<$fh>) {
        push @lines, $1 if ( $_ =~ /^log:\s+(.+)/i );
        if ( $_ =~ /log: Finished copying/i ) {
            $status = 1;
            last;
        }
        elsif ( $_ =~ /log: Added genome id(\d+)/i ) {
            $gid = $1;
        }
        elsif ( $_ =~ /log: error/i ) {
            $status = -1;
            last;
        }
    }
    close($fh);

    return encode_json(
        { status => $status, genome_id => $gid, log => \@lines } );
}

sub export_fasta_irods {
    my %opts    = @_;
    my $gid = $opts{gid};
    #print STDERR "export_fasta_irods $gid\n";

    my $genome = $coge->resultset('Genome')->find($gid);
	return 0 unless ($USER->has_access_to_genome($genome));

	my $src = $genome->file_path;
	my $dest = get_irods_path() . "/genome_$gid.faa";
	
	unless ($src and $dest) {
		print STDERR "GenomeInfo:export_fasta_irods: error, undef src or dest\n";
		return;
	}
	
	# Send to iPlant Data Store using iput
	CoGe::Accessory::IRODS::irods_iput($src, $dest);
	#TODO need to check rc of iput and abort if failure occurred

	# Set IRODS metadata for object
	my %meta = (
		    'Imported From' => "CoGe: http://genomevolution.org",
		    'CoGe OrganismView Link' => "http://genomevolution.org/CoGe/OrganismView.pl?gid=".$genome->id,
		    'CoGe GenomeInfo Link'=> "http://genomevolution.org/CoGe/GenomeInfo.pl?gid=".$genome->id,
		    'CoGe Genome ID'   => $genome->id,
		    'Organism Name'    => $genome->organism->name,
		    'Organism Taxonomy'    => $genome->organism->description,
		    'Version'     => $genome->version,
		    'Type'        => $genome->type->info,
		   );
    my $i = 1;
    my @sources = $genome->source;
    foreach my $item (@sources)
      {
	my $source = $item->name;
	$source.= ": ".$item->description if $item->description;
	my $met_name = "Source";
	$met_name .= $i if scalar @sources > 1;
	$meta{$met_name}= $source;
	$meta{$met_name." Link"} = $item->link if $item->link;
	$i++;
      }


	$meta{'Genome Link'} = $genome->link if ($genome->link);
	$meta{'Addition Info'} = $genome->message if ($genome->message);
	$meta{'Genome Name'} = $genome->name if ($genome->name);
	$meta{'Genome Description'} = $genome->description if ($genome->description);

	CoGe::Accessory::IRODS::irods_imeta($dest, \%meta);
}

sub get_irods_path {
	my $username = $USER->user_name;
	my $dest = $P->{IRODSDIR};
    $dest =~ s/\<USER\>/$username/;
    return $dest;
}

sub generate_html {
    my $name = $USER->user_name;
    $name = $USER->first_name if $USER->first_name;
    $name .= ' ' . $USER->last_name
      if ( $USER->first_name && $USER->last_name );

    my $template =
      HTML::Template->new( filename => $P->{TMPLDIR} . 'generic_page.tmpl' );
    $template->param(
        PAGE_TITLE => $PAGE_TITLE,
        PAGE_LINK  => $LINK,
        HELP       => '/wiki/index.php?title=' . $PAGE_TITLE . '.pl',
        USER       => $name,
        LOGO_PNG   => $PAGE_TITLE . "-logo.png",
        BODY       => generate_body(),
        ADJUST_BOX => 1,
        LOGON      => ( $USER->user_name ne "public" )
    );

    return $template->output;
}

sub generate_body {
    my $template =
      HTML::Template->new( filename => $P->{TMPLDIR} . $PAGE_TITLE . '.tmpl' );
    $template->param( MAIN => 1, PAGE_NAME => "$PAGE_TITLE.pl" );

    my $gid = $FORM->param('gid');
    return "No genome specified" unless $gid;

    my $genome = $coge->resultset('Genome')->find($gid);
    return "Genome id$gid not found" unless ($genome);
    return "Access denied" unless $USER->has_access_to_genome($genome);

    my $user_can_edit = $USER->is_admin || $USER->is_owner_editor( dsg => $gid );
    my $user_can_delete = $USER->is_admin || $USER->is_owner( dsg => $gid );

    $template->param(
    	LOAD_ID         => $LOAD_ID,
        GID             => $gid,
        GENOME_INFO     => get_genome_info( genome => $genome ),
        GENOME_DATA     => get_genome_data( genome => $genome ),
        EXPERIMENTS     => get_experiments( genome => $genome ),
        ANNOTATION      => get_datasets( genome => $genome, exclude_seq => 1 ),
        USER_CAN_EDIT   => $user_can_edit,
        USER_CAN_ADD    => ( !$genome->restricted or $user_can_edit ),
        USER_CAN_DELETE => $user_can_delete,
        DELETED         => $genome->deleted,
        IRODS_HOME      => get_irods_path()
    );

    if ( $USER->is_admin ) {
        $template->param(
            ADMIN_AREA => 1,
        );
    }

    return $template->output;
}

# FIXME this routine is duplicated elsewhere
sub experimentcmp {
    no warnings 'uninitialized';    # disable warnings for undef values in sort
    versioncmp( $b->version, $a->version )
      || $a->name cmp $b->name
      || $b->id cmp $a->id;
}
