#! /usr/bin/perl -w
use strict;
use CGI;
use CGI::Carp 'fatalsToBrowser';
use CoGeX;
use CoGe::Accessory::Web;
use CoGe::Accessory::Utils qw( commify );
use CoGe::Core::Genome qw(genomecmp genomecmp2);
use HTML::Template;
use Data::Dumper;
#use CGI::Ajax;
use JSON::XS;
use Benchmark;
use File::Path;
use Benchmark qw(:all);
use Statistics::Basic::Mean;
use POSIX;
use Sort::Versions;

no warnings 'redefine';

use vars qw(
    $P $PAGE_NAME $PAGE_TITLE $LINK
    $TEMPDIR $TEMPURL $USER $FORM $coge $HISTOGRAM
    %FUNCTION $P $SERVER $MAX_NUM_ORGANISM_RESULTS $MAX_NUM_CHROMOSOME_RESULTS
    $MAX_DS_LENGTH $EMBED
);

$| = 1;    # turn off buffering

$PAGE_TITLE = 'OrganismView';
$PAGE_NAME  = "$PAGE_TITLE.pl";

$FORM = new CGI;
( $coge, $USER, $P, $LINK ) = CoGe::Accessory::Web->init(
    cgi => $FORM,
    page_title => $PAGE_TITLE
);

$TEMPDIR   = $P->{TEMPDIR} . "/$PAGE_TITLE";
$TEMPURL   = $P->{TEMPURL} . "/$PAGE_TITLE";
$SERVER    = $P->{SERVER};
$HISTOGRAM = $P->{HISTOGRAM};

$MAX_NUM_ORGANISM_RESULTS = 500; #5000; # mdb changed 9/25/14, COGE-504
$MAX_NUM_CHROMOSOME_RESULTS = 50; # mdb changed 2/29/16 from 20
$MAX_DS_LENGTH = 10000000;

%FUNCTION = (
    get_genomes             => \&get_genomes,
    get_genome_info         => \&get_genome_info,
    get_datasets            => \&get_datasets,
    get_dataset_info        => \&get_dataset_info,
    get_chr_info            => \&get_chr_info,
    gen_data                => \&gen_data,
    get_orgs                => \&get_orgs,
    get_org_info            => \&get_org_info,
    get_start_stop          => \&get_start_stop,
    get_gc_for_chromosome   => \&get_gc_for_chromosome,
    get_gc_for_noncoding    => \&get_gc_for_noncoding,
    get_gc_for_feature_type => \&get_gc_for_feature_type,
    get_aa_usage            => \&get_aa_usage,
    get_wobble_gc           => \&get_wobble_gc,
    get_wobble_gc_diff      => \&get_wobble_gc_diff,
    get_total_length_for_ds => \&get_total_length_for_ds,
    get_chr_length_hist     => \&get_chr_length_hist,
    get_genome_name         => \&get_genome_name,
    add_to_irods            => \&add_to_irods,
    make_genome_public      => \&make_genome_public,
    make_genome_private     => \&make_genome_private,
);

CoGe::Accessory::Web->dispatch( $FORM, \%FUNCTION, \&gen_html );

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
    #    else{
    #  print $FORM->header, gen_html();
    #    }
}

sub gen_html {
    my $template;

    $EMBED = $FORM->param('embed') || 0;
    if ($EMBED) {
        $template = HTML::Template->new( filename => $P->{TMPLDIR} . 'embedded_page.tmpl' );
    }
    else {
        $template = HTML::Template->new( filename => $P->{TMPLDIR} . 'generic_page.tmpl' );
        $template->param( PAGE_TITLE => 'OrganismView',
            TITLE                    => 'OrganismView: Search Organisms and Genomes',
            PAGE_LINK                => $LINK,
            SUPPORT_EMAIL            => $P->{SUPPORT_EMAIL},
            HEAD                     => qq{},
            HOME                     => $P->{SERVER},
            HELP                     => 'OrganismView',
            WIKI_URL                 => $P->{WIKI_URL} || '',
            ADMIN_ONLY               => $USER->is_admin,
            USER                     => $USER->display_name || '',
            CAS_URL                  => $P->{CAS_URL} || '',
            COOKIE_NAME              => $P->{COOKIE_NAME} || ''
        );
        $template->param( LOGON => 1 ) unless ($USER->user_name eq "public");
    }

    my ( $body, $seq_names, $seqs ) = gen_body();
    $template->param( BODY => $body );

    return $template->output;
}

sub gen_body {
    my $form = shift || $FORM;
    my $org_name = $form->param('org_name');
    my $oid      = $form->param('oid');
    my $dsname   = $form->param('dsname');
    my $dsid     = $form->param('dsid');
    my $gid      = $form->param('dsgid');
       $gid      = $form->param('gid') if ( $form->param('gid') );

    # Initialize template
    my $template = HTML::Template->new( filename => $P->{TMPLDIR} . 'OrganismView.tmpl' );

    $template->param( EMBED => $EMBED );
    $template->param( ORG_SEARCH  => $org_name ) if $org_name;
    
    # Get organisms and insert into template
    my ( $org_list, $org_count, $selected_oid ) = get_orgs( name => $org_name, oid => $oid, gid => $gid, dsid => $dsid, dsname => $dsname, output => 'html' );
    if ($org_count) {
        $template->param( 
            ORG_LIST  => $org_list,
            ORG_COUNT => $org_count 
        );
        
        # Get organism info
        my $org_info = get_org_info(oid => $selected_oid, output => 'html');
        $template->param( ORG_INFO => $org_info ) if $org_info;
    }
    else {
        $template->param( NO_RESULTS  => $org_list );
        return $template->output;
    }
    
    my ( $genome_list, $genome_count, $selected_gid ) = get_genomes( oid => $selected_oid, gid => $gid, output => 'html' );
    $template->param( GENOME_LIST  => $genome_list );
    $template->param( GENOME_COUNT => $genome_count );
    
    my $genome_info = get_genome_info( gid => $selected_gid, output => 'html' );
    $template->param( GENOME_INFO => $genome_info );
    
    # Get datasets and insert into template
    my ( $dslist, $dscount, $selected_dsid ) = get_datasets( gid => $selected_gid, dsid => $dsid, dsname => $dsname, output => 'html' );
    $template->param( DS_LIST  => $dslist )  if $dslist;
    $template->param( DS_COUNT => $dscount ) if $dscount;

    my ( $ds_info, $chr_list, $chr_count, $selected_chr ) = get_dataset_info( dsid => $selected_dsid, output => 'html' );
    $template->param( DS_INFO  => $ds_info );
    $template->param( CHR_LIST => $chr_list );

    # Get chromosome info and insert into template
    my ( $chr_info, $viewer, $seqview ) = get_chr_info( dsid => $selected_dsid, chr => $selected_chr, output => 'html' );
    $template->param( CHR_INFO => $chr_info ) if ($chr_info);
    $template->param( VIEWER   => $viewer ) if ($viewer);
    $template->param( GET_SEQ  => $seqview ) if ($seqview);

    # Finish template
    $template->param( ORGANISM_ID => $selected_oid)  if ($selected_oid);
    $template->param( GENOME_ID   => $selected_gid ) if ($selected_gid);
    $template->param( DATASET_ID  => $dsid)  if ($dsid);
    $template->param( CHR_ID      => $selected_chr) if ($selected_chr);
    
    return $template->output;
}

sub make_genome_public {
    my %opts  = @_;
    my $dsgid = $opts{dsgid};
    return "No DSGID specified" unless $dsgid;
    return "Permission denied."
      unless $USER->is_admin || $USER->is_owner( dsg => $dsgid );
    my $dsg = $coge->resultset('Genome')->find($dsgid);
    $dsg->restricted(0);
    $dsg->update;
    foreach my $ds ( $dsg->datasets ) {
        $ds->restricted(0);
        $ds->update;
    }
    return 1;
}

sub make_genome_private {
    my %opts  = @_;
    my $dsgid = $opts{dsgid};
    return "No DSGID specified" unless $dsgid;
    return "Permission denied."
      unless $USER->is_admin || $USER->is_owner( dsg => $dsgid );
    my $dsg = $coge->resultset('Genome')->find($dsgid);
    $dsg->restricted(1);
    $dsg->update;
    foreach my $ds ( $dsg->datasets ) {
        $ds->restricted(1);
        $ds->update;
    }
    return 1;
}

sub get_orgs {
    my %opts = @_;
    my $name    = $opts{name};   # organism name
    my $desc    = $opts{desc};   # organism description
    my $oid     = $opts{oid};    # organism id
    my $gid     = $opts{gid};    # genome id
    my $dsid    = $opts{dsid};   # dataset id
    my $dsname  = $opts{dsname}; # dataset name
    my $output  = $opts{output} || 'json';
    #print STDERR "get_orgs\n", Dumper \%opts, "\n";

    my (%organisms, %genomes, %datasets);
    
    # Set organism by specified genome id
    my $genome;
    if ($gid) {
        my $genome = $coge->resultset('Genome')->find($gid);
        if ($USER->has_access_to_genome($genome)) {
            $genomes{$genome->id} = $genome;
            $organisms{$genome->organism_id} = $genome->organism;
        }
    }
    elsif ($dsid || $dsname) {
        my @datasets;
        
        if ($dsid) {
            my $ds = $coge->resultset('Dataset')->find($dsid);
            push @datasets, $ds if ($ds);
        }
        if ($dsname && length($dsname) > 5) {
            my $search_term = '%' . $dsname . '%';
            my @ds = $coge->resultset('Dataset')->search(\[
                    'dataset_id = ? OR name LIKE ? OR description LIKE ?',
                    [ 'dataset_id', ($dsid ? $dsid : 0) ], [ 'name', $search_term ], [ 'description', $search_term ]
                ]);
            push @datasets, @ds if (@ds);
        }
        
        foreach my $dataset (@datasets) {
            my @ds_genomes = $dataset->genomes;
            foreach my $genome (@ds_genomes) {
                if ($USER->has_access_to_genome($genome)) {
                    $genomes{$genome->id} = $genome;
                    $organisms{$genome->organism_id} = $genome->organism };
            }
        }
    }
    elsif ($oid || $name || $desc) {
        my @orgs;
        if ($oid) {
            my $org = $coge->resultset("Organism")->find($oid);
            push @orgs, $org if ($org);
        }
        else {
            my $search_term = '%' . $name . '%';
            @orgs = $coge->resultset("Organism")->search(
                \[
                    'name LIKE ? OR description LIKE ?',
                    [ 'name', $search_term ], [ 'description', $search_term ]
                ]
            );
        }
        
        # Limit results
        if (@orgs > $MAX_NUM_ORGANISM_RESULTS) {
            my $msg = qq{Too many results to show (}.@orgs.qq{), please refine your search.};
            return $msg if ($output eq 'html');
            return encode_json({ error => $msg });
        }
        elsif (@orgs) {
            %organisms = map { $_->id => $_ } @orgs;
        }
    }
    else {
        my $msg = 'Please enter a search term.';
        return $msg, 0 if ($output eq 'html');
        return encode_json({ error => $msg });
    }
    
    unless (keys %organisms) {
        my $msg = 'No matching results were found.';
        return $msg, 0 if ($output eq 'html');
        return encode_json({ error => $msg });
    }
    
    # Build html options list of organisms
    my @opts;
    foreach my $org ( sort { uc( $a->name ) cmp uc( $b->name ) } values %organisms ) {
        my $this_id = $org->id;
        my $this_name = $org->name;

        # Set selected option if first or specified by user - FIXME there is a bettery way
        my $selected = '';
        if ( ($oid && $this_id == $oid) || ($gid && $this_id == $genomes{$gid}->organism->id) ) {
            $selected = "selected";
            $oid = $this_id;
        }
        elsif ($dsid) {
            $selected = "selected";
            $oid = $this_id;
        }
        elsif (scalar(@opts) == 0) { # first in list
            $selected = "selected"; # index
            $oid = $this_id;    
        }
        
        $name = substr( $this_name, 0, 50 ) . "..." if ( length($this_name) > 50 );
        my $text = "$this_name (id$this_id)";
        my $option = qq{<option $selected title="$text" value="$this_id">$text</option>};
        push @opts, $option;
    }
    unless (@opts) {
        return ('', 0) if ($output eq 'html');
        return encode_json({ organisms => '', count => 0 });
    }
    
    my $html = join('\n', @opts);
#    print STDERR "html:\n$html\n";
    
    return $html, scalar(@opts), $oid if ($output eq 'html');
    return encode_json({ organisms => $html, count => scalar(@opts), selected_id => $oid });
}

sub get_genome_name {
    my %opts = @_;
    my $gid  = $opts{gid};
    return encode_json({ error => 'Missing gid parameter' }) unless $gid;
    
    my $genome = $coge->resultset("Genome")->find($gid);
    return encode_json({ error => 'Could not find the genome' }) unless $genome;
    return encode_json({ gid => $gid, info => $genome->info });
}

sub get_org_info {
    my %opts = @_;
    my $oid = $opts{oid};
    my $output = $opts{output} || 'json';
#    print STDERR "get_org_info\n", Dumper \%opts, "\n";
    
    return encode_json({error => "An organism was not specified"}) unless $oid;

    # Get organisms
    my $organism = $coge->resultset("Organism")->find($oid);
    unless ($organism) {
        return encode_json({ organism => "Unable to find an organism for id: $oid\n" });
    }

    my $html = 
          qq{<table class='small annotation_table' style="margin:3px;">}
        . qq{<tr><td>Name:</td>}
        . qq{<td>} . $organism->name . qq{</td>}
        . qq{</tr>};

    if ( $organism->description ) {
        $html .= qq{<tr><td class='top'>Description:<td>};

        my $desc_len;
        foreach my $item ( split( /;/, $organism->description ) ) {
            $item =~ s/^\s+//;
            $item =~ s/\s+$//;
            $html .= qq{<span class='link' onclick='window.open("OrganismView.pl?org_desc=$item");'>$item</span>;};
            $desc_len += length($item);
            if ( $desc_len > 100 ) {
                $html .= '<br>';
                $desc_len = 0;
            }
        }
    }
    
    $html .= qq{<tr>}
        . qq{<td>Tools:</td><td><span class="link" onclick="window.open('OrganismView.pl?oid=$oid', 'target=_new');">OrganismView</span>}
        . qq{&nbsp|&nbsp}
        . qq{<span class="link" onclick="window.open('CodeOn.pl?oid=$oid', 'target=_new')">CodeOn</span></td>}
        . qq{<tr><td>Search:</td><td>};
    my $search_term = $organism->name;
    $html .= qq{<img onclick="window.open('http://www.ncbi.nlm.nih.gov/taxonomy?term=$search_term')" src="picts/other/NCBI-icon.png" title="NCBI" class="link">&nbsp}
          .  qq{<img onclick="window.open('http://en.wikipedia.org/w/index.php?title=Special%3ASearch&search=$search_term')" src="picts/other/wikipedia-icon.png" title="Wikipedia" class="link">&nbsp};
    $search_term =~ s/\s+/\+/g;
    $html .= qq{<img onclick="window.open('http://www.google.com/search?q=$search_term')" src="picts/other/google-icon.png" title="Google" class="link">}
          .  qq{</table>};

    #print STDERR "html:\n$html";
    return $html if ($output eq 'html');
    return encode_json({ organism => $html });
}

sub get_genomes {
    my %opts  = @_;
    my $oid   = $opts{oid};
    my $gid = $opts{gid};
    my $output = $opts{output} || 'json';
#    print STDERR "get_genomes oid=",($oid ? $oid : '')," gid=",($gid ? $gid : ''), "\n";
    
    my $org = $coge->resultset("Organism")->find($oid);
    return unless $org;
    
    my $selected_gid;
    my @opts;
    
    if ($org) {
        my $favorites = CoGe::Core::Favorites->new(user => $USER);
        
        no warnings 'uninitialized'; # disable warnings for undef values in sort
        foreach my $genome (sort { genomecmp2($a, $b, $favorites) } $org->genomes) {
            next if $genome->deleted;
            next unless $USER->has_access_to_genome($genome);
            
            # Set option to selected if specified or first in list - FIXME rewrite this
            my $selected = '';
            if ($gid) { 
                if ($gid == $genome->id ) {
                    $selected = 'selected';
                    $selected_gid = $genome->id;
                }
            }
            elsif (scalar(@opts) == 0) {
                $selected = 'selected';
                $selected_gid = $genome->id;
            }
            
            my $option;
            $option .= "&#11088; "  if ($favorites->is_favorite($genome));
            $option .= "&#x2705; "  if $genome->certified;
            $option .= "&#x1f512; " if $genome->restricted;
            my $info = $genome->info( hideRestrictedSymbol => 1 );
            $option .= $info;
            push @opts, qq{<OPTION title="$info" value="}.$genome->id.qq{" $selected>$option</OPTION>};
        }
    }
    my $html;
    if (@opts) {
        $html .= join( "\n", @opts );
        $html =~ s/OPTION/OPTION SELECTED/ unless $html =~ /SELECTED/i;
    }

    my $message = "No genomes found";

    unless (scalar(@opts) > 0) {
        return $message, 0 if ($output eq 'html');
        return encode_json({ error => $message });
    }
    
    #print STDERR Dumper { genomes => $html, count => scalar(@opts), selected_id => $selected_gid }, "\n";
    return $html, scalar(@opts), $selected_gid if ($output eq 'html');
    return encode_json({ genomes => $html, count => scalar(@opts), selected_id => $selected_gid });
}

sub add_to_irods {
    my %opts             = @_;
    my $dsgid            = $opts{dsgid};
    my $dsg              = $coge->resultset('Genome')->find($dsgid);
    my $add_to_irods_bin = $P->{BINDIR} . "/irods/add_to_irods.pl";
    my $cmd              = $add_to_irods_bin . " -file " . $dsg->file_path;
    my $new_name         = $dsg->organism->name . " " . $dsg->id . ".faa";
    $cmd .= " -new_name '$new_name'";
    $cmd .= " -dir collections";
    $cmd .= " -tag 'organism=" . $dsg->organism->name . "'";
    $cmd .= " -tag version=" . $dsg->version;
    $cmd .= " -tag 'sequence_type=" . $dsg->sequence_type->name . "'";
    my ($ds) = $dsg->datasets;
    $cmd .= " -tag 'source_name=" . $ds->name . "'";
    $cmd .= " -tag 'source_link=" . $ds->link . "'" if $ds->link;
    $cmd .= " -tag 'imported_from=CoGe: http://genomevolution.org/CoGe/OrganismView.pl?dsgid=$dsgid'";
    system($cmd);
#    print STDERR $cmd;

    #return $cmd;
    return "Complete!";
}

sub get_genome_info {
    my %opts  = @_;
    my $gid = $opts{gid};
    my $output = $opts{output} || 'json';
#    print STDERR "get_genome_info\n", Dumper \%opts, "\n";
    
    unless ($gid) {
        return "A genome was not specified" if ($output eq 'html');
        return encode_json({ error  => "A genome was not specified" });
    }
    
    my $dsg = $coge->resultset("Genome")->find($gid);
    unless ($gid) {
        return "Unable to get genome object for id: $gid" if ($output eq 'html');
        return encode_json({ error => "Unable to get genome object for id: $gid" });
    }    

    my $total_length = $dsg->length;
    my $chr_num = $dsg->chromosome_count();
    
    my $html;
    $html .= qq{<table>};
    $html .= qq{<table class='small annotation_table'>};
    $html .= qq{<tr><td>Name:</td><td>} . $dsg->name . qq{</td></tr>} if $dsg->name;
    $html .= qq{<tr><td>Description:</td><td>} . $dsg->description . qq{</td></tr>} if $dsg->description;

    my $gstid    = $dsg->genomic_sequence_type->id;
    my $gst_name = $dsg->genomic_sequence_type->name;
    $gst_name .= ": " . $dsg->type->description if $dsg->type->description;

    my $owner = $dsg->owner;
    my $owner_field;
    if ($owner) {
        my $username = $owner->display_name;
        $owner_field = qq{<tr><td>Owner:</td><td>$username</td></tr>};
    }

    $html .=
        qq{<tr><td>Genome ID: </td><td>$gid</td>}
      . qq{<tr><td>Sequence type: </td>}
      . qq{<td>$gst_name (gstid$gstid) </td>}
      . qq{</tr>}
      . $owner_field
      . qq{<tr><td>Length: </td>}
      . qq{<td class='l'> }. commify($total_length) . qq{ bp</td>}
      . qq{</tr>}
      . qq{</div>};

    # mdb removed 7/31/13 issue 77
    #    my $seq_file = $dsg->file_path;
    #    my $cogedir  = $P->{COGEDIR};
    #    my $cogeurl  = $P->{URL};
    #    $seq_file =~ s/$cogedir/$cogeurl/i;
    #my $seq_url = api_url_for("genomes/$gid/sequence"); #"api/v1/legacy/sequence/$gid"; # mdb changed 2/12/16 for hypnotoad

    $html .= "<tr><td>Tools:</td>"
     . qq{<td>}
     . qq{<span class="link" onclick="window.open('GenomeInfo.pl?gid=$gid', '_blank');">GenomeInfo</span>}
     . qq{&nbsp|&nbsp}
     . qq{<span class="link" onclick="window.open('OrganismView.pl?dsgid=$gid', '_blank');">OrganismView</span>}
     . qq{&nbsp|&nbsp}
     . qq{<span class="link" onclick="window.open('CodeOn.pl?dsgid=$gid', '_blank');">CodeOn</span>}
     . qq{&nbsp|&nbsp}
     . qq{<span class='link' onclick="window.open('SynMap.pl?dsgid1=$gid;dsgid2=$gid', '_blank');">SynMap</span>}
     . qq{&nbsp|&nbsp}
     . qq{<span class='link' onclick="window.open('CoGeBlast.pl?dsgid=$gid', '_blank');">CoGeBlast</span>}
     . qq{&nbsp|&nbsp}
     . qq{<span class="link" onClick="add_genome_to_list($gid);">Add to GenomeList</span>}
     . qq{</td></tr>}
     #. qq{</table></td>}
     #. qq{<td id="dsg_features"></td>}
     . qq{</table>};

#temporarily removed until this is connected correctly for individual users
#    $html .= qq{&nbsp|&nbsp};
#    $html .= qq{<span id=irods class='link' onclick="gen_data(['args__loading...'],['irods']);add_to_irods(['args__dsgid','args__$gid'],['irods']);">Send To CyVerse Data Store</span>};

    return $html if ($output eq 'html');
    return encode_json({ genome => $html });
}

sub get_datasets {
    my %opts   = @_;
    my $gid    = $opts{gid};
    my $dsid   = $opts{dsid};
    my $dsname = $opts{dsname};
    my $output = $opts{output} || 'json';
#    print STDERR "get_datasets\n", Dumper \%opts, "\n";

    unless ($dsid || $dsname || $gid) {
        return "No datasets found" if ($output eq 'html');
        return encode_json({
            error => "No datasets found",
            counts => 0
        });
    }
    
    my @datasets;

    if ($gid) {
        my $genome = $coge->resultset("Genome")->find($gid);
        my @ds = $genome->datasets;
        push @datasets, @ds;
    }
    elsif ($dsid) {
        my $dataset = $coge->resultset("Dataset")->find($dsid);
        push @datasets, $dataset if ($dataset);
    }
    elsif ($dsname) {
        # Search datasets in db
        my $search_term = '%' . $dsname . '%' if ($dsname);
        @datasets = $coge->resultset("Dataset")->search(
            \[
                'name LIKE ? OR description LIKE ?',
                [ 'name', $search_term ], [ 'description', $search_term ]
            ]
        );
    }
    
    # Limit results
    if (@datasets > $MAX_NUM_ORGANISM_RESULTS) {
        my $results = qq{<option class="small alert">Too many results to show (}.scalar(@datasets).qq{), please refine your search</option>};
        return ( $results, scalar(@datasets) ) if ($output eq 'html');
        return encode_json({
            datasets => $results,
            count => scalar(@datasets)
        });
    }
    
    my @opts;
    my $selected_id;
    foreach my $item (sort { $b->version <=> $a->version || uc( $a->name ) cmp uc( $b->name ) } @datasets) {
        #next unless $USER->has_access_to_dataset($item);
        my $selected = '';
        if ($dsid && $dsid == $item->id) {
            $selected = 'selected'; # index
            $selected_id = $item->id;
        }
        elsif ( scalar(@opts) == 0 ) {
            $selected = 'selected'; # index
            $selected_id = $item->id;
        }
        push @opts,
            qq{<option $selected value="}
          . $item->id . qq{">}
          . $item->name . "(v"
          . $item->version . ", id"
          . $item->id
          . qq{)</option>};
    }

    my $html .= join( "\n", @opts );

    return $html, scalar(@opts), $selected_id if ($output eq 'html');
    return encode_json({ datasets => $html, count => scalar(@opts), selected_id => $selected_id });
}

sub get_dataset_info {
    my %opts = @_;
    my $dsid = $opts{dsid};
    my $output = $opts{output} || 'json';
#    print STDERR "get_dataset_info $dsid\n";
    
    # Get dataset
    my $ds;
    if ($dsid) {
        $ds = $coge->resultset("Dataset")->find($dsid);
    }
    unless ($ds and $dsid) {
        my $error = "Dataset dsid" . ( $dsid ? $dsid : '[undef]' ) . " not found";
        return $error,0 if ($output eq 'html');
        return encode_json({
            error => $error,
            counts => 0
        });
    }

    #
    # Get dataset info
    #
    my $ds_name  = $ds->name;
    $ds_name .= ": " . $ds->description if $ds->description;
    $ds_name = " <a href=\"" . $ds->link . "\" target=_new\>" . $ds_name . "</a>" if $ds->link;
    my $source_name = $ds->data_source->name;
    $source_name .= ": " . $ds->data_source->description
      if $ds->data_source->description;
    my $link = $ds->data_source->link;

    $link = "http://" . $link if ( $link && $link !~ /http/ );
    $source_name = qq{<a href="$link" target=_new>$source_name</a>} if $ds->data_source->link;
    my $html_ds = 
          qq{<table><tr valign='top'><td><table class="small annotation_table">}
        . qq{<tr>}
        . qq{<td>Name: <td>$ds_name (id} . $ds->id . qq{)</td>}
        . qq{<tr><td>Source: <td>$source_name (id} . $ds->data_source->id . qq{)}
        . qq{<tr><td>Version: <td>} . $ds->version
        . qq{<tr><td>Organism:<td class="link"><a href="OrganismView.pl?oid=}
        . $ds->organism->id
        . qq{" target=_new>}
        . $ds->organism->name 
        . qq{</a>}
        . qq{</tr>}
        . qq{<tr><td>Created: </td><td>} . $ds->date . qq{</td></tr>};
    #$html_ds .= qq{<div class='alert small'>Restricted dataset</div>} if $ds->restricted;

    #
    # Get chromosome list - FIXME why is this here and not in separate routine? mdb, 9/26/14
    #
    my $html_chr;
    my $total_length = $ds->total_length( ftid => 4 ); # FIXME hardcoded ftid value
    my $chr_num = $ds->chromosome_count( ftid => 4 );
    my %chr;
    my $tmp_count = 0;
    foreach my $c ( $ds->get_chromosomes( ftid => 4, length => 1, limit => $MAX_NUM_CHROMOSOME_RESULTS, max => 1000 ) ) {
        if (ref($c) =~ /CoGeX/ ) { $chr{$c->chromosome} = { length => $c->stop }; }
        else { $chr{ $c } = { length => 0 }; }
    }
    # Build chromosome options list
    my $selected_chr;
    if (keys %chr) {
        my $count = 100000;
        foreach my $item ( sort keys %chr ) {
            my ($num) = $item =~ /(\d+)/;
            $num = $count unless $num;
            $chr{$item}{num} = $num;
            $count++;
        }
        
        my @chr =
          $chr_num > $MAX_NUM_CHROMOSOME_RESULTS
          ? sort { $chr{$b}{length} <=> $chr{$a}{length} } keys %chr
          : sort { $chr{$a}{num} <=> $chr{$b}{num} || $a cmp $b } keys %chr;
      
        my @opts;
        foreach (@chr) {
            my $selected = '';
            if (@opts == 0) {
                $selected_chr = $_;
                $selected = 'selected';
            }
            
            push @opts, 
                qq{<option $selected value="$_">$_ (}.commify( $chr{$_}{length} ).qq{ bp)</option>};
        }
        $html_chr .= join('', @opts);
    }
    elsif ($chr_num) {
        $html_chr .= "<option>Too many chromosomes (".commify($chr_num).") to query</option>";
    }
    else {
        $html_chr .= "<option>No chromosomes</option>";
    }
    $html_ds .= "<tr><td>Chromosomes:</td><td><div style=\"float: left;\">" . commify($chr_num)
          .  "<tr><td>Total length:<td><div style=\"float: left;\">" . commify($total_length) . " bp ";
    my $gc = 
      $total_length && $total_length < $MAX_DS_LENGTH && $chr_num && $chr_num < $MAX_NUM_CHROMOSOME_RESULTS
      ? get_gc_for_chromosome( dsid => $ds->id )
      : 0;
    $html_ds .= $gc if $gc;
    $html_ds .= qq{<tr><td>Tools:</td>}
          . "<td>"
          . "<a href='OrganismView.pl?dsid=$dsid' target=_new>OrganismView</a>"
          . qq{</td></tr>};

    $html_ds .= qq{</table></td>}
          . qq{<td id="ds_features"></td>}
          . qq{</table>};
    my $chr_count = $chr_num;
    $chr_count .= " <span class='note'> (only $MAX_NUM_CHROMOSOME_RESULTS largest listed)</span>"
      if ( $chr_count > $MAX_NUM_CHROMOSOME_RESULTS );
    
#    print STDERR Dumper {
#            dataset => $html_ds,
#            chromosomes => $html_chr,
#            count => $chr_count,
#            selected_chr => $selected_chr
#        }, "\n";
    return $html_ds, $html_chr, $chr_count, $selected_chr if ($output eq 'html');
    return encode_json({
        dataset => $html_ds,
        chromosomes => $html_chr,
        count => scalar(keys %chr),
        selected_chr => $selected_chr
    });
}

sub get_chr_info {
    my %opts = @_;
    my $dsid  = $opts{dsid};
    my $chr   = $opts{chr};
    my $gid   = $opts{gid};
    my $output  = $opts{output} || 'json';
#    print STDERR "get_chr_info $dsid ",($chr ? $chr : ''),"\n";

    $dsid  = 0 unless $dsid;
    unless ( $dsid && defined $chr && $chr ne '' ) # error flag for empty dataset
    {
        return '', '', '' if ($output eq 'html');
        return encode_json({
            chr_info => '',
            viewer => '',
            seqview => ''
        });
    }
    
    # Get dataset from db
    my $ds = $coge->resultset("Dataset")->find($dsid);
    unless ($ds) {
        return '', '', '' if ($output eq 'html');
        return encode_json({
            chr_info => '',
            viewer => '',
            seqview => ''
        });
    }

    my $start = "'start'";
    my $stop  = "'stop'";
    my $html = "<table class=\"small annotation_table\">";
    
    my $length = 0;
    $length = $ds->last_chromosome_position($chr) if defined $chr;
    my $gc =
      $length < $MAX_DS_LENGTH
      ? get_gc_for_chromosome( dsid => $ds->id, chr => $chr )
      : 0;
    $length = commify($length) . " bp ";
    $html .= qq{<tr><td>Chromosome ID:</td><td>$chr</td></tr>}
          .  qq{<tr><td>Nucleotides:</td><td>$length</td>}
          #.  qq{<td>$gc</td>}
          .  qq{</tr>};
    $html .= qq{</table>};
    
    # Generate html for launch & seq retrieval buttons
    $gid = 0 unless defined $gid;
    my ($viewer, $seq_grab);
    if ( defined $chr ) {
        $viewer .=
           qq{<span class='coge-button coge-button-sm' style="width:11em;" onClick="launch_viewer()">Launch Genome Viewer</span>}
         . qq{&nbsp;<span class="small text">Start: <input type="small text" size=10 value="20000" id="x"></span>};
         #. qq{Zoom level: <input type = "text" size=10 value ="6" id = "z">};

        $seq_grab .=
           qq{<span class='coge-button coge-button-sm' style="width:11em;" onClick="launch_seqview()">Get Sequence</span>}
         . qq{&nbsp;<span class="small text">Start: <input type="text" size=10 value="1" id="start"></span>}
         . qq{&nbsp;<span class="small text">End: <input type="text" size=10 value="100000" id="stop"></span>};
    }

    return $html, $viewer, $seq_grab if ($output eq 'html');
    return encode_json({
        chr_info => $html,
        viewer => $viewer,
        seqview => $seq_grab
    });
}

sub gen_data {
    my $message = shift;
    return qq{<font class="small alert">$message</font>};
}

sub get_gc_for_feature_type {
    my %opts   = @_;
    my $dsid   = $opts{dsid};
    my $dsgid  = $opts{dsgid};
    my $chr    = $opts{chr};
    my $typeid = $opts{typeid};
    my $gstid  = $opts{gstid};  #genomic sequence type id
    my $min    = $opts{min};    #limit results with gc values greater than $min;
    my $max    = $opts{max};    #limit results with gc values smaller than $max;
    my $hist_type = $opts{hist_type};
    $hist_type = "counts" unless $hist_type;
    $min       = undef if $min       && $min       eq "undefined";
    $max       = undef if $max       && $max       eq "undefined";
    $chr       = undef if $chr       && $chr       eq "undefined";
    $dsid      = undef if $dsid      && $dsid      eq "undefined";
    $hist_type = undef if $hist_type && $hist_type eq "undefined";
    $typeid = 1 if $typeid eq "undefined";
    return unless $dsid || $dsgid;
    my $gc   = 0;
    my $at   = 0;
    my $n    = 0;
    my $type = $coge->resultset('FeatureType')->find($typeid);
    my @data;
    my @fids;    #storage for fids that passed.  To be sent to FeatList

	my (@items, @datasets);
	if ($dsid) {
		my $ds = $coge->resultset('Dataset')->find($dsid);
		return "unable to find dataset id$dsid\n" unless $ds;
		push @items, $ds;
		push @datasets, $ds;
	}
    if ($dsgid) {
        my $dsg = $coge->resultset('Genome')->find($dsgid);
        return "unable to find genome id$dsgid\n" unless $dsgid;
        $gstid = $dsg->type->id;
        push @items, $dsg;
        push @datasets, $dsg->datasets;
    }

    my %seqs; # prefetch the sequences with one call to genomic_sequence (slow for many seqs)
    foreach my $item (@items) {
        map {
        	$seqs{$_} = $item->get_genomic_sequence( chr => $_, seq_type => $gstid )
        } (defined $chr ? ($chr) : $item->chromosomes);
    }

    my $search = { "feature_type_id" => $typeid };
    $search->{"me.chromosome"} = $chr if defined $chr;

    foreach my $ds (@datasets) {
        my @feats = $ds->features(
            $search,
            {
                join => [
                    'locations',
                    { 'dataset' => { 'dataset_connectors' => 'genome' } }
                ],
                prefetch => [
                    'locations',
                    { 'dataset' => { 'dataset_connectors' => 'genome' } }
                ],
            }
        );
        foreach my $feat (@feats) {
            my $seq = substr(
                $seqs{ $feat->chromosome },
                $feat->start - 1,
                $feat->stop - $feat->start + 1
            );

            $feat->genomic_sequence( seq => $seq );
            my @gc = $feat->gc_content( counts => 1 );

            $gc += $gc[0] if $gc[0] =~ /^\d+$/;
            $at += $gc[1] if $gc[1] =~ /^\d+$/;
            $n  += $gc[2] if $gc[2] =~ /^\d+$/;
            my $total = 0;
            $total += $gc[0] if $gc[0];
            $total += $gc[1] if $gc[1];
            $total += $gc[2] if $gc[2];
            my $perc_gc = 100 * $gc[0] / $total if $total;
            next unless $perc_gc;    #skip if no values
            next
              if defined $min
                  && $min =~ /\d+/
                  && $perc_gc < $min;    #check for limits
            next
              if defined $max
                  && $max =~ /\d+/
                  && $perc_gc > $max;    #check for limits
            push @data, sprintf( "%.2f", $perc_gc );
            push @fids, $feat->id . "_" . $gstid;
        }
    }
    my $total = $gc + $at + $n;
    return "error" unless $total;

	my @dsids = map { $_->id } @datasets;
    my $file = $TEMPDIR . "/" . join( "_", @dsids );

    #perl -T flag
    ($min) = $min =~ /(.*)/ if defined $min;
    ($max) = $max =~ /(.*)/ if defined $max;
    ($chr) = $chr =~ /(.*)/ if defined $chr;
    $file .= "_" . $chr . "_" if defined $chr;
    $file .= "_min" . $min    if defined $min;
    $file .= "_max" . $max    if defined $max;
    $file .= "_$hist_type"    if $hist_type;
    $file .= "_" . $type->name . "_gc.txt";
    open( OUT, ">" . $file );
    print OUT "#wobble gc for dataset ids: " . join( " ", @dsids ), "\n";
    print OUT join( "\n", @data ), "\n";
    close OUT;
    my $cmd = $HISTOGRAM;
    $cmd .= " -f $file";
    my $out = $file;
    $out =~ s/txt$/png/;
    $cmd .= " -o $out";
    $cmd .= " -t \"" . $type->name . " gc content\"";
    $cmd .= " -min 0";
    $cmd .= " -max 100";
    $cmd .= " -ht $hist_type" if $hist_type;
    `$cmd`;

    $min = 0   unless defined $min && $min =~ /\d+/;
    $max = 100 unless defined $max && $max =~ /\d+/;
    my $info;
    $info .= qq{<div class="small">
Min: <input type="text" size="3" id="feat_gc_min" value="$min">
Max: <input type="text" size="3" id="feat_gc_max" value="$max">
Type: <select id="feat_hist_type">
<option value="counts">Counts</option>
<option value="percentage">Percentage</option>
</select>
};
    $info =~ s/>Per/ selected>Per/ if $hist_type =~ /per/;
    my $gc_args;
    $gc_args = "chr: '$chr'," if defined $chr;
    $gc_args .= "dsid: $dsid," if $dsid; #set a var so that histograms are only calculated for the dataset and not hte genome
    $gc_args .= "typeid: '$typeid'";
    $info .= qq{<span class="link" onclick="get_feat_gc({$gc_args})">Regenerate histogram</span>};
    $info .= "</div>";
    $info .=
        "<div class = small>Total length: "
      . commify($total)
      . " bp, GC: "
      . sprintf( "%.2f", 100 * $gc / ($total) )
      . "%  AT: "
      . sprintf( "%.2f", 100 * $at / ($total) )
      . "%  N: "
      . sprintf( "%.2f", 100 * ($n) / ($total) )
      . "%</div>";

    if ( $min || $max ) {
        $min = 0   unless defined $min;
        $max = 100 unless defined $max;
        $info .= qq{<div class=small style="color: red;">Limits set:  MIN: $min  MAX: $max</div>};
    }
    my $stuff = join "::", @fids;
    $info .= qq{<div class="link small" onclick="window.open('FeatList.pl?fid=$stuff')">Open FeatList of Features</div>};

    $out =~ s/$TEMPDIR/$TEMPURL/;
    $info .= "<br><img src=\"$out\">";
    return $info;
}

sub get_gc_for_chromosome {
    my %opts  = @_;
    my $dsid  = $opts{dsid};
    my $chr   = $opts{chr};
    my $gstid = $opts{gstid};
    my $dsgid = $opts{dsgid};
    my @ds;
    if ($dsid) {
        my $ds = $coge->resultset('Dataset')->find($dsid);
        push @ds, $ds if $ds;
    }
    if ($dsgid) {
        my $dsg = $coge->resultset('Genome')->find($dsgid);
        $gstid = $dsg->type->id;
        map { push @ds, $_ } $dsg->datasets;
    }
    return unless @ds;
    my ( $gc, $at, $n, $x ) = ( 0, 0, 0, 0 );
    my %chr;
    foreach my $ds (@ds) {
        if ( defined $chr ) {
            $chr{$chr} = 1;
        }
        else {
            map { $chr{$_} = 1 } $ds->chromosomes;
        }
        foreach my $chr ( keys %chr ) {
            my @gc =
              $ds->percent_gc( chr => $chr, seq_type => $gstid, count => 1 );
            $gc += $gc[0] if $gc[0];
            $at += $gc[1] if $gc[1];
            $n  += $gc[2] if $gc[2];
            $x  += $gc[3] if $gc[3];
        }
    }
    my $total = $gc + $at + $n + $x;
    return "error" unless $total;
    my $results =
        "&nbsp(GC: "
      . sprintf( "%.2f", 100 * $gc / $total )
      . "%  AT: "
      . sprintf( "%.2f", 100 * $at / $total )
      . "%  N: "
      . sprintf( "%.2f", 100 * $n / $total )
      . "%  X: "
      . sprintf( "%.2f", 100 * $x / $total ) . "%)"
      if $total;
    return $results;
}

sub get_gc_for_noncoding {
    my %opts  = @_;
    my $dsid  = $opts{dsid};
    my $dsgid = $opts{dsgid};
    my $chr   = $opts{chr};
    my $gstid = $opts{gstid};    #genomic sequence type id
    return "error" unless $dsid || $dsgid;
    my $gc = 0;
    my $at = 0;
    my $n  = 0;
    my $x  = 0;
    my $search = { "feature_type_id" => 3 };
    $search->{"me.chromosome"} = $chr if defined $chr;
    my @data;

	my (@items, @datasets);
	if ($dsid) {
		my $ds = $coge->resultset('Dataset')->find($dsid);
		return "unable to find dataset id$dsid\n" unless $ds;
		push @items, $ds;
		push @datasets, $ds;
	}
    if ($dsgid) {
        my $dsg = $coge->resultset('Genome')->find($dsgid);
        return "unable to find genome id$dsgid\n" unless $dsgid;
        $gstid = $dsg->type->id;
        push @items, $dsg;
        push @datasets, $dsg->datasets;
    }

    my %seqs; # prefetch the sequences with one call to genomic_sequence (slow for many seqs)
    foreach my $item (@items) {
        map {
        	$seqs{$_} = $item->get_genomic_sequence( chr => $_, seq_type => $gstid )
        } (defined $chr ? ($chr) : $item->chromosomes);
    }

    foreach my $ds (@datasets) {
        foreach my $feat (
            $ds->features(
                $search,
                {
                    join => [
                        'locations',
                        { 'dataset' => { 'dataset_connectors' => 'genome' } }
                    ],
                    prefetch => [
                        'locations',
                        { 'dataset' => { 'dataset_connectors' => 'genome' } }
                    ]
                }
            )
          )
        {
            foreach my $loc ( $feat->locs ) {
                if ( $loc->stop > length( $seqs{ $feat->chromosome } ) ) {
                    print STDERR "feature "
                      . $feat->id
                      . " stop exceeds sequence length: "
                      . $loc->stop . " :: "
                      . length( $seqs{ $feat->chromosome } ), "\n";
                }
                substr(
                    $seqs{ $feat->chromosome },
                    $loc->start - 1,
                    ( $loc->stop - $loc->start + 1 )
                ) = "-" x ( $loc->stop - $loc->start + 1 );
            }

            #push @data, sprintf("%.2f",100*$gc[0]/$total) if $total;
        }
    }
    foreach my $seq ( values %seqs ) {
        $gc += $seq =~ tr/GCgc/GCgc/;
        $at += $seq =~ tr/ATat/ATat/;
        $n  += $seq =~ tr/nN/nN/;
        $x  += $seq =~ tr/xX/xX/;
    }
    my $total = $gc + $at + $n + $x;
    return "error" unless $total;
    return
        commify($total) . " bp"
      . "&nbsp(GC: "
      . sprintf( "%.2f", 100 * $gc / ($total) )
      . "%  AT: "
      . sprintf( "%.2f", 100 * $at / ($total) ) . "% N: "
      . sprintf( "%.2f", 100 * $n /  ($total) ) . "% X: "
      . sprintf( "%.2f", 100 * $x /  ($total) ) . "%)";

	my @dsids = map { $_->id } @datasets;
    my $file = $TEMPDIR . "/" . join( "_", @dsids ) . "_wobble_gc.txt";
    open( OUT, ">" . $file );
    print OUT "#wobble gc for dataset ids: " . join( " ", @dsids ), "\n";
    print OUT join( "\n", @data ), "\n";
    close OUT;
    my $cmd = $HISTOGRAM;
    $cmd .= " -f $file";
    my $out = $file;
    $out =~ s/txt$/png/;
    $cmd .= " -o $out";
    $cmd .= " -t \"CDS wobble gc content\"";
    $cmd .= " -min 0";
    $cmd .= " -max 100";
    `$cmd`;
    my $info =
        "<div class = small>Total: "
      . commify($total)
      . " codons.  Mean GC: "
      . sprintf( "%.2f", 100 * $gc / ($total) )
      . "%  AT: "
      . sprintf( "%.2f", 100 * $at / ($total) )
      . "%  N: "
      . sprintf( "%.2f", 100 * ($n) / ($total) )
      . "%</div>";
    $out =~ s/$TEMPDIR/$TEMPURL/;
    my $hist_img = "<img src=\"$out\">";

    return $info, $hist_img;
}

sub get_aa_usage {
    my %opts  = @_;
    my $dsid  = $opts{dsid};
    my $chr   = $opts{chr};
    my $dsgid = $opts{dsgid};
    my $gstid = $opts{gstid};
    return unless $dsid || $dsgid;

    my $search;
    $search = { "feature_type.name" => "CDS" };
    $search->{"me.chromosome"} = $chr if defined $chr;

	my (@items, @datasets);
	if ($dsid) {
		my $ds = $coge->resultset('Dataset')->find($dsid);
		return "unable to find dataset id$dsid\n" unless $ds;
		push @items, $ds;
		push @datasets, $ds;
	}
    if ($dsgid) {
        my $dsg = $coge->resultset('Genome')->find($dsgid);
        return "unable to find genome id$dsgid\n" unless $dsgid;
        $gstid = $dsg->type->id;
        push @items, $dsg;
        push @datasets, $dsg->datasets;
    }

    my %seqs; # prefetch the sequences with one call to genomic_sequence (slow for many seqs)
    foreach my $item (@items) {
        map {
        	$seqs{$_} = $item->get_genomic_sequence( chr => $_, seq_type => $gstid )
        } (defined $chr ? ($chr) : $item->chromosomes);
    }

    my %codons;
    my $codon_total = 0;
    my %aa;
    my $aa_total   = 0;
    my $feat_count = 0;
    my ( $code, $code_type );

    foreach my $ds (@datasets) {
        foreach my $feat (
            $ds->features(
                $search,
                {
                    join => [
                        "feature_type", 'locations',
                        { 'dataset' => { 'dataset_connectors' => 'genome' } }
                    ],
                    prefetch => [
                        'locations',
                        { 'dataset' => { 'dataset_connectors' => 'genome' } }
                    ]
                }
            )
          )
        {
            my $seq = substr(
                $seqs{ $feat->chromosome },
                $feat->start - 1,
                $feat->stop - $feat->start + 1
            );
            $feat->genomic_sequence( seq => $seq );
            $feat_count++;
            ( $code, $code_type ) = $feat->genetic_code() unless $code;
            my ($codon) = $feat->codon_frequency( counts => 1 );
            grep { $codon_total += $_ } values %$codon;
            grep { $codons{$_}  += $codon->{$_} } keys %$codon;
            foreach my $tri ( keys %$code ) {
                $aa{ $code->{$tri} } += $codon->{$tri};
                $aa_total += $codon->{$tri};
            }
            print STDERR ".($feat_count)" if !$feat_count % 10;
        }
    }
    %codons = map { $_, $codons{$_} / $codon_total } keys %codons;

    %aa = map { $_, $aa{$_} / $aa_total } keys %aa;

#    my $html1 = "Codon Usage: $code_type";
#    $html1 .= CoGe::Accessory::genetic_code->html_code_table(data=>\%codons, code=>$code);
    my $html2 .= "<div class='coge-table-header'>Predicted Amino Acid Usage: $code_type</div>"
        . CoGe::Accessory::genetic_code->html_aa_new( data => \%aa )
        . '<br>';
    return $html2; #return $html1, $html2;
}

sub get_wobble_gc {
    my %opts  = @_;
    my $dsid  = $opts{dsid};
    my $dsgid = $opts{dsgid};
    my $chr   = $opts{chr};
    my $gstid = $opts{gstid};   #genomic sequence type id
    my $min   = $opts{min};     #limit results with gc values greater than $min;
    my $max   = $opts{max};     #limit results with gc values smaller than $max;
    my $hist_type = $opts{hist_type};
    return "error" unless $dsid || $dsgid;
    my $gc = 0;
    my $at = 0;
    my $n  = 0;
    my $search;
    $search = { "feature_type_id" => 3 };
    $search->{"me.chromosome"} = $chr if defined $chr;
    my @data;
    my @fids;
    my @dsids;
    push @dsids, $dsid if $dsid;

    if ($dsgid) {
        my $dsg = $coge->resultset('Genome')->find($dsgid);
        unless ($dsg) {
            my $error = "unable to create genome object using id $dsgid\n";
            return $error;
        }
        $gstid = $dsg->type->id;
        foreach my $ds ( $dsg->datasets() ) {
            push @dsids, $ds->id;
        }
    }
    foreach my $dsidt (@dsids) {
        my $ds = $coge->resultset('Dataset')->find($dsidt);
        unless ($ds) {
            warn "no dataset object found for id $dsidt\n";
            next;
        }
        foreach my $feat (
            $ds->features(
                $search,
                {
                    join => [
                        'locations',
                        { 'dataset' => { 'dataset_connectors' => 'genome' } }
                    ],
                    prefetch => [
                        'locations',
                        { 'dataset' => { 'dataset_connectors' => 'genome' } }
                    ],
                }
            )
          )
        {
            my @gc = $feat->wobble_content( counts => 1 );
            $gc += $gc[0] if $gc[0] && $gc[0] =~ /^\d+$/;
            $at += $gc[1] if $gc[1] && $gc[1] =~ /^\d+$/;
            $n  += $gc[2] if $gc[2] && $gc[2] =~ /^\d+$/;
            my $total = 0;
            $total += $gc[0] if $gc[0];
            $total += $gc[1] if $gc[1];
            $total += $gc[2] if $gc[2];
            my $perc_gc = 100 * $gc[0] / $total if $total;
            next unless $perc_gc;    #skip if no values
            next
              if defined $min
                  && $min =~ /\d+/
                  && $perc_gc < $min;    #check for limits
            next
              if defined $max
                  && $max =~ /\d+/
                  && $perc_gc > $max;    #check for limits
            push @data, sprintf( "%.2f", $perc_gc );
            push @fids, $feat->id . "_" . $gstid;

            #push @data, sprintf("%.2f",100*$gc[0]/$total) if $total;
        }
    }
    my $total = $gc + $at + $n;
    return "error" unless $total;

    my $file = $TEMPDIR . "/" . join( "_", @dsids );    #."_wobble_gc.txt";
    ($min) = $min =~ /(.*)/ if defined $min;
    ($max) = $max =~ /(.*)/ if defined $max;
    ($chr) = $chr =~ /(.*)/ if defined $chr;
    $file .= "_" . $chr . "_" if defined $chr;
    $file .= "_min" . $min    if defined $min;
    $file .= "_max" . $max    if defined $max;
    $file .= "_$hist_type"    if $hist_type;
    $file .= "_wobble_gc.txt";
    open( OUT, ">" . $file );
    print OUT "#wobble gc for dataset ids: " . join( " ", @dsids ), "\n";
    print OUT join( "\n", @data ), "\n";
    close OUT;
    my $cmd = $HISTOGRAM;
    $cmd .= " -f $file";
    my $out = $file;
    $out =~ s/txt$/png/;
    $cmd .= " -o $out"
         . " -t \"CDS wobble gc content\""
         . " -min 0"
         . " -max 100";
    $cmd .= " -ht $hist_type" if $hist_type;
    `$cmd`;
    $min = 0   unless defined $min && $min =~ /\d+/;
    $max = 100 unless defined $max && $max =~ /\d+/;
    my $info;
    $info .= qq{<div class="coge-table-header">Wobble Codon GC Content</div><div class="small">
Min: <input type="text" size="3" id="wobble_gc_min" value="$min">
Max: <input type="text" size="3" id="wobble_gc_max" value="$max">
Type: <select id="wobble_hist_type">
<option value="counts">Counts</option>
<option value="percentage">Percentage</option>
</select>
};
    $info =~ s/>Per/ selected>Per/ if $hist_type =~ /per/;
    my $args;
    $args .= "'args__dsid','ds_id',"   if $dsid;
    $args .= "'args__dsgid','dsg_id'," if $dsgid;
    $args .= "'args__chr','chr',"      if defined $chr;
    $args .= "'args__min','wobble_gc_min',";
    $args .= "'args__max','wobble_gc_max',";
    $args .= "'args__max','wobble_gc_max',";
    $args .= "'args__hist_type', 'wobble_hist_type',";
    $info .= qq{<span class="link" onclick="get_wobble_gc([$args],['wobble_gc_histogram']).html('loading...');">Regenerate histogram</span>};
    $info .= "</div>";

    $info .=
        "<div class = small>Total: "
      . commify($total)
      . " codons.  Mean GC: "
      . sprintf( "%.2f", 100 * $gc / ($total) )
      . "%  AT: "
      . sprintf( "%.2f", 100 * $at / ($total) )
      . "%  N: "
      . sprintf( "%.2f", 100 * ($n) / ($total) )
      . "%</div>";

    if ( $min || $max ) {
        $min = 0   unless defined $min;
        $max = 100 unless defined $max;
        $info .= qq{<div class=small style="color: red;">Limits set:  MIN: $min  MAX: $max</div>};
    }
    my $stuff = join "::", @fids;
    $info .= qq{<div class="link small" onclick="window.open('FeatList.pl?fid=$stuff')">Open FeatList of Features</div>};
    $out =~ s/$TEMPDIR/$TEMPURL/;
    my $hist_img = "<img src=\"$out\">";
    return $info . "<br>" . $hist_img . '<br>';
}

sub get_wobble_gc_diff {
    my %opts  = @_;
    my $dsid  = $opts{dsid};
    my $dsgid = $opts{dsgid};
    my $chr   = $opts{chr};
    my $gstid = $opts{gstid};    #genomic sequence type id
    return "error", " " unless $dsid || $dsgid;
    my $search;
    $search = { "feature_type_id" => 3 };
    $search->{"me.chromosome"} = $chr if defined $chr;
    my @data;
    my @dsids;
    push @dsids, $dsid if $dsid;

    if ($dsgid) {
        my $dsg = $coge->resultset('Genome')->find($dsgid);
        unless ($dsg) {
            my $error = "unable to create genome object using id $dsgid\n";
            return $error;
        }
        $gstid = $dsg->type->id;
        foreach my $ds ( $dsg->datasets() ) {
            push @dsids, $ds->id;
        }
    }
    foreach my $dsidt (@dsids) {
        my $ds = $coge->resultset('Dataset')->find($dsidt);
        unless ($ds) {
            warn "no dataset object found for id $dsidt\n";
            next;
        }
        foreach my $feat (
            $ds->features(
                $search,
                {
                    join => [
                        'locations',
                        { 'dataset' => { 'dataset_connectors' => 'genome' } }
                    ],
                    prefetch => [
                        'locations',
                        { 'dataset' => { 'dataset_connectors' => 'genome' } }
                    ]
                }
            )
          )
        {
            my @wgc  = $feat->wobble_content();
            my @gc   = $feat->gc_content();
            my $diff = $gc[0] - $wgc[0] if defined $gc[0] && defined $wgc[0];
            push @data, sprintf( "%.2f", 100 * $diff ) if $diff;
        }
    }
    return "error", " " unless @data;
    my $file = $TEMPDIR . "/" . join( "_", @dsids ) . "_wobble_gc_diff.txt";
    open( OUT, ">" . $file );
    print OUT "#wobble gc for dataset ids: " . join( " ", @dsids ), "\n";
    print OUT join( "\n", @data ), "\n";
    close OUT;
    my $cmd = $HISTOGRAM;
    $cmd .= " -f $file";
    my $out = $file;
    $out =~ s/txt$/png/;
    $cmd .= " -o $out"
         . " -t \"CDS GC - wobble gc content\"";
    `$cmd`;
    my $sum = 0;
    map { $sum += $_ } @data;
    my $mean = sprintf( "%.2f", $sum / scalar @data );
    my $info = "<div class='coge-table-header'>CDS GC vs. Wobble Codon GC</div>"
        . "<div class='small'>Mean $mean%"
        . " (";
    $info .= $mean > 0 ? "CDS" : "wobble";
    $info .= " is more GC rich)</div>";
    $out =~ s/$TEMPDIR/$TEMPURL/;
    my $hist_img = "<img src=\"$out\">";
    return $info . "<div class='padded'>" . $hist_img . "</div>";
}

sub get_chr_length_hist {
    my %opts  = @_;
    my $dsgid = $opts{dsgid};
    return "error", " " unless $dsgid;
    my ($dsg) = $coge->resultset('Genome')->find($dsgid);
    unless ($dsg) {
        my $error = "unable to create genome object using id $dsgid\n";
        return $error;
    }
    my @data = $dsg->chromosome_lengths;
    return "error", " " unless @data;
    @data = sort { $a <=> $b } @data;
    my $mid  = floor( scalar(@data) / 2 );
    my $mode = $data[$mid];
    my $file = $TEMPDIR . "/" . join( "_", $dsgid ) . "_chr_length.txt";
    open( OUT, ">" . $file );
    print OUT "#chromosome/contig lenghts for $dsgid\n";
    print OUT join( "\n", @data ), "\n";
    close OUT;
    my $cmd = $HISTOGRAM;
    $cmd .= " -f $file";
    my $out = $file;
    $out =~ s/txt$/png/;
    $cmd .= " -o $out";
    $cmd .=
        " -t \"Chromosome length for "
      . $dsg->organism->name . " (v"
      . $dsg->version . ")\"";
    `$cmd`;
    my $sum = 0;
    map { $sum += $_ } @data;
    my $n50;

    foreach my $val (@data) {
        $n50 += $val;
        if ( $n50 >= $sum / 2 ) {
            $n50 = $val;
            last;
        }
    }
    my $mean = sprintf( "%.0f", $sum / scalar @data );
    my $info =
        "<table class=small><TR><td>Count:<TD>"
      . commify( $dsg->chromosome_count )
      . " chromosomes (contigs, scaffolds, etc.)<tr><Td>Mean:<td>"
      . commify($mean)
      . " nt<tr><Td>Mode:<td>"
      . commify($mode)
      . " nt<tr><td>N50:<td>"
      . commify($n50)
      . " nt</table>";
    $out =~ s/$TEMPDIR/$TEMPURL/;
    my $hist_img = "<img src=\"$out\">";
    return $info . "<div class='padded'>" . $hist_img . "</div>";
}

sub get_total_length_for_ds {
    my %opts   = @_;
    my $dsid   = $opts{dsid};
    my $ds     = $coge->resultset('Dataset')->find($dsid);
    my $length = 0;
    map { $length += $ds->last_chromosome_position($_) } $ds->get_chromosomes();
    return commify($length);
}

1;
