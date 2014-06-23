#! /usr/bin/perl -w
use strict;
use CGI;
use CGI::Carp 'fatalsToBrowser';
use CoGeX;
use CoGe::Accessory::Web;
use CoGe::Accessory::Utils qw( commify );
use HTML::Template;
use Data::Dumper;
use CGI::Ajax;
use Benchmark;
use File::Path;
use Benchmark qw(:all);
use Statistics::Basic::Mean;
use POSIX;
use Sort::Versions;

no warnings 'redefine';

use vars qw($P $PAGE_NAME $PAGE_TITLE $LINK
  $TEMPDIR $TEMPURL $USER $FORM $coge $HISTOGRAM
  %FUNCTION $P $SERVER $MAX_NUM_ORGANISM_RESULTS);

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

$MAX_NUM_ORGANISM_RESULTS = 5000;

%FUNCTION = (
    get_genomes             => \&get_genomes,
    get_genome_info         => \&get_genome_info,
    get_dataset             => \&get_dataset,
    get_dataset_info        => \&get_dataset_info,
    get_dataset_chr_info    => \&get_dataset_chr_info,
    gen_data                => \&gen_data,
    get_orgs                => \&get_orgs,
    get_org_info            => \&get_org_info,
    get_recent_orgs         => \&get_recent_orgs,
    get_start_stop          => \&get_start_stop,
    get_feature_counts      => \&get_feature_counts,
    get_gc_for_chromosome   => \&get_gc_for_chromosome,
    get_gc_for_noncoding    => \&get_gc_for_noncoding,
    get_gc_for_feature_type => \&get_gc_for_feature_type,
    get_codon_usage         => \&get_codon_usage,
    get_aa_usage            => \&get_aa_usage,
    get_wobble_gc           => \&get_wobble_gc,
    get_wobble_gc_diff      => \&get_wobble_gc_diff,
    get_total_length_for_ds => \&get_total_length_for_ds,
    get_chr_length_hist     => \&get_chr_length_hist,
    update_genomelist       => \&update_genomelist,
    parse_for_GenoList      => \&parse_for_GenoList,
    get_genome_list_for_org => \&get_genome_list_for_org,
    add_to_irods            => \&add_to_irods,
    make_genome_public      => \&make_genome_public,
    make_genome_private     => \&make_genome_private,
    edit_genome_info        => \&edit_genome_info,
    update_genome_info      => \&update_genome_info,
);
my $pj = new CGI::Ajax(%FUNCTION);
$pj->JSDEBUG(0);
$pj->DEBUG(0);

if ( $FORM->param('jquery_ajax') ) {
    dispatch();
}
else {
    print $pj->build_html( $FORM, \&gen_html );
}

#print "Content-Type: text/html\n\n";print gen_html($FORM);

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
    #	print $FORM->header, gen_html();
    #    }
}

sub parse_for_GenoList {
    my $genomelist = shift;
    my $url        = "GenomeList.pl?dsgid=$genomelist";
    return $url;
}

sub gen_html {
    my $html;
    my ( $body, $seq_names, $seqs ) = gen_body();
    my $template =
      HTML::Template->new( filename => $P->{TMPLDIR} . 'generic_page.tmpl' );
    $template->param( PAGE_TITLE => 'OrganismView',
    				  PAGE_LINK  => $LINK,
    				  HEAD       => qq{},
    				  HELP       => "/wiki/index.php?title=OrganismView" );
    my $name = $USER->user_name;
    $name = $USER->first_name if $USER->first_name;
    $name .= " " . $USER->last_name if $USER->first_name && $USER->last_name;
    $template->param( USER     => $name );
    $template->param( LOGON    => 1 ) unless $USER->user_name eq "public";
    $template->param( LOGO_PNG => "OrganismView-logo.png" );
    $template->param( BODY     => $body );

    $html .= $template->output;
    return $html;
}

sub gen_body {
    my $form = shift || $FORM;
    my $template =
      HTML::Template->new( filename => $P->{TMPLDIR} . 'OrganismView.tmpl' );
    my $org_name = $form->param('org_name');
    my $desc     = $form->param('org_desc');
    my $oid      = $form->param('oid');
    my $org      = $coge->resultset('Organism')->resolve($oid) if $oid;
    my $dsname   = $form->param('dsname');
    my $dsid     = $form->param('dsid');
    my $dsgid    = $form->param('dsgid');
    $dsgid = $form->param('gid') if ( $form->param('gid') );
    my ($dsg) = $coge->resultset('Genome')->find($dsgid) if $dsgid;
    $template->param( SHOW_RESULTS => 1 ) if ($dsgid or $dsid or $desc or $oid or $org or $dsname);

    my $link = "http://" . $P->{SERVER} . $PAGE_NAME . "?";
    my @params;
    push @params, "org_name=" . $org_name if $org_name;
    push @params, "org_desc=" . $desc     if $desc;
    push @params, "oid=" . $oid           if $oid;
    push @params, "dsname=" . $dsname     if $dsname;
    push @params, "dsid=" . $dsid         if $dsid;
    push @params, "dsgid=" . $dsgid       if $dsgid;
    $link .= join ";", @params;
    $link = CoGe::Accessory::Web::get_tiny_link(
        db      => $coge,
        user_id => $USER->id,
        page    => $PAGE_NAME,
        url     => $link
    ) if @params;

    $template->param( SERVER => $SERVER );
    $org      = $dsg->organism if $dsg;
    $org_name = $org->name     if $org;
    $org_name .= $desc if $desc; # mdb added 4/22/14 combined org name & desc searches into single input
    $template->param( ORG_SEARCH => $org_name ) if $org_name;

    my ( $org_list, $org_count ) =
      get_orgs( name => $org_name, oid => $oid, dsgid => $dsgid );
    $template->param( ORG_LIST  => $org_list );
    $template->param( ORG_COUNT => $org_count );

    #$template->param(RECENT=>get_recent_orgs());
    my ($ds) = $coge->resultset('Dataset')->resolve($dsid) if $dsid;
    $dsname = $ds->name if $ds;
    $dsname = "Search" unless $dsname;
    $template->param( DS_NAME => $dsname );
    $dsname = "" if $dsname =~ /Search/;
    my ( $dslist, $dscount ) = get_dataset( dsname => $dsname, dsid => $dsid )
      if $dsname;
    $template->param( DS_LIST  => $dslist )  if $dslist;
    $template->param( DS_COUNT => $dscount ) if $dscount;
    my $dsginfo = "<input type='hidden' id='gstid'>";
    $dsginfo .=
      $dsgid
      ? "<input type='hidden' id='dsg_id' value='$dsgid'>"
      : "<input type='hidden' id='dsg_id'>";
    $template->param( DSG_INFO => $dsginfo );
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

sub update_genome_info {
    my %opts  = @_;
    my $dsgid = $opts{dsgid};
    return "No DSGID specified" unless $dsgid;
    return "Permission denied."
      unless $USER->is_admin || $USER->is_owner( dsg => $dsgid );
    my $dsg     = $coge->resultset('Genome')->find($dsgid);
    my $name    = $opts{name};
    my $desc    = $opts{desc};
    my $ver     = $opts{ver};
    my $message = $opts{message};
    my $link    = $opts{link};
    $dsg->name($name);
    $dsg->description($desc);
    $dsg->version($ver);
    $dsg->message($message);
    $dsg->link($link);
    $dsg->update;
    return 1;
}

sub edit_genome_info {
    my %opts  = @_;
    my $dsgid = $opts{dsgid};
    return "No DSGID specified" unless $dsgid;
    return "Permission denied."
      unless $USER->is_admin || $USER->is_owner( dsg => $dsgid );
    my $dsg     = $coge->resultset('Genome')->find($dsgid);
    my $name    = $dsg->name;
    my $desc    = $dsg->description;
    my $ver     = $dsg->version;
    my $message = $dsg->message;
    my $link    = $dsg->link;
    my $html = qq{

<table class=small>
 <tr>
  <td>Name:</td>
  <td><input type=text id=dsg_name size=50 value="$name"></td>
 </tr>
 <tr>
  <td>Description:</td>
  <td><input type=text id=dsg_desc size=50 value="$desc"></td>
 </tr>
 <tr>
  <td>Version:</td>
  <td><input type=text id=dsg_ver size=5 value="$ver"></td>
 </tr>
 <tr>
  <td>Message:</td>
  <td><textarea id=dsg_message cols=50 rows=5>$message</textarea></td>
 </tr>
 <tr>
  <td>Link:</td>
  <td><input type=text id=dsg_link size=50 value="$link"></td>
 </tr>
</table>
<span class="ui-button ui-corner-all ui-button-go" onClick="update_genome_info('$dsgid')">Update</span>
};
    return $html;
}

sub get_recent_orgs {
    my %opts  = @_;
    my $limit = $opts{limit} || 100;
    my @db    = $coge->resultset("Dataset")->search(
        { restricted => 0 },
        {
            distinct => "organism.name",
            join     => "organism",
            order_by => "me.date desc",
            rows     => $limit
        }
    );

    my $i = 0;
    my @opts;
    my %org_names;
    foreach my $item (@db) {
        $i++;
        my $date = $item->date;
        $date =~ s/\s.*//;

        #next if $USER->user_name =~ /public/i && $item->organism->restricted;
        next if $org_names{ $item->organism->name };
        $org_names{ $item->organism->name } = 1;
        push @opts,
            "<OPTION value=\""
          . $item->organism->id . "\">"
          . $date . " "
          . $item->organism->name . " (id"
          . $item->organism->id . ") "
          . "</OPTION>";
    }
    
    my $html;
    #$html .= qq{<FONT CLASS ="small">Organism count: }.scalar @opts.qq{</FONT>\n<BR>\n};
    unless (@opts) {
        $html .= qq{<input type="hidden" name="org_id" id="org_id">};
        return $html;
    }
    #print STDERR $i . '+++++++++++++\n';
    $html .= qq{<SELECT class="ui-widget-content ui-corner-all" id="recent_org_id" SIZE="5" MULTIPLE onChange="recent_dataset_chain()" >\n}
          . join( "\n", @opts )
          . "\n</SELECT>\n";
    $html =~ s/OPTION/OPTION SELECTED/;
    return $html;
}

sub get_orgs {
    my %opts  = @_;
    my $name  = $opts{name};
    my $desc  = $opts{desc};
    my $oid   = $opts{oid};
    my $dsgid = $opts{dsgid};

    return ('', 0) unless ($name || $desc);

    my $dsg   = $coge->resultset('Genome')->find($dsgid) if $dsgid;
    $dsg = undef unless $USER->has_access_to_genome($dsg);

    my $search_term = '%' . $name . '%';
    my @db = $coge->resultset("Organism")->search(
        \[
            'name LIKE ? OR description LIKE ?',
            [ 'name', $search_term ], [ 'description', $search_term ]
        ]
    );
    
    if (@db > $MAX_NUM_ORGANISM_RESULTS) {
        return (
            qq{<input type="hidden" name="org_id" id="org_id"><span class="small alert">Please refine your search</span>},
            scalar(@db)
        );
    }

    my @opts;
    foreach my $item ( sort { uc( $a->name ) cmp uc( $b->name ) } @db ) {
        my $option = "<OPTION value=\"" . $item->id . "\"";
        $option .= " SELECTED" if $oid && $item->id == $oid;
        $option .= " SELECTED" if $dsg && $item->id == $dsg->organism->id;
        my $name = $item->name;
        if ( length($name) > 50 ) {
            $name = substr( $name, 0, 50 ) . "...";
        }
        $option .= ">" . $name . " (id" . $item->id . ")</OPTION>";
        push @opts, $option;
    }
    my $html;
    if ( ( $name || $desc ) && !@opts ) {
        $html .= qq{<input type="hidden" name="org_id" id="org_id"><span class="small alert">No organisms found</span>};
        return $html, 0;
    }
    $html .= qq{<SELECT class="ui-widget-content ui-corner-all" id="org_id" size="5" MULTIPLE onChange="get_org_info(['args__oid','org_id'],[genome_chain])" >\n};
    $html .= join( "\n", @opts );
    $html .= "\n</SELECT>\n";
    $html =~ s/OPTION/OPTION SELECTED/ unless $html =~ /SELECTED/;
    my $opts = "?";
    $opts .= "name=$name;"   if $name;
    $opts .= "desc=$desc;"   if $desc;
    $opts .= "oid=$oid;"     if $oid;
    $opts .= "dsgid=$dsgid;" if $dsgid;
    $html .= qq{<!--<div class="padded"><span class='ui-button ui-corner-all' onclick="window.open('get_org_list.pl$opts');">Download Organism List</span></div>-->};
    return $html, scalar(@db);
}

sub update_genomelist {
    my %opts      = @_;
    my $genome_id = $opts{genomeid};
    return unless $genome_id;
    my $dsg = $coge->resultset("Genome")->find($genome_id);
    my $genome_name;
    $genome_name = $dsg->name;
    $genome_name = $dsg->organism->name unless $genome_name;
    $genome_name .= " (v" . $dsg->version . ")";
    return $genome_name, $genome_id;
}

sub get_org_info {
    my %opts = @_;

    my $oid = $opts{oid};
    return " " unless $oid;
    my $org = $coge->resultset("Organism")->find($oid);
    return "Unable to find an organism for id: $oid\n" unless $org;
    my $html;    # = qq{<div class="backbox small">};
    $html .=
      "<span class=alert>Private Organism!  Authorized Use Only!</span><br>"
      if $org->restricted;
    $html .= qq{<table class='small annotation_table'>};
    $html .= qq{<tr><td>Name:};
    $html .= qq{<td>} . $org->name;

    if ( $org->description ) {
        $html .= qq{<tr><td valign='top'>Description:<td>};

        my $desc_len;
        foreach my $item ( split( /;/, $org->description ) ) {
            $item =~ s/^\s+//;
            $item =~ s/\s+$//;
            $html .= "<a href=OrganismView.pl?org_desc=$item>$item</a>;";
            $desc_len += length($item);
            if ( $desc_len > 100 ) {
                $html .= '<br>';
                $desc_len = 0;
            }
        }
    }
    $html .=
"<tr><td>Links:<td><a href='OrganismView.pl?oid=$oid' target=_new>OrganismView</a>&nbsp|&nbsp<a href='CodeOn.pl?oid=$oid' target=_new>CodeOn</a>";
    $html .= "<tr><td>Search:<td>";
    my $search_term = $org->name;
    $html .=
qq{<img onclick="window.open('http://www.ncbi.nlm.nih.gov/taxonomy?term=$search_term')" src = "picts/other/NCBI-icon.png" title="NCBI" class=link>&nbsp};
    $html .=
qq{<img onclick="window.open('http://en.wikipedia.org/w/index.php?title=Special%3ASearch&search=$search_term')" src = "picts/other/wikipedia-icon.png" title="Wikipedia" class=link>&nbsp};
    $search_term =~ s/\s+/\+/g;
    $html .=
qq{<img onclick="window.open('http://www.google.com/search?q=$search_term')" src="picts/other/google-icon.png" title="Google" class=link>};
    $html .= "</table>";

    #    $html .= "</div>";
    return $html;
}

sub get_genome_list_for_org {
    my %opts = @_;
    my $oid  = $opts{oid};
    my $org  = $coge->resultset("Organism")->find($oid);
    my @opts;
    if ($org) {
        my @dsg;
        foreach my $dsg ( $org->genomes ) {
            next if $dsg->deleted;
            next unless $USER->has_access_to_genome($dsg);
            $dsg->name( $org->name ) unless $dsg->name;
            push @dsg, $dsg;
        }
        @opts = map {
                $_->id . "%%"
              . $_->name . " (v"
              . $_->version
              . ", dsgid"
              . $_->id . "): "
              . $_->genomic_sequence_type->name
          } sort {
            versioncmp( $b->version, $a->version )
              || $a->type->id <=> $b->type->id
              || $a->name cmp $b->name
              || $b->id cmp $a->id
          } @dsg;
    }

    my $res = join( "&&", @opts );
}

sub get_genomes {
    my %opts  = @_;
    my $oid   = $opts{oid};
    my $dsgid = $opts{dsgid};
    my $org   = $coge->resultset("Organism")->find($oid);
    my @opts;
    my %selected;
    $selected{$dsgid} = "SELECTED" if $dsgid;
    if ($org) {
        my @dsg;
        foreach my $dsg ( $org->genomes ) {
            next if $dsg->deleted;
            next unless $USER->has_access_to_genome($dsg);
            push @dsg, $dsg;
        }
        foreach my $dsg (@dsg) {
            $selected{ $dsg->id } = " " unless $selected{ $dsg->id };
        }
        no warnings 'uninitialized'; # disable warnings for undef values in sort
        foreach my $item (
            sort {
                     versioncmp( $b->version, $a->version )
                  || $a->type->id <=> $b->type->id
                  || $a->name cmp $b->name
                  || $b->id cmp $a->id
            } @dsg
          )
        {
            my $string =
                "<OPTION value=\""
              . $item->id . "\" "
              . $selected{ $item->id } . ">"
              . $item->info
              . "</OPTION>";
            push @opts, $string;
        }
    }
    my $html;
    if (@opts) {

#	$html = qq{<FONT CLASS ="small">Dataset group count: }.scalar (@opts).qq{</FONT>\n<BR>\n};
        $html .= qq{<SELECT class="ui-widget-content ui-corner-all" style="max-width:500px;" id="dsg_id" size="5" MULTIPLE onChange="get_genome_info(['args__dsgid','dsg_id'],[dataset_chain]);" >\n};
        $html .= join( "\n", @opts );
        $html =~ s/OPTION/OPTION SELECTED/ unless $html =~ /SELECTED/i;
    }
    else {
        $html .= qq{<input type = hidden name="dsg_id" id="dsg_id">};
    }
    return $html, scalar @opts;
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
    $cmd .=
" -tag 'imported_from=CoGe: http://genomevolution.org/CoGe/OrganismView.pl?dsgid=$dsgid'";
    system($cmd);
    print STDERR $cmd;

    #return $cmd;
    return "Complete!";
}

sub get_genome_info {
    my %opts  = @_;
    my $dsgid = $opts{dsgid};
    return " " unless $dsgid;
    my $dsg = $coge->resultset("Genome")->find($dsgid);
    return "Unable to get genome object for id: $dsgid" unless $dsg;

    my $html;
    my $total_length = $dsg->length;

    my $chr_num = $dsg->chromosome_count();#$dsg->genomic_sequences->count();
    $html .= qq{<table>} . "<tr valign=top><td><table class='small annotation_table'>";
    $html .= qq{<tr><td>Name:</td><td>} . $dsg->name . qq{</td></tr>}
      if $dsg->name;
    $html .=
      qq{<tr><td>Description:</td><td>} . $dsg->description . qq{</td></tr>}
      if $dsg->description;

    my $gstid    = $dsg->genomic_sequence_type->id;
    my $gst_name = $dsg->genomic_sequence_type->name;
    $gst_name .= ": " . $dsg->type->description if $dsg->type->description;
    $html .=
        qq{<tr><td>Sequence type <a href="SeqType.pl">?</a>: <td>}
      . $gst_name
      . qq{ (gstid$gstid)<input type=hidden id=gstid value=}
      . $gstid
      . qq{></td></tr>}
      . qq{<tr><td>Length: </td>}
      . qq{<td><div style="float: left;"> }
      . commify($total_length)
      . " bp </div>";

    # mdb removed 7/31/13 issue 77
    #    my $seq_file = $dsg->file_path;
    #    my $cogedir  = $P->{COGEDIR};
    #    my $cogeurl  = $P->{URL};
    #    $seq_file =~ s/$cogedir/$cogeurl/i;
    my $seq_url = "services/service.pl/sequence/$dsgid";#my $seq_url = "services/JBrowse/service.pl/sequence/$dsgid"; # mdb added 7/31/13 issue 77 # mdb changed 8/9/14

    $html .= "<tr><td>Links:</td>"
     . qq{<td>}
     . qq{<a href="GenomeInfo.pl?gid=$dsgid"><strong>Click here for more info</strong></a>}
     . qq{&nbsp|&nbsp}
     . "<a href='OrganismView.pl?dsgid=$dsgid' target=_new>OrganismView</a>&nbsp|&nbsp<a href='CodeOn.pl?dsgid=$dsgid' target=_new>CodeOn</a>"
     . qq{&nbsp|&nbsp}
     . qq{<span class='link' onclick="window.open('SynMap.pl?dsgid1=$dsgid;dsgid2=$dsgid');">SynMap</span>}
     . qq{&nbsp|&nbsp}
     . qq{<span class='link' onclick="window.open('CoGeBlast.pl?dsgid=$dsgid');">CoGeBlast</span>};

#temporarily removed until this is connected correctly for individual users
#    $html .= qq{&nbsp|&nbsp};
#    $html .= qq{<span id=irods class='link' onclick="gen_data(['args__loading...'],['irods']);add_to_irods(['args__dsgid','args__$dsgid'],['irods']);">Send To iPlant Data Store</span>};
    $html .= "</td></tr>"
          . qq{<tr><td colspan=2><div class="padded"><span class="ui-button ui-corner-all" onClick="update_genomelist(['args__genomeid','args__$dsgid'],[add_to_genomelist]);\$('#geno_list').dialog('option', 'width', 500).dialog('open');">Add to Genome List</span>}
          . qq{</div></td></tr>}
          . "</table></td>"
          . qq{<td id=dsg_features></td>}
          . "</table>";

    return $html;
}

sub get_dataset {
    my %opts   = @_;
    my $dsgid  = $opts{dsgid};
    my $dsname = $opts{dsname};
    my $dsid   = $opts{dsid};
    
    return qq{<input type="hidden" name="ds_id" id="ds_id">}, 0 unless ($dsid || $dsname || $dsgid);
    
    if ($dsid) {
        my ($ds) = $coge->resultset('Dataset')->resolve($dsid);
        $dsname = $ds->name;
    }
    
    my $html;
    my @opts;
    if ($dsgid) {
        my $dsg = $coge->resultset("Genome")->find($dsgid);
        @opts = map {
                "<OPTION value=\""
              . $_->id . "\">"
              . $_->name . " (v"
              . $_->version
              . ", dsid"
              . $_->id
              . ")</OPTION>"
          } sort { versioncmp ($b->version, $a->version) || $a->name cmp $b->name }
          $dsg->datasets
          if $dsg;

    }
    elsif ($dsname) {
        my @ds =
          $coge->resultset("Dataset")
          ->search( { name => { like => "%" . $dsname . "%" } } );

        my %orgs;
        foreach my $item (
            sort {
                $b->version <=> $a->version
                  || uc( $a->name ) cmp uc( $b->name )
            } @ds
          )
        {
            next unless $USER->has_access_to_dataset($item);
            my $option =
                "<OPTION value=\""
              . $item->id . "\">"
              . $item->name . "(v"
              . $item->version . ", id"
              . $item->id
              . ")</OPTION>";
            if ( $dsid && $dsid == $item->id ) {
                $option =~ s/(<OPTION)/$1 selected/;
            }
            push @opts, $option;
            $orgs{ $item->organism->id } = $item->organism;
        }
    }
    if (@opts) {
        $html .=
qq{<SELECT class="ui-widget-content ui-corner-all" id="ds_id" SIZE="5" MULTIPLE onChange="dataset_info_chain()" >\n};
        $html .= join( "\n", @opts );
        $html .= "\n</SELECT>\n";
        $html =~ s/OPTION/OPTION SELECTED/ unless $dsid;
    }
    else {
        $html .= qq{<input type="hidden" name="ds_id" id="ds_id">};
    }
    return $html, scalar @opts;
}

sub get_dataset_info {
    my $dsid           = shift;
    return qq{<input type="hidden" id="chr" value="">}, " ", 0
      unless ($dsid); # error flag for empty dataset

    my $ds = $coge->resultset("Dataset")->find($dsid);
    unless ($ds) {
        print STDERR "get_dataset_info: unable to find dataset object for id: $dsid";
        return qq{<input type="hidden" id="chr" value="">}, " ", 0;
    }

	my $chr_num_limit = 20;    
    my $html = "";
    $html .= "<span class='alert small'>Restricted dataset</span><br>"
      if $ds->restricted;
    $html .= "<table>";
    $html .= "<tr valign=top><td><table class=\"small annotation_table\">";
    my $dataset = $ds->name;
    $dataset .= ": " . $ds->description if $ds->description;
    $dataset =
      " <a href=\"" . $ds->link . "\" target=_new\>" . $dataset . "</a>"
      if $ds->link;
    my $source_name = $ds->data_source->name;
    $source_name .= ": " . $ds->data_source->description
      if $ds->data_source->description;
    my $link = $ds->data_source->link;

    $link = "http://" . $link if ( $link && $link !~ /http/ );
    $source_name =
      "<a href =\"" . $link . "\" target=_new\>" . $source_name . "</a>"
      if $ds->data_source->link;
    $html .= qq{<tr><td>Name: <td>$dataset} . "\n"
          . qq{<TR><TD><span class="link" onclick="window.open('Sources.pl')">Data Source:</span> <TD>$source_name (id}
          . $ds->data_source->id . qq{)} . "\n"
          . qq{<tr><td>Version: <td>} . $ds->version . "\n"
          . qq{<tr><td>Organism:<td class="link"><a href="OrganismView.pl?oid=}
          . $ds->organism->id
          . qq{" target=_new>}
          . $ds->organism->name
          . "</a>\n"
          . qq{<tr><td>Date deposited: <td>} . $ds->date . "\n";

    my $html2;
    my $total_length = $ds->total_length( ftid => 4 );
    my $chr_num = $ds->chromosome_count( ftid => 4 );

    my %chr;
    
# mdb removed 5/29/14 - doesn't work when chr name strings are returned from get_chromosomes() instead of db objects
#    map { $chr{ $_->chromosome } = { length => $_->stop } }
#      ( $ds->get_chromosomes( ftid => 4, length => 1, limit => $chr_num_limit )
#      );    #the chromosome feature type in coge is 301
      
    # mdb added 5/29/14
    my $tmp_count = 0;
    foreach my $c ( $ds->get_chromosomes( ftid => 4, length => 1, limit => $chr_num_limit ) ) {
        if (ref($c) =~ /CoGeX/ ) {
            $chr{ $c->chromosome } = { length => $c->stop }
        }
        else {
            $chr{ $c } = { length => 0 };
        }
    }
       
    my $count = 100000;
    foreach my $item ( sort keys %chr ) {
        my ($num) = $item =~ /(\d+)/;
        $num = $count unless $num;
        $chr{$item}{num} = $num;
        $count++;
    }
    my @chr =
      $chr_num > $chr_num_limit
      ? sort { $chr{$b}{length} <=> $chr{$a}{length} } keys %chr
      : sort { $chr{$a}{num} <=> $chr{$b}{num} || $a cmp $b } keys %chr;
    if (@chr) {
        my $size = scalar @chr;
        $size = 5 if $size > 5;
        my $select;
        $select .= qq{<SELECT class="ui-widget-content ui-corner-all" id="chr" size="$size" onChange="dataset_chr_info_chain()" >\n};
        $select .= join(
            "\n",
            map {
                    "<OPTION value=\"$_\">" 
                  . $_ . " ("
                  . commify( $chr{$_}{length} )
                  . " bp)</OPTION>"
              } @chr
        ) . "\n";
        $select =~ s/OPTION/OPTION SELECTED/;
        $select .= "\n</SELECT>\n";

        $html2 .= $select;
    }
    else {
        $html2 .= qq{<input type="hidden" id="chr" value="">};
        $html2 .= "<tr><td>No chromosomes";
    }
    $html .= "<tr><td>Chromosome count:<td><div style=\"float: left;\">"
      . commify($chr_num);
    $html .= "<tr><td>Total length:<td><div style=\"float: left;\">"
      . commify($total_length) . " bp ";
    my $gc =
      $total_length < 10000000 && $chr_num < $chr_num_limit
      ? get_gc_for_chromosome( dsid => $ds->id )
      : 0;
    $gc =
        $gc
      ? $gc
      : qq{  </div><div style="float: left; text-indent: 1em;" id=dataset_gc class="link" onclick="gen_data(['args__loading...'],['dataset_gc']);\$('#dataset_gc').removeClass('link'); get_gc_for_chromosome(['args__dsid','ds_id','args__gstid', 'gstid'],['dataset_gc']);">  Click for percent GC content</div>}
      if $total_length;
    $html .= $gc if $gc;
    $html .= qq{<tr><td>Links:</td>};
    $html .= "<td>";
    $html .= "<a href='OrganismView.pl?dsid=$dsid' target=_new>OrganismView</a>";
    $html .= qq{</td></tr>};
    my $feat_string = qq{
<tr><td><div id=ds_feature_count class="small link" onclick="gen_data(['args__loading...'],['ds_features']);get_feature_counts(['args__dsid','ds_id','args__gstid', 'gstid'],['ds_features']);" >Click for Features</div></td></tr>};
    $html .= $feat_string;

    $html .= qq{</table></td>};
    $html .= qq{<td id=ds_features></td>};
    $html .= qq{</table>};

    my $chr_count = $chr_num;
    $chr_count .= " <span class='note'> (only $chr_num_limit largest listed)</span>"
      if ( $chr_count > $chr_num_limit );
    return $html, $html2, $chr_count;
}

sub get_dataset_chr_info {
    my $dsid  = shift;
    my $chr   = shift;
    my $dsgid = shift;
    $dsgid = 0 unless defined $dsgid;
    $dsid  = 0 unless $dsid;
    unless ( $dsid && defined $chr )    # error flag for empty dataset
    {
        return "", "", "";
    }
    my $start = "'start'";
    my $stop  = "'stop'";
    my $html .= "<table>";
    $html .= "<tr valign=top><td><table class=\"small annotation_table\">";
    my $ds = $coge->resultset("Dataset")->find($dsid);
    return $html unless $ds;
    my $length = 0;
    $length = $ds->last_chromosome_position($chr) if defined $chr;
    my $gc =
      $length < 10000000
      ? get_gc_for_chromosome( dsid => $ds->id, chr => $chr )
      : 0;
    $gc =
        $gc
      ? $gc
      : qq{<div style="float: left; text-indent: 1em;" id=chromosome_gc class="link" onclick="\$('#chromosome_gc').removeClass('link'); get_gc_for_chromosome(['args__dsid','ds_id','args__chr','chr','args__gstid', 'gstid'],['chromosome_gc']);">Click for percent GC content</div>};
    $length = commify($length) . " bp ";
    $html .= qq{
<tr><td>Chromosome:</td><td>$chr</td></tr>
<tr><td>Nucleotides:</td><td>$length</td><td>$gc</td></tr>
};

    $html .= qq{
<tr><td>Noncoding sequence:<td colspan=2><div id=noncoding_gc class="link" onclick = "gen_data(['args__loading...'],['noncoding_gc']);\$('#noncoding_gc').removeClass('link');  get_gc_for_noncoding(['args__dsid','ds_id','args__chr','chr','args__gstid', 'gstid'],['noncoding_gc']);">Click for percent GC content</div>
} if $length;

    my $feat_string = qq{
<tr><td><div class=small id=feature_count onclick="gen_data(['args__loading...'],['chr_features']);get_feature_counts(['args__dsid','ds_id','args__chr','chr','args__gstid', 'gstid'],['chr_features']);" >Click for Features</div></td></tr>};

    $html .= $feat_string;
    $html .= "</table></td>";
    $html .= qq{<td id=chr_features></td>};
    $html .= qq{</table>};
    my $viewer;
    if ( defined $chr ) {
        $viewer .= qq{<div class="coge-table-header _orgviewresult">Genome Viewer</div>}
         . "<table class=\"small ui-corner-all ui-widget-content _orgviewresult\">"
         . "<tr><td nowrap>Starting location: "
         . qq{<td><input type="text" size=10 value="20000" id="x">}
         . qq{<tr><td >Zoom level:<td><input type = "text" size=10 value ="6" id = "z">}
         . qq{<tr><td colspan=2><span style="font-size:1em" class='ui-button ui-button-icon-left ui-corner-all' onClick="launch_viewer('$dsgid', '$chr')"><span class="ui-icon ui-icon-newwin"></span>Launch Genome Viewer</span>}
         . "</table>";
    }
    my $seq_grab;
    if ( defined $chr ) {
        $seq_grab .= qq{<div class="coge-table-header _orgviewresult">Genomic Sequence Retrieval</div>}
         . qq{<table class=\"small ui-corner-all ui-widget-content _orgviewresult padded\">}
         . "<tr><td>Start position: "
         . qq{<td><input type="text" size=10 value="1" id="start">}
         . "<tr><td>End position: "
         . qq{<td><input type="text" size=10 value="100000" id="stop">}
         . qq{<tr><td colspan=2><span style="font-size:1em" class='ui-button ui-button-icon-left ui-corner-all' onClick="launch_seqview('$dsgid', '$chr','$dsid')"><span class="ui-icon ui-icon-newwin"></span>Get Sequence</span>}
         . qq{</table>};
    }
    return $html, $viewer, $seq_grab;
}

sub get_feature_counts {
    my %opts  = @_;
    my $dsid  = $opts{dsid};
    my $dsgid = $opts{dsgid};
    my $gstid = $opts{gstid};
    my $chr   = $opts{chr};
    my $query;
    my $name;
    if ($dsid) {
        my $ds = $coge->resultset('Dataset')->find($dsid);
        $name  = "dataset " . $ds->name;
        $query = qq{
SELECT count(distinct(feature_id)), ft.name, ft.feature_type_id
  FROM feature
  JOIN feature_type ft using (feature_type_id)
 WHERE dataset_id = $dsid
};
        $query .= qq{AND chromosome = '$chr'} if defined $chr;
        $query .= qq{
  GROUP BY ft.name
};
        $name .= " chromosome $chr" if defined $chr;
    }
    elsif ($dsgid) {
        my $dsg = $coge->resultset('Genome')->find($dsgid);
        $name = "dataset group ";
        $name .= $dsg->name ? $dsg->name : $dsg->organism->name;
        $query = qq{
SELECT count(distinct(feature_id)), ft.name, ft.feature_type_id
  FROM feature
  JOIN feature_type ft using (feature_type_id)
  JOIN dataset_connector dc using (dataset_id)
 WHERE genome_id = $dsgid
  GROUP BY ft.name

};
    }

    my $dbh = $coge->storage->dbh;  #DBI->connect( $connstr, $DBUSER, $DBPASS );
    my $sth = $dbh->prepare($query);
    $sth->execute;
    my $feats = {};
    while ( my $row = $sth->fetchrow_arrayref ) {
        my $name = $row->[1];
        $name =~ s/\s+/_/g;
        $feats->{$name} = {
            count => $row->[0],
            id    => $row->[2],
            name  => $row->[1],
        };
    }
    my $gc_args;
    $gc_args = "chr: '$chr'," if defined $chr;
    $gc_args .= "dsid: $dsid,"
      if $dsid
    ; #set a var so that histograms are only calculated for the dataset and not hte genome
    $gc_args .= "typeid: ";
    my $feat_list_string = $dsid ? "dsid=$dsid" : "dsgid=$dsgid";
    $feat_list_string .= ";chr=$chr" if defined $chr;
    my $feat_string;    # .= qq{<div>Features for $name</div>};
    $feat_string .= qq{<div class = " ui-corner-all ui-widget-content small">};
    $feat_string .= qq{<table class=small>};

    foreach my $type ( sort { $a cmp $b } keys %$feats ) {
        $feat_string .= "<tr valign=top>";
        $feat_string .=
            "<td valign=top><div id=$type  >"
          . $feats->{$type}{name}
          . " (ftid"
          . $feats->{$type}{id}
          . ")</div>";
        $feat_string .=
          "<td valign=top align=right>" . commify( $feats->{$type}{count} );
        $feat_string .= "<td><div id=" . $type . "_type class=\"link small\"
  onclick=\"
  \$('#gc_histogram').dialog('option','title', 'Histogram of GC content for "
          . $feats->{$type}{name} . "s');
  \$('#gc_histogram').dialog('open');"
          . "get_feat_gc({$gc_args"
          . $feats->{$type}{id} . "})\">"
          . '%GC Hist</div>';
        $feat_string .= "<td>|</td>";
        $feat_string .=
"<td class='small link' onclick=\"window.open('FeatList.pl?$feat_list_string"
          . "&ftid="
          . $feats->{$type}{id}
          . ";gstid=$gstid')\">FeatList";
        $feat_string .= "<td>|</td>";
        $feat_string .=
"<td class='small link' onclick=\"window.open('bin/get_seqs_for_feattype_for_genome.pl?ftid="
          . $feats->{$type}{id} . ";";
        $feat_string .= "dsgid=$dsgid;" if $dsgid;
        $feat_string .= "dsid=$dsid;"   if $dsid;
        $feat_string .= "')\">DNA Seqs";

        if ( $feats->{$type}{name} eq "CDS" ) {
            $feat_string .= "<td>|</td>";
            $feat_string .=
"<td class='small link' onclick=\"window.open('bin/get_seqs_for_feattype_for_genome.pl?p=1;ftid="
              . $feats->{$type}{id};
            $feat_string .= ";dsgid=$dsgid" if $dsgid;
            $feat_string .= ";dsid=$dsid"   if $dsid;
            $feat_string .= "')\">Prot Seqs";
        }
    }
    $feat_string .= "</table>";

    if ( $feats->{CDS} ) {
        my $args;
        $args .= "'args__dsid','ds_id',"   if $dsid;
        $args .= "'args__dsgid','dsg_id'," if $dsgid;
        $args .= "'args__chr','chr',"      if defined $chr;
        $feat_string .=
            "<div class=\"small link\" id=wobble_gc onclick=\"\$('#wobble_gc_histogram').html('loading...').show().dialog('open');"
          . "get_wobble_gc([$args],['wobble_gc_histogram']);\">"
          . "Histogram of wobble codon GC content"
          . "</div>";
        $feat_string .=
            "<div class=\"small link\" id=wobble_gc_diff onclick=\"\$('#wobble_gc_diff_histogram').html('loading...').show().dialog('open');"
          . "get_wobble_gc_diff([$args],['wobble_gc_diff_histogram']);\">"
          . "Histogram of diff(CDS GC vs. wobble codon GC) content"
          . "</div>";
        $feat_string .= "<div class=\"small link \" id='codon_usage' onclick=\"
        \$('#codon_usage_table').html('loading...').show().dialog('open');
        get_codon_usage([$args],['codon_usage_table']); 
        \">" . "Codon usage table" . "</div>";
        $feat_string .= "<div class=\"small link\" id=aa_usage onclick=\"
        \$('#aa_usage_table').html('loading...').show().dialog('open');
        get_aa_usage([$args],[open_aa_usage_table]); 
        \">" . "Amino acid usage table" . "</div>";

    }
    $feat_string .= "</div>";
    $feat_string .= "None" unless keys %$feats;
    return $feat_string;
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
Max: <input type=text size=3 id=feat_gc_max value=$max>
Type: <select id=feat_hist_type>
<option value ="counts">Counts</option>
<option value = "percentage">Percentage</option>
</select>
};
    $info =~ s/>Per/ selected>Per/ if $hist_type =~ /per/;
    my $gc_args;
    $gc_args = "chr: '$chr'," if defined $chr;
    $gc_args .= "dsid: $dsid,"
      if $dsid
    ; #set a var so that histograms are only calculated for the dataset and not hte genome
    $gc_args .= "typeid: '$typeid'";
    $info .=
qq{<span class="link" onclick="get_feat_gc({$gc_args})">Regenerate histogram</span>};
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
        $info .=
qq{<div class=small style="color: red;">Limits set:  MIN: $min  MAX: $max</div>
};
    }
    my $stuff = join "::", @fids;
    $info .=
qq{<div class="link small" onclick="window.open('FeatList.pl?fid=$stuff')">Open FeatList of Features</div>};

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
            foreach my $loc ( $feat->locations ) {
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

sub get_codon_usage {
    my %opts  = @_;
    my $dsid  = $opts{dsid};
    my $chr   = $opts{chr};
    my $dsgid = $opts{dsgid};
    my $gstid = $opts{gstid};
    return unless $dsid || $dsgid;

    my $search = { "feature_type.name" => "CDS" };
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
    my $feat_count  = 0;
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
            print STDERR ".($feat_count)" if !$feat_count % 10;
        }
    }
    %codons = map { $_, $codons{$_} / $codon_total } keys %codons;

    my $html = "<div class='coge-table-header'>Codon Usage: $code_type</div>" .
     	CoGe::Accessory::genetic_code->html_code_table(
        	data => \%codons,
        	code => $code
    	) . '<br>';
    return $html;
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
    my @data;
    my ($dsg) = $coge->resultset('Genome')->find($dsgid);
    unless ($dsg) {
        my $error = "unable to create genome object using id $dsgid\n";
        return $error;
    }
    foreach my $gs ( $dsg->genomic_sequences ) {
        push @data, $gs->sequence_length;
    }
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
      . commify( $dsg->genomic_sequences->count )
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
