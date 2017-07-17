#! /usr/bin/perl -w

use strict;
use CGI;
use CoGeX;
use CoGeX::Result::Genome qw(ERROR LOADING);
use CoGe::JEX::Jex;
use CoGe::Accessory::Web;
use CoGe::Accessory::Utils qw(sanitize_name get_unique_id commify execute);
use CoGe::Accessory::IRODS qw(irods_iput irods_imeta_add);
use CoGe::Core::Chromosomes;
use CoGe::Core::Genome;
use CoGe::Core::Metadata;
use CoGe::Core::Experiment qw(experimentcmp);
use CoGe::Core::Storage;
use CoGe::Core::Favorites;
use CoGe::Builder::CommonTasks;

use HTML::Template;
use JSON::XS;
use Sort::Versions;
use File::Basename qw(basename);
use File::Path qw(mkpath);
use File::Spec::Functions;
use File::Copy;
use POSIX qw(floor);
use Data::Dumper;

no warnings 'redefine';

use vars qw(
  $P $PAGE_TITLE $TEMPDIR $SECTEMPDIR $LOAD_ID $USER $config $DB $FORM %FUNCTION
  $MAX_SEARCH_RESULTS $LINK $node_types $ERROR $HISTOGRAM $TEMPURL $SERVER $JEX
  $JOB_ID $EMBED
);

$PAGE_TITLE = 'GenomeInfo';

$node_types = CoGeX::node_types();

use constant MAX_TITLE_LENGTH => 150;

$FORM = new CGI;
( $DB, $USER, $config, $LINK ) = CoGe::Accessory::Web->init(
    cgi => $FORM,
    page_title => $PAGE_TITLE
);

$JEX = CoGe::JEX::Jex->new( host => $config->{JOBSERVER}, port => $config->{JOBPORT} );
$JOB_ID  = $FORM->Vars->{'job_id'};
$LOAD_ID = ( defined $FORM->Vars->{'load_id'} ? $FORM->Vars->{'load_id'} : get_unique_id() );
$SECTEMPDIR = $config->{SECTEMPDIR} . $PAGE_TITLE . '/' . $USER->name . '/' . $LOAD_ID . '/';
$TEMPDIR   = $config->{TEMPDIR} . "/$PAGE_TITLE";
$TEMPURL   = $config->{TEMPURL} . "/$PAGE_TITLE";
$HISTOGRAM = $config->{HISTOGRAM};
$SERVER    = $config->{SERVER};

$MAX_SEARCH_RESULTS = 100;
$ERROR = encode_json({ error => 1 });

my %ajax = CoGe::Accessory::Web::ajax_func();

%FUNCTION = (
    get_genome_info            => \&get_genome_info,
    get_genome_data            => \&get_genome_data,
    edit_genome_info           => \&edit_genome_info,
    update_genome_info         => \&update_genome_info,
    update_owner               => \&update_owner,
    update_certified           => \&update_certified,
    toggle_favorite            => \&toggle_favorite,
    search_organisms           => \&search_organisms,
    search_users               => \&search_users,
    delete_genome              => \&delete_genome,
#    delete_dataset             => \&delete_dataset,
    check_login                => \&check_login,
    copy_genome                => \&copy_genome,
    export_fasta		       => \&export_fasta,
    export_file_chr		       => \&export_file_chr,
    annotate                   => \&annotate,
    get_datasets               => \&get_datasets,
    get_bed                    => \&get_bed,
    get_gff                    => \&get_gff,
    get_tbl                    => \&get_tbl,
    get_load_log               => \&get_load_log,
    get_progress_log           => \&get_progress_log,
    export_bed                 => \&export_bed,
    export_gff                 => \&export_gff,
    export_tbl                 => \&export_tbl,
    export_features            => \&export_features,
    get_genome_info_details    => \&get_genome_info_details,
    get_features               => \&get_feature_counts,
    gen_data                   => \&gen_data,
    get_gc_for_genome          => \&get_gc_for_genome,
    get_gc_for_noncoding       => \&get_gc_for_noncoding,
    get_gc_for_feature_type    => \&get_gc_for_feature_type,
    get_chr_length_hist        => \&get_chr_length_hist,
    get_chromosomes            => \&get_chromosomes,
    cache_chr_fasta			   => \&cache_chr_fasta,
    get_aa_usage               => \&get_aa_usage,
    get_wobble_gc              => \&get_wobble_gc,
    get_wobble_gc_diff         => \&get_wobble_gc_diff,
    get_codon_usage            => \&get_codon_usage,
    get_experiments            => \&get_experiments,
    recommend_certification    => \&recommend_certification,
    %ajax
);

CoGe::Accessory::Web->dispatch( $FORM, \%FUNCTION, \&generate_html );

sub get_genome_info_details {
    my %opts  = @_;
    my $dsgid = $opts{dsgid};
    return " " unless $dsgid;
    my $dsg = $DB->resultset("Genome")->find($dsgid);
    return "Unable to get genome object for id: $dsgid" unless $dsg;
    my $html;

    #TABLE
    $html .= qq{<div class="left coge-table-header">Statistics</div>};
    $html .= qq{<table style="padding: 2px; margin-bottom: 5px;" class="border-top">};
    my $total_length = $dsg->length;

    # Count
    my $chr_num = $dsg->chromosome_count();
    $html .= qq{<tr><td class="title5">Chromosomes/contigs:</td</td>}
        . qq{<td class="data5">} . commify($chr_num);

    # Histogram
    $html .= qq{&nbsp;&nbsp;&nbsp;<span class="link" onclick="chr_hist($dsgid);">Histogram</span>};
  	# chromosome list
    $html .= qq{ | <span class="link" onclick="chr_list();">List</span>};
    $html .= qq{</td></tr>};
    
 
    my $gstid    = ($dsg->genomic_sequence_type ? $dsg->genomic_sequence_type->id : '');
    my $gst_name = ($dsg->genomic_sequence_type ? $dsg->genomic_sequence_type->name : '');
    $gst_name .= ": " . $dsg->type->description if ( $dsg->type && $dsg->type->description );

    # Sequence Type
    $html .= qq{<tr><td class="title5">Sequence type:<td class="data5" title="gstid$gstid">}
      . $gst_name
      . qq{ <input type=hidden id=gstid value=}
      . $gstid
      . qq{></td></tr>};
    $html .= qq{<tr><td class="title5">Length: </td>};
    $html .=
        qq{<td class="data5"><div style="float: left;"> }
      . commify($total_length)
      . " bp </div>";

    my ($gc, $noncoding);

    eval {
        $gc = (has_statistic($dsg, "gc") or ($total_length < 10000000
        && $chr_num < 20)) ? get_gc_for_genome( dsgid => $dsgid ) : 0;

        $noncoding = (has_statistic($dsg, "noncoding_gc") and $total_length)
        ? get_gc_for_noncoding(dsgid => $dsgid) : 0;
    };

    $gc =
        $gc
      ? $gc :
      qq{  <div style="float: left; text-indent: 1em;" id="dsg_gc" class="link" onclick="get_gc_content('#dsg_gc', 'get_gc_for_genome');">%GC</div><br/>};
    $html .= "$gc</td></tr>";

    if ($noncoding) {
        # Non-coding Sequence
        $html .= qq{
        <tr><td class="title5">Noncoding sequence:<td class="data5">$noncoding</td></tr>
        } if $total_length;

    } else {
        # Non-coding Sequence
        $html .= qq{
    <tr><td class="title5">Noncoding sequence:<td><span id=dsg_noncoding_gc class="link small" onclick="get_gc_content('#dsg_noncoding_gc', 'get_gc_for_noncoding');">%GC</div></td></tr>
    } if $total_length;
    }

#temporarily removed until this is connected correctly for individual users
#    $html .= qq{&nbsp|&nbsp};
#    $html .= qq{<span id=irods class='link' onclick="gen_data(['args__loading...'],['irods']);add_to_irods(['args__dsgid','args__$dsgid'],['irods']);">Send To CyVerse Data Store</span>};
    $html .= "</td></tr>";
#    if ( my $exp_count = $dsg->experiments->count( { deleted => 0 } ) ) {
#        $html .= qq{<tr><td class="title5">Experiment count:</td>};
#        $html .=
#            qq{<td class="data5"><span class=link onclick=window.open('ExperimentList.pl?dsgid=}
#          . $dsg->id . "')>"
#          . $exp_count
#          . "</span></td></tr>";
#    }
    $html .= "</table>";

    $html .= qq{<div class="left coge-table-header">Features</div>}
          .  qq{<div id="genome_features" style="margin-bottom: 5px;" class="small padded link border-top" onclick="get_features('#genome_features');" >Click for Features</div>};

    return $html;
}

sub gen_data {
    my $message = shift;
    return qq{<font class="small alert">$message</font>};
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
        my $ds = $DB->resultset('Dataset')->find($dsid);
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
        my $dsg = $DB->resultset('Genome')->find($dsgid);
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

    my $dbh = $DB->storage->dbh;  #DBI->connect( $connstr, $DBUSER, $DBPASS );
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
    $gc_args .= "dsid: $dsid," if $dsid; #set a var so that histograms are only calculated for the dataset and not hte genome
    $gc_args .= "typeid: ";
    my $feat_list_string = $dsid ? "dsid=$dsid" : "dsgid=$dsgid";
    $feat_list_string .= ";chr=$chr" if defined $chr;
    my $feat_string;
    $feat_string .= qq{<table style="padding: 2px; margin-bottom: 5px;">};

    foreach my $type ( sort { $a cmp $b } keys %$feats ) {
        $feat_string .= "<tr valign=top>";
        $feat_string .=
            qq{<td valign=top class="title5"><div id="$type" }
          . 'title="ftid'
          . $feats->{$type}{id}
          . '">'
          . $feats->{$type}{name}
          . "</div>";
        $feat_string .=
          qq{<td class="data5"valign=top align=right>} . commify( $feats->{$type}{count} );

        $feat_string .= "<td><div id=" . $type . "_type class=\"link small\" onclick=\"\$('#gc_histogram').dialog('option','title', 'Histogram of GC content for "
          . $feats->{$type}{name} . "s');"
          . "get_feat_gc({$gc_args"
          . $feats->{$type}{id} . "})\">"
          . '%GC Hist</div>';

        $feat_string .= "<td>|</td>";
        $feat_string .= "<td class='small link' onclick=\"window.open('FeatList.pl?$feat_list_string"
          . "&ftid="
          . $feats->{$type}{id}
          . ";gstid=$gstid')\">FeatList";
        $feat_string .= "<td>|</td>";
        $feat_string .= "<td class='small link' onclick=\"window.open('get_seqs_for_feattype_for_genome.pl?ftid="
          . $feats->{$type}{id} . ";";
        $feat_string .= "dsgid=$dsgid;" if $dsgid;
        $feat_string .= "dsid=$dsid;"   if $dsid;
        $feat_string .= "')\">Nuc Seqs</td>";

        my $fid = $feats->{$type}{id};

        if ($dsgid) {
            $feat_string .= qq{<td>|</td><td class="small link" onclick="export_features_to_irods($dsgid, $fid, true, 0);">Export Nuc Seqs</td>};
        }
        else {
            $feat_string .= qq{<td>|</td><td class="small link" onclick="export_features_to_irods($dsid, $fid, false, 0);">Export Nuc Seqs</td>};
        }

        if ( $feats->{$type}{name} eq "CDS" ) {
            $feat_string .= "<td>|</td>";
            $feat_string .= "<td class='small link' onclick=\"window.open('get_seqs_for_feattype_for_genome.pl?p=1;ftid="
              . $feats->{$type}{id};
            $feat_string .= ";dsgid=$dsgid" if $dsgid;
            $feat_string .= ";dsid=$dsid"   if $dsid;
            $feat_string .= "')\">Prot Seqs";

            if ($dsgid) {
                $feat_string .= qq{<td>|</td><td class="small link" onclick="export_features_to_irods($dsgid, $fid, true, 1);">Export Prot Seqs</td>};
            }
            else {
                $feat_string .= qq{<td>|</td><td class="small link" onclick="export_features_to_irods($dsid, $fid, false, 1);">Export Prot Seqs</td>};
            }
        }

    }

    if ( $feats->{CDS} ) {
        my $param = defined $chr ? $chr : "";

        # Wobble codon
        $feat_string .=qq{<tr><td colspan="13" class="small link" id="wobble_gc" onclick="get_content_dialog('#wobble_gc_histogram', 'get_wobble_gc', '$param');">Histogram of wobble codon GC content</td></tr>};

        # Diff content
        $feat_string .= qq{<tr><td colspan="13" class="small link" id="wobble_gc_diff" onclick="get_content_dialog('#wobble_gc_diff_histogram','get_wobble_gc_diff', '$param');">Histogram of diff(CDS GC vs. wobble codon GC) content</td></tr>};

        #Codon usage tables
        $feat_string .= qq{<tr><td colspan="13" class="small link" id="codon_usage" onclick="get_content_dialog('#codon_usage_table', 'get_codon_usage', '$param');">Codon usage table</td></tr>};

        #Amino acid usage table
        $feat_string .= qq{<tr><td colspan="13" class="small link" id="aa_usage" onclick="open_aa_usage_table('$param');">Amino acid usage table</td></tr>};
    }

    $feat_string .= "</table>";
    $feat_string .= "None" unless keys %$feats;
    return $feat_string;
}

sub get_codon_usage {
    my %opts  = @_;
    my $chr   = $opts{chr};
    my $dsgid = $opts{dsgid};
    return unless $dsgid;

    my $genome = $DB->resultset('Genome')->find($dsgid);
    return "unable to find genome id $dsgid\n" unless $dsgid;
    return $genome->get_codon_usage($chr);
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
    my $search;
    $search = { "feature_type_id" => 3 };
    $search->{"me.chromosome"} = $chr if defined $chr;

    my $genome = $DB->resultset('Genome')->find($dsgid);
    my $raw = get_wobble_histogram($genome);

    my @dsids;
    push @dsids, map { $_->id } $genome->datasets();

    my @fids = keys $raw;
    my @data = map { $_->{percent_gc} } values $raw;

    my ($gc, $at, $n) = (0, 0, 0);

    for my $item (values $raw) {
        $gc += $item->{gc} if $item->{gc};
        $at += $item->{at} if $item->{at};
        $n  += $item->{n}  if $item->{n};
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
    $cmd .= " -o $out";
    $cmd .= " -t \"CDS wobble gc content\"";
    $cmd .= " -min 0";
    $cmd .= " -max 100";
    $cmd .= " -ht $hist_type" if $hist_type;
    `$cmd`;
    $min = 0   unless defined $min && $min =~ /\d+/;
    $max = 100 unless defined $max && $max =~ /\d+/;
    my $info;
    $info .= qq{<div class="small">
Min: <input type="text" size="3" id="wobble_gc_min" value="$min">
Max: <input type=text size=3 id=wobble_gc_max value=$max>
Type: <select id=wobble_hist_type>
<option value ="counts">Counts</option>
<option value = "percentage">Percentage</option>
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
    $info .=
qq{<span class="link" onclick="get_wobble_gc([$args],['wobble_gc_histogram']);\$('#wobble_gc_histogram').html('loading...');">Regenerate histogram</span>};
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
        $info .=
qq{<div class=small style="color: red;">Limits set:  MIN: $min  MAX: $max</div>
};
    }
    my $stuff = join "::", @fids;
    $info .=
qq{<div class="link small" onclick="window.open('FeatList.pl?fid=$stuff')">Open FeatList of Features</div>};
    $out =~ s/$TEMPDIR/$TEMPURL/;
    my $hist_img = "<img src=\"$out\">";
    return $info . "<br>" . $hist_img;
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

    my $genome = $DB->resultset('Genome')->find($dsgid);
    my $data = get_wobble_gc_diff_histogram($genome);
    foreach my $ds ( $genome->datasets() ) {
        push @dsids, $ds->id;
    }

    return "error", " " unless @$data;
    my $file = $TEMPDIR . "/" . join( "_", @dsids ) . "_wobble_gc_diff.txt";
    open( OUT, ">" . $file );
    print OUT "#wobble gc for dataset ids: " . join( " ", @dsids ), "\n";
    print OUT join( "\n", @$data ), "\n";
    close OUT;
    my $cmd = $HISTOGRAM;
    $cmd .= " -f $file";
    my $out = $file;
    $out =~ s/txt$/png/;
    $cmd .= " -o $out";
    $cmd .= " -t \"CDS GC - wobble gc content\"";
    `$cmd`;
    my $sum = 0;
    map { $sum += $_ } @$data;
    my $mean = sprintf( "%.2f", $sum / scalar @$data );
    my $info = "Mean $mean%";
    $info .= " (";
    $info .= $mean > 0 ? "CDS" : "wobble";
    $info .= " is more GC rich)";
    $out =~ s/$TEMPDIR/$TEMPURL/;
    my $hist_img = "<img src=\"$out\">";
    return $info . "<br>" . $hist_img;
}

sub export_features {
    my %args = @_;
    my $genome = $DB->resultset('Genome')->find($args{gid});
    my $ft = $DB->resultset('FeatureType')->find($args{fid});

    return 1 unless ($USER->has_access_to_genome($genome));

    my $workflow = $JEX->create_workflow(name => "Export features");

    my $basename = sanitize_name($genome->organism->name . "-ft-" . $ft->name);

    $args{basename} = $basename;

    my ($output, %task) = generate_features(%args);
    $workflow->add_job(\%task);

    my $response = $JEX->submit_workflow($workflow);
    say STDERR "RESPONSE ID: " . $response->{id};
    $JEX->wait_for_completion($response->{id});

    my %json;
    $json{file} = basename($output);

    my %meta = (
        'Imported From' => "CoGe: http://genomevolution.org",
        'CoGe OrganismView Link' => "http://genomevolution.org/CoGe/OrganismView.pl?gid=".$genome->id,
        'CoGe GenomeInfo Link'=> "http://genomevolution.org/CoGe/GenomeInfo.pl?gid=".$genome->id,
        'CoGe Genome ID'   => $genome->id,
        'Organism Name'    => $genome->organism->name,
        'Organism Taxonomy'    => $genome->organism->description,
        'Version'     => $genome->version,
        'Type'        => $genome->type->info,
        'Feature Type' => $ft->name,
        'Data Type'    => "FASTA",
    );
    $meta{'Feature Description'} = $ft->description if $ft->description;

    $json{error} = export_to_irods( file => $output, meta => \%meta );

    return encode_json(\%json);
}

sub get_gc_for_genome {
    my %opts  = @_;
    my $dsid  = $opts{dsid};
    my $chr   = $opts{chr};
    my $gstid = $opts{gstid};
    my $dsgid = $opts{dsgid};

    my $genome = $DB->resultset('Genome')->find($dsgid);
    my $data = get_stats($genome, 'gc', \&calc_gc);

    # Skip if no data
    return unless $data;

    my $results =
        "&nbsp(GC: "
      . sprintf( "%.2f", 100 * $data->{gc})
      . "%  AT: "
      . sprintf( "%.2f", 100 * $data->{at})
      . "%  N: "
      . sprintf( "%.2f", 100 * $data->{n})
      . "%  X: "
      . sprintf( "%.2f", 100 * $data->{x}) . "%)";

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

    my $genome = $DB->resultset('Genome')->find($dsgid);
    my $data = get_stats($genome, 'noncoding_gc', \&calc_noncoding_gc);

    return
        commify($data->{total}) . " bp"
      . "&nbsp(GC: "
      . sprintf( "%.2f", 100 * $data->{gc})
      . "%  AT: "
      . sprintf( "%.2f", 100 * $data->{at}) . "% N: "
      . sprintf( "%.2f", 100 * $data->{n}) . "% X: "
      . sprintf( "%.2f", 100 * $data->{x}) . "%)";
}

sub get_gc_for_feature_type {
    my %opts   = @_;
    my $dsid   = $opts{dsid};
    my $dsgid  = $opts{dsgid};
    my $chr    = $opts{chr}; my $typeid = $opts{typeid};
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

    my $genome = $DB->resultset('Genome')->find($dsgid);
    my $raw = get_feature_type_gc_histogram($genome, $typeid);

    my (@items, @datasets);

    my @dsids;
    push @dsids, map { $_->id } $genome->datasets();

    #storage for fids that passed.  To be sent to FeatList
    my @fids = keys $raw;
    my @data = map { $_->{percent_gc} } values $raw;

    my ($gc, $at, $n) = (0) x 3;

    for my $item (values $raw) {
        $gc += $item->{gc} if $item->{gc};
        $at += $item->{at} if $item->{at};
        $n  += $item->{n}  if $item->{n};
    }

    my $total = $gc + $at + $n;
    return "error" unless $total;

    my $type = $DB->resultset('FeatureType')->find($typeid);
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

sub get_chr_length_hist { #TODO use API Genome Fetch
    my %opts  = @_;
    my $dsgid = $opts{dsgid};
    return "error", " " unless $dsgid;
    my ($dsg) = $DB->resultset('Genome')->find($dsgid);
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
    print OUT "#chromosome/contig lengths for $dsgid\n";
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
        "<table><TR><td>Count:<TD>"
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
    return $info . "<br>" . $hist_img;
}

# major hack, rewrite to get json from api, don't send chr_len to update_chromosome_list_plot_button
sub get_chromosomes { #TODO use API Genome Fetch
    my %opts  = @_;
    my $gid = $opts{gid};
	my $c = CoGe::Core::Chromosomes->new($gid);
    my $html = '[';
    my $first = 1;
	while ($c->next) {
		if ($first) {
			$first = 0;
		} else {
			$html .= ",";
		}
		$html .= '["' . $c->name . ' <a href=\\"GenomeView.pl?gid=' . $gid . '&loc=' . $c->name . '%3A1..' . ($c->length - 1) .
            '\\" target=\\"_blank\\"><span class=\\"glyphicon glyphicon-eye-open\\" style=\\"color:black;padding-left:20px;\\" title=\\"Browse\\"></span></a>","' .
            $c->length . '","<input type=\\"radio\\" name=\\"chr\\" id=\\"f' . $c->name . '\\" onchange=\\"update_chromosome_list_plot_button(' .
            $c->length . ')\\" /> FASTA","<input type=\\"radio\\" name=\\"chr\\" id=\\"g' . $c->name .
            '\\" onchange=\\"update_chromosome_list_plot_button(' .
            $c->length . ')\\" /> GFF","<input type=\\"radio\\" name=\\"chr\\" id=\\"n' . $c->name . '\\" onchange=\\"update_chromosome_list_plot_button(' .
            $c->length . ')\\" /> %GC/AT"]';
  	}
	$html .= ']';
	return $html;	
}

sub cache_chr_fasta { #FIXME remove this, use API genome/sequence
    my %opts  = @_;
    my $gid = $opts{gid};
    my $chr = $opts{chr};

    my $path = catfile($config->{SECTEMPDIR}, "downloads/genome", $gid);
	mkpath( $path, 0, 0777 ) unless -d $path;
    my $file = catfile($path, $gid . "_" . $chr . ".faa");
	my %json;
	$json{file} = $file;
    if (-e $file) {
    	return encode_json(\%json);
    }

    my $genome = $DB->resultset('Genome')->find($gid);
    return "unable to create genome object using id $gid" unless ($genome);

	# Get sequence from file
	my $seq = $genome->get_genomic_sequence(
		chr   => $chr,
		start => 1,
		stop  => $genome->get_chromosome_length($chr)
	);
	open(my $fh, '>', $file) or return "Could not open file '$file' $!";
    print $fh '>chromosome ' . $chr . "\n";
	for (my $i=0; $i<length($seq); $i+=70) {
		print $fh substr($seq, $i, 70);
		print $fh "\n";
	}
	close $fh;
    return encode_json(\%json);
}

sub get_genome_info {
    my %opts   = @_;
    my $gid    = $opts{gid};
    my $genome = $opts{genome};
    return unless ( $gid or $genome );

    unless ($genome) {
        $genome = $DB->resultset('Genome')->find($gid);
        return unless ($genome);
    }

    my $template = HTML::Template->new( filename => $config->{TMPLDIR} . $PAGE_TITLE . '.tmpl' );
    my $users_with_access = $genome->restricted ? join(', ', sort map { $_->display_name } $USER->users_with_access($genome)) : 'Everyone';

    $template->param(
        DO_GENOME_INFO => 1,
        ORGANISM       => $genome->organism->name,
        VERSION        => $genome->version,
        TYPE           => ($genome->type ? $genome->type->info : ''),
        SOURCE         => get_genome_sources($genome) // '',
        LINK           => $genome->link,
        USER_CAN_EDIT   => $genome->is_editable($USER),
        USER_CAN_DELETE => $genome->is_deletable($USER),
        RESTRICTED     => $genome->restricted,
        USERS_WITH_ACCESS => $users_with_access,
        NAME           => $genome->name,
        DESCRIPTION    => $genome->description,
        DELETED        => $genome->deleted,
        CERTIFIED       => $genome->certified,
        LOGON           => ( $USER->user_name ne "public" ),
    );

    my $owner = $genome->owner;
    my $creator = ($genome->creator_id ? $genome->creator->display_name : '');
    my $creation_date = $genome->get_date();
    my $groups = ($genome->restricted ? join(', ', sort map { $_->name } $USER->groups_with_access($genome)) : undef);
    $template->param( groups_with_access => $groups) if $groups;
    $template->param( OWNER => $owner->display_name ) if $owner;
    $template->param( CREATOR => $creator ) if $creator;
    $template->param( CREATION_DATE => $creation_date ) if $creation_date;
    $template->param( GID => $genome->id );

    return $template->output;
}

sub edit_genome_info {
    my %opts = @_;
    my $gid  = $opts{gid};
    return unless ($gid);

    my $genome = $DB->resultset('Genome')->find($gid);
    return unless ($genome);

    my $template =
      HTML::Template->new( filename => $config->{TMPLDIR} . $PAGE_TITLE . '.tmpl' );
    $template->param(
        EDIT_GENOME_INFO => 1,
        ORGANISM         => $genome->organism->name,
        VERSION          => $genome->version,
        TYPE             => $genome->type->name,
        SOURCE           => get_genome_sources($genome),
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

    #print STDERR "gid=$gid organism=$org_name version=$version source=$source_name\n";
    return "Error: missing params."
      unless ( $gid and $org_name and $version and $source_name );

    my $genome = $DB->resultset('Genome')->find($gid);
    return "Error: can't find genome." unless ($genome);

    my $organism = $DB->resultset('Organism')->find( { name => $org_name } );
    return "Error: can't find organism." unless ($organism);

    my $source = $DB->resultset('DataSource')->find( { name => $source_name } );
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

    #print STDERR "$search_term $timestamp\n";
    return unless $search_term;

    # Perform search
    $search_term = '%' . $search_term . '%';
    my @organisms = $DB->resultset("Organism")->search(
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
    my @users = $DB->resultset("User")->search(
        \[
            'user_name LIKE ? OR first_name LIKE ? OR last_name LIKE ?',
            [ 'user_name',  $search_term ],
            [ 'first_name', $search_term ],
            [ 'last_name',  $search_term ]
        ]
    );

    # Limit number of results displayed
    # if (@users > $MAX_SEARCH_RESULTS) {
    #   return encode_json({timestamp => $timestamp, items => undef});
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

    # Get user to become owner
    my $user = $DB->resultset('User')->find( { user_name => $user_name } );
    unless ($user) {
        return "error finding user '$user_name'\n";
    }

    # Remove existing owner
    foreach my $conn ( # should only be one owner but loop just in case
        $DB->resultset('UserConnector')->search(
            {
                parent_type => $node_types->{user},
                child_id    => $gid,
                child_type  => $node_types->{genome},
                role_id     => 2                        # FIXME hardcoded
            }
        )) 
    {
        $conn->delete;
    }
    
    # Remove existing user connection (should only be one, but loop just in case)
    foreach my $conn (
        $DB->resultset('UserConnector')->search(
            {
                parent_id   => $user->id,
                parent_type => $node_types->{user},
                child_id    => $gid,
                child_type  => $node_types->{genome},
                role_id     => {'!=' => 2}           # FIXME hardcoded
            }
        ))
    {
        $conn->delete;
    }
    
    # Make user owner of genome
    my $conn = $DB->resultset('UserConnector')->find_or_create(
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

    return;
}

sub update_certified {
    my %opts      = @_;
    my $gid       = $opts{gid};
    my $certified = $opts{certified};
    return unless ($gid and defined $certified);

    # Poweruser-only function (includes admins)
    return unless $USER->is_poweruser;

    # Get genome
    my $genome = $DB->resultset('Genome')->find($gid);
    return "Error: can't find genome." unless ($genome);
    
    # Update "certified" field for given genome
    $genome->certified($certified);
    $genome->update();
    
    # Record in log
    CoGe::Accessory::Web::log_history(
        db          => $DB,
        user_id     => $USER->id,
        page        => $PAGE_TITLE,
        description => ($certified ? 'Certify' : 'Uncertify') . ' genome ' . $genome->info_html,
        parent_id   => $gid,
        parent_type => 2 #FIXME magic number
    );
    
    return;
}

sub toggle_favorite {
    my %opts   = @_;
    my $gid = $opts{gid};
    return unless $gid;
    return if ($USER->is_public); # must be logged in
    
    # Get genome
    my $genome = $DB->resultset('Genome')->find($gid);
    return unless $genome;
    
    # Toggle favorite
    my $favorites = CoGe::Core::Favorites->new( user => $USER );
    my $is_favorited = $favorites->toggle($genome);
    
    # Record in log
    CoGe::Accessory::Web::log_history(
        db          => $DB,
        user_id     => $USER->id,
        page        => $PAGE_TITLE,
        description => ($is_favorited ? 'Favorited' : 'Unfavorited') . ' genome ' . $genome->info_html,
        parent_id   => $gid,
        parent_type => 2 #FIXME magic number
    );
    
    return $is_favorited;
}

sub get_genome_sources {
    my $genome = shift;
    return unless ($genome && $genome->first_dataset);
    
    # mdb removed 7/8/15
    #my %sources = map { $_->name => 1 } $genome->source;
    #return join( ',', sort keys %sources);
    
    # mdb added 7/8/15 - limit to one source
    return $genome->first_dataset->source->name;
}

sub get_genome_data {
    my %opts   = @_;
    my $gid    = $opts{gid};
    my $genome = $opts{genome};
    return unless ( $gid or $genome );

    $gid = $genome->id if $genome;
    unless ($genome) {
        $genome = $DB->resultset('Genome')->find($gid);
        return unless ($genome);
    }

    my $template =
      HTML::Template->new( filename => $config->{TMPLDIR} . $PAGE_TITLE . '.tmpl' );
    $template->param(
        #CHROMOSOME_COUNT => commify( $genome->chromosome_count() ),
        #LENGTH           => commify( $genome->length ),
    );

    return $template->output;
}

sub get_genome_download_links {
    my $genome = shift;

}

# Leave this commented out until issue 212 resolved
#sub delete_dataset {
#    my %opts = @_;
#    my $gid  = $opts{gid};
#    my $dsid  = $opts{dsid};
#    print STDERR "delete_dataset $gid $dsid\n";
#    return 0 unless ($gid and $dsid);
#
#    my $genome = $DB->resultset('Genome')->find($gid);
#    return 0 unless $genome;
#    return 0 unless ( $USER->is_admin or $USER->is_owner( dsg => $gid ) );
#
#
#    my $dataset = $DB->resultset('Dataset')->find($dsid);
#    return 0 unless $dataset;
#    return unless ($genome->dataset_connectors({ dataset_id => $dsid })); # make sure genome has specified dataset
#
#    print STDERR "delete_dataset: deleted dataset id$dsid\n";
#    $dataset->deleted(1);
#    $dataset->update;
#
#    # Record in log
#    CoGe::Accessory::Web::log_history(
#          db          => $DB,
#          user_id     => $USER->id,
#          page        => $PAGE_TITLE,
#          description => "delete dataset id $dsid in genome id $gid"
#      );
#
#    return 1;
#}

sub get_sequence_types {
    my $type_id = shift;

    my $html;
    foreach my $type ( sort { $a->info cmp $b->info }
        $DB->resultset('GenomicSequenceType')->all() )
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
    my %unique;
    foreach ( $DB->resultset('DataSource')->all() ) {
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
        $genome = $DB->resultset('Genome')->find($gid);
        return unless ($genome);
    }

    my @rows;
    foreach my $exp ( sort experimentcmp $genome->experiments ) {
        next if ( $exp->deleted );
        next if ( $exp->restricted && !$USER->has_access_to_experiment($exp) );
        my $info = $exp->info;
        $info =~ tr{\n}{ };
        push @rows, '["<span class=\\"link\\" onclick=\\"window.open(\'ExperimentView.pl?eid=' . $exp->id . '\')\\">' . $info . '</span>"]';
    }
    return '[' . join(',', @rows) . ']';
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
        $genome = $DB->resultset('Genome')->find($gid);
        return unless ($genome);
    }

    my @rows;
    foreach my $ds ( sort { $a->id <=> $b->id } $genome->datasets ) {
        #next if ($exclude_seq && filter_dataset($ds)); #FIXME add dataset "type" field instead?
        my $info = $ds->info;
        $info =~ tr{\n}{ };
        push @rows, '["<span>' . $info . '</span>' . ($ds->link ? '&nbsp;&nbsp;&nbsp;&nbsp;<a href=\\"' . $ds->link . '\\" target=_new>Link</a>' : '') . '"]';
    }
    return '[' . join(',', @rows) . ']';
}

sub delete_genome {
    my %opts = @_;
    my $gid  = $opts{gid};
    print STDERR "delete_genome $gid\n";
    return 0 unless $gid;

    my $genome = $DB->resultset('Genome')->find($gid);
    return 0 unless $genome;
    return 0 unless ( $USER->is_admin or $USER->is_owner( dsg => $gid ) );
    my $delete_or_undelete = ($genome->deleted ? 'undeleted' : 'deleted');
    #print STDERR "delete_genome " . $genome->deleted . "\n";
    $genome->deleted( !$genome->deleted ); # do undelete if already deleted
    $genome->update;

    # Record in log
    CoGe::Accessory::Web::log_history(
        db          => $DB,
        user_id     => $USER->id,
        page        => $PAGE_TITLE,
        description => "$delete_or_undelete genome " . $genome->info_html,
        parent_id   => $gid,
        parent_type => 2 #FIXME magic number
    );

    return 1;
}

sub check_login {
    #print STDERR $USER->user_name . ' ' . int($USER->is_public) . "\n";
    return ($USER && !$USER->is_public);
}

sub copy_genome {
    my %args = @_;
    my $gid  = $args{gid};
    my $mask = $args{mask};
    my $seq_only = $args{seq_only};

    print STDERR "GenomeInfo::copy_and_mask_genome: gid=$gid mask=$mask seq_only=$seq_only\n";

    if ($USER->is_public) {
        return 'Not logged in';
    }

    my $desc = $mask ? "Copying and masking genome" : "Copying genome";
    $desc .= " (no annotations)" if $seq_only;

    my $workflow = $JEX->create_workflow(name => $desc, init => 1);

    my ($staging_dir, $result_dir) = get_workflow_paths($USER->name, $workflow->id);
    $workflow->logfile( catfile($result_dir, 'debug.log') );

    $args{uid} = $USER->id;
    $args{wid} = $workflow->id;
    $args{staging_dir} = $staging_dir;

    my %task = copy_and_mask(%args);
    $workflow->add_job(\%task);

    my $response = $JEX->submit_workflow($workflow);

    unless ($JEX->is_successful($response)) {
        return encode_json({error => "The job could not be scheduled"});
    }

    # Get tiny link
    my $tiny_link = CoGe::Accessory::Web::get_tiny_link(
        url => $SERVER . "$PAGE_TITLE.pl?gid=" . $gid . ";job_id=" . $response->{id}
    );

    # Record in log
    CoGe::Accessory::Web::log_history(
        db          => $DB,
        user_id     => $USER->id,
        parent_id   => $workflow->id,
        parent_type => 7, #FIXME magic number
        page        => $PAGE_TITLE,
        description => $desc,
        link => $tiny_link
    );

    return encode_json({
        job_id => $response->{id},
        link   => $tiny_link
    });
}

sub get_load_log {
    my %opts = @_;
    my $gid = $opts{gid};
    my $getEverything = $opts{get_everything};
    return unless $gid;
    
    my $log = get_log( item_id => $gid, item_type => 'genome', getEverything => $getEverything, html => 1 );
    
    return $log;
}

sub get_progress_log {
    my %opts         = @_;
    my $workflow_id = $opts{workflow_id};
    return unless $workflow_id;
    
    my $results = get_workflow_results($USER->name, $workflow_id);
    return unless $results;

    my $genome_id;
    foreach (@$results) {
        if ($_->{type} eq 'genome') {
            $genome_id = $_->{id};
            last;
        }
    }

    return encode_json({ genome_id => $genome_id });
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
    $search->{"me.chromosome"} = $chr if defined $chr and $chr;

    my (@items, @datasets);
    if ($dsid) {
        my $ds = $DB->resultset('Dataset')->find($dsid);
        return "unable to find dataset id$dsid\n" unless $ds;
        push @items, $ds;
        push @datasets, $ds;
    }
    if ($dsgid) {
        my $dsg = $DB->resultset('Genome')->find($dsgid);
        return "unable to find genome id $dsgid\n" unless $dsgid;
        $gstid = $dsg->type->id;
        push @items, $dsg;
        push @datasets, $dsg->datasets;
    }

    my %seqs; # prefetch the sequences with one call to genomic_sequence (slow for many seqs)
    foreach my $item (@items) {
        my @chrs = (defined $chr and $chr) ? ($chr) : $item->chromosomes;

        for my $chr (@chrs) {
            $seqs{$chr} = $item->get_genomic_sequence(chr => $chr,
                                                      seq_type => $gstid);
        }
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
    my $html2 .= "Predicted amino acid usage using $code_type";
    $html2 .= CoGe::Accessory::genetic_code->html_aa_new( data => \%aa );
    return $html2; #return $html1, $html2;
}

sub export_fasta { #TODO migrate to API
    my %opts    = @_;
    my $gid = $opts{gid};
    #print STDERR "export_fasta $gid\n";

    my $genome = $DB->resultset('Genome')->find($gid);
    return unless ($USER->has_access_to_genome($genome));

    my $src = $genome->file_path;
    my $dest_filename = "genome_$gid.faa";
    my $dest = get_irods_path() . '/' . $dest_filename;
    unless ($src and $dest) {
        print STDERR "GenomeInfo:export_fasta: error, undef src or dest\n";
        return;
    }

    # Send to CyVerse Data Store using iput
    CoGe::Accessory::IRODS::irods_iput($src, $dest);
    #TODO need to check rc of iput and abort if failure occurred

    # Set IRODS metadata for object
    my $md = get_irods_metadata($genome);
    CoGe::Accessory::IRODS::irods_imeta_add($dest, $md);

	my %json;
	$json{file} = $dest_filename;
    return encode_json(\%json);
}

sub export_file_chr { #TODO migrate to API
    my %args = @_;
    my $dsg = $DB->resultset('Genome')->find($args{gid});
    my $file = $args{file};
#    print STDERR Dumper \%args;

    # # ensure user is logged in
    # return $ERROR if $USER->is_public;

    # ensure user has permission
    return $ERROR unless $USER->has_access_to_genome($dsg);

    my $md = get_irods_metadata($dsg);

    my %json;
    $json{file} = basename($file);
    $json{error} = export_to_irods( file => $file, meta => $md );
    return encode_json(\%json);
}

sub get_irods_path {
    my $username = $USER->user_name;
    my $dest = $config->{IRODSDIR};
    $dest =~ s/\<USER\>/$username/;
    return $dest;
}

sub linkify { # FIXME: this routine should be moved into the Utils module
    my ( $link, $desc ) = @_;
    return "<span class='link' onclick=\"window.open('$link')\">$desc</span>";
}

#
# TBL FILE
#

sub get_tbl {
    my %args = @_;
    my $dsg = $DB->resultset('Genome')->find($args{gid});

    # ensure user has permission
    return $ERROR unless $USER->has_access_to_genome($dsg);

    $args{basename} = sanitize_name($dsg->organism->name);

    my $workflow = $JEX->create_workflow(name => "Export Tbl");
    my ($output, $task) = generate_tbl(%args);
    $workflow->add_job($task);

    my $response = $JEX->submit_workflow($workflow);
    say STDERR "RESPONSE ID: " . $response->{id};
    $JEX->wait_for_completion($response->{id});

    my %json;
    $json{file} = basename($output);
    $json{files} = [ download_url_for(gid => $args{gid}, file => $output) ];

    return encode_json(\%json);
}

sub export_tbl { #TODO migrate to API
    return export_file_to_irods("Tbl", \&generate_tbl, @_);
}

#
# BED FILE
#

sub get_bed { #TODO migrate to API
    my %args = @_;
    my $dsg = $DB->resultset('Genome')->find($args{gid});

    # ensure user has permission
    return $ERROR unless $USER->has_access_to_genome($dsg);

    $args{basename} = sanitize_name($dsg->organism->name);

    my $workflow = $JEX->create_workflow(name => "Export bed file");
    my ($output, $task) = generate_bed(%args);
    $workflow->add_job($task);

    my $response = $JEX->submit_workflow($workflow);
    say STDERR "RESPONSE ID: " . $response->{id};
    $JEX->wait_for_completion($response->{id});

    my %json;
    $json{file} = basename($output);
    $json{files} = [ download_url_for(gid => $args{gid}, file => $output) ];

    return encode_json(\%json);
}

sub export_bed { #TODO migrate to API
    return export_file_to_irods("bed file", \&generate_bed, @_);
}

sub get_gff { #TODO migrate to API
    my %args = @_;
    
    # Get genome and check permission
    my $dsg = $DB->resultset('Genome')->find($args{gid});
    return $ERROR unless $USER->has_access_to_genome($dsg);

    # Create workflow
    my $workflow = $JEX->create_workflow( name => "Export gff" );
    $args{basename} = sanitize_name($dsg->organism->name);
    my ($output, $task) = generate_gff(%args);
    $workflow->add_job($task);

    # Submit workflow and block until completion
    my $response = $JEX->submit_workflow($workflow);
    say STDERR "GenomeInfo::get_gff: wid=" . $response->{id};
    $JEX->wait_for_completion($response->{id});

    return encode_json({
        file => basename($output),
        files => [ download_url_for(gid => $args{gid}, file => $output) ]
    });
}

sub export_gff { #TODO migrate to API
    return export_file_to_irods("gff", \&generate_gff, @_);
}

# %args:
# - gid
# - file_type, used to name workflow
# - generate_func, function to call to generate file to export
# %args is also passed to generate_func
sub export_file_to_irods { #TODO migrate to API
	my $file_type = shift;
	my $generate_func = shift;
    my %args = @_;
    my $dsg = $DB->resultset('Genome')->find($args{gid});

    # ensure user is logged in
    return $ERROR if $USER->is_public;

    # ensure user has permission
    return $ERROR unless $USER->has_access_to_genome($dsg);

    my $md = get_irods_metadata($dsg);

    $args{basename} = sanitize_name($dsg->organism->name);

    my $workflow = $JEX->create_workflow(name => "Export " . $file_type);
    my ($output, $task) = $generate_func->(%args);
    $workflow->add_job($task);

    my $response = $JEX->submit_workflow($workflow);
    say STDERR "RESPONSE ID: " . $response->{id};
    my $success = $JEX->wait_for_completion($response->{id});

    my %json;
    $json{file} = basename($output);
    if($success) {
        $json{error} = export_to_irods( file => $output, meta => $md );
    } else {
        $json{error} = 1;
    }

    return encode_json(\%json);
}

#XXX: Add error checking
sub export_to_irods { #TODO migrate to API
    my %args = @_;
    my $file = $args{file};
    my $meta = $args{meta};

    say STDERR "IFILE: $file";

    #Exit if the file does not exist
    return 1 unless -r $file;

    my $ipath = get_irods_path();
    my $ifile = catfile($ipath, basename($file));

    CoGe::Accessory::IRODS::irods_iput($file, $ifile);
    CoGe::Accessory::IRODS::irods_imeta_add($ifile, $meta);

    return 0;
}

sub annotate { #TODO create API "annotate" job instead
    my %args = @_;
    my $gid = $args{gid};
    
    # Get genome and check permission
    my $genome = $DB->resultset('Genome')->find($gid);
    return $ERROR unless $USER->has_access_to_genome($genome);

    # Create workflow
    my $workflow = $JEX->create_workflow( name => "Annotate genome", init => 1 );
    my ($staging_dir) = get_workflow_paths($USER->name, $workflow->id);
    mkpath($staging_dir);
    
    my $fasta_file = get_genome_file($gid);
    my $staged_fasta_file = catfile($staging_dir, basename($fasta_file));
    copy($fasta_file, $staged_fasta_file);
    
    # Add TransDecoder tasks
    my $task = create_transdecoder_longorfs_job($staged_fasta_file);
    $workflow->add_job($task);
    my $done_file = $task->{outputs}->[1];
    
    $task = create_transdecoder_predict_job($staged_fasta_file, $staging_dir, $done_file);
    $workflow->add_job($task);
    my $gff_file = $task->{outputs}->[0];
    
    $task = create_load_annotation_job(
        user => $USER,
        staging_dir => $staging_dir,
        wid => $workflow->id,
        gid => $gid,
        input_file => $gff_file,
        metadata => {
            name => 'Putative features',
            description => 'Generated by CoGe using TransDecoder',
            version => '1',
            source => 'CoGe'
        }
    );
    $workflow->add_job($task);

    # Submit workflow and block until completion
    my $response = $JEX->submit_workflow($workflow);
    say STDERR "GenomeInfo::annotate: wid=" . $response->{id};
    $JEX->wait_for_completion($response->{id});

    return encode_json({
        status => 'success'
    });
}

sub recommend_certification {
    my %args = @_;
    my $gid = $args{gid};
    my $reason = $args{reason};

    my $body = 'User ' . $USER->name . ' (' . $USER->email . ') has recommended genome <a href="https://genomevolution.org/coge/GenomeInfo.pl?gid=' . $gid . '">' . $gid . '</a> for certification in CoGe.<br><br>' . $reason;
    my $email = $config->{SUPPORT_EMAIL};
    CoGe::Accessory::Web::send_email(
        from => $email,
        to => $email,
        subject => "genome certification recommendation\nContent-Type: text/html; charset=ISO-8859-1",
        body => '<html><head></head><body>' . $body . '</body></html>'
    );
}

sub generate_html {
    my $template;

    $EMBED = $FORM->param('embed') || 0;
    if ($EMBED) {
        $template =
          HTML::Template->new(
            filename => $config->{TMPLDIR} . 'embedded_page.tmpl' );
    }
    else {
        $template =
          HTML::Template->new( filename => $config->{TMPLDIR} . 'generic_page.tmpl' );
        $template->param(
            PAGE_TITLE => $PAGE_TITLE,
	        TITLE      => 'GenomeInfo',
            PAGE_LINK  => $LINK,
	        HOME       => $config->{SERVER},
	        HELP       => 'GenomeInfo',
	        WIKI_URL   => $config->{WIKI_URL} || '',
            USER       => $USER->display_name || '',
            LOGON      => ( $USER->user_name ne "public" ),
            ADMIN_ONLY => $USER->is_admin,
            CAS_URL    => $config->{CAS_URL} || '',
            COOKIE_NAME => $config->{COOKIE_NAME} || ''
        );
    }

    $template->param( BODY => generate_body() );

    return $template->output;
}

sub generate_body {
    my $template = HTML::Template->new( filename => $config->{TMPLDIR} . $PAGE_TITLE . '.tmpl' );
    $template->param( MAIN => 1, PAGE_NAME => "$PAGE_TITLE.pl" );

    my $gid = $FORM->param('gid');
    return "No genome specified" unless $gid;

    my $genome = $DB->resultset('Genome')->find($gid);
    return "Genome id $gid not found" unless ($genome);
    return "Access denied" unless $USER->has_access_to_genome($genome);
    return "There was an error while loading this genome" if $genome->status && $genome->status == ERROR;
    return "This genome is still being loaded" if $genome->status && $genome->status == LOADING;

    $template->param( OID => $genome->organism->id );

    my $title = $genome->info;
    $title = substr($title, 0, MAX_TITLE_LENGTH) . '...' if (length($title) > MAX_TITLE_LENGTH);
    
    my $favorites = CoGe::Core::Favorites->new(user => $USER);

    $template->param(
        EMBED           => $EMBED,
        LOAD_ID         => $LOAD_ID,
        JOB_ID          => $JOB_ID,
        GID             => $gid,
        GENOME_TITLE    => $title,
        GENOME_INFO     => get_genome_info( genome => $genome ) || undef,
        GENOME_DATA     => get_genome_info_details( dsgid => $genome->id) || undef,
        LOGON           => ( $USER->user_name ne "public" ),
        DEFAULT_TYPE    => 'note', # default annotation type
        USER_CAN_EDIT   => $genome->is_editable($USER),
        USER_CAN_ADD    => ( (!$genome->restricted and !$USER->is_public) or $genome->is_editable($USER) ),
        USER_CAN_DELETE => $genome->is_deletable($USER),
        DELETED         => $genome->deleted,
        RESTRICTED      => $genome->restricted,
        CERTIFIED       => $genome->certified,
        FAVORITED       => int($favorites->is_favorite($genome)),
        TRANSCRIPTOME   => ($genome->type->name =~ /transcriptome/i) ? 1 : 0,
        IRODS_HOME      => get_irods_path(),
        USER            => $USER->user_name,
        API_BASE_URL => 'api/v1/',
        DOWNLOAD_URL    => url_for(api_url_for("genomes/$gid/sequence")) #$config->{SERVER}."api/v1/legacy/sequence/$gid" # mdb changed 2/12/16 for hypnotoad
    );

    if ( $USER->is_admin ) {
        $template->param( ADMIN_USER => 1 );
    }
    if ( $USER->is_poweruser ) {
        $template->param( POWER_USER => 1 );
    }

    return $template->output;
}
