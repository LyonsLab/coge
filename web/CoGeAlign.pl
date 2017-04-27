#! /usr/bin/perl -w
use strict;
use CGI;
use CGI::Ajax;
use CoGeX;
use CoGeX::Result::Feature;
use Data::Dumper;
use CoGe::Accessory::Web;
use HTML::Template;
use URI::Escape;
use Spreadsheet::WriteExcel;
use LWP::UserAgent;

#use LWP::Simple;
#use LWP::Simple::Post qw(post post_xml);
use CoGe::Accessory::genetic_code;
use File::Path;
no warnings 'redefine';

use vars
  qw($P $DBNAME $DBHOST $DBPORT $DBUSER $DBPASS $connstr $TEMPDIR $TEMPURL $USER $DATE $CLUSTAL $BASEFILE $coge $cogeweb $FORM $NEWICKTOPS $CONVERT $MAX_PROC $COOKIE_NAME);
$P            = CoGe::Accessory::Web::get_defaults();
$ENV{PATH}    = $P->{COGEDIR};
$MAX_PROC     = $P->{MAX_PROC};
$ENV{THREADS} = $MAX_PROC;

$DATE = sprintf(
    "%04d-%02d-%02d %02d:%02d:%02d",
    sub { ( $_[5] + 1900, $_[4] + 1, $_[3] ), $_[2], $_[1], $_[0] }
      ->(localtime)
);

$TEMPDIR = $P->{TEMPDIR} . "CoGeAlign";
$TEMPURL = $P->{TEMPURL} . "CoGeAlign";
mkpath( $TEMPDIR, 0, 0777 ) unless -d $TEMPDIR;

$FORM   = new CGI;
$DBNAME = $P->{DBNAME};
$DBHOST = $P->{DBHOST};
$DBPORT = $P->{DBPORT};
$DBUSER = $P->{DBUSER};
$DBPASS = $P->{DBPASS};
$connstr =
  "dbi:mysql:dbname=" . $DBNAME . ";host=" . $DBHOST . ";port=" . $DBPORT;
$coge = CoGeX->connect( $connstr, $DBUSER, $DBPASS );

$COOKIE_NAME = $P->{COOKIE_NAME};

my ($cas_ticket) = $FORM->param('ticket');
$USER = undef;
($USER) = CoGe::Accessory::Web->login_cas(
    cookie_name => $COOKIE_NAME,
    ticket      => $cas_ticket,
    coge        => $coge,
    this_url    => $FORM->url()
) if ($cas_ticket);
($USER) = CoGe::Accessory::Web->get_user(
    cookie_name => $COOKIE_NAME,
    coge        => $coge
) unless $USER;

$CLUSTAL    = get_command_path('CLUSTALW', 'clustalw2');
$NEWICKTOPS = get_command_path('NEWICKTOPS'); # from njplot package
$CONVERT    = get_command_path('CONVERT'); # from ImageMagic package

#$CLUSTAL = "/usr/bin/clustalw-mtv";

my $pj = new CGI::Ajax(
    gen_html          => \&gen_html,
    refresh_seq       => \&refresh_seq,
    create_tree_image => \&create_tree_image,
    loading           => \&loading,
    run               => \&run,
);
$pj->js_encode_function('escape');
print $pj->build_html( $FORM, \&gen_html );

#print $FORM->header; gen_html();

sub gen_html {
    my $html;
    my ($body) = gen_body();
    my $template =
      HTML::Template->new( filename => $P->{TMPLDIR} . 'generic_page.tmpl' );
    $template->param( TITLE      => 'ClustalW2 Alignments',
                      PAGE_TITLE => 'Align',
                      SUPPORT_EMAIL => $CONF->{SUPPORT_EMAIL},
                      HOME       => $P->{SERVER},
                      HELP       => 'CoGeAlign',
                      WIKI_URL   => $P->{WIKI_URL} || '',
                      USER       => $USER->display_name || '' );
    $template->param( LOGON      => 1 ) unless $USER->user_name eq "public";
    $template->param( DATE       => $DATE );
    $template->param( BODY       => $body );
    $template->param( ADMIN_ONLY => $USER->is_admin );
    $template->param( CAS_URL    => $P->{CAS_URL} || '' );
    my $prebox =
      HTML::Template->new( filename => $P->{TMPLDIR} . 'CoGeAlign.tmpl' );
    $prebox->param( RESULTS_DIV => 1 );
    $template->param( PREBOX => $prebox->output );
    $html .= $template->output;
}

sub gen_body {
    my $template =
      HTML::Template->new( filename => $P->{TMPLDIR} . 'CoGeAlign.tmpl' );
    $template->param( JAVASCRIPT => 1 );
    $template->param( MAIN       => 1 );
    my $form = $FORM;
    my $no_values;
    $BASEFILE = $form->param('basename');
    my $feat_list = [];
    $feat_list = read_file() if $BASEFILE;    #: $opts{feature_list};

    foreach my $item ( $form->param('fid') ) {
        push @$feat_list, $item;              # if $item =~ /^\d+$/;
    }

    my $dsid = $form->param('dsid') if $form->param('dsid');
    my $chr  = $form->param('chr')  if $form->param('chr');
    my $align  = $form->param('align');
    my $autogo = $form->param('autogo');
    $autogo = 0 unless $autogo;
    if ( $align eq "codon" ) {
        $template->param( 'CODON_ALIGN_YES' => "selected" );
        $template->param( 'PROTEIN_ALIGN'   => "checked" );
    }
    elsif ( $align eq "protein" ) {
        $template->param( 'PROTEIN_ALIGN'  => "checked" );
        $template->param( 'CODON_ALIGN_NO' => "selected" );
    }
    else {
        $template->param( 'DNA_ALIGN'      => "checked" );
        $template->param( 'CODON_ALIGN_NO' => "selected" );
    }
    $template->param( 'AUTOGO' => $autogo );
    push @$feat_list, @{ get_fids_from_dataset( dsid => $dsid, chr => $chr ) }
      if $dsid;

    my $seqs = generate_sequence( feature_list => $feat_list );
    my $fid_string;
    foreach my $fid (@$feat_list) {
        $fid_string .= $fid . ":";
    }
    $fid_string =~ s/:$//;

    $seqs =~ s/(\sGenomic.+-\d+)?//g;
    if ($seqs) {
        my $num_seqs = $seqs =~ tr/>/>/;
        $template->param( COGE_SEQS => 1 );
        $template->param( SEQUENCE  => $seqs );
        $template->param( FIDS      => $fid_string );
    }
    else {
        $template->param( SEQUENCE => "Enter Fasta Sequences here" );
    }

    return $template->output;
}

sub generate_sequence {
    my %opts        = @_;
    my $featid_list = $opts{feature_list};
    my $prot        = $opts{prot} || 0;
    return unless @$featid_list;
    my $seqs;
    foreach my $featid (@$featid_list) {
        my ( $fid, $gstidt );
        if ( $featid =~ /_/ ) {
            ( $fid, $gstidt ) = split /_/, $featid;
        }
        else {
            $fid = $featid;
        }
        my ($feat) = $coge->resultset('Feature')->find($fid);
        next unless $feat;
        my $seq = $feat->fasta(
            col       => 0,
            prot      => $prot,
            name_only => 1,
            gstid     => $gstidt
        );
        $seq =~ s/>/>fid:$fid /g;
        $seqs .= $seq;
    }
    return $seqs;
}

sub loading {
    my $message = shift || "Generating results. . .";
    return qq{<div class="dna"><div id="loading">$message</div></div>};
}

sub refresh_seq {
    my %opts     = @_;
    my $fids     = $opts{fids};
    my $seq_type = $opts{seq_type};
    my $protein  = $seq_type =~ /dna/i ? 0 : 1;
    my $seqs;

    foreach my $fid ( split /:/, $fids ) {
        my ( $fidt, $gstidt );
        if ( $fid =~ /_/ ) {
            ( $fidt, $gstidt ) = split /_/, $fid;
        }
        else {
            $fidt = $fid;
        }
        my ($feat) = $coge->resultset("Feature")->find($fidt);
        next unless $feat;
        my $seq = $feat->fasta(
            col       => 0,
            prot      => $protein,
            name_only => 1,
            gstid     => $gstidt
        );
        $seq =~ s/>/>fid:$fid /g;
        $seqs .= $seq;

        #seqs .= $feat->fasta(col=>0,prot=>$protein, name_only=>1, add_fid=>1);
    }
    $seqs =~ s/(\sGenomic.+-\d+)?//g;
    return $seqs;
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

sub run {
    my %opts   = @_;
    my $inseqs = $opts{seq};

    my $seq_type   = $opts{seq_type};
    my $matrix     = $opts{matrix};
    my $gap_open   = $opts{gap_open};
    my $gap_ext    = $opts{gap_ext};
    my $iteration  = $opts{iteration};
    my $format     = "clustal";           #$opts{format}; #set by default now
    my $gen_matrix = 1;                   #$opts{gen_matrix}; #on by default now
    my $codon      = $opts{codon_align};
    $seq_type = $seq_type =~ /prot/i ? "PROTEIN" : "DNA";

    #my $num_interations = $opts{num_iterations};

    my $num_seqs = $inseqs =~ tr/>/>/;

    my %file_format =
      ( NEXUS => 'nxs', PHYLIP => 'phy', GDE => 'gde', PIR => 'pir' );

    $cogeweb = CoGe::Accessory::Web::initialize_basefile( tempdir => $TEMPDIR );
    my $seq_file = $cogeweb->basefile . "_clustalw.infile";

##Check to make sure no spaces in fasta header
    # $inseqs =~ s/\s+/_/g;

    open( NEW, "> $seq_file" );
    print NEW $inseqs;
    close NEW;
    my $suffix = $format =~ /(jalview|clustal)/ ? 'aln' : $file_format{$format};
    my $outfile     = $cogeweb->basefile . "_clustalw." . $suffix;
    my $tree_out    = $TEMPURL . "/" . $cogeweb->basefilename . "_clustalw.dnd";
    my $phylip_file = $TEMPURL . "/" . $cogeweb->basefilename . "_clustalw.ph";
    my $pre_command = "-INFILE=$seq_file";

    $pre_command .= " -TYPE=$seq_type";
    if ( $format && $format !~ /jalview/i && $format !~ /clustal/ii ) {
        $pre_command .= " -OUTPUT=$format";
    }
    $pre_command .=
      $seq_type =~ /dna/i ? " -DNAMATRIX=$matrix" : " -MATRIX=$matrix";
    $pre_command .= " -GAPOPEN=$gap_open";
    $pre_command .= " -GAPEXT=$gap_ext";
    $pre_command .= " -ITERATION=$iteration";
    $pre_command .= " -OUTPUTTREE=phylip";

    #$pre_command .= " -NUMITER=$num_interations";

    my $x;
    ( $x, $pre_command ) = CoGe::Accessory::Web::check_taint($pre_command);

    my $command = "$CLUSTAL $pre_command";

    print STDERR $command, "\n";

    `$command`;

    $command .= " -tree";
    `$command`;

    my $output;

    open( IN, $outfile ) || die "$! can't open $outfile for reading";
    while (<IN>) {

        $output .= $_;
    }
    close IN;

    my $outfile_jalview = $outfile;
    $outfile_jalview =~ s/$TEMPDIR/$TEMPURL/;
    my $name_conversion =
      convert_phylip_names( file => $phylip_file, seqs => $inseqs );
    my ( $header_html, $seq_html, $seqs, $codon_alignment, $clustal_alignment,
        $name_order )
      = parse_results(
        clustal      => $output,
        num_seqs     => $num_seqs,
        codon_align  => $codon,
        seq_type     => $seq_type,
        name_convert => $name_conversion
      );
    my $box_template =
      HTML::Template->new( filename => $P->{TMPLDIR} . 'box.tmpl' );

    #print STDERR $html,"\n";
    my $html;
    my $phylip_file_jalview = $phylip_file;
    $phylip_file_jalview =~ s/http.+edu//;

    #jalview screws up jquery for some reason
    $html .=
qq{<applet width="140" height="35" code="jalview.bin.JalviewLite" archive="/CoGe/bin/JalView/jalviewApplet.jar"><param name="file" value="$outfile_jalview"><param name="tree" value="$phylip_file_jalview"><param name="showbutton" value="true"><param name="defaultColour" value="Clustal"></applet><br><a class=small href="http://www.jalview.org" target=_new>Information on JalView</a>};
    $clustal_alignment =~ s/\n/<br\/>/g;
    my $tree = create_tree_image($phylip_file);
    $html .= qq{
<br>
<span style="display:none" id=show_tree class='ui-button ui-corner-all' onClick="\$('#hide_tree').show();\$('#show_tree').hide();\$('#tree_box').show();">Show Tree</span>
<span id=hide_tree class='ui-button ui-corner-all' onClick="\$('#hide_tree').hide();\$('#show_tree').show();\$('#tree_box').hide();">Hide Tree</span>
};

    $html .= "<div id=tree_box>";
    $html .=
qq{<span class='ui-button ui-corner-all ' onclick="\$('#select_feats').dialog('open')">Open Feature Selection Box</span><br>}
      if keys %$name_conversion;
    $html .= "<img src=$tree></div>";

    #let's make a series of checkboxes to send sequence to featlist if possible
    if ( keys %$name_conversion ) {
        my @select_feats;
        foreach my $id (@$name_order) {
            my $name = $name_conversion->{$id};
            next unless $name;
            $id =~ s/fid_//;
            push @select_feats, { ID => $id, NAME => $name . " ($id)" };
        }
        if (@select_feats) {
            my $template =
              HTML::Template->new(
                filename => $P->{TMPLDIR} . 'CoGeAlign.tmpl' );
            $template->param( FEATURE_SELECT => 1 );
            $template->param( FEATS          => \@select_feats );
            $html .= $template->output;
        }
    }

    $html .= qq{
<br>
<span style="display:none" id=show_alignment class='ui-button ui-corner-all' onClick="\$('#hide_alignment').show();\$('#show_alignment').hide();\$('#alignment_box').show();">Show Alignment</span>
<span  id=hide_alignment class='ui-button ui-corner-all' onClick="\$('#hide_alignment').hide();\$('#show_alignment').show();\$('#alignment_box').hide();">Hide Alignment</span>
};
    $html .= qq{<div id=alignment_box align=left class=resultborder "><br/>
<table><tr valign=top>
<td><pre>$header_html</pre>
<td><pre  style= "overflow:auto;width:700px;max-height:700px;">$seq_html</pre>
</table>
</div>};
    $html .= qq{
<br>
<span id=show_clustalw_alignment class='ui-button ui-corner-all' onClick="\$('#hide_clustalw_alignment').show();\$('#show_clustalw_alignment').hide();\$('#clustalw_alignment_box').show();">Show ClustalW Alignment</span>
<span  style="display:none" id=hide_clustalw_alignment class='ui-button ui-corner-all' onClick="\$('#hide_clustalw_alignment').hide();\$('#show_clustalw_alignment').show();\$('#clustalw_alignment_box').hide();">Hide ClustalW Alignment</span>
};

    $html .=
qq{<div id=clustalw_alignment_box align=left class=resultborder style="display:none;overflow:auto;width:700px;max-height:700px;"><br/><pre>$clustal_alignment</pre></div>};

    $html .= qq{
<br>
<span id=show_fasta_alignment class='ui-button ui-corner-all' onClick="\$('#hide_fasta_alignment').show();\$('#show_fasta_alignment').hide();\$('#fasta_alignment_box').show();">Show Fasta Alignment</span>
<span  style="display:none" id=hide_fasta_alignment class='ui-button ui-corner-all' onClick="\$('#hide_fasta_alignment').hide();\$('#show_fasta_alignment').show();\$('#fasta_alignment_box').hide();">Hide Fasta Alignment</span>
};

    my $fasta;
    foreach my $name ( keys %$seqs ) {
        next if $name eq "alignment_stuff";
        $fasta .= ">" . $name . "\n";
        $fasta .= $seqs->{$name} . "\n";
    }
    $html .=
qq{<div id=fasta_alignment_box align=left class=resultborder style="display:none;overflow:auto;width:700px;max-height:700px;"><br/><pre>$fasta</pre></div>};

#  my $select_menu = qq{<select id=file_types><option value=GCG>GCG<option value=GDE>GDE<option value=PHYLIP>Phylip<option value=PIR>PIR<option value=NEXUS>NEXUS</select>};

#</tr><tr><td>Select Alt. Output File for Download:</tr><tr><td>$select_menu<input type=button value="Go">
#    $html .= "</table>";
    if ($gen_matrix) {
        $html .= qq{
<br>
<span id=show_matrix class='ui-button ui-corner-all' onClick="\$('#hide_matrix').show();\$('#show_matrix').hide();\$('#matrix_box').show();">Show Scoring Matrix</span>
<span  style="display:none" id=hide_matrix class='ui-button ui-corner-all' onClick="\$('#hide_matrix').hide();\$('#show_matrix').show();\$('#matrix_box').hide();">Hide Scoring Matrix</span>
<div id="matrix_box" style="display:none" >
<span class="small">(Using the method describe in Henikoff and Henikoff, 1992)</span>
};

        if ($codon) {
            my $matrix = gen_matrix( seqs => $seqs, type => 'p' );
            $html .= $matrix if $matrix;
            $matrix = gen_matrix( seqs => $codon_alignment, type => 'c' );
            $html .= "<br>" . $matrix if $matrix;
        }
        else {
            my $matrix = gen_matrix( seqs => $seqs, type => $seq_type );
            $html .= $matrix if $matrix;
        }
        $html .= "</div>";
    }

    $html .= qq{<div class=small>Downloads:
<br><a href="$outfile" target="_blank">ClustalW Alignment File</a>
<br><a href="$tree_out" target="_blank">ClustalW Tree File</a>
<br><a href="$phylip_file" target="_blank">ClustalW Phylip File</a>
</div>
};

    $box_template->param( BODY       => $html );
    my $outhtml = $box_template->output;
    return $outhtml;
}

sub parse_results {
    my %opts        = @_;
    my $clustal     = $opts{clustal};
    my $num_seqs    = $opts{num_seqs};
    my $codon_align = $opts{codon_align};
    my $matrix      = $opts{gen_matrix};
    my $seq_type    = $opts{seq_type};
    my $names       = $opts{name_convert};
    my @lines;
    my @headers;
    my @seqs;
    my $clustal_title;
    my %seqs;
    my $pad = 0;

    foreach my $line ( split /\n/, $clustal ) {
        next unless $line;
        next if $line =~ /^CLUSTAL/;
        if ( $line =~ /^\s/ ) {
            $line = substr( $line, $pad );
            $seqs{alignment_stuff} .= $line;
            next;
        }
        unless ($pad) {
            my ($tmp) = $line =~ /(^\S+\s+)/;
            $pad = length $tmp;
        }
        my ( $name, $seq ) = split /\s+/, $line, 2;

        $seq = uc($seq);
        push @headers, $name unless $seqs{$name};
        $seqs{$name} .= $seq;
    }

    my ( $header_output, $seq_output );
    my %codon_alignment;
    foreach my $name (@headers) {
        my $seq = uc( $seqs{$name} );
        if ($codon_align) {
            my ($fid) = $name =~ /fid_(\d+)/;
            if ($fid) {
                my $feat = $coge->resultset('Feature')->find($fid);
                my $tseq = $feat->genomic_sequence;
                my $pos  = 0;
                my $dna;
                my $prot;
                foreach my $chr ( split //, $seqs{$name} ) {
                    if ( $chr eq "-" ) {
                        $dna .= "-" x 3 . " ";
                    }
                    else {
                        $dna .= substr( $tseq, $pos, 3 ) . " ";
                        $pos += 3;
                    }
                    $prot .= "  $chr ";
                }
                $seq = $prot . "<br>" . $dna;
                $dna =~ s/\s//g;
                $codon_alignment{$name} = $dna;
            }

        }
        if ( !$codon_align ) {
            if ( $seq_type eq "DNA" ) {
                $seq =~
s/(T+)/<span style='background-color:deepskyblue'>$1<\/span>/g;
                $seq =~ s/(A+)/<span style='background-color:red'>$1<\/span>/g;
                $seq =~
                  s/(G+)/<span style='background-color:yellow'>$1<\/span>/g;
                $seq =~
                  s/(C+)/<span style='background-color:green'>$1<\/span>/g;
            }
            else {
                $seq =~
s/([HKRDESTNQ])/<span style='background-color:deepskyblue'>$1<\/span>/g;
                $seq =~
                  s/([CPG])/<span style='background-color:red'>$1<\/span>/g;
                $seq =~
s/([AILMFWYV])/<span style='background-color:green'>$1<\/span>/g;
            }
        }
        my $outname = $names->{$name}
          if $names->{$name};    #convert name if possible
        $header_output .= $outname . "<br/>";
        $header_output .= "<br>" if $codon_align;    #needs and extra <br>
        $seq_output    .= $seq . "<br/>";

    }
    $header_output .= "<br/>";
    $codon_alignment{alignment_stuff} =
      "  " . join( "  ", ( split //, $seqs{alignment_stuff} ) )
      if $codon_align;
    $seq_output .=
      $codon_align
      ? "<pre>" . "  "
      . join( "  ", map { $_ . " " } ( split //, $seqs{alignment_stuff} ) )
      . "</pre><br/>"
      : "<pre>" . $seqs{alignment_stuff} . "</pre><br/>";

    #convert names in clustal file
    #1 find length of logest sequence name
    my $name_len = 0;
    foreach my $name ( values %$names ) {
        $name_len = length($name) if length($name) > $name_len;
    }
    my $name_space = 0;
    while ( my ( $k, $v ) = each %$names ) {
        my $name = $v;
        $name .= " " x ( $name_len - length($v) + 1 );
        $clustal =~ s/($k\s+)/$name/g;
        $name_space = length($1) if $1 && !$name_space;
    }
    my $space = "\n" . " " x ( $name_len + 1 );
    $clustal =~ s/\n {$name_space}/$space/g;
    return $header_output, $seq_output, \%seqs, \%codon_alignment, $clustal,
      \@headers;
}

sub gen_matrix {
    my %opts = @_;
    my $seqs = $opts{seqs};
    my $type = $opts{type};
    my %sub_count;
    my $aln = $seqs->{alignment_stuff};
    delete $seqs->{alignment_stuff};
    my $length = length $aln;
    my %freq;
    my %total_chrs;

#log-odds score based on Henikoff and Henikoff 1992 paper:  Amino acid substitution matrices from protein blocks
    my $add = $type =~ /^c/ ? 3 : 1;
    for ( my $i = 0 ; $i < $length ; $i += $add ) {
        my %chrs;
        if ( $type =~ /^c/ ) {
            map { $chrs{$_}++ } map { substr( $_, $i, 3 ) } values %$seqs;
        }
        else {
            map { $chrs{$_}++ } map { substr( $_, $i, 1 ) } values %$seqs;
        }
        my @chrs = keys %chrs;
        foreach my $c1 (@chrs) {
            next if $c1 eq "-" || $c1 eq "---";    #skip gaps
                    #skip non ATCG DNA characters in codons/DNA
            next if ( $type =~ /^c/i || $type =~ /^d/i ) && $c1 =~ /[^ATCG]/;

            #calculate pair count for identical match and add to total
            my $x = $chrs{$c1} - 1;
            $freq{$c1}{$c1} += $x * ( $x + 1 ) / 2;

            #add count to total character count
            $total_chrs{$c1} += $chrs{$c1};
            foreach my $c2 (@chrs) {
                next if $c2 eq "-";
                next if $c1 eq $c2;    #already done
                      #calculate pair count for non-indentical matches
                $freq{$c1}{$c2} += $chrs{$c2} * $chrs{$c1}
                  unless $freq{$c2}{$c1};

            }
        }
    }
    my $total_pairs = 0;
    foreach my $c1 ( keys %total_chrs ) {
        foreach my $c2 ( keys %total_chrs ) {
            $total_pairs += $freq{$c1}{$c2} if $freq{$c1}{$c2};

        }
    }
    my %q;
    foreach my $c1 ( keys %total_chrs ) {
        foreach my $c2 ( keys %total_chrs ) {
            $q{$c1}{$c2} = $freq{$c1}{$c2} / $total_pairs if $freq{$c1}{$c2};
        }
    }
    my %p;
    foreach my $c1 ( keys %total_chrs ) {
        $p{$c1} = $q{$c1}{$c1};
        foreach my $c2 ( keys %total_chrs ) {
            next if $c1 eq $c2;
            $p{$c1} += $q{$c1}{$c2} if $q{$c1}{$c2};
        }
    }

    #    return "<pre>".Dumper (\%q)."</pre>";
    my %matrix;
    foreach my $c1 ( keys %total_chrs ) {
        foreach my $c2 ( keys %total_chrs ) {
            my $q = $q{$c1}{$c2} ? $q{$c1}{$c2} : $q{$c2}{$c1};
            next unless $q;
            my $denom = ( $p{$c1} * $p{$c2} );
            my $val = sprintf( "%.0f", 2 * log( $q / $denom ) ) if $denom;
            $val = 0 if $val eq "-0";
            $matrix{ uc($c1) }{ uc($c2) } = $val;
        }
    }
    $seqs->{alignment_stuff} = $aln;
    my $html = gen_matrix_output( matrix => \%matrix, type => $type );

    #    my $html = "<pre>".Dumper (\%matrix)."</pre>";
    return $html;
}

sub gen_matrix_output {
    my %opts = @_;
    my $data = $opts{matrix};
    my $type = $opts{type};
    my $code = CoGe::Accessory::genetic_code->code;
    $code = $code->{code};
    my $html = "<table>";
    if ( $type =~ /^p/i )    #proteins
    {
        my $aa_sort = CoGe::Accessory::genetic_code->sort_aa_by_gc();
        foreach my $c1 ( keys %$aa_sort ) {
            foreach my $c2 ( keys %$aa_sort ) {
                $data->{$c1}{$c2} = "-50" unless defined $data->{$c1}{$c2};
            }
        }
        $html .= "<tr><th><th>"
          . join( "<th>",
            sort { $aa_sort->{$b} <=> $aa_sort->{$a} || $a cmp $b }
              keys %$aa_sort );
        $html .= "<th>Total:";
        $html .= "<tr>";
        foreach my $aa1 (
            sort { $aa_sort->{$b} <=> $aa_sort->{$a} || $a cmp $b }
            keys %$aa_sort
          )
        {
            $html .= "<th>$aa1";
            my %vals;
            map { $vals{$_}++ } map { $data->{$aa1}{$_} } keys %$aa_sort;
            my ($max) = sort { $b <=> $a } keys %vals;
            my ( $min1, $min2 ) = sort { $a <=> $b } keys %vals;
            $min2 = $min1 unless $min2;
            my $total = 0;
            foreach my $aa2 (
                sort { $aa_sort->{$b} <=> $aa_sort->{$a} || $a cmp $b }
                keys %$aa_sort
              )
            {
                my $val = $data->{$aa1}{$aa2};
                $total += $val if $val =~ /^\d+$/;
                my $color;
                if ( $val > 0 ) {
                    $color = color_by_usage( $max, $val );
                    $html .=
                      "<td style=\"background-color: rgb($color,255,$color)\">"
                      . $val;
                }
                else {
                    $color = color_by_usage( abs($min2), $val * -1 );
                    $html .=
                      "<td style=\"background-color: rgb(255,$color,$color)\">"
                      . $val;
                }
            }
            $html .= "<td>$total<tr>";
        }
    }
    elsif ( $type =~ /^c/i )    #codons
    {
        my @dna = qw(A T C G);
        my @codons;
        foreach my $c1 (@dna) {
            foreach my $c2 (@dna) {
                foreach my $c3 (@dna) {
                    push @codons, $c1 . $c2 . $c3;
                }
            }
        }
        foreach my $c1 (@codons) {
            foreach my $c2 (@codons) {
                $data->{$c1}{$c2} = "-50" unless defined $data->{$c1}{$c2};
            }
        }
        $html .= "<tr><th><th>" . join(
            "<th>",
            map { $_ . "(" . $code->{$_} . ")" } sort {
                sort_nt1( substr( $a, 0, 1 ) )
                  <=> sort_nt1( substr( $b, 0, 1 ) )
                  || sort_nt2( substr( $a, 1, 1 ) )
                  <=> sort_nt2( substr( $b, 1, 1 ) )
                  || sort_nt3( substr( $a, 2, 1 ) )
                  <=> sort_nt3( substr( $b, 2, 1 ) )
              } @codons
        );
        $html .= "<th>Total";
        $html .= "<tr>";
        foreach my $aa1 (
            sort {
                sort_nt1( substr( $a, 0, 1 ) )
                  <=> sort_nt1( substr( $b, 0, 1 ) )
                  || sort_nt2( substr( $a, 1, 1 ) )
                  <=> sort_nt2( substr( $b, 1, 1 ) )
                  || sort_nt3( substr( $a, 2, 1 ) )
                  <=> sort_nt3( substr( $b, 2, 1 ) )
            } @codons
          )
        {
            $html .= "<th>$aa1" . "(" . $code->{$aa1} . ")";
            my %vals;
            map   { $vals{$_}++ }
              map { $data->{$aa1}{$_} } keys %{ $data->{$aa1} };
            my ($max) = sort { $b <=> $a } keys %vals;
            my ( $min1, $min2 ) = sort { $a <=> $b } keys %vals;
            $min2 = $min1 unless $min2;
            my $total = 0;

            foreach my $aa2 (
                sort {
                    sort_nt1( substr( $a, 0, 1 ) )
                      <=> sort_nt1( substr( $b, 0, 1 ) )
                      || sort_nt2( substr( $a, 1, 1 ) )
                      <=> sort_nt2( substr( $b, 1, 1 ) )
                      || sort_nt3( substr( $a, 2, 1 ) )
                      <=> sort_nt3( substr( $b, 2, 1 ) )
                } @codons
              )    #	    foreach my $aa2 (sort keys %$data)
            {
                my $val = $data->{$aa1}{$aa2};
                $total += $val;
                my $color;
                if ( $val == $min1 ) {
                    $html .=
                        "<td style=\"background-color: rgb(125,125,125)\">"
                      . $val . " "
                      . $code->{$aa1} . "-"
                      . $code->{$aa2};
                }
                elsif ( $val > 0 ) {
                    $color = color_by_usage( $max, $val );
                    $html .=
                      "<td style=\"background-color: rgb($color,255,$color)\">"
                      . $val . " "
                      . $code->{$aa1} . "-"
                      . $code->{$aa2};
                }
                else {
                    $color = color_by_usage( abs($min2), $val * -1 );
                    $html .=
                      "<td style=\"background-color: rgb(255,$color,$color)\">"
                      . $val . " "
                      . $code->{$aa1} . "-"
                      . $code->{$aa2};
                }
            }
            $html .= "<td>$total<tr>";
        }
    }
    else    #dna
    {
        my @dna = qw(A T C G);
        $html .= "<tr><th><th>" . join( "<th>", @dna );
        $html .= "<th>Total";
        $html .= "<tr>";
        foreach my $c1 (@dna) {
            $html .= "<th>$c1";
            my ($max) =
              sort { $b <=> $a }
              map  { $data->{$c1}{$_} } keys %{ $data->{$c1} };
            my ($min) =
              sort { $a <=> $b }
              map  { $data->{$c1}{$_} } keys %{ $data->{$c1} };
            my $total = 0;
            foreach my $c2 (@dna) {
                my $val = $data->{$c1}{$c2};
                $total += $val;
                my $color;
                if ( $val > 0 ) {
                    $color = color_by_usage( $max, $val );
                    $html .=
                      "<td style=\"background-color: rgb($color,255,$color)\">"
                      . $val;
                }
                else {
                    $color = color_by_usage( abs($min), $val * -1 );
                    $html .=
                      "<td style=\"background-color: rgb(255,$color,$color)\">"
                      . $val;
                }
            }
            $html .= "<td>$total<tr>";
        }
    }
    $html .= "</table>";
    return $html;
}

sub color_by_usage {
    my ( $max, $value, $opt ) = @_;
    $opt = 255 unless $opt;
    return $opt unless $max;
    my $g = $opt * ( ( $max - $value ) / $max );
    return int( $g + .5 );
}

sub sort_nt1 {
    my $chr = uc(shift);

    $chr = substr( $chr, -1, 1 ) if length($chr) > 1;
    my $val = 0;
    if ( $chr eq "G" ) {
        $val = 1;
    }
    elsif ( $chr eq "A" ) {
        $val = 2;
    }
    elsif ( $chr eq "T" ) {
        $val = 3;
    }
    return $val;
}

sub sort_nt2 {
    my $chr = uc(shift);

    $chr = substr( $chr, -1, 1 ) if length($chr) > 1;
    my $val = 0;
    if ( $chr eq "G" ) {
        $val = 1;
    }
    elsif ( $chr eq "A" ) {
        $val = 2;
    }
    elsif ( $chr eq "T" ) {
        $val = 3;
    }
    return $val;
}

sub sort_nt3 {
    my $chr = uc(shift);

    $chr = substr( $chr, -1, 1 ) if length($chr) > 1;
    my $val = 0;
    if ( $chr eq "G" ) {
        $val = 1;
    }
    elsif ( $chr eq "T" ) {
        $val = 2;
    }
    elsif ( $chr eq "C" ) {
        $val = 3;
    }
    return $val;
}

sub create_tree_image {
    my $treefile = shift;
    $treefile =~ s/$TEMPURL/$TEMPDIR/;
    my $treebase = $treefile;
    $treebase =~ s/\.ph$//;
    my $treeps  = $treebase . ".ps";
    my $treepng = $treebase . ".png";
    unless ( -e $treepng ) {
        `$NEWICKTOPS $treefile -us -notitle`;
        print STDERR "       $NEWICKTOPS $treefile -us\n";
        `$CONVERT $treeps $treepng`;
        $treepng =~ s/$TEMPDIR/$TEMPURL/;
        $treepng =~ s/CoGeAlignCoGeAlign/CoGeAlign/;
    }
    return $treepng;
}

sub convert_phylip_names {
    my %opts = @_;
    my $file = $opts{file};
    my $seqs = $opts{seqs};
    $file =~ s/$TEMPURL/$TEMPDIR/;
    my %names;
    foreach my $item ( split /\n/, $seqs ) {
        next unless $item =~ /^>(fid:\S+)/;
        my $id = $1;
        $item =~ s/$id\s*//;
        $id   =~ s/:/_/;
        $item =~ s/^>//;
        $names{$id} = $item;
    }
    my $output;
    open( IN, $file ) || warn "Can't open $file for reading: $!";
    while (<IN>) {
        $output .= $_;
    }
    close IN;
    while ( my ( $k, $v ) = each %names ) {
        $output =~ s/$k/$v/xsg;
    }
    open( OUT, ">" . $file );
    print OUT $output;
    close OUT;
    return \%names;
}
