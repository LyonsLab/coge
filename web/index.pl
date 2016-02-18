#!/usr/bin/perl -w

use AuthCAS;
use strict;
use CGI;
use CGI::Cookie;
use CGI::Carp 'fatalsToBrowser';
use HTML::Template;
use Data::Dumper;
use CGI::Log;
use CoGeX;
use CoGe::Accessory::Web qw(url_for);
use CoGe::Accessory::Utils qw( units commify sanitize_name );
use JSON qw(encode_json);
use POSIX 'ceil';

no warnings 'redefine';
use vars qw($CONF $USER $FORM $DB $LINK);

$FORM = new CGI;
( $DB, $USER, $CONF, $LINK ) = CoGe::Accessory::Web->init( cgi => $FORM );

# Logout is only called through this program!  All logouts from other pages are redirected to this page.
# mdb changed 2/24/14, issue 329 - added confirmation for CoGe-only or CyVerse-all logout
if ($FORM->param('logout_coge')) {
    CoGe::Accessory::Web->logout_coge(
        cookie_name => $CONF->{COOKIE_NAME},
        coge        => $DB,
        user        => $USER,
        form        => $FORM,
        url         => url_for('index.pl')
    );
}
elsif ($FORM->param('logout_all')) {
    CoGe::Accessory::Web->logout_cas(
        cookie_name => $CONF->{COOKIE_NAME},
        coge        => $DB,
        user        => $USER,
        form        => $FORM,
        url         => url_for('index.pl')
    );
}

my %FUNCTION = ( get_latest_genomes => \&get_latest_genomes );

CoGe::Accessory::Web->dispatch( $FORM, \%FUNCTION, \&generate_html );

sub generate_html {
    my $template = HTML::Template->new( filename => $CONF->{TMPLDIR} . 'generic_page.tmpl' );
    $template->param(
        TITLE      => 'Accelerating <span style="color: #119911">Co</span>mparative <span style="color: #119911">Ge</span>nomics',
        PAGE_TITLE => 'Comparative Genomics',
        PAGE_LINK  => $LINK,
        HOME       => $CONF->{SERVER},
        HELP       => '',
        WIKI_URL   => $CONF->{WIKI_URL} || '',
        USER       => $USER->display_name || undef,
        BODY       => generate_body(),
        ADMIN_ONLY => $USER->is_admin,
        CAS_URL    => $CONF->{CAS_URL} || ''
    );

    $template->param( LOGON => 1 ) unless $USER->user_name eq "public";
    $template->param( ADMIN_ONLY => $USER->is_admin );

    return $template->output;
}

sub generate_body {
    my $tmpl = HTML::Template->new( filename => $CONF->{TMPLDIR} . 'index.tmpl' );

    $tmpl->param(
        ACTIONS => [
            map {
                {
                    ACTION => $_->{NAME},
                    DESC   => $_->{DESC},
                    LINK   => $_->{LINK},
                    LOGO   => $_->{LOGO}
                }
            } sort { $a->{ID} <=> $b->{ID} } @{ actions() }
        ]
    );
    $tmpl->param(
        INTRO => 1,
        ORG_COUNT => commify( $DB->resultset('Organism')->count() ),
        GEN_COUNT => commify( $DB->resultset('Genome')->search( { deleted => 0 } )->count() ),
        FEAT_COUNT => commify( $DB->resultset('Feature')->count() ),
        ANNOT_COUNT => commify( $DB->resultset('FeatureAnnotation')->count() ),
        EXP_COUNT => commify( $DB->resultset('Experiment')->search( { deleted => 0 } )->count() ),
        QUANT_COUNT => commify(
            units(
                $DB->resultset('Experiment')->search( { deleted => 0 } )->get_column('row_count')->sum
            )
        )
    );

    $tmpl->param( wikifeed => $CONF->{WIKI_URL}."/CoGepedia:Current_events" ) if $CONF->{WIKI_URL};

    return $tmpl->output;
}

sub actions {
    my @actions = (
        {
            ID         => 5,
            LOGO       => "picts/GEvo.svg",
            ACTION     => qq{<a href="./GEvo.pl">GEvo</a>},
            LINK       => qq{./GEvo.pl},
            SCREENSHOT => "picts/preview/GEvo.png",
            NAME       => "GEvo",
            DESC       => qq{Compare sequences and genomic regions to discover patterns of genome evolution.<br><a href ="GEvo.pl?prog=blastz;accn1=at1g07300;fid1=4091274;dsid1=556;chr1=1;dr1up=20000;dr1down=20000;gbstart1=1;gblength1=772;accn2=at2g29640;fid2=4113333;dsid2=557;chr2=2;dr2up=20000;dr2down=20000;gbstart2=1;rev2=1;num_seqs=2;autogo=1">Example</a>},
        },
        {
            NAME       => "OrganismView",
            ID         => 1,
            LOGO       => "picts/OrganismView.svg",
            ACTION     => qq{<a href="./OrganismView.pl">OrganismView</a>},
            LINK       => qq{./OrganismView.pl},
            SCREENSHOT => "picts/preview/OrganismView.png",
            DESC       => qq{Search for organisms, get an overview of their genomic make-up, and visualize them using a dynamic, interactive genome browser.<br><a href="OrganismView.pl?org_name=W3110">Example</a>},
        },
        {
            NAME       => "CoGeBlast",
            ID         => 2,
            LOGO       => "picts/CoGeBlast.svg",
            ACTION     => qq{<a href="./CoGeBlast.pl">CoGeBlast</a>},
            LINK       => qq{./CoGeBlast.pl},
            SCREENSHOT => "picts/preview/Blast.png",
            DESC       => qq{Blast sequences against any number of organisms in CoGe.<br><a href="CoGeBlast.pl?dsgid=3068;fid=40603528">Example</a>},
        },
        {
            NAME       => "SynMap",
            ID         => 3,
            LOGO       => "picts/SynMap-inverse.svg",
            ACTION     => qq{<a href="./SynMap.pl">SynMap</a>},
            LINK       => qq{./SynMap.pl},
            SCREENSHOT => "picts/preview/SynMap.png",
            DESC       => qq{Compare any two genomes to identify regions of synteny.<br><a href="SynMap.pl?dsgid1=3068;dsgid2=8;D=20;g=10;A=5;w=0;b=1;ft1=1;ft2=1;dt=geneorder;ks=1;autogo=1">Example</a>},
        },
        {
            NAME       => "SynFind",
            ID         => 4,
            LOGO       => "picts/SynFind.svg",
            ACTION     => qq{<a href="./SynFind.pl">SynFind</a>},
            LINK       => qq{./SynFind.pl},
            SCREENSHOT => "picts/preview/SynMap.png",
            DESC       => qq{Search CoGe's annotation database for homologs.<br><a href="SynFind.pl?dsgid=3068;fid=40603528;run=1">Example</a>},
        },
    	{
    	    NAME       => "FeatView",
    	    ID         => 6,
    	    LOGO       => "picts/FeatView.svg",
    	    ACTION     => qq{<a href="./FeatView.pl">FeatView</a>},
    	    LINK       => qq{./FeatView.pl},
    	    SCREENSHOT => "picts/preview/FeatView.png",
    	    DESC       => qq{Search for a gene by name across all genomes in CoGe.<br><a href="FeatView.pl?fid=306206343&gstid=1">Example</a>},
    	},
    );
    return \@actions;
}

sub get_latest_genomes {
    my %opts = @_;
    my @latest = $DB->resultset("Genome")->get_recent_public($opts{limit});
    my @genomes;

    foreach my $dsg (@latest) {
        my $name = $dsg->organism->name;
        my $orgview_link = "GenomeInfo.pl?gid=" . $dsg->id;

        my @datetime = split " ", $dsg->get_date;

        push @genomes, { organism => $name, added => $datetime[0], url => $orgview_link };
    }

    return encode_json(\@genomes);
}

sub get_latest_old {
    my %opts = @_;
    my $limit = $opts{limit} || 20;

    my @db = $DB->resultset("Genome")->search(
        {},
        {
            distinct => "organism.name",
            join     => "organism",
            prefetch => "organism",
            order_by => "genome_id desc",
            rows     => $limit * 10,
        }
    );

    my $html = "<table class='small'>";
    $html .= "<tr><th>"
      . join( "<th>", qw( Organism &nbsp Length&nbsp(nt) &nbsp Related Link ) );
    my @opts;
    my %org_names;
    my $genome_count = 0;
    foreach my $dsg (@db) {
        next unless $USER->has_access_to_genome($dsg);
        next if $org_names{ $dsg->organism->name };
        last if $genome_count >= $limit;
        $org_names{ $dsg->organism->name } = 1;
        my $orgview_link = "OrganismView.pl?oid=" . $dsg->organism->id;
        my $entry        = qq{<tr>};
        $entry .= qq{<td><span class="link" onclick=window.open('$orgview_link')>};
        my $name = $dsg->organism->name;
        $name = substr( $name, 0, 40 ) . "..." if length($name) > 40;
        $entry .= $name;
        $entry .= qq{</span>};
        $entry .= "<td>(v" . $dsg->version . ")&nbsp";
        $entry .= "<td align=right>" . commify( $dsg->length ) . "<td>";
        my @desc = split( /;/, $dsg->organism->description );
        while ( $desc[0] && !$desc[-1] ) { pop @desc; }
        $desc[-1] =~ s/^\s+// if $desc[-1];
        $desc[-1] =~ s/\s+$// if $desc[-1];
        my $orgview_search = "OrganismView.pl?org_desc=" . $desc[-1];
        $entry .=
qq{<td><span class="link" onclick="window.open('$orgview_search')">Search</span>};
        $entry .= qq{<td>};
        $entry .=
qq{<img onClick="window.open('$orgview_link')" src="picts/other/CoGe-icon.png" title="CoGe" class="link">};

        my $search_term = $dsg->organism->name;
        $entry .=
qq{<img onclick="window.open('http://www.ncbi.nlm.nih.gov/taxonomy?term=$search_term')" src="picts/other/NCBI-icon.png" title="NCBI" class="link">};
        $entry .=
qq{<img onclick="window.open('http://en.wikipedia.org/w/index.php?title=Special%3ASearch&search=$search_term')" src="picts/other/wikipedia-icon.png" title="Wikipedia" class="link">};
        $search_term =~ s/\s+/\+/g;
        $entry .=
qq{<img onclick="window.open('http://www.google.com/search?q=$search_term')" src="picts/other/google-icon.png" title="Google" class="link">};
        $entry .= qq{</tr>};
        push @opts, $entry
          ; #, "<OPTION value=\"".$item->organism->id."\">".$date." ".$item->organism->name." (id".$item->organism->id.") "."</OPTION>";
        $genome_count++;
    }
    $html .= join "\n", @opts;
    $html .= "</table>";
    return $html;
}
