#! /usr/bin/perl -w
use strict;
use CoGeX;
use CoGe::Accessory::Web;
use CGI;
use CGI::Ajax;
use HTML::Template;
use Text::Wrap qw($columns &wrap);
use Data::Dumper;
use URI::Escape;
use POSIX;
use Digest::MD5 qw(md5_hex);
use DBIxProfiler;
use File::Temp;
use File::Basename;
use Spreadsheet::WriteExcel;
use Mail::Mailer;
use File::Path;
no warnings 'redefine';

use vars
  qw($P $DBNAME $DBHOST $DBPORT $DBUSER $DBPASS $connstr $ARAGORN $FASTADIR $DATADIR $TEMPDIR $TEMPURL $DATE $USER $FORM $coge $cogeweb $COOKIE_NAME);

$P = CoGe::Accessory::Web::get_defaults();
$ENV{PATH} = $P->{COGEDIR};

$TEMPDIR = $P->{TEMPDIR} . "tRNA";
mkpath( $TEMPDIR, 0, 0777 ) unless -d $TEMPDIR;
$TEMPURL  = $P->{TEMPURL} . "tRNA";
$DATADIR  = $P->{DATADIR};
$FASTADIR = $P->{FASTADIR};

$ARAGORN = get_command_path('ARAGORN');

$DATE = sprintf(
    "%04d-%02d-%02d %02d:%02d:%02d",
    sub { ( $_[5] + 1900, $_[4] + 1, $_[3] ), $_[2], $_[1], $_[0] }
      ->(localtime)
);
$FORM = new CGI;
my %ajax = CoGe::Accessory::Web::ajax_func();

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

my $pj = new CGI::Ajax(
    gen_html               => \&gen_html,
    get_orgs               => \&get_orgs,
    run_aragorn            => \&run_aragorn,
    check_address_validity => \&check_address_validity,
);

$pj->js_encode_function('escape');
print $pj->build_html( $FORM, \&gen_html );
#print $FORM->header;
#print gen_html();

sub gen_html {
    my $html;
    my ($body) = gen_body();
    my $template =
      HTML::Template->new( filename => $P->{TMPLDIR} . 'generic_page.tmpl' );
    $template->param( PAGE_TITLE => "tRNAView",
                      TITLE      => 'CoGe tRNA and tmRNA Search Tool',
                      SUPPORT_EMAIL => $P->{SUPPORT_EMAIL},
                      HOME       => $P->{SERVER},
                      HELP       => 'tRNAView',
                      WIKI_URL   => $P->{WIKI_URL} || '',
                      USER       => $USER->display_name || '' );
    $template->param( LOGON => 1 ) unless $USER->user_name eq "public";
    $template->param( DATE => $DATE );
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
      HTML::Template->new( filename => $P->{TMPLDIR} . 'tRNAView.tmpl' );
    my $form   = $FORM;
    my $featid = join( ",", $form->param('featid'), $form->param('fid') ) || 0;
    my $chr    = $form->param('chr') || 0;
    my $upstream   = $form->param('upstream') || 0;
    my $downstream = $form->param('downstream') || 0;
    my $dsid       = $form->param('dsid') || 0;
    my $rc         = $form->param('rc') || 0;
    my $tab_opt    = $featid || $dsid ? 1 : 0;
    $template->param( JAVASCRIPT => 1 );
    $template->param( TAB_OPT    => $tab_opt );
    $template->param( ORG_LIST   => get_orgs() );
    $template->param(
        SEQ => get_sequence(
            fid        => $featid,
            chr        => $chr,
            start      => $upstream,
            start      => $upstream,
            stop       => $downstream,
            downstream => $downstream,
            rc         => $rc,
            dsid       => $dsid
        )
    );
    my $html = $template->output;
    return $html;
}

sub get_sequence {
    my %opts       = @_;
    my $fids       = $opts{fid};
    my $dsid       = $opts{dsid};
    my $chr        = $opts{chr};
    my $start      = $opts{start};
    my $stop       = $opts{stop};
    my $upstream   = $opts{upstream};
    my $downstream = $opts{downstream};
    my $rc         = $opts{rc};
    my $fasta;

    if ($fids) {
        foreach my $fid ( split(/,/, $fids) ) {
            my $feat = $coge->resultset('Feature')->find($fid);
            next unless $USER->has_access_to_dataset( $feat->dataset );
            $fasta .=
              ref($feat) =~ /Feature/i
              ? $feat->fasta(
                rc         => $rc,
                upstream   => $upstream,
                downstream => $downstream,
              )
              : ">Unable to retrieve Feature object for id: $fid\n";
        }
    }
    else {
        my $ds = $coge->resultset('Dataset')->find($dsid);
        $fasta =
          ref($ds) =~ /dataset/i
          ? $ds->fasta(
            start => $start,
            stop  => $stop,
            chr   => $chr,
            rc    => $rc,
          )
          : ">Unable to retrieve dataset object for id: $dsid";
    }
    return $fasta;
}

sub get_orgs {
    my $name = shift;
    my $html;
    if ($name) {
        my @db =
          $coge->resultset('Organism')
          ->search( { name => { like => "%" . $name . "%" } } );
        my @opts;
        foreach my $item ( sort { uc( $a->name ) cmp uc( $b->name ) } @db ) {
            push @opts,
                "<OPTION value=\""
              . $item->id
              . "\" id=\"o"
              . $item->id . "\">"
              . $item->name
              . "</OPTION>";
        }

        $html .=
            qq{<FONT CLASS ="small" id="org_count">Organism count: }
          . scalar @opts
          . qq{</FONT>\n<BR>\n};
        unless (@opts) {
            $html .= qq{<input type = hidden name="org_id" id="org_id"><br>};
            $html .= "No results";
            return $html;
        }

        $html .=
qq{<SELECT id="org_id" SIZE="8" MULTIPLE onClick="\$('#remove').hide(0);\$('#add').show(0);" ondblclick="get_from_id(['org_id'],[add_to_list]);">\n};
        $html .= join( "\n", @opts );
        $html .= "\n</SELECT>\n";
        $html =~ s/OPTION/OPTION SELECTED/;
    }
    else {
        $html .=
            qq{<FONT CLASS ="small" id="org_count">Organism count: }
          . $coge->resultset('Organism')->count()
          . qq{</FONT>\n<BR>\n};
    }
    return $html;
}

sub run_aragorn {
    my %opts = @_;
    #print STDERR Dumper \%opts;
    #put options here
    my $search_type = $opts{type};
    my $gcode       = $opts{gcode};
    my $dna_top     = $opts{dna_top};
    my $strand      = $opts{strand};
    my $pseudo      = $opts{pseudo};
    my $format      = $opts{format};
    my $orgid       = $opts{orgid};
    my $email       = $opts{email};
    my $seq         = $opts{seq};

    $cogeweb = CoGe::Accessory::Web::initialize_basefile( tempdir => $TEMPDIR );
    my $outfile = $cogeweb->basefile . "_aragorn.output";

    my $fasta_file =
      $seq !~ /undefined/ ? make_fasta_file($seq) : get_genome_fasta($orgid);

    my $precommand = "$search_type $gcode $dna_top $strand ";

    $precommand .= "$format " unless $format =~ /none/;
    $precommand .= "$pseudo " unless $pseudo =~ /undefined/;
    $precommand .= "-i -seq -br ";

    $precommand .= "-o $outfile $fasta_file";
    my $x;
    ( $x, $precommand ) = CoGe::Accessory::Web::check_taint($precommand);

    my $command = "$ARAGORN $precommand";

    #print STDERR $command,"\n";

    `$command`;

    my $output = read_file($outfile);

    my ( $output_parse, $excel_file ) =
      $format !~ /none/ ? parse_deliminated($output) : ( "", "" );

    $excel_file =~ s/$TEMPDIR/$TEMPURL/ if $excel_file;

    my $outfile_parse = save_file($output_parse) if $output_parse;

    my $html =
qq{<div align=center><pre><div align=left class=resultborder style="overflow:auto;width:700px;max-height:700px;"><br/>$output</div></pre>};

    $outfile =~ s/\/opt\/apache//;

    $html .= qq!Download output: <a href="$outfile">Aragorn Output File</a>!;
    $html .=
qq!<br>Download tab deliminated output: <a href="$outfile_parse">Tab Deliminated Output</a>!
      if $outfile_parse;
    $html .=
qq!<br>Download Output as Excel File: <a href="$excel_file">Excel Formatted Output</a>!
      if $excel_file;
    $html .= "</div><br><br>";

    if ( $email !~ /undefined/ ) {
        email_results(
            email  => $email,
            raw    => $outfile,
            tabbed => $outfile_parse,
            excel  => $excel_file
        );
    }

    return $html;

}

sub create_excel_file {
    my $results  = shift;
    my $filename = $cogeweb->basefile . "_aragorn.xls";
    my $workbook = Spreadsheet::WriteExcel->new("$filename");
    $workbook->set_tempdir("$TEMPDIR");
    my $worksheet = $workbook->add_worksheet();
    my $i         = 0;
    foreach my $result (@$results) {
        $worksheet->write( $i, 0, $result->{name} );
        $i++;
        $worksheet->write( $i, 0, $result->{gene_no} );
        $i++;

        $worksheet->write( $i, 0, "No." );
        $worksheet->write( $i, 1, "Type" );
        $worksheet->write( $i, 2, "Location" );
        $worksheet->write( $i, 3, "Anticodon" );
        $worksheet->write( $i, 4, "Anticodon Location" );
        $worksheet->write( $i, 5, "Intron Location" );
        $worksheet->write( $i, 6, "Sequence" );
        $i++;

        foreach my $trna ( @{ $result->{results} } ) {
            $worksheet->write( $i, 0, $trna->{num} );
            $worksheet->write( $i, 1, $trna->{type} );
            $worksheet->write( $i, 2, $trna->{loc} );
            $worksheet->write( $i, 3, $trna->{anticodon} );
            $worksheet->write( $i, 4, $trna->{antiloc} );
            $worksheet->write( $i, 5, $trna->{intron} );
            $worksheet->write( $i, 6, $trna->{prim_seq} );
            $i++;
            $worksheet->write( $i, 6, $trna->{sec_seq} );
            $i++;
        }

        $i++;
    }
    $workbook->close() or die "Error closing file: $!";
    return $filename;
}

sub parse_deliminated {
    my $output = shift;
    my @ds = split( />/, $output );
    #print STDERR Dumper \@ds;
    my @results;

    foreach my $ds (@ds) {
        next unless $ds;
        my @lines     = split( /\n/, $ds );
        my $ds_name   = $lines[0];
        my $gene_no   = $lines[1];
        my $trna_data = [];
        for ( my $i = 2 ; $i < scalar @lines ; $i += 3 ) {
            my ( $num, $type, $loc, $antiloc, $anticodon, $intron ) =
              $lines[$i] =~
/(\d+)\s+(tRNA-\w+)\s+(c?\[\d+,\d+\])\s+(\d+)\s+(\(\w{3}\))(i\(\d+,\d+\))?/;
            $intron = "none" unless $intron;
            my $prim_seq   = $lines[ $i + 1 ];
            my $second_seq = $lines[ $i + 2 ];
            push @$trna_data,
              {
                num       => $num,
                type      => $type,
                loc       => $loc,
                antiloc   => $antiloc,
                anticodon => $anticodon,
                intron    => $intron,
                prim_seq  => $prim_seq,
                sec_seq   => $second_seq
              };
        }
        push @results,
          { name => $ds_name, gene_no => $gene_no, results => $trna_data };
    }
    my $data;

    foreach my $result (@results) {
        $data .= $result->{name} . "\n";
        $data .= $result->{gene_no} . "\n";
        $data .=
"\nNo.\tGene Type\tLocation\tLocation of Anticodon\tAnticodon\tLocation of Intron (if exists)\tSequence (Primary & Secondary Prediction)\n";
        #print STDERR $data,"\n";
        foreach my $trna ( @{ $result->{results} } ) {
            $data .=
                $trna->{num} . "\t"
              . $trna->{type} . "\t"
              . $trna->{loc} . "\t"
              . $trna->{antiloc} . "\t"
              . $trna->{anticodon} . "\t"
              . $trna->{intron} . "\t"
              . $trna->{prim_seq} . "\n";
            $data .= "\t\t\t\t\t\t" . $trna->{sec_seq} . "\n";
        }
    }

    my $excel_file = create_excel_file( \@results );

    return $data, $excel_file;

#		if ($line =~ />/ && $i > 0) {$i++;}
#		if ($line =~/>/ || $line =~/gene/)
#		{
#			$line .= "\n";
#			$header_flag++;
#		}
#		else{
#			if($line =~ /tRNA-/)
#			{
#				$line =~ s/\s+/\t/g;
#				if($line=~/\)i\(/)
#				{
#					$line =~ s/\)i\(/\)\ti\(/;
#				}
#				else{
#					$line .= "\t";
#				}
#
#				$line =~ s/\n//;
#				$line .= "\t";
#			}
#			else
#			{
#				if($flag)
#				{
#					$line = $line."\n";
#				}
#				else{
#					$line = "\t\t\t\t\t\t".$line."\n";
#				}
#			}
#
#			$flag = $line =~ /tRNA-/ ? 1 : 0;
#			$line = "No.\tGene Type\tLocation\tLocation of Anticodon\tAnticodon\tLocation of Intron (if exists)\tSequence (Primary & Secondary Prediction)\n".$line if $header_flag > 1;
#			$header_flag = 0;
#		}
#$str .= $line;

}

sub read_file {
    my $file = shift;
    my $tmp;
    open( IN, $file ) || die "can't open $file for reading: $!";
    while (<IN>) {
        $tmp .= $_;
    }
    close IN;
    return $tmp;
}

sub save_file {
    my $output = shift;
    my $file   = $cogeweb->basefile . "_aragorn.tabbed";
    #print STDERR "lalalalalalalala\n\n";
    #print STDERR $file, "\n";
    open( NEW, "> $file" ) || die "Cannot Save $!\n";
    print NEW $output;
    close NEW;
    $file =~ s/$TEMPDIR/$TEMPURL/;
    return $file;

}

sub make_fasta_file {
    my $seq  = shift;
    my $file = $cogeweb->basefile . "fasta";
    open( NEW, "> $file" ) || die "Cannot Save $!\n";
    print NEW $seq;
    close NEW;
    return $file;
}

sub get_genome_fasta {
    my $orgid = shift;
    my $org   = $coge->resultset('Organism')->find($orgid);
    my @ds    = $org->current_datasets();
    @ds = sort { $a->id <=> $b->id } @ds;
    return unless @ds;
    my $org_name = $ds[0]->organism->name;
    my $title    = $org_name . " (";
    my %vers     = map { $_->version, 1 } @ds;
    if ( keys %vers > 1 ) {
        my @chrs;
        foreach my $ds (@ds) {
            push @chrs,
              join( ", ",
                map { "chr:" . $_ . " v:" . $ds->version . " ds:" . $ds->id }
                  $ds->get_chromosomes );
        }
        $title .= join( ", ", @chrs );
    }
    else {
        $title .= "v" . join( "", keys %vers ) . " ds:" . $ds[0]->id();
    }
    $title .= ")";
    $title =~ s/(`|')//g;
    my $md5  = md5_hex($title);
    my $file = $FASTADIR . "/$md5.fasta";
    my $res;
    if ( -r $file ) {
        CoGe::Accessory::Web::write_log( "*$org_name* fasta file ($md5) exists",
            $cogeweb->logfile );
        $res = 1;
    }
    else {
        $res = generate_fasta( dslist => \@ds, file => $file ) unless -r $file;
    }
    return $file if $res;
    return 0;

}

sub generate_fasta {
    my %opts   = @_;
    my $dslist = $opts{dslist};
    my $file   = $opts{file};
    $file = $FASTADIR . "/$file" unless $file =~ /$FASTADIR/;
    CoGe::Accessory::Web::write_log( "creating fasta file.",
        $cogeweb->logfile );
    open( OUT, ">$file" ) || die "Can't open $file for writing: $!";
    foreach my $ds (@$dslist) {
        next unless $USER->has_access_to_dataset($ds);
        foreach my $chr ( sort $ds->get_chromosomes ) {
            my $title =
                $ds->organism->name . " (v"
              . $ds->version . ") "
              . "chromosome: $chr"
              . ", CoGe database id: "
              . $ds->id;
            $title =~ s/^>+/>/;
            CoGe::Accessory::Web::write_log( "adding sequence $title",
                $cogeweb->logfile );
            print OUT ">" . $title . "\n";
            print OUT $ds->get_genomic_sequence( chr => $chr ), "\n";
        }
    }
    close OUT;
    CoGe::Accessory::Web::write_log( "Completed fasta creation",
        $cogeweb->logfile );
    return 1 if -r $file;
    CoGe::Accessory::Web::write_log( "Error with fasta file creation",
        $cogeweb->logfile );
    return 0;
}

sub email_results {
    my %opts          = @_;
    my $email_address = $opts{email};
    my $excel         = $opts{excel};
    my $tabbed        = $opts{tabbed};
    my $raw_file      = $opts{raw};
    #print STDERR Dumper \%ENV;
    #print STDERR "my server is $server\n";
    return
      unless $email_address =~
/^[_a-zA-Z0-9-]+(\.[_a-zA-Z0-9-]+)*@[a-zA-Z0-9-]+(\.[a-zA-Z0-9-]+)*\.(([0-9]{1,3})|([a-zA-Z]{2,3})|(aero|coop|info|museum|name))$/;
    my $server = $ENV{HTTP_HOST};

    my $raw_url = "http://" . $server . $raw_file;
    my $tabbed_url;
    $tabbed_url = "http://" . $server . $tabbed if $tabbed;
    my $excel_url;
    $excel_url = "http://" . $server . $excel if $excel;

    my $url = "Raw Aragorn Output: $raw_url";
    $url .= "\nTab-deliminated Aragorn Output: $tabbed_url" if $tabbed_url;
    $url .= "\nExcel File Aragorn Output: $excel_url"       if $excel_url;

    my $mailer = Mail::Mailer->new("sendmail");
    $mailer->open(
        {
            From    => 'CoGE <coge_results@genomevolution.org>',
            To      => $email_address,
            Subject => 'Aragorn Search Results',
        }
    ) or die "Can't open: $!\n";
    my $username = $USER->user_name;
    $username = $USER->first_name if $USER->first_name;
    $username .= " " . $USER->last_name
      if $USER->first_name && $USER->last_name;
    my $body = qq{Dear $username,

Thank you for using CoGe and Aragorn! The results from your latest analysis are ready, and can be viewed here:

$url

These results will remain on our servers for approximately 24 hours; please save them to your own computer before they are removed.

The following is a message from the creators of Aragorn:
--------------------------------------------------------
Please reference the following papers if you use this
program as part of any published research.

Laslett, D. and Canback, B. (2004) ARAGORN, a
program for the detection of transfer RNA and transfer-messenger
RNA genes in nucleotide sequences
Nucleic Acids Research, 32;11-16

Laslett, D. and Canback, B. (2008) ARWEN: a
program to detect tRNA genes in metazoan mitochondrial
nucleotide sequences
Bioinformatics, 24(2); 172-175.
--------------------------------------------------------

Thank you for using the CoGe Software Package.

- The CoGe Team
};

    print $mailer $body;
    $mailer->close();
}

sub check_address_validity {
    my $address = shift;
    return 'valid' unless $address;
    my $validity =
      $address =~
/^[_a-zA-Z0-9-]+(\.[_a-zA-Z0-9-]+)*@[a-zA-Z0-9-]+(\.[a-zA-Z0-9-]+)*\.(([0-9]{1,3})|([a-zA-Z]{2,3})|(aero|coop|info|museum|name))$/
      ? 'valid'
      : 'invalid';
    return $validity;
}
