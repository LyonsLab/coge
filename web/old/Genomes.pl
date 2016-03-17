#! /usr/bin/perl -w

use strict;
use CGI;
use CoGeX;
use CoGe::Accessory::Web;
use CoGe::Core::Genome qw(genomecmp);
use HTML::Template;
use Sort::Versions;
no warnings 'redefine';

use vars qw( $P $PAGE_TITLE $USER $coge $FORM %FUNCTION $LINK );

$PAGE_TITLE = 'Genomes';

$FORM = new CGI;
( $coge, $USER, $P, $LINK ) = CoGe::Accessory::Web->init(
	cgi => $FORM,
	page_title => $PAGE_TITLE
);

%FUNCTION = (
    generate_html        => \&generate_html,
    create_genome        => \&create_genome,
    delete_genome        => \&delete_genome,
    get_genomes_for_user => \&get_genomes_for_user,
);

CoGe::Accessory::Web->dispatch( $FORM, \%FUNCTION, \&generate_html );

sub delete_genome {
    my %opts = @_;
    my $gid  = $opts{gid};
    return "Must have valid genome id\n" unless ($gid);
    print STDERR "delete_genome $gid\n";

    # Check permission
    return 0 unless ( $USER->is_admin or $USER->is_owner( dsg => $gid ) );

    # Delete the genome and associated connectors & datasets
    my $genome = $coge->resultset('Genome')->find($gid);
    return 0 unless $genome;
    $genome->deleted(1);
    $genome->update;

    CoGe::Accessory::Web::log_history(
        db          => $coge,
        user_id     => $USER->id,
        page        => "$PAGE_TITLE",
        description => 'delete genome id' . $genome->id
    );

    return 1;
}

sub get_genomes_for_user {

    #my %opts = @_;

    my @genomes;

    #	if ( $USER->is_admin ) {
    #		@genomes = $coge->resultset('Genome')->all();
    #	}
    #	else {
    push @genomes, $USER->genomes;

    #	}

    my @genome_info;
    foreach my $g (@genomes) {    #( sort genomecmp @genomes ) {
        push @genome_info,
          {
            NAME =>
              qq{<span class="link" onclick='window.open("GenomeInfo.pl?gid=}
              . $g->id
              . qq{")'>}
              . $g->info
              . "</span>",
            VERSION => $g->version,
            DATE    => $g->get_date,
            EDIT_BUTTON =>
"<span class='link ui-icon ui-icon-gear' onclick=\"window.open('GenomeInfo.pl?gid="
              . $g->id
              . "')\"></span>",
            DELETE_BUTTON =>
"<span class='link ui-icon ui-icon-trash' onclick=\"dialog_delete_genome({gid: '"
              . $g->id
              . "'});\"></span>"
          };
    }

    my $template =
      HTML::Template->new( filename => $P->{TMPLDIR} . $PAGE_TITLE . '.tmpl' );
    $template->param( DO_GENOME_TABLE => 1 );
    $template->param( GENOME_LOOP     => \@genome_info );

    #	$template->param( TYPE_LOOP => get_list_types() );

    return $template->output;
}

sub generate_html {
    my $template =
      HTML::Template->new( filename => $P->{TMPLDIR} . 'generic_page.tmpl' );
    $template->param( PAGE_TITLE => $PAGE_TITLE,
				      TITLE      => 'Genomes',
    				  PAGE_LINK  => $LINK,
    				  HOME       => $P->{SERVER},
                      HELP       => 'Genomes',
                      WIKI_URL   => $P->{WIKI_URL} || '' );
    my $name = $USER->user_name;
    $name = $USER->first_name if $USER->first_name;
    $name .= ' ' . $USER->last_name
      if ( $USER->first_name && $USER->last_name );
    $template->param( USER     => $name );
    $template->param( LOGON    => 1 ) unless $USER->user_name eq "public";
    my $link = "http://" . $ENV{SERVER_NAME} . $ENV{REQUEST_URI};
    $link = CoGe::Accessory::Web::get_tiny_link( url => $link );

    $template->param( BODY       => generate_body() );
    $template->param( ADMIN_ONLY => $USER->is_admin );
    $template->param( CAS_URL    => $P->{CAS_URL} || '' );

    return $template->output;
}

sub generate_body {
    my $template =
      HTML::Template->new( filename => $P->{TMPLDIR} . $PAGE_TITLE . '.tmpl' );
    $template->param( MAIN         => 1 );
    $template->param( PAGE_NAME    => $PAGE_TITLE );
    $template->param( GENOME_TABLE => get_genomes_for_user() );

    return $template->output;
}
