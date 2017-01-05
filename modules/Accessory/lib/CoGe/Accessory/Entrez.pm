package CoGe::Accessory::Entrez;

=head1 NAME

CoGe::Accessory::SRA

=head1 SYNOPSIS

Interface to NCBI-SRA

=head1 DESCRIPTION

Provide simple API for accessing NCBI-SRA via Entrez

=head1 AUTHOR

Matt Bomhoff

=head1 COPYRIGHT

The full text of the license can be found in the
LICENSE file included with this module.

=head1 SEE ALSO

=cut

use strict;
use warnings;

use LWP::Simple qw($ua get);
use XML::Simple qw(XMLin);
use Data::Dumper;

$ua->timeout(10); # Set LWP timeout to 10 seconds

BEGIN {
    use vars qw ($VERSION @ISA @EXPORT $BASE_URL $ESEARCH_URL $ESUMMARY_URL $EFETCH_URL $RETMAX );
    require Exporter;

    $BASE_URL     = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/';
    $ESEARCH_URL  = $BASE_URL . 'esearch.fcgi';
    $ESUMMARY_URL = $BASE_URL . 'esummary.fcgi';
    $EFETCH_URL   = $BASE_URL . 'efetch.fcgi';
    $RETMAX       = 100_000; # maximum number of search results (Entrez max value)

    $VERSION = 0.1;
    @ISA     = qw( Exporter );
    @EXPORT  = qw( esearch esummary efetch );
}

sub esearch {
    my ($db, $term) = @_;
    $db = lc($db);

    my $result = get($ESEARCH_URL . '?db=' . $db . '&term=' . $term . '&retmax=' . $RETMAX);
    my $record = XMLin($result, ForceArray => ['Id']);
    #warn Dumper $record;
    my $id_list = $record->{IdList}->{Id};
    #warn Dumper $id_list;

    return $id_list; # arrayref of IDs
}

sub esummary {
    my ($db, $ids) = @_;
    $db = lc($db);

    my $xmlResponse = get($ESUMMARY_URL . '?db=' . $db . '&id=' . ( ref($ids) eq 'ARRAY' ? join(',', @$ids) : $ids ) );
    #warn Dumper $xmlResponse;

    if ($db eq 'sra') {
        my $response = XMLin($xmlResponse, KeyAttr => {Item => 'Name'}, ForceArray => ['Item'], ContentKey => '-content');
        #warn Dumper $response;
        unless ($response) {
            warn "esummary: No response";
            return;
        }

        my %result;
        foreach my $d (@{$response->{DocSum}}) {
            # Parse name
            my $exp = '<Exp>'.$d->{Item}->{ExpXml}->{content}.'</Exp>'; # Kludge for NCBI, see COGE-778
            if ($exp) {
                my $record = XMLin($exp);
                @result{keys %$record} = values %$record;
            }

            # Parse runs
            my $runs = '<Runs>'.$d->{Item}->{Runs}->{content}.'</Runs>'; # Kludge for NCBI, see COGE-778
            if ($runs) {
                my $record = XMLin($runs, KeyAttr => { Run => 'acc' });
                @result{keys %$record} = values %$record;
            }
        }

        #warn Dumper \%result;
        return \%result;
    }

    return $xmlResponse;
}

sub efetch {
    my ($db, $id) = @_;
    $db = lc($db);

    my $result = get($EFETCH_URL . '?db=' . $db . '&id=' . $id);
    #warn $result;

    if ($db eq 'sra') {
        my $record = XMLin($result);
        #warn Dumper $record;
    }

    return $result;
}

1;
