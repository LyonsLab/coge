#!/usr/bin/perl -w

use strict;
use Data::Dumper;
use CoGeX;
use CGI;
use CGI::Carp 'fatalsToBrowser';
use CoGe::Accessory::Web;
use File::Path;

use vars
  qw($P $DBNAME $DBHOST $DBPORT $DBUSER $DBPASS $connstr $PAGE_NAME $TEMPDIR $TEMPURL $USER $FORM $DATE $DEBUG $SERVER $coge $COOKIE_NAME);

$P                 = CoGe::Accessory::Web::get_defaults();
$ENV{PATH}         = $P->{COGEDIR};
$ENV{irodsEnvFile} = "/var/www/.irods/.irodsEnv";

# set this to 1 to print verbose messages to logs
$DEBUG   = 0;
$TEMPDIR = $P->{TEMPDIR} . "CoGe2Bed";
$TEMPURL = $P->{TEMPURL} . "CoGe2Bed";
$SERVER  = $P->{SERVER};
mkpath( $TEMPDIR, 0, 0777 ) unless -d $TEMPDIR;

$|    = 1;         # turn off buffering
$DATE = sprintf(
    "%04d-%02d-%02d %02d:%02d:%02d",
    sub { ( $_[5] + 1900, $_[4] + 1, $_[3] ), $_[2], $_[1], $_[0] }
      ->(localtime)
);

$FORM = new CGI;

$PAGE_NAME = 'CoGe2Bed.pl';
$DBNAME    = $P->{DBNAME};
$DBHOST    = $P->{DBHOST};
$DBPORT    = $P->{DBPORT};
$DBUSER    = $P->{DBUSER};
$DBPASS    = $P->{DBPASS};
$connstr =
  "dbi:mysql:dbname=" . $DBNAME . ";host=" . $DBHOST . ";port=" . $DBPORT;
$coge = CoGeX->connect( $connstr, $DBUSER, $DBPASS );
#$coge->storage->debugobj(new DBIxProfiler());
#$coge->storage->debug(1);

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

my $gid = $FORM->param('gid');
$gid = $FORM->param('dsgid') unless $gid;

unless ($gid) {
    print $FORM->header;
    print "Error: no genome id\n";
    exit;
}

my $genome = $coge->resultset('Genome')->find($gid);
unless ($genome) {
    print $FORM->header;
    print "Error: no genome for id $gid\n";
    exit;
}
if ( $genome && $genome->restricted ) {
    if ( !$USER->has_access_to_genome($genome) ) {
        print $FORM->header;
        print "Error: permission denied\n";
        exit;
    }
}
my @datasets = map { $_->dataset_id } $genome->datasets;

my $org    = $genome->organism->name;
$org =~ s/\///g;
$org =~ s/\s+/_/g;
$org =~ s/\(//g;
$org =~ s/\)//g;
$org =~ s/://g;
$org =~ s/;//g;
$org =~ s/#/_/g;
$org =~ s/'//g;
$org =~ s/"//g;
my $header = "Content-disposition: attachement; filename=";    #test.gff\n\n";
$header .= $org;
$header .= "_gid$gid";
$header .= ".bed\n\n";
print $header;

my $chrs = get_locs( \@datasets );

sub get_sizes {
    my $locs_ref = shift;
    my $n        = scalar( @{$locs_ref} );
    my $i        = 0;
    my $s        = "";
    for ( $i = 0 ; $i < $n ; $i += 2 ) {
        $s .= ( $locs_ref->[ $i + 1 ] - $locs_ref->[$i] + 1 ) . ",";
    }
    chop $s;
    return $s;
}

sub get_starts {
    my $locs_ref = shift;
    my $start    = shift;
    my $s        = "";
    for ( my $i = 0 ; $i < scalar( @{$locs_ref} ) ; $i += 2 ) {
        $s .= ( $locs_ref->[$i] - $start ) . ",";
    }
    chop $s;
    return $s;
}

sub get_locs {

    my $datasets = shift;
    my $SEP      = "\t";
    my %chrs;
    foreach my $ds (@$datasets) {
        my $dso = $coge->resultset('Dataset')->resolve($ds);
        foreach my $chr ( $dso->get_chromosomes ) {
            $chrs{$chr} = $dso->last_chromosome_position($chr);
        }
    }
    my @chrs  = sort { $a cmp $b } keys %chrs;
    my %names = ();
    my %fids  = ();
    my $index = 0;

    foreach my $chr (@chrs) {
        my %seen    = ();
        my $gene_rs = $coge->resultset('Feature')->search(
            {
                'me.dataset_id' => { 'IN' => $datasets },
                'me.chromosome' => $chr,
                # NOTE: should probably check for pseudogenes as well!!
                'feature_type.name' => {
                    'IN' => [
                        'gene',                 'pseudogene',
                        'transposon',           'Transposable Element',
                        'transposable_element', 'transposable_element_gene'
                    ]
                }
            },
            {
                'prefetch' => [ 'feature_type', 'feature_names' ],
                'order_by' => ['me.start']
            }
        );

        my %chrs;
#        print STDERR "dataset_ids: " . join(",", @$datasets) . ";  chr: $chr\n";
        while ( my $g = $gene_rs->next() ) {
            if ( $fids{ $g->feature_id } ) { next; }
            $fids{ $g->feature_id } = 1;
            my @gene_names;
            @gene_names = $g->feature_names();
            my $gene_name;
            foreach my $name (@gene_names) {
                $gene_name = $name->name;
                if ( !$names{$gene_name} ) { last; }
                $gene_name = "";
            }
            if ( !$gene_name || $names{$gene_name} ) { next; }
            #if ($index > 100) {print STDERR "breaking at 100\n";return ; }

            $index++;
            my $cds_rs = $coge->resultset('Feature')->search(
                {
                    'me.dataset_id'      => { 'IN' => $datasets },
                    'me.chromosome'      => $chr,
                    'feature_names.name' => $gene_name,
                    'feature_type.name'  => 'CDS'
                },
                {
                    'join'     => 'feature_names',
                    'prefetch' => [ 'feature_type', 'locations' ],
                    'order_by' => [ 'me.start', 'locations.start' ]
                }
            );

            my $strand = $g->strand == 1 ? '+' : '-';
            my $clean_name = $gene_name;
            $names{$gene_name} = 1;
            $clean_name =~ s/\s+/_/g;

            my $parent  = $clean_name;
            my $has_cds = 0;
            my @locs;
            while ( my $f = $cds_rs->next() ) {
                if ( $fids{ $f->feature_id } ) { next; }
                $fids{ $f->feature_id } = 1;
                $has_cds = 1;
                foreach
                  my $loc ( $f->locations( {}, { order_by => ['me.start'] } ) )
                {
                    @locs = @{ add_sorted_locs( \@locs, $loc ) };
                }
            }

            my @data;
            $chrs{ $g->chr } = 1;
            if ( scalar(@locs) != 0 ) {
                @data = (
                    $chr, $g->start - 1,
                    $g->stop, $clean_name, $g->stop - $g->start,
                    $strand, ".", ".", ".", scalar(@locs) / 2
                );
            }
            else {
                my $sub_rs = $coge->resultset('Feature')->search(
                    {
                        'me.dataset_id'      => { 'IN' => $datasets },
                        'feature_names.name' => $gene_name,
                        'feature_type.name' =>
                          { 'NOT IN' => [ 'gene', 'CDS', 'intron' ] }
                    },
                    {
                        'join'     => ['feature_names'],
                        'prefetch' => [ 'feature_type', 'locations' ],
                        'order_by' => [ 'me.chromosome', 'me.start' ]
                    }
                );

                undef @locs;
                my $ftype;
                while ( my $f = $sub_rs->next() ) {
                    if ( $fids{ $f->feature_id } ) { next; }
                    $fids{ $f->feature_id } = 1;
                    $ftype = $f->type->name;
                    foreach my $loc (
                        $f->locations( {}, { order_by => ['me.start'] } ) )
                    {
                        @locs = @{ add_sorted_locs( \@locs, $loc ) };
                    }
                }

                @data = (
                    $chr, $g->start - 1,
                    $g->stop, $clean_name, $g->stop - $g->start,
                    $strand, ".", ".", "."
                );    #, scalar(@locs)/2);
                push( @data, $ftype ? scalar(@locs) / 2 : 1 );
            }
            if ( scalar(@locs) == 0 ) {
                @locs = ( $g->start, $g->stop );
            }
            my $start = $g->start < $locs[0] ? $g->start : $locs[0];
            $data[1] = $start - 1;
            push( @data, get_sizes( \@locs ) );
            push( @data, get_starts( \@locs, $start ) );
            my $gstr = join( "\t", @data, $g->type->name );
            print $gstr . "\n";
        }
    }
    return \%chrs;
}

sub add_sorted_locs {
    my $loc_ref = shift;
    my @locs    = @{$loc_ref};
    my $loc     = shift;
    my $l       = scalar(@locs);
    # dont add exons repeatedly.
    if (   $l > 0
        && $locs[ $l - 2 ] == $loc->start
        && $locs[ $l - 1 ] == $loc->stop )
    {
        return \@locs;
    }
    # merge overlapping / alternative splicings.
    my $added = 0;
    for ( my $i = 1 ; $i < $l ; $i += 2 ) {
        if ( $loc->start <= $locs[$i] && $loc->stop >= $locs[$i] ) {
            $locs[$i] = $loc->stop;
            $added = 1;
        }
        if ( $loc->start <= $locs[ $i - 1 ] && $loc->stop >= $locs[ $i - 1 ] ) {
            $locs[ $i - 1 ] = $loc->start;
            $added = 1;
        }
        if ($added) { last; }
    }
    if ( !$added ) {
        push( @locs, $loc->start );
        push( @locs, $loc->stop );
    }
    @locs = sort { $a <=> $b } @locs;
    return \@locs;
}
