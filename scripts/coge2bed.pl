#!/usr/bin/perl -w

use strict;
use CoGeX;
use Getopt::Long;
use File::Copy;
use File::Path;
use File::Spec::Functions;
use URI::Escape::JavaScript qw(unescape);

our ($DEBUG, $dbname, $user, $pass, $gid, $config, $host, $port, $P,
     $filename, $download_dir, $db);

GetOptions(
    "debug=s"           => \$DEBUG,
    "gid=i"             => \$gid,
    "download_dir=s"    => \$download_dir,
    "filename|f=s"      => \$filename,

    # Database params
    "host|h=s"          => \$host,
    "port|p=s"          => \$port,
    "database|db=s"     => \$dbname,
    "user|u=s"          => \$user,
    "password|pw=s"     => \$pass,

    # Or use config file
    "config=s"          => \$config,
);

$| = 1;
sub main {
    $download_dir //= ".";
    mkpath($download_dir, 0, 0777) unless -r $download_dir;

    my $logfile = catfile($download_dir, "$filename.log");
    open (my $logh, ">", $logfile) or die "Error opening log file";

    if ($config) {
        $P     = CoGe::Accessory::Web::get_defaults($config);
        $dbname= $P->{DBNAME};
        $host = $P->{DBHOST};
        $port = $P->{DBPORT};
        $user = $P->{DBUSER};
        $pass = $P->{DBPASS};
    }

    # Verify parameters
    $gid = unescape($gid) if $gid;
    $filename = unescape($filename) if $filename;

    if (not $gid) {
        say $logh "log: error: genome not specified use gid";
        exit(-1);
    }

    if (not $filename) {
        say $logh "log: error: output file not specified use output";
        exit(-1);
    }

    my $connstr = "dbi:mysql:db=$dbname;host=$host;port=$port;";
    $db = CoGeX->connect( $connstr, $user, $pass );
    #$db->storage->debugobj(new DBIxProfiler());
    #$db->storage->debug(1);

    unless ($db) {
        say $logh "log: error: couldn't connect to database";
        exit(-1);
    }

    my $file = catfile($download_dir, $filename);
    my $file_temp = $file . ".tmp";

    # Exit if file already exists
    return if -r $file;

    my $genome = $db->resultset('Genome')->find($gid);
    my @datasets = map { $_->dataset_id } $genome->datasets;

    get_locs(datasets => \@datasets, file => $file_temp);

    exit 1 unless move($file_temp, $file);
}

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
    my %args = @_;
    my $datasets = $args{datasets};
    my $file = $args{file};

    my $SEP      = "\t";
    my %chrs;
    foreach my $ds (@$datasets) {
        my $dso = $db->resultset('Dataset')->resolve($ds);
        foreach my $chr ( $dso->get_chromosomes ) {
            $chrs{$chr} = $dso->last_chromosome_position($chr);
        }
    }
    my @chrs  = sort { $a cmp $b } keys %chrs;
    my %names = ();
    my %fids  = ();
    my $index = 0;

    open(my $fh, ">", $file) or die "Error creating bed file";

    foreach my $chr (@chrs) {
        my %seen    = ();
        my $gene_rs = $db->resultset('Feature')->search(
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
            my $cds_rs = $db->resultset('Feature')->search(
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
                my $sub_rs = $db->resultset('Feature')->search(
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
            print $fh $gstr . "\n";
        }
    }
    close($fh);

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

main;
