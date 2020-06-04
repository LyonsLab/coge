#!/usr/bin/perl -w

use strict;
use CoGeX;
use Getopt::Long;
use File::Copy;
use File::Path;
use File::Spec::Functions;
use URI::Escape::JavaScript qw(unescape);
use Data::Dumper;

our ($DEBUG, $dbname, $user, $pass, $gid, $config, $host, $port, $P,
     $filename, $download_dir, $db, $feature_type_id, $feature_type);

GetOptions(
    "debug=s"           => \$DEBUG,
    "gid=i"             => \$gid,
    "download_dir=s"    => \$download_dir,
    "filename|f=s"      => \$filename,
    "feature_type|ft=s"      => \$feature_type,
    "feature_type_id|ftid=i"   => \$feature_type_id,
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

    if ($feature_type) {
        my $rs_feature_type = $db->resultset('FeatureType')->search({'name'=> $feature_type}); 
                $feature_type_id = $rs_feature_type->next->id;
    }

    get_locs(datasets => \@datasets, file => $file_temp, ftid => $feature_type_id);

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
    #print "In get_locs\n";
    my %args = @_;
    my $datasets = $args{datasets};
    my $file = $args{file};
    my $ftid = $args{ftid};

    my $SEP      = "\t";
    my %chrs;
    foreach my $ds (@$datasets) {
        my $dso = $db->resultset('Dataset')->resolve($ds);
        foreach my $chr ( $dso->get_chromosomes ) {
            #print $chr,"\n";
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
        my $gene_rs;
        if ($ftid) {
           $gene_rs = $db->resultset('Feature')->search(
            {
                'me.dataset_id' => { 'IN' => $datasets },
                'me.chromosome' => $chr,
                'me.feature_type_id' => $ftid
            },
            {
                'prefetch' => [ 'feature_names' ],
                'order_by' => ['me.start']
            }
          );
        } else {
            $gene_rs = $db->resultset('Feature')->search(
            {
                'me.dataset_id' => { 'IN' => $datasets },
                'me.chromosome' => $chr,    
            },
            {
                'prefetch' => [ 'feature_names' ],
                'order_by' => ['me.start']
            });
        }

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
            #print $index,"\t", $g->id,"\t",$gene_name,"\n";

            my $strand = $g->strand == 1 ? '+' : '-';
            my $clean_name = $gene_name;
            $names{$gene_name} = 1;
            $clean_name =~ s/\s+/_/g;

            my @locs;
            foreach  my $loc ( $g->locations( {}, { order_by => ['me.start'] } ) )
                {
                    @locs = @{ add_sorted_locs( \@locs, $loc ) };
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
