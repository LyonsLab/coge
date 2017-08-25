#!/usr/bin/perl -w

use strict;
use Data::Dumper;
use CoGeX;
use Getopt::Long;
use File::Path;
use File::Basename;
use File::Spec::Functions qw( catdir catfile );
use File::Touch;
use URI::Escape;
use URI::Escape::JavaScript qw(unescape);
use CoGe::Accessory::Web qw(get_defaults get_command_path);
use CoGe::Accessory::Utils qw( commify units );
use CoGe::Core::Genome qw(fix_chromosome_id);
use CoGe::Core::Storage;
use List::Util qw( min max );
use Benchmark;

my $t1 = new Benchmark;

my $GO          = 1;
#my $DEBUG       = 0;
my $DB_BATCH_SZ = 10 * 1000;
use vars qw($staging_dir $data_file $wid
  $name $description $link $version $restricted $creator_id
  $gid $source_name $user_name $user_id $config $allow_all_chr
  $host $port $db $user $pass $P $GUNZIP);

GetOptions(
    "staging_dir=s" => \$staging_dir,
#    "result_dir=s"  => \$result_dir,     # results path
    "wid=s"         => \$wid,            # workflow id
    "data_file=s"   => \$data_file,      # data file (JS escape)
    "name=s"        => \$name,           # experiment name (JS escaped)
    "desc=s"        => \$description,    # experiment description (JS escaped)
    "link=s"        => \$link,           # experiment description (JS escaped)
    "version=s"     => \$version,        # experiment version (JS escaped)
    "restricted=i"  => \$restricted,     # experiment restricted flag
    "gid=s"         => \$gid,            # genome id
    "source_name=s" => \$source_name,    # data source name (JS escaped)
    "user_id=i"     => \$user_id,        # user ID to assign dataset
    "user_name=s"   => \$user_name,      # user name to assign dataset (alternative to user_id)
    "creator_id=i"  => \$creator_id,     # user ID to set as dataset creator
    "config=s"      => \$config,         # configuration file

    # Optional Flags
    "allow_all_chr=i" => \$allow_all_chr # Allow non-existent chromosomes
);
$allow_all_chr = 1 unless defined $allow_all_chr; # mdb added 8/3/16 -- set default, many users are loading GFFs with contigs not present in the genome

$| = 1;
print STDOUT "Starting $0 (pid $$)\n", qx/ps -o args $$/;

# Setup staging path
unless ($staging_dir) {
    print STDOUT "log: error: staging_dir argument is missing\n";
    exit(-1);
}
mkpath($staging_dir, 0, 0777) unless -r $staging_dir;

# Prevent loading again (issue #417)
my $logdonefile = "$staging_dir/log.done";
if (-e $logdonefile) {
    print STDOUT "log: error: done file already exists: $logdonefile\n";
    exit(-1);
}

# Open log file for detailed feature info -- everything else goes to STDOUT for JEX
my $logfile = "$staging_dir/load_annotation.log";
open( my $log, ">$logfile" ) or die "Error opening log file $logfile";
$log->autoflush(1);

# Process and verify parameters
$data_file   = unescape($data_file);
$name        = unescape($name);
$description = unescape($description);
$link        = unescape($link);
$version     = unescape($version);
$source_name = unescape($source_name);
$restricted  = '0' if ( not defined $restricted ); # doesn't really matter, datasets inherit permissions from genome

if (not defined $user_id and not defined $user_name) {
    print STDOUT "log: error: user not specified, use user_id or user_name\n";
    exit(-1);
}

if ((defined $user_name and $user_name eq 'public') || (defined $user_id and $user_id eq '0')) {
    print STDOUT "log: error: not logged in\n";
    exit(-1);
}

# Load config file
unless ($config) {
    print STDOUT "log: error: can't find config file\n";
    exit(-1);
}
$P    = CoGe::Accessory::Web::get_defaults($config);
$db   = $P->{DBNAME};
$host = $P->{DBHOST};
$port = $P->{DBPORT};
$user = $P->{DBUSER};
$pass = $P->{DBPASS};

$GUNZIP = get_command_path('GUNZIP');

# Validate the data file
print STDOUT "log: Validating data file ...\n";
unless ( -e $data_file ) {
    print STDOUT "log: can't find data file\n";
    exit(-1);
}

# Connect to database
my $connstr = "dbi:mysql:dbname=$db;host=$host;port=$port;";
my $coge = CoGeX->connect( $connstr, $user, $pass );
unless ($coge) {
  print STDOUT "log: couldn't connect to database\n";
    exit(-1);
}

# Retrieve user
my $user;
if ($user_id) {
    $user = $coge->resultset('User')->find($user_id);
}
elsif ($user_name) {
    $user = $coge->resultset('User')->find( { user_name => $user_name } );
}
else {
    print STDOUT "log: error user not specified, see user_id or user_name\n";
    exit(-1);
}

unless ($user) {
    print STDOUT "log: error finding user ", ($user_name ? $user_name : $user_id) , "\n";
    exit(-1);
}

# Retrieve creator
my $creator;
if ($creator_id) {
    $creator = $coge->resultset('User')->find($creator_id);
    unless ($creator) {
        print STDOUT "log: error finding creator $creator_id\n";
        exit(-1);
    }
}
$creator = $user unless $creator;

# Retrieve genome
my $genome = $coge->resultset('Genome')->find( { genome_id => $gid } );
unless ($genome) {
    print STDOUT "log: error finding genome id$gid\n";
    exit(-1);
}

# Get list of chromosomes for genome
my %valid_chrs = map { $_ => 1 } $genome->chromosomes;

# Some defaults to check for in names and annotations
my @check_names = (
    "ID",           "name",   "Name",     "Alias",
    "gene",         "Parent", "Locus_id", "ID_converter",
    "Gene_symbols", "gene_id",
);
my %check_names = map { $_ => 1 } @check_names;

my @skip_attr =
  ( "Link_to", "References", "Sequence_download", "transcript_id", );
my %skip_attr = map { $_ => 1 } @skip_attr;

my @anno_names = (
    "Source",
    "Note",
    "NIAS_FLcDNA",
    "Comment",
    "GO",
    "ORF_evidence",    # can we link to SGD?
    "Transcript_evidence",
    "Status",
    "InterPro",
    "Pfam",
    "Description",
    "Function",
    "Derives_from",
);
my %anno_names = map { $_ => 1 } @anno_names;

my @skip_names_re = qw(
  :five_prime
  :three_prime
  :exon
  \.exon
  :utr
  \.utr
  _utr
  :cds
  \.cds
  cds\.
  :hsp
  \.hsp

  intron
  _E\d
  ^CDS$
  ^exon$
  ^tss$
  ^tts$
  ^intron$
);

# Decompress file if necessary
if ( $data_file =~ /\.gz$/ ) {
    print STDOUT "log: Decompressing '$data_file'\n";
    $data_file =~ s/\.gz$//;
    execute( $GUNZIP . ' -c ' . $data_file . '.gz' . ' > ' . $data_file );
}

# Load GFF file into %data
my (%data, %annos, %seen_types, %seen_attr);
#TODO copy gff file into staging directory to read from instead of upload directory
unless ( process_gff_file() ) {
    print STDOUT "log: error: no annotations found, perhaps your file is missing required information, <span class='bold alert'>please check the <a href='http://genomevolution.org/wiki/index.php/GFF_ingestion'>GFF documentation</a></span>\n";
    exit(-1);
}

# Create gene annotations if none present in GFF file
unless ( $seen_types{gene} ) {
    print STDOUT "log: Creating gene entities\n";
    foreach my $chr_loc ( keys %data ) {
      name: foreach my $name ( keys %{ $data{$chr_loc} } ) {
            my ( $chr, $start, $stop, $strand );
            my %names;
            my $name = $data{$chr_loc}{$name};
            foreach my $type ( keys %$name ) {
                map { $names{$_} = 1 } keys %{ $name->{$type}{names} };
                foreach my $loc ( @{ $name->{$type}{loc} } ) {
                    next name if $type eq "gene";
                    $start = $loc->{start} unless $start;
                    $start = $loc->{start} if $loc->{start} < $start;
                    $stop = $loc->{stop} unless $stop;
                    $stop   = $loc->{stop}   if $loc->{stop} > $stop;
                    $strand = $loc->{strand} if $loc->{strand};
                    $strand = 1 unless (defined $strand); # mdb added 11/7/13 issue 248 - set default strand to '+'
                    $chr    = $loc->{chr};
                }
                foreach my $loc ( @{ $name->{$type}{loc} } ) {
                    $loc->{strand} = $strand;
                }
            }
            $name->{gene}{loc} = [
                {
                    start  => $start,
                    stop   => $stop,
                    strand => $strand,
                    chr    => $chr,
                }
            ] if (defined $start and defined $stop and defined $strand);
            $name->{gene}{names} = \%names;
            $seen_types{gene}++;
        }
    }
}

print STDOUT "log: Annotation types:\n", join(
    "\n",
    map {
        "log: &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;" . $_ . ": "
          . commify( $seen_types{$_} )
      } sort keys %seen_types
  ),
  "\n";

print STDOUT "log: Data types:\n", join(
    "\n",
    map {
        "log: &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;" . $_ . ": "
          . commify( $seen_attr{$_} )
      } sort keys %seen_attr
  ),
  "\n";

my $t2 = new Benchmark;

################################################################################
# If we've made it this far without error then we can feel confident about
# the input data.  Now we can go ahead and create the db entities.
################################################################################

# Create datasource
my $datasource = $coge->resultset('DataSource')->find_or_create( { name => $source_name, description => "" } );
unless ($datasource) {
    print STDOUT "log: error creating data source\n";
    exit(-1);
  }

# Create dataset
my $dataset = $coge->resultset('Dataset')->create(
    {
        data_source_id => $datasource->id,
        name           => basename($data_file),
        description    => $description,
        link           => $link,
        version        => $version,
        restricted     => $restricted,
        creator_id     => $creator->id,
    }
);
unless ($dataset) {
    print STDOUT "log: error creating dataset\n";
    exit(-1);
}

#TODO set link field if loaded from FTP
print STDOUT "dataset id: " . $dataset->id . "\n";

# Add dataset to the genome
my $dsconn = $coge->resultset('DatasetConnector')->find_or_create( { dataset_id => $dataset->id, genome_id => $genome->id } );
unless ($dsconn) {
    print STDOUT "log: error creating dataset connector\n";
    exit(-1);
}

# Count total annotations to load -- mdb added 1/8/14, issue 260
my $total_annot = 0;
foreach my $chr_loc ( keys %data ) {
    foreach my $name ( keys %{ $data{$chr_loc} } ) {
        foreach my $feat_type ( keys %{ $data{$chr_loc}{$name} } ) {
        	foreach ( @{$data{$chr_loc}{$name}{$feat_type}{loc}} ) {
        		$total_annot++;
        	}
        }
    }
}

# Populate all feature-related tables in DB
print STDOUT "log: Loading database ...\n";
my %anno_types;    # hash to store annotation type objects
my %feat_types;    # store feature type objects
my $loaded_annot = 0;
my @loc_buffer;     # buffer for bulk inserts into Location table
my @anno_buffer;    # buffer for bulk inserts into FeatureAnnotation table
my @name_buffer;    # buffer for bulk inserts into FeatureName table
foreach my $chr_loc ( sort { $a cmp $b } keys %data ) {
    foreach my $name ( sort { $a cmp $b } keys %{ $data{$chr_loc} } ) {
      	my $pctLoaded = int( 100 * $loaded_annot / $total_annot );
        print STDOUT "log: Loaded ", commify($loaded_annot), " annotations (", ( $pctLoaded ? $pctLoaded : '<1' ), "%)\n\n"
          if ( $loaded_annot and ( $loaded_annot % 1000 ) == 0 );

        foreach my $feat_type ( sort { $a cmp $b } keys %{ $data{$chr_loc}{$name} } ) {
	        print $log "\n" if $log;

            my ($start, $stop, $strand, $chr);
            my $loc = $data{$chr_loc}{$name}{$feat_type}{loc};
            if (@$loc) {
                $start    = min map { $_->{start} } @$loc;
                $stop     = max map { $_->{stop}  } @$loc;
                ($strand) = map { $_->{strand} } @$loc;
                ($chr)    = map { $_->{chr}    } @$loc;
            }
            else { # mdb added else 4/8/14 issue 358 - no locations (e.g. tRNA w/o parent)
                my $coords = $data{$chr_loc}{$name}{$feat_type}{coords};
                $start  = $coords->{start};
                $stop   = $coords->{stop};
                $strand = $coords->{strand};
                $chr    = $coords->{chr};
            }

            # Add feature type to DB
            $feat_types{$feat_type} = $coge->resultset('FeatureType')->find_or_create( { name => $feat_type } )
              if $GO && !$feat_types{$feat_type};
            my $feat_type_obj = $feat_types{$feat_type};

            print $log "Creating feature of type $feat_type\n" if $log;

            # mdb added check 4/8/14 issue 358
            unless (defined $start and defined $stop and defined $chr) {
                print $log "warning: feature '", (defined $name ? $name : ''), "' (type '$feat_type') missing coordinates", "\n" if $log;
                #print STDOUT Dumper $data{$chr_loc}{$name}{$feat_type}, "\n";
                next; #exit(-1);
            }

            # Add feature
            #TODO this could be batched by nesting location & other inserts,
            # see http://search.cpan.org/~abraxxa/DBIx-Class-0.08209/lib/DBIx/Class/ResultSet.pm#populate
            my $feat = $dataset->add_to_features(
                {
                    feature_type_id => $feat_type_obj->id,
                    start           => $start,
                    stop            => $stop,
                    chromosome      => $chr,
                    strand          => $strand
                }
            ) if $GO;
            my $featid = $feat ? $feat->id : "no_go";
            
            # Add locations
            my %seen_locs;
            my $loc_count = 0;
            foreach my $loc ( sort { $a->{start} <=> $b->{start} } @$loc ) {
                my ($start, $stop) = ($loc->{start}, $loc->{stop});
                $loc_count++;
                next if $feat_type eq "gene" && $loc_count > 1; #only use the first one as this will be the full length of the gene.  Stupid hack
                next if $seen_locs{$start}{$stop};
                $seen_locs{$start}{$stop} = 1;
                print $log "Adding location $chr:(" . $start . "-" . $stop . ", $strand)\n" if $log;
                $loaded_annot++;
                batch_add(
                    \@loc_buffer,
                    'location',
                    {
                        feature_id => $feat->id,
                        chromosome => $loc->{chr},
                        start      => $loc->{start},
                        stop       => $loc->{stop},
                        strand     => $loc->{strand}
                    }
                ) if $GO;
            }

            # Add feature names and annotations
            my %names = map { $_ => 1 } keys %{ $data{$chr_loc}{$name}{$feat_type}{names} };
            my %seen_annos; #hash to store annotations so duplicates aren't added
            master_names: foreach my $tmp ( keys %names ) {
                foreach my $re (@skip_names_re) {
                    next master_names if $tmp =~ /$re/i;
                }
                my $master = 0;
                $master = 1 if $tmp eq $name;
                print $log "Adding name $tmp to feature ", $featid,
                  ( $master ? " (MASTER)" : '' ), "\n"
                  if $log;

                batch_add(
                    \@name_buffer,
                    'feature_name',
                    {    
                        feature_id   => $feat->id,
                        name         => $tmp,
                        primary_name => $master
                    }
                ) if $GO;

                if ( $annos{$tmp} ) {
                    foreach my $anno ( keys %{ $annos{$tmp} } ) {
                        next unless $anno;
                        next if $seen_annos{$anno};
                        $seen_annos{$anno} = 1;
                        
                        # Add annotation type to DB
                        my $type_name = $annos{$tmp}{$anno}{type} || "Note";
                        my ($anno_type) = $anno_types{$type_name};
                        unless ($anno_type) {
                            ($anno_type) = $coge->resultset('AnnotationType')->find_or_create( { name => $type_name } );
                            $anno_types{$type_name} = $anno_type;
                        }
                        
                        # Add feature annotation to DB
                        my $link = $annos{$tmp}{$anno}{link};
                        print $log "Adding annotation ($type_name): $anno\n" . ( $link ? "\tlink: $link" : '' ) . "\n" if $log;
                        batch_add(
                            \@anno_buffer,
                            'feature_annotation',
                            {
                                feature_id         => $feat->id,
                                annotation_type_id => $anno_type->id,
                                annotation         => $anno,
                                link               => $link
                            }
                        ) if $GO && $anno;
                    }
                }
            }
        }
    }
}

# Flush DB insertion buffers
batch_add( \@loc_buffer,  'location' );
batch_add( \@name_buffer, 'feature_name' );
batch_add( \@anno_buffer, 'feature_annotation' );
print STDOUT "log: " . commify($loaded_annot) . " annotations loaded\n";

# Print times to parse and load
my $t3 = new Benchmark;
print STDOUT "Time to parse: "
  . timestr( timediff( $t2, $t1 ) )
  . ", Time to load: "
  . timestr( timediff( $t3, $t2 ) ) . "\n";

# Save result
unless (add_workflow_result($user_name, $wid, 
        {
            type           => 'dataset',
            id             => int($dataset->id),
            genome_id      => int($genome->id),
            info           => $genome->info
        })
    )
{
    print STDOUT "log: error: could not add workflow result\n";
    exit(-1);
}

# Create "log.done" file to indicate completion to JEX
touch($logdonefile);
print STDOUT "All done!\n";

exit;

#-------------------------------------------------------------------------------
sub batch_add {
    my $buffer     = shift;
    my $table_name = shift;
    my $item       = shift;

    if ( defined $buffer ) {
        push @$buffer, $item if ( defined $item );
        if ( @$buffer >= $DB_BATCH_SZ or not defined $item ) {
            print STDOUT "Populate $table_name " . @$buffer . "\n";
            my $startTime = time;
#            $coge->resultset($table_name)->populate($buffer) if (@$buffer); # mdb removed 1/7/14 -- defaulting to single-insert due to missing primary key value, see http://search.cpan.org/~ribasushi/DBIx-Class-0.082810/lib/DBIx/Class/ResultSet.pm#populate
            my $dbh = $coge->storage->dbh;
            my @columns = keys %{$buffer->[0]};
            my $stmt = "INSERT $table_name ( " . join(',', @columns) . " ) VALUES " . # TODO use prepared statement
                join(',', map { '(' . join(',', map { $dbh->quote($_) } @{$_}{@columns}) . ')' } @$buffer) . ';'; # mdb changed 4/28/15 -- fix ordering of values
#            print $log $stmt, "\n";
            unless ($dbh->do($stmt)) {
                print STDOUT "log: error: database batch insertion for table '$table_name' failed: $stmt\n";
                exit(-1);
            }
            print STDOUT "Time: ", (time - $startTime), "\n";
            @$buffer = ();
        }
    }
}

#sub batch_add_async {
#    my $buffer     = shift;
#    my $table_name = shift;
#    my $item       = shift;
#
#    if ( defined $buffer ) {
#        push @$buffer, $item if ( defined $item );
#        if ( @$buffer >= $DB_BATCH_SZ or not defined $item ) {
#            print STDOUT "Async populate $table_name " . @$buffer . "\n";
#            if ( !defined( my $child_pid = fork() ) ) {
#                print STDOUT "Cannot fork: $!";
#                batch_add(@_);
#                return;
#            }
#    	    elsif ( $child_pid == 0 ) {
#    	       print STDOUT "child running to populate $table_name\n";
#    	       $coge->resultset($table_name)->populate($buffer) if (@$buffer);
#    	       exit;
#    	    }
#            @$buffer = ();
#        }
#    }
#}

sub process_gff_file {
    print STDOUT "process_gff_file: $data_file\n";

    open( my $in, $data_file ) || die "can't open $data_file for reading: $!";

    my $line_num   = 0;
    my $gene_count = 0;
    my $last_RNA   = "mRNA"; #storage for the last RNA type seen.  For converting exons to appropriate RNA type.
    while ( my $line = <$in> ) {
        $line_num++;
        last if $line =~ /^##FASTA/; # mdb added 9/28/15 COGE-656
        next if $line =~ /^#/;
        next if $line =~ /^Error/;
        chomp $line;
        next unless $line;
        print STDOUT "log: Processed " . commify($line_num) . " lines\n"
          unless $line_num % 100000;

        my @line = split( /\t/, $line );
        if ( @line != 9 ) {
            log_line("Incorrect format (wrong number of columns, expecting 9)", $line_num, $line);
            return 0;
        }
        my ($chr, $type, $start, $stop, $strand, $attr) = ($line[0], $line[2], $line[3], $line[4], $line[6], $line[8]);

        # Ignore these types
        #next if $type eq "";
        next if $type eq "clone";
        next if $type eq "intron";
        next if $type eq "chromosome";
        next if $type eq "start_codon";
        next if $type eq "stop_codon";
#        next if $type eq "transcript";
        next if $type eq "protein";

        # Process and check chromosomes
        ($chr) = split(/\s+/, $chr);
        $chr = fix_chromosome_id($chr, \%valid_chrs);
        unless ( defined $chr && $valid_chrs{$chr} ) {
            log_line("Chromosome '$chr' does not exist in the dataset", $line_num, $line);
            next if ($allow_all_chr);
            return 0;
        }

        $type = "mRNA" if $type eq "transcript";
        # In many GFF files, the mRNA is what CoGe calls a Gene (the full extent
        # of the transcribed sequence including introns and exons.  Instead,
        # what the GFF calls an exon is really the transcribed mRNA.  In this
        # cases, we want to hold the mRNA information to link to Parents and
        # whatever annotation it contains, but don't want to actually add the
        # location.  We will change the feature type to something weird that
        # can be handled downstream correctly -- specifically the locations
        if ( $type =~ /([m]RNA)/ ) { # mdb changed from /(.*RNA.*)/, 9/3/13 issue 198 # mdb changed from [mt]RNA 7/14/14
            $last_RNA = $type;
            $type     = "$1_no_locs";
        }
        $type = $last_RNA if $type eq "exon";
        $type = $last_RNA if $type eq "five_prime_UTR";
        $type = $last_RNA if $type eq "three_prime_UTR";

        $seen_types{$type}++;

        my %names;
        my $name;
        my ( $parent, $id );
        foreach my $item ( split( /;/, $attr ) ) {
            $item =~ s/"//g;
            $item =~ s/^\s+//;
            $item =~ s/\s+$//;
            next unless $item;

            my ( $key, $value ) = ( split( /[\s=]/, $item, 2 ) );
            $seen_attr{$key}++;
            $parent = $value if $key eq "Parent";
            $id     = $value if $key eq "ID";
            next if $skip_attr{$key};

            if ( $check_names{$key} ) {
              outer:
                foreach my $item ( split( /,/, $value ) ) {
                    $names{$item} = 1;

                    # these nexts will skip from using the primary name as the ID to the Parent name
                    foreach my $re (@skip_names_re) {
                        next outer if $item =~ /$re/i;
                    }
                    $name = $item unless $name;
                    if ( $item =~ /^LOC_/ ) {
                        my $tmp = $item;
                        $tmp =~ s/^LOC_//;
                        $names{$tmp} = 1;
                    }
		            $name =~ s/_cds_\d+$//i;
		            $name =~ s/_exon_\d+$//i;
                }
            }
            next unless $name; # no name, don't know what to do!
            $value = uri_unescape($value); # remove URL formatting
            $annos{$name}{$value} = { type => $key } if $anno_names{$key};
        }
        next unless $name; # no name, don't know what to do!

        if    ( $strand =~ /-/ )  { $strand = -1; }
        elsif ( $strand =~ /\./ ) { $strand = 0;  }
        else                      { $strand = 1;  } # mdb changed 11/7/13 issue 248 - made '+' strand the default

        my @types = ($type);
        #push @types, "CDS" if $add_cds && $type eq "mRNA";
        # phytozome replications of CDS to mRNA
        #push @types, "mRNA" if $type eq "CDS";
        #push @types, "mRNA" if $type =~ /UTR/;
        # replicate mRNA to gene
        #push @types, "gene" if $type eq "mRNA";

        foreach my $tmp (@types) {
            my $tmp_name = $name;
            my $type     = $tmp;
            $type =~ s/_no_locs//;
            $tmp_name = $parent if $type eq "gene" && $parent; # ugly hack
            #print join ("\t", $name, $tmp_name),"\n";

            # initialize data structure
            $data{$chr}{$tmp_name}{$type} = {}
              unless $data{$chr}{$tmp_name}{$type};
            foreach my $n ( keys %names ) {
                $data{$chr}{$tmp_name}{$type}{names}{$n} = 1;
            }

            # mdb added 4/8/14 issue 358 - save location for later
#            if ($tmp =~ /_no_locs/) {
#                $data{$chr}{$tmp_name}{$type}{coords} =
#                  {
#                    start  => $start,
#                    stop   => $stop,
#                    strand => $strand,
#                    chr    => $chr
#                };
#            }

            next if ($tmp =~ /_no_locs/); # skip adding locations for things like mRNA
            push @{ $data{$chr}{$tmp_name}{$type}{loc} }, 
              {
                start  => $start,
                stop   => $stop,
                strand => $strand,
                chr    => $chr
              };
        }
        $total_annot++;
    }

    close($in);
    print STDOUT "log: Processed " . commify($line_num) . " total lines\n";

    return $total_annot;
}

sub log_line {
    my ( $msg, $line_num, $line ) = @_;
    print STDOUT "log: error at line $line_num: $msg\n", "log: ", substr($line, 0, 100), "\n";    
}

sub execute { # FIXME move into Util.pm
    my $cmd = shift;
    print STDOUT "$cmd\n";
    my @cmdOut    = qx{$cmd};
    my $cmdStatus = $?;
    if ( $cmdStatus != 0 ) {
        print STDOUT "log: error: command failed with rc=$cmdStatus: $cmd\n";
        exit(-1);
    }
}
