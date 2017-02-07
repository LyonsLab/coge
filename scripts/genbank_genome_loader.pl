#!/usr/bin/perl -w

use strict;
use LWP::Simple;
use Data::Dumper;
use Getopt::Long;
use CoGeX;
use File::Path;
use File::Touch;
use File::Spec::Functions qw( catdir catfile );
use CoGe::Accessory::Web;
use CoGe::Core::Storage;
use CoGe::Accessory::GenBank;
use POSIX qw(ceil);

my (
    $DEBUG,         $GO,              $autoupdate,  $autoskip,
    @accns,         $staging_dir,     $help,        $user_chr,
    $ds_link,       $delete_src_file, $test,        $auto_increment_chr,
    $base_chr_name, $accn_file,       $max_entries, $user_id,
    $user_name,     $config,          $host,        $port,
    $db,            $dbuser,          $pass,        $install_dir,
    $server,        $force,           $result_dir,  $wid
);

GetOptions(
    "debug"                     => \$DEBUG,
    "go"                        => \$GO,
    "accn|a=s"                  => \@accns,
    "staging_dir=s"             => \$staging_dir,
    "help|h"                    => \$help,
    "user_chr|chr=s"            => \$user_chr,
    "base_chr_name|basechr=s"   => \$base_chr_name,
    "dataset_link=s"            => \$ds_link,
    "test"                      => \$test,          #to add the name test to dataset name for testing purposes
    "autoupdate"                => \$autoupdate,    #automatically say 'yes' to any question about proceeding
    "autoskip"                  => \$autoskip,      #automatically say 'no' to any question about proceeding
    "delete_src_file"           => \$delete_src_file,
    "auto_increment_chr"        => \$auto_increment_chr,
    "accn_file|af=s"            => \$accn_file,
    "max_entries=i"             => \$max_entries,
    "user_id|uid=i"             => \$user_id,       #CoGe user id to which genome is associated
    "user_name=s"               => \$user_name,     #or CoGe user name to which genome is associated
    "wid=s"                     => \$wid,            # workflow id
    "install_dir=s"             => \$install_dir,   #optional install path for sequences
#    "result_dir=s"              => \$result_dir,    #optional result path for JSON info
    "force|f"                   => \$force,         #force the install of a the genome
    "config=s"                  => \$config
);

$| = 1;
print STDOUT "Starting $0 (pid $$)\n", qx/ps -o args $$/;

help() if $help;

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

$user_chr = 1 unless $user_chr;
$DEBUG    = 0 unless defined $DEBUG; # set to 1 to enable debug printing
$GO       = 0 unless defined $GO;    # set to 1 to actually make db calls (instead of testing)
$autoupdate = 0 unless $autoupdate;
$autoskip   = 0 unless $autoskip;

if ($config) {
    my $P = CoGe::Accessory::Web::get_defaults($config);

    #database
    $db   = $P->{DBNAME};
    $host = $P->{DBHOST};
    $port = $P->{DBPORT};
    $dbuser = $P->{DBUSER};
    $pass = $P->{DBPASS};

    #other stuff
    $install_dir = $P->{SEQDIR};
    $server      = $P->{SERVER};
}

# Connect to database
my $connstr = "dbi:mysql:dbname=$db;host=$host;port=$port;";
my $coge = CoGeX->connect( $connstr, $dbuser, $pass );
unless ($coge) {
    print STDOUT "log: couldn't connect to database\n";
    exit(-1);
}

# Retreive User
my $user;
if ($user_id) {
    $user = $coge->resultset('User')->find($user_id);
}
else {
    $user = $coge->resultset('User')->find( { user_name => $user_name } );
}
unless ($user) {
    print STDOUT "log: error finding user '$user_name'\n";
    exit(-1);
}

print STDOUT "Go = $GO.  Will be adding to the databaes. \n" if $DEBUG;
print STDOUT "Force is on: will not be checking if genome has been previously loaded.\n"
  if $force;
my $data_source = get_data_source();    #for NCBI

# loop through all the accessions
my %previous_datasets;    #storage for previously loaded accessions
my $genome;               #storage for coge genome_obj
my %chromosomes;          #storage for used chromosome names;

# Load optional accn file
if ( $accn_file && -r $accn_file ) {
    open( IN, $accn_file );
    while (<IN>) {
        next if /#^/;
        s/'//g;
        chomp;
        foreach my $item ( split(/\s|,/) ) {
            next unless $item;
            push @accns, $item;
        }
    }
}

my $fasta_output;       # storage of main fasta output
my @links_to_existing;  # links to previously loaded datasets

accn: foreach my $accn (@accns) {
    print STDOUT "log: Working on $accn...\n";
    unless ($force) {
        my $previous = check_accn($accn);
        foreach my $item (@$previous) {
            if ( !$item->{version_diff} && !$item->{length_diff} ) {
                #push @previous_datasets, $item->{ds};
                $previous_datasets{ $item->{ds}->id } = $item->{ds};
                my $link = $server . "OrganismView.pl?dsid=" . $item->{ds}->id;
                print STDOUT "log: Dataset previously loaded: <a href='$link'>$link</a>.  Skipping ...\n";
                push @links_to_existing, $link;
                next accn;
            }
            elsif ( !$item->{version_diff} && $item->{length_diff} ) {
                print STDOUT
                  "Detected a difference in total genomic length between CoGe ("
                  . $item->{coge_length}
                  . ") and NCBI("
                  . $item->{ncbi_length}
                  . ").  Including new dataset\n";
            }
        }
    }

    my $genbank = new CoGe::Accessory::GenBank();
    $genbank->debug(1);
    $genbank->logfile(*STDOUT);
    $genbank->max_entries($max_entries) if $max_entries;
    $genbank->get_genbank_from_ncbi(
        file => "$staging_dir/$accn.gbk",
        accn => $accn
    );
    if ( !$genbank->sequence && !@{ $genbank->wgs_data } ) {
        print STDOUT "Skipping sequence.  No sequence or wgs data\n";
        next;
    }
    my $EXIT = 0;
    my $chromosome;
    $chromosome = $user_chr if $user_chr;
    $user_chr++ if $auto_increment_chr;
    if ( $genbank->chromosome ) {
        $chromosome = $genbank->chromosome;
        print STDOUT "#" x 20, "\n";
        print STDOUT $genbank->data_source();
        print STDOUT "GBChr:  $chromosome\n";
        print STDOUT "#" x 20, "\n";
    }
    $chromosome =~ s/chromosome//i;
    $chromosome =~ s/chr//i;
    $chromosome =~ s/^\s+//;
    $chromosome =~ s/\s+$//;
    $chromosome =~ s/\s+/_/g;
    $chromosome =~ s/\.0$//;
    $chromosome = $base_chr_name . $chromosome if $base_chr_name;

    if ( $chromosomes{$chromosome} ) {
        print STDOUT "log: previously seen $chromosome.  Updating name.\n";
        while ( $chromosomes{$chromosome} ) {
            $chromosome .= ".1";
        }
    }
    $chromosomes{$chromosome} = 1;
    print STDOUT "\tchromosome: $chromosome", "\n";
    my ( $organism, $dataset );

    unless ($organism) {
        ($organism) = get_organism($genbank);
        if ($organism) {
            print STDOUT "Organism info:";
            print STDOUT "\t", $organism->id,          ": ";
            print STDOUT "\t", $organism->name,        "\n";
            print STDOUT "\t", $organism->description, "\n";
        }
        else {
            print STDOUT "WARNING:  Unable to retrieve an organism object for $accn.  Probably a problem with the genbank object file parsing\n"
              unless !$GO;
            next if $GO;
        }
    }
    my $dataset_desc = "LOCUS: " . $genbank->locus();
    $dataset_desc .= ", ACCESSION: " . $genbank->accession();
    $dataset_desc .= ", VERSION: " . $genbank->version();

    # testing to see if already in database
    my $version = $genbank->version();
    my $name    = $genbank->accession;
    $name .= ".gbk" unless $name =~ /gbk$/;
    $name .= ".test" if $test;

    # actually create a dataset object
    $ds_link = "http://www.ncbi.nlm.nih.gov/entrez/viewer.fcgi?db=nucleotide&dopt=gbwithparts&list_uids=" . $accn unless $ds_link;
    $dataset = $coge->resultset('Dataset')->create(
        {
            name           => $name,
            description    => $dataset_desc,
            link           => $ds_link,
            data_source_id => $data_source->id(),
            version        => $genbank->version,
        }
      )
      if $GO;

    ### Main Feature Processing Loop ###
    my @gbs;
    if ( @{ $genbank->wgs_data } ) {
      entry: foreach my $entry ( @{ $genbank->wgs_data } ) {
            my $accn = $entry->accession;
            my $try  = 1;
            while ( !$accn ) {

                $entry->get_genbank_from_ncbi( reload => 1 );
                print STDOUT "#" x 20, "\n";
                print STDOUT "Warning.  Didn't not retrieve a valid accession for entry.\n";
                print STDOUT "Requested id: " . $entry->requested_id . "\n";
                print STDOUT "Trying to retrieve valid entry.\n";
                print STDOUT "#" x 20, "\n";
                $accn = $entry->accession;
                $try++;

                if ( $try >= 10 ) {
                    print STDOUT "Giving up!  Skipping "
                      . $entry->requested_id
                      . " after $try trys.\n";
                    next entry;
                }
            }
            print STDOUT "\tChecking WGS $accn...\n";
            unless ($force) {
                my $previous = check_accn($accn);
                foreach my $item (@$previous) {
                    if ( !$item->{version_diff} && !$item->{length_diff} ) {
                        $previous_datasets{ $item->{ds}->id } = $item->{ds};
                        my $link = $server . "OrganismView.pl?dsid=" . $item->{ds}->id;
                        print STDOUT "log: Dataset previously loaded: $link\n";
                        next entry;
                    }
                    elsif ( !$item->{version_diff} && $item->{length_diff} ) {
                        print STDOUT "Detected a difference in total genomic length between CoGe ("
                          . $item->{coge_length}
                          . ") and NCBI("
                          . $item->{ncbi_length}
                          . ").  Including new dataset.\n";
                    }
                    else {
                        print STDOUT "not present.  Will be loaded.\n";
                    }
                }
                push @gbs, $entry;
            }
        }
        #@gbs = @{$genbank->wgs_data};
    }
    else {
        push @gbs, $genbank;
    }

    foreach my $entry (@gbs) {
        $chromosome = "contig_" . $entry->accession if @{ $genbank->wgs_data };
        print STDOUT "Processing features for " . $entry->accession . "...\n"
          unless $EXIT;
        foreach my $feature ( @{ $entry->features() } ) {
            unless ( $feature->type() ) {
                print STDOUT
                  "Feature has no feature type name \$feature->type():\n";
                print STDOUT Dumper $feature;
                next;
            }
            if ( $feature->type() =~ /source/i ) {
                #change source to chromosome
                $feature->type('chromosome');
                #next;
            }
            my $feat_type =
              $coge->resultset('FeatureType')->find_or_create( { name => $feature->type() } )
              if $GO;

           # create a db_feature for to link this feature with the dataset table
            my ( $start, $stop, $strand ) = get_feature_location($feature);
            my $db_feature = $coge->resultset('Feature')->create(
                {
                    feature_type_id => $feat_type->id,
                    dataset_id      => $dataset->id,
                    chromosome      => $chromosome,
                    strand          => $strand,
                    start           => $start,
                    stop            => $stop,
                }
              ) if $GO;

            # expect first feature to be the source feature!
            if ( $feature->type() =~ /chromosome/i ) {
                # generate name based on chromosome
                my $feat_name = $coge->resultset('FeatureName')->create(
                    {
                        name => $chromosome,
                        #description => "Chromosome " . $chromosome,
                        feature_id => $db_feature->id
                    }
                  ) if $GO;

                # generate name for accession
                $feat_name = $coge->resultset('FeatureName')->create(
                    {
                        name       => $entry->accession,
                        feature_id => $db_feature->id()
                    }
                  ) if $GO;

                # generate name for version
                $feat_name = $coge->resultset('FeatureName')->create(
                    {
                        name       => $entry->accession . "." . $entry->version,
                        feature_id => $db_feature->id
                    }
                  ) if $GO;

                # generate name for GI
                if ($entry->gi) { # mdb added if condition 2/7/17
                    $feat_name = $coge->resultset('FeatureName')->create(
                        {
                            name        => $entry->gi,
                            description => "GI number",
                            feature_id  => $db_feature->id
                        }
                    ) if $GO;
                }
                else {
                    print STDOUT "Missing gi for entry ", $entry->accession, "\n";
                }
            }

            # add a location entry
            my $loc_string = $feature->location;
            $loc_string =~ s/complement//g;
            $loc_string =~ s/join//;
            $loc_string =~ s/order//;
            $loc_string =~ s/\(|\)//g;

            # loop through the locatons
            foreach my $loc ( split(/,/, $loc_string) ) {
                $loc =~ s/<|>//g;
                my ( $start, $stop ) = split /\.\./, $loc;
                $stop = $start unless $stop;
                $start =~ s/\^.*//;
                $stop  =~ s/\^.*//;

                die "problem with $accn start $start or stop $stop\n"
                  unless $start =~ /^\d+$/ && $stop =~ /^\d+$/;
                my $location = $db_feature->add_to_locations(
                    {
                        start      => $start,
                        stop       => $stop,
                        strand     => $feature->strand,
                        chromosome => $chromosome
                    }
                  ) if $GO;
            }

            # now work through the qualifiers for this feature
            # start by getting the hashref of qualifiers
            my $annot = $feature->qualifiers();
            my %names;

            #print STDOUT Dumper $annot;
            foreach my $anno ( keys %{$annot} ) {
                my $stuff = $annot->{$anno};

                # deal with db_xref: (taxon:3702) (GeneID:821318) (GI:18379324)
                if ( $anno =~ /xref/i ) {
                    my $anno_type_group =
                      $coge->resultset('AnnotationTypeGroup')->find_or_create( { name => $anno } )
                      if $GO;

                    # go through each of the entries in the db_xref qualifier values and split on ':', then add entries individually
                    foreach my $xref ( @{$stuff} ) {
                        my @inner = split( /:/, $xref );

                        # first add the annot_type_obj
                        my ($anno_type) =
                          $coge->resultset('AnnotationType')->search(
                            {
                                name                     => $inner[0],
                                annotation_type_group_id =>
                                  $anno_type_group->id()
                            }
                          );
                        ($anno_type) =
                          $coge->resultset('AnnotationType')->create(
                            {
                                name                     => $inner[0],
                                annotation_type_group_id =>
                                  $anno_type_group->id()
                            }
                          ) if $GO && !$anno_type;

                        # now create the row for the data value of the xref
                        my $sub_anno = $db_feature->add_to_feature_annotations(
                            {
                                annotation         => $inner[1],
                                annotation_type_id => $anno_type->id()
                            }
                          ) if $GO;
                    }
                }
                elsif (
                       $anno =~ /locus_tag/i
                    || $anno =~ /transcript_id/i
                    || $anno =~ /protein_id/i
                    || $anno =~ /gene/i
                    || $anno =~ /standard_name/i
                    || $anno =~
                    /synonym/i # synonyms are embedded in the /note= tag! these are names
                    || $anno eq "names"
                  )
                {
                    my $master = 1;    #make first one master
                    foreach my $item ( @{$stuff} ) {
                        foreach my $thing ( split( /;/, $item ) ) {
                            $thing =~ s/^\s+//;
                            $thing =~ s/\s+$//;
                            $names{$thing} = 0 unless defined $names{$thing};
                            $names{$thing} = 1
                              if $anno =~ /locus_tag/i;    #primary_name;
                        }
                    }
                }
                elsif ( $anno =~ /translation/i
                  )    # this needs to be entered into the sequence table
                {
                    next
                      ; #skip this.  Protein sequences are translated on the fly from DNA sequence
                    my $seq_type =
                      $coge->resultset('SequenceType')->find_or_create(
                        {
                            name        => "protein",
                            description => "translation"
                        }
                      )
                      if $GO;
                    foreach my $item ( @{$stuff} ) {
                        $item =~ s/\s+//g;
                        my $sequence = $db_feature->add_to_sequences(
                            {
                                sequence_type_id => $seq_type->id(),
                                sequence_data    => $item,
                            }
                          )
                          if $GO;
                    }
                }
                elsif ( $anno eq "note" ) {

                    # if go annot are present, they'll be in the note qualifier,
                    # so process is specifically
                    foreach my $item ( @{$stuff} ) {
                        my $leftover = "";
                        my @temp = split( /;/, $item );
                        foreach my $go_raw (@temp) {
                            if ( $go_raw =~ /go_/ ) {
                                while ( $go_raw =~
                                    /(go_.*?):\s+(.*?)\[goid G?O?:?(.*?)\]/g )
                                {

                              # example:
                              # go_function: nucleic acid binding [goid 0003676]
                                    my $anno_type_group =
                                      $coge->resultset('AnnoationTypeGroup')
                                      ->find_or_create( { name => $1 } )
                                      if $GO;

                                    # $1 should be "go_function"
                                    my ($anno_type) =
                                      $coge->resultset('AnnotationType')
                                      ->search(
                                        {
                                            name                     => $3,
                                            annotation_type_group_id =>
                                              $anno_type_group->id()
                                        }
                                      );
                                    ($anno_type) =
                                      $coge->resultset('AnnotationType')
                                      ->create(
                                        {
                                            name                     => $3,
                                            annotation_type_group_id =>
                                              $anno_type_group->id()
                                        }
                                      )
                                      if $GO && !$anno_type;
                                    my $sub_anno =
                                      $db_feature->add_to_feature_annotations(
                                        {
                                            annotation => $2
                                            , #this should be "nucleic acid binding"
                                            annotation_type_id => $anno_type->id
                                        }
                                      )
                                      if $GO;
                                }
                            }
                            else {
                                $leftover .= " " . $go_raw if $go_raw;
                            }

                            # now just add the note remainder
                            $leftover =~ s/^\s+//;
                            $leftover =~ s/\s+$//;
                            if ($leftover) {
                                my ($anno_type) =
                                  $coge->resultset('AnnotationType')
                                  ->search( { name => $anno } );
                                ($anno_type) =
                                  $coge->resultset('AnnotationType')
                                  ->create( { name => $anno } )
                                  if $GO && !$anno_type;

                                my $sub_anno =
                                  $db_feature->add_to_feature_annotations(
                                    {
                                        annotation         => $leftover,
                                        annotation_type_id => $anno_type->id(),
                                    }
                                  )
                                  if $GO;
                            }
                        }
                    }
                }
                else    ##everything else
                {
                    foreach my $item ( @{$stuff} ) {
                        my ($anno_type) =
                          $coge->resultset('AnnotationType')
                          ->search( { name => $anno } );
                        ($anno_type) =
                          $coge->resultset('AnnotationType')
                          ->create( { name => $anno } )
                          if $GO && !$anno_type;
                        my $sub_anno = $db_feature->add_to_feature_annotations(
                            {
                                annotation         => $item,
                                annotation_type_id => $anno_type->id(),
                            }
                          )
                          if $GO;
                    }
                }
            }
            foreach my $name ( keys %names ) {
                my $master = $names{$name} ? 1 : 0;
                $name =~ s/\s+$//g;
                $name =~ s/^\s+//g;
                my $feat_name = $db_feature->add_to_feature_names(
                    {
                        name         => $name,
                        primary_name => $master,
                    }
                  )
                  if $GO;
            }
        }

        #initialize genome object if needed
        if ( $organism && !$genome )    #gst_id 1 is for unmasked sequence data
        {
            print STDOUT "Creating Genome Object ...\n";
            $genome = generate_genome(
                version => $version,
                org_id  => $organism->id,
                user_id => $user->id,
                gst_id  => 1,
            );
            unless ($genome) {
                print STDOUT "Error adding genome to database\n";
                exit(-1);
            }

            print STDOUT "log: Added genome id", $genome->id, "\n"; # !!!! don't change, gets parsed by calling code
            print STDOUT "Creating install path for sequences ...\n";
            $install_dir = catdir($install_dir, CoGe::Core::Storage::get_tiered_path( $genome->id ));
            if ($GO) {
                mkpath($install_dir);
                unless ( -d $install_dir ) {
                    print STDOUT "log: error in mkpath $install_dir\n";
                    exit(-1);
                }
                print STDOUT "completed: $install_dir\n";
            }
        }

        #fasta_format sequence for the chromosome
        $fasta_output .= fasta_genomic_sequence(
            genome => $genome,
            seq    => $entry->sequence,
            chr    => $chromosome
        );

        if ($GO) {
            my $load = 1;
            foreach my $dsc ( $genome->dataset_connectors
              ) #check to see if there is a prior link to the dataset -- this will happen when loading whole genome shotgun sequence
            {
                $load = 0 if $dsc->dataset_id == $dataset->id;
            }
            $genome->add_to_dataset_connectors( { dataset_id => $dataset->id } )
              if $load;
            if ( $dataset->version > $genome->version ) {
                $genome->version( $dataset->version );
                $genome->update;
            }
        }
    }
    print STDOUT "completed parsing for $accn!\n";    # if $DEBUG;

    if ($delete_src_file) {
        print STDOUT "Deleting genbank src file: " . $genbank->srcfile . "\n";
        my $cmd = "rm " . $genbank->srcfile;
        `$cmd`;
        #`rm /tmp/gb/*`;
    }
}

#need to add previous datasets if new dataset was added to a genome
if ($GO) {
    if ( $genome && keys %previous_datasets ) {
        my $ver;    #need a higher version number than previous
        foreach my $ds ( values %previous_datasets ) {
            my ($test) =
              $genome->dataset_connectors( { dataset_id => $ds->id } );
            if ($test) {
                my $name = $ds->name;
                print STDOUT "$name has been previously added to this dataset group.  Skipping\n";
                next;
            }
            foreach my $item ( $ds->genomes ) {
                $ver = $item->version unless $ver;
                $ver = $item->version if $item->version > $ver;
            }
            foreach my $chr ( $ds->chromosomes ) {
                $fasta_output .= fasta_genomic_sequence(
                    genome => $genome,
                    seq    => $ds->get_genomic_sequence( chr => $chr ),
                    chr    => $chr
                );
            }
            $genome->add_to_dataset_connectors( { dataset_id => $ds->id } );
        }

        #incement and update genome version if new version number is higher
        $ver++;
        if ( $ver > $genome->version ) {
            $genome->version($ver);
            $genome->update();
        }
    }
}
my $output_file = catfile($install_dir, "genome.faa");
print STDOUT "log: Creating output sequence file and indexing\n";
print STDOUT "file: $output_file\n";
add_and_index_sequence( fasta => $fasta_output, file => $output_file ) if $GO;

# Make user owner of new genome
if ( $GO && $user ) {
    # mdb removed 6/17/14 issue 394 - don't create a user connector for NCBI loaded genomes
    #my $node_types = CoGeX::node_types();
    #my $conn       = $coge->resultset('UserConnector')->create(
    #   {
    #       parent_id   => $user->id,
    #       parent_type => $node_types->{user},
    #       child_id    => $genome->id,
    #       child_type  => $node_types->{genome},
    #       role_id     => 2                       # FIXME hardcoded creator role
    #   }
    #);
    #unless ($conn) {
    #   print STDOUT "log: error creating user connector\n";
    #   exit(-1);
    #}
}

# Save result
unless (add_workflow_result($user_name, $wid, 
        {
            type           => 'genome',
            id             => int($genome->id),
            name           => $genome->name,
            description    => $genome->description,
            info           => $genome->info,
            version        => $genome->version,
            link           => $genome->link,
            links          => \@links_to_existing,
            organism_id    => $genome->organism_id,
            restricted     => $genome->restricted
        })
    )
{
    print STDOUT "log: error: could not add workflow result\n";
    exit(-1);
}

# Add genome ID to log - mdb added 7/8/15, needed after log output was moved to STDOUT for jex
my $logtxtfile = "$staging_dir/log.txt";
open(my $logh, '>', $logtxtfile);
print $logh "genome id: " . $genome->id . "\n";
close($logh);

# Save workflow_id in genome data path -- #TODO move into own routine in Storage.pm
CoGe::Accessory::TDS::write(
    catfile($install_dir, 'metadata.json'),
    {
        workflow_id => int($wid)
    }
);

# Create "log.done" file to indicate completion to JEX
touch($logdonefile);

exit;

#-------------------------------------------------------------------------------

sub add_and_index_sequence {
    my %opts     = @_;
    my $fasta    = $opts{fasta};
    my $file     = $opts{file};
    my $compress = $opts{compress};
    unless ($fasta) {
        print STDOUT
          "log: No data to add to file $file.  Not creating fasta file\n";
        exit;
    }
    if ( -r $file ) {
        print STDOUT "log: error:  $file already exists.  Will not overwrite existing sequence.  Fatal error!\n";
        exit(-1);
    }
    open( OUT, ">$file" ) || die "Died: can't open $file for writing: !$\n";
    print OUT $fasta;
    close OUT;
    print STDOUT "Indexing genome file\n";
    my $rc = CoGe::Core::Storage::index_genome_file(
        file_path => $file,
        compress  => $compress
    );
    if ( $rc != 0 ) {
        print STDOUT "log: error: couldn't index fasta file\n";
        exit(-1);
    }
}

sub fasta_genomic_sequence {
    my %opts   = @_;
    my $seq    = $opts{seq};
    my $chr    = $opts{chr};
    my $genome = $opts{genome};

    #    my $file = $opts{file};
    return unless $genome && $seq && $chr;

    $seq =~ s/\s//g;
    $seq =~ s/\n//g;

    my $seqlen = length $seq;
#    if ( my ($item) = $genome->get_genomic_sequence( chromosome => $chr ) ) {
#        my $prev_length = $item->sequence_length;
    my $prev_length = $genome->get_chromosome_length($chr);
    if ( $prev_length ) {
        print STDOUT "$chr has previously been added to this genome.  Previous length: $prev_length.  Currently length: $seqlen.  Skipping.\n";
        return;
    }
    print STDOUT "Loading genomic sequence ($seqlen nt)\n";

# no longer add to db, we use the index file instead. not sure if the rest of fasta_genomic_sequence is necessary
#    $genome->add_to_genomic_sequences(
#        {
#            sequence_length => $seqlen,
#            chromosome      => $chr,
#        }
#      )
#      if $GO;

    my $head = $chr =~ /^\d+$/ ? ">gi" : ">lcl";
    $head .= "|" . $chr;

    #    open( OUT, ">>" . $file );
    #    print OUT "$head\n$seq\n";
    #    close OUT;
    return "$head\n$seq\n";
}

sub get_feature_location {
    my $feat       = shift;
    my $loc_string = $feat->location;
    my $strand     = $feat->location =~ /complement/ ? "-1" : 1;
    $loc_string =~ s/complement//g;
    $loc_string =~ s/join//;
    $loc_string =~ s/order//;
    $loc_string =~ s/\(|\)//g;
    my ( $rstart, $rstop );

    foreach my $loc ( split /,/, $loc_string ) {
        my ( $start, $stop ) = split /\.\./, $loc;
        $stop = $start unless $stop;
        ($start) = $start =~ /(\d+)/;
        ($stop)  = $stop  =~ /(\d+)/;
        $start =~ s/\^.*//;
        $stop  =~ s/\^.*//;
        $rstart = $start unless $rstart;
        $rstop  = $stop  unless $rstop;
        $rstart = $start if $start < $rstart;
        $rstop  = $stop  if $stop > $rstop;
    }
    return ( $rstart, $rstop, $strand );
}

sub get_organism {
    my ($entry) = shift;
    my $name = $entry->data_source();
    if ( $entry->strain ) {
        my $strain = $entry->strain;
        $strain =~ s/\(/;/g;
        $strain =~ s/\)//g;
        $strain =~ s/strain://g;
        $strain =~ s/\\//g;
        $strain =~ s/=/;/g;
        $name   =~ s/strain//;
        $name   =~ s/\sstr\.?\s/ /;
        my @strains = split /;/, $strain;
        my @parse_strains;

        foreach ( sort @strains ) {
            s/str\.\s//;
            s/^\s+//;
            s/\s+$//;
            next unless $_;
            $name =~ s/$_//;
            push @parse_strains, $_;
        }
        $name =~ s/\s\s+/ /g;
        $name =~ s/\s+$//;
        my $add = join( "; ", sort @parse_strains );
        $name .= " strain" unless $add =~ /strain/;
        $name .= " $add";
    }
    if ( $entry->substrain ) {
        my $sstrain = $entry->substrain;
        $sstrain =~ s/\(/\\\(/g;
        $sstrain =~ s/\)/\\\)/g;
        $sstrain =~ s/substrain://g;
        $name    =~ s/\ssubstr\.?\s/ /;
        $sstrain =~ s/\\//g;
        $sstrain =~ s/=/;/g;
        my @sstrains = split( /;/, $sstrain );
        my @parse_sstrains;

        foreach ( sort @sstrains ) {
            s/substr\.\s//;
            s/^\s+//;
            s/\s+$//;
            next unless $_;
            $name =~ s/$_//;
            push @parse_sstrains, $_;
        }
        $name =~ s/\s\s+/ /g;
        $name =~ s/\s+$//;
        my $add = join( "; ", sort @parse_sstrains );
        $name .= " substrain" unless $add =~ /substrain/;
        $name .= " $add";
    }
    return unless $name;
    $name =~ s/'//g;
    $name =~ s/\(\s*\)//g;
    $name =~ s/\s\s+/ /g;
    $name =~ s/^\s+//;
    $name =~ s/\s+$//;

    #  print STDOUT $name,"\n";
    #  print STDOUT $entry->organism,"\n";
    my $desc = $entry->organism();
    $desc =~ s/^.*?::\s*//;
    $desc =~ s/\.$//;
    print STDOUT qq{
Organism Information from Genbank Entry:
  $name
  $desc
} if $DEBUG;
    my $org =
      $coge->resultset('Organism')
      ->find( { name => $name, description => $desc } );
    unless ($org) {

        unless ($name) {
            print STDOUT "WARNING: ", $entry->accession,
              " has no organism name\n";
            return;
        }
        $org = $coge->resultset('Organism')->find_or_create(
            {
                name        => $name,
                description => $desc,
            }
          )
          if $GO && $name;
    }
    return $org;
}

sub get_data_source {
    return $coge->resultset('DataSource')->find_or_create(
        {
            name        => 'NCBI',
            description => "National Center for Biotechnology Information",
            link        => 'www.ncbi.nih.gov'
        }
    );
}

sub generate_genome {
    my %opts      = @_;
    my $name      = $opts{name};
    my $desc      = $opts{desc};
    my $version   = $opts{version};
    my $org_id    = $opts{org_id};
    my $gst_id    = $opts{gst_id};
    my $genome_id = $opts{genome_id};
    my $user_id   = $opts{user_id};
    my $genome    = $genome_id
      ? $coge->resultset('Genome')->find($genome_id)
      : $coge->resultset('Genome')->create(
        {
            name                     => $name,
            description              => $desc,
            version                  => $version,
            organism_id              => $org_id,
            creator_id               => $user_id,
            genomic_sequence_type_id => $gst_id,
        }
      )
      if $GO;
    return unless $genome;

    return $genome;
}

sub check_accn {
    my $accn = shift;
    my $gi   = get_gi($accn);

    #    print STDOUT "gi|".$gi."...";
    my $summary = get_gi_summary($gi);
    my ($version) = $summary =~ /ref\|.*?\.(\d+)\|/i;

    #    print STDOUT "version: $version\n";
    $version = 1 unless $version;
    my ($length) =
      $summary =~ /<Item Name="Length" Type="Integer">(\d+)<\/Item>/i;
    my ($taxaid) =
      $summary =~ /<Item Name="TaxId" Type="Integer">(\d+)<\/Item>/i;

    my @results;
    my %tmp;
    foreach my $ds ( $coge->resultset('Dataset')
        ->search( { name => $accn . ".gbk", version => $version } ) )
    {
        my $version_diff = $ds->version eq $version ? 0 : 1;
        my $length_diff;
        my $cogelength;

        foreach my $feat (
            $ds->features(
                { "feature_type.name" => "chromosome" },
                { "join"              => "feature_type" },
            )
          )
        {
            $cogelength += $feat->length;
        }
        $length_diff = $cogelength eq $length ? 0 : 1;
        push @results, {
            ds           => $ds,
            version_diff => $version_diff,
            length_diff  => $length_diff,
            coge_length  => $cogelength,
            ncbi_length  => $length,

        };
    }
    return \@results;
}

###NCBI eutils stuff

sub get_gi {
    my $accn    = shift;
    my $esearch = "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=nucleotide&term=";
    my $result = get( $esearch . "$accn" );
    my ($id) = $result =~ /<id>(.*?)<\/id>/i;
    return $id;
}

sub get_gi_summary {
    my $gi       = shift;
    my $esummary = "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=nucleotide&rettype=gb&retmode=text&complexity=0&id=";
    my $result = get( $esummary . $gi );
    return $result;
}

sub help {
    print qq
    {
Welcome to $0!  This program loads genbank entries into the CoGe genomes database

Example usage: genbank_genome_loader.pl -config /opt/apache/coge/web/coge.conf -go -accn <accn_1> -accn <accn_2>

Options:
 -help|h       :      Print this message

 -go (flag)    :      Make database inserts happen.  Otherwise, it is a fake load.  No options.

 -accn         :      Specify the genbank/refseq accession(s) to load.  It is recommended that a version
                      number is not used.  Multiple accessions for a single genome may be specified by:
                      -accn <accn_1> -accn <accn_2> -accn <accn_3>

 -accn_file    :      Specify a file of genbank/refseq accessions to use

 -config       :      Optional:  CoGe conf file for determining database connection stuff, temp_dir, server
                      for links, and sequence installation directory for fasta files.

 -temp_dir|staging_dir|td:  Directory to which genbank files are downloaded for processing.
                            Default: /tmp/ncbi/".ceil(rand(9999999999))

 -install_dir  :      Directory for installing sequences.  Includes indexing them.

 -force (flag) :      For the loading of the specified genbank accessions.  Otherwise, program will
                      check to see if they have been previously loaded for that version and not load
                      that file.  Also does a check to make sure that the specified length of the sequence
                      in the genbank file matches what has been previously loaded.  If they do not match,
                      the genbank file is loaded as a new version. No options.

 -base_chr_name|basechr : Specifies a string to add to the begining of a chromosome name

 -test (flag)  :      Add the extension ".test" to the dataset name in CoGe.  Used for testing stuff.

 Database connection stuff
 -host         :      Host for database connections
 -port         :      Port for database connections
 -database     :      Database name for inserts
 -user         :      User name to log into database
 -password     :      Password to log into database

    };
}
