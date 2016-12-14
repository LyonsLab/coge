#!/usr/bin/perl -w
use strict;
use CoGeX;
use CoGe::Accessory::Web qw(get_defaults get_command_path);
use CoGe::Accessory::Utils qw( trim execute );
use CoGe::Core::Notebook qw(add_items_to_notebook);
use File::Path;
use File::Basename;
use File::Spec::Functions qw(catdir catfile);
use File::Touch;
use File::Copy qw(copy);
use Getopt::Long;
use URI::Escape::JavaScript qw(escape unescape);
use Data::Dumper;

use vars qw(
    $staging_dir $file_str $notebook_name $notebook_desc $gid $nid $user_name
    $config $log_file $user $genome $notebook $result_dir $wid @failed_experiments
);

my $DELIMITER = '\t';

GetOptions(
    "staging_dir=s" => \$staging_dir,    # temporary staging path
#    "result_dir=s"  => \$result_dir,     # results path
    "files=s"       => \$file_str,       # input data file (JS escape)
    "name=s"        => \$notebook_name,  # notebook name (JS escaped) if not nid
    "desc=s"        => \$notebook_desc,  # notebook description (JS escaped)
    "gid=s"         => \$gid,            # genome id
    "nid=s"         => \$nid,            # optional notebook id, otherwise new notebook created
    "wid=s"         => \$wid,            # workflow id
    "user_name=s"   => \$user_name,      # user name
    "config=s"      => \$config,         # CoGe config file
);

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

# Process and verify parameters
$file_str      = unescape($file_str);
$notebook_name = unescape($notebook_name);
$notebook_desc = unescape($notebook_desc);

if ($user_name eq 'public') {
    print "log: error: not logged in\n";
    exit(-1);
}

unless ($gid) {
    print "log: error: genome id parameter 'gid' not specified\n";
    exit(-1);
}

# Load config file
my $P;
if ($config) {
    $P = CoGe::Accessory::Web::get_defaults($config);
}
else {
    $P = CoGe::Accessory::Web::get_defaults();
}

# Connect to database
my $connstr = "dbi:mysql:dbname=".$P->{DBNAME}.";host=".$P->{DBHOST}.";port=".$P->{DBPORT}.";";
my $coge = CoGeX->connect( $connstr, $P->{DBUSER}, $P->{DBPASS} );
unless ($coge) {
    print "log: error: couldn't connect to database\n";
    exit(-1);
}

# Retrieve genome
if ($gid) {
    $genome = $coge->resultset('Genome')->find( { genome_id => $gid } );
    unless ($genome) {
        print "log: error finding genome id$gid\n";
        exit(-1);
    }
}

# Retrieve notebook
if ($nid) {
    $notebook = $coge->resultset('List')->find( { list_id => $nid } );
    unless ($notebook) {
        print "log: error finding notebook id$nid\n";
        exit(-1);
    }
}

# Retrieve user
$user = $coge->resultset('User')->find( { user_name => $user_name } );
unless ($user) {
    print "log: error finding user '$user_name'\n";
    exit(-1);
}

# Process each file into staging area
my $data_dir = catdir($staging_dir, 'data');
mkpath($data_dir);
my @files = split( ',', $file_str );
foreach my $file (@files) {
    my $filename = basename($file);

    # Decompress file (if necessary) - removed, let the load_genome/experiment do this
#    if ( $file =~ /\.gz$/ ) {
#        print "log: Decompressing '$filename'\n";
#        $filename =~ s/\.gz$//;
#        run( $P->{GUNZIP}.' -c '.$file.' > '.catfile($data_dir, $filename) );
#        next;
#    }

    # Untar file (if necessary) - TODO: do this before gunzip if tar.gz file
    if ( $file =~ /\.tar$/ || $file =~ /\.tar\.gz$/ || $file =~ /\.tgz$/ ) {
        print "log: Extracting files\n";
        run( get_command_path('TAR').' -xf '.$file.' --directory '.$data_dir );
        next;
    }
    
    # Otherwise just copy the file as is
    run("cp $file $data_dir/$filename");
#    unless ( copy( $file, catfile($data_dir, $filename) ) ) {
#        print "log: error copying file:\n";
#        print "log: $!\n";
#        exit(-1);
#   }
}

# Find metadata file
my ($metadata_file) = glob("$data_dir/*.txt");
unless ($metadata_file) { # require a metadata file for now but could be optional in the future
    print "log: error: cannot find .txt metadata file\n";
    exit(-1);
}

# Load metadata file
my $metadata;
print "log: Loading metadata file '".fileparse($metadata_file)."'\n";
$metadata = get_metadata($metadata_file);
unless ($metadata) {
    print "log: error: no metadata loaded\n";
    exit(-1);
}
print "METADATA:\n", Dumper($metadata), "\n";

# Load each data file
($notebook, my $items) = process_dir($data_dir, $metadata, $notebook);

# Save result
unless (add_workflow_result($user_name, $wid, 
        {
            type        => 'notebook',
            id          => int($notebook->id),
            name        => $notebook->name,
            description => $notebook->description,
            genome_id   => int($genome->id)
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
# Format defined here:  http://genomevolution.org/wiki/index.php/Experiment_Metadata
sub get_metadata {
    my $file = shift;
    my $metadata = shift;
    my (@header, %data);
    my $lineNum = 0;

    open(my $fh, $file);
    while (my $line = <$fh>) {
        $lineNum++;
        chomp $line;
        next unless ($line && $line =~ /\S+/);
        next if $line =~ /^#/;

        # First line of file (after any comment lines) should be the header
        unless (@header) {
            @header = map { trim($_) } split($DELIMITER, $line);
            if (!@header) {
                print "log: error: empty header line\n";
                exit(-1);
            }
            next;
        }

        # Parse data line
        my @tok = map { trim($_) } split($DELIMITER, $line);
        if (@tok < @header) {
            print "log: error: missing fields (", scalar(@tok), "<", scalar(@header), ") on line $lineNum\n";
            exit(-1);
        }
        my $i = 0;
        my %fields = map { $_ => $tok[$i++] } @header;

        # Make sure required fields are present
        my $filename = $fields{Filename};
        if (!$filename or !$fields{Name}) {
            print "log: error: missing required column at line $lineNum: $line\n";
            exit(-1);
        }
        if ($data{$filename}) {
            print "log: error: duplicate filename '$filename'\n";
            exit(-1);
        }

        $data{$filename} = \%fields;
    }
    close($fh);
    return \%data;
}

sub process_dir {
    my $dir = shift;
    my $metadata = shift;
    my $notebook = shift;

    print "process_dir: $dir\n";

    # Load all data files in directory
    my @loaded;
    my $total_count = 0;
    my $load_count = 0;
    opendir( my $fh, $dir ) or die;
    my @contents = sort readdir($fh);
    foreach my $item ( @contents ) {
        print "process_dir: $dir/$item\n";
        next if ($item =~ /^\./); # skip hidden files
        next if ($item =~ /\.txt$/); # skip metadata file
        
        # Check for subdirectories in tarball
        if (-d "$dir/$item") {
            print "log: warning: '$item' is a directory, skipping ...\n";
            next;
        }
        
        # Get metadata entry for file
        my $md = ();
        if ($metadata->{$item}) {
            $md = $metadata->{$item};
        }
        else {
            print "log: warning: no metadata for file $item, skipping file ...\n";
            next;
        }
        
        # Call appropriate loader based on file type
        my $path = catdir($staging_dir, ++$total_count);
        my ($item_id, $item_type);
        if ( $item =~ /\.csv|\.tsv|.bed|\.gff|\.vcf|\.bam|\.fastq/i ) { # experiment
            $item_id = process_experiment( file => "$dir/$item", metadata => $md, path => $path );
            $item_type = 3; # FIXME hardcoded type
        }
        elsif ( $item =~ /\.fa|\.faa|\.fasta/i ) { # genome
            $item_id = process_genome( file => "$dir/$item", metadata => $md, path => $path );
            $item_type = 2; # FIXME hardcoded type
        }
        else {
            print "log: warning: unknown file type for $item, skipping file ...\n";
            next;
        }
        
        next unless ($item_id); # skip if load error
        $load_count++;
        push @loaded, { item_type => int($item_type), item_id => int($item_id) };
    }
    closedir($fh);

    unless ($load_count) {
        print "log: error: no data files can be loaded, check that file extensions and metadata fields are correct ",
              "<a href='https://genomevolution.org/wiki/index.php/LoadBatch'>here</a>.", "\n";
        exit(-1);
    }
    
    unless (@loaded) {
        print "log: error: none of the data files were loaded successfully\n";
        exit(-1);
    }
    
    print "log: ----------------------------------------------------\n";
    
    # Print summary of files loaded
    my %itemsTypeCount; map { $itemsTypeCount{ CoGeX->node_type_name($_->{item_type}) }++ } @loaded;
    print "log: Number of items loaded: $load_count (", ($total_count - $load_count) ," failed or skipped)\n";
    foreach my $type (keys %itemsTypeCount) {
        print 'log: ', $type, ': ', $itemsTypeCount{$type};
    }
    
    # Use specified notebook
    if ($notebook) {
        my @item_list = map { [ $_->{item_id}, $_->{item_type} ] } @loaded;
        my $error = add_items_to_notebook(db => $coge, user => $user, notebook => $notebook, item_list => \@item_list); # FIXME change item_list to hash param
        if ($error) {
            print "log: error: $error\n";
            exit(-1);
        }
        print "log: Added items to existing notebook '", $notebook->name, "' id", $notebook->id, "\n";
    }
    # Otherwise, create notebook of experiments
    else {
        $notebook = create_notebook(name => $notebook_name, desc => $notebook_desc, item_list => \@loaded); #FIXME use CoGe::Core::Notebook::create_notebook instead
        unless ($notebook) {
            print "log: error: failed to create notebook '$notebook_name'\n";
            exit(-1);
        }
        print "log: Created notebook '", $notebook->name, "' id", $notebook->id, "\n";
    }
    
    return ($notebook, \@loaded);
}


sub process_genome { #TODO merge with process_experiment?
    my %opts = @_;
    my $file = $opts{file};
    my $md   = $opts{metadata};
    my $path = $opts{path};
    die unless ($file && $md && $path);

    # Check params and set defaults - FIXME better way to handle capitalized keys as one-liner
    my $name;
    $name = $md->{Name} if ($md and $md->{Name});
    $name = $md->{name} if ($md and $md->{name});
    my $description = '';
    $description = $md->{Description} if ($md and $md->{Description});
    $description = $md->{description} if ($md and $md->{description});
    my $source = $user->display_name;
    $source      = $md->{Source} if ($md and $md->{Source});
    $source      = $md->{source} if ($md and $md->{source});
    my $version = 1;
    $version     = $md->{Version} if ($md and $md->{Version});
    $version     = $md->{version} if ($md and $md->{version});
    my $restricted = 1;
    $restricted = ($md->{Restricted} eq 'yes') if ($md and $md->{Restricted});
    $restricted = ($md->{restricted} eq 'yes') if ($md and $md->{restricted});
    my $organism_id;
    $organism_id = $md->{Organism} if ($md and $md->{Organism});
    $organism_id = $md->{organism} if ($md and $md->{organism});
    unless ($name && $organism_id) {
        print "log: error: missing required metadata field for file ", basename($file), ", skipping file ...\n";
        return;
    }

    # Run load script
    print "log: ----------------------------------------------------\n";
    print "log: Loading genome '$name'\n";
    $file = escape($file);
    my $install_dir = $P->{SEQDIR};
    my $cmd = catfile($P->{SCRIPTDIR}, 'load_genome.pl') . ' ' .
        "-config $config -user_name '".$user->user_name."' -restricted 1 -name '$name' -desc '$description' " .
        "-version '$version' -source_name '$source' -organism_id $organism_id " .
        "-staging_dir $path -install_dir $install_dir -fasta_file '$file' ";
    print "Running: ", $cmd, "\n";
    my $output = qx{ $cmd }; # TODO: use run() here instead?
    print $output;
    if ( $? != 0 ) {
        print "load_genome.pl failed with rc=$?\n",
              "log: error: Genome '$name' was not loaded due to an error\n";
        return;
    }

    # Extract genome id from output
    my ($id) = $output =~ /genome id: (\d+)/;
    unless ($id) {
        print "log: error: unable to retrieve genome id from result\n";
        exit(-1);
    }
    print "log: Genome '$name' id$id loaded successfully\n";

    # Add user-defined metadata fields to genome
    my $g = $coge->resultset('Genome')->find($id);
    die unless $g;

    foreach my $column_name (keys %$md) {
        # Skip built-in fields
        next if (grep { /$column_name/i } ('filename', 'name', 'description', 'source', 'version', 'restricted', 'genome', 'organism'));

        # Add annotation
        my $type = $coge->resultset('AnnotationType')->find_or_create( { name => $column_name } );
        die unless ($type);
        $g->add_to_experiment_annotations(
            {
                annotation_type_id => $type->id,
                annotation         => $md->{$column_name}
            }
        );
    }

    # TODO add support for "_link" option

    return $id;
}

sub process_experiment {
    my %opts = @_;
    my $file = $opts{file};
    my $md   = $opts{metadata};
    my $path = $opts{path};
    die unless ($file && $md && $path);

    # Check params and set defaults
    my $name;
    $name = $md->{Name} if ($md and $md->{Name});
    $name = $md->{name} if ($md and $md->{name});
    my $description = '';
    $description = $md->{Description} if ($md and $md->{Description});
    $description = $md->{description} if ($md and $md->{description});
    my $source = $user->display_name;
    $source      = $md->{Source} if ($md and $md->{Source});
    $source      = $md->{source} if ($md and $md->{source});
    my $version = 1;
    $version     = $md->{Version} if ($md and $md->{Version});
    $version     = $md->{version} if ($md and $md->{version});
    my $restricted = 1;
    $restricted = ($md->{Restricted} eq 'yes') if ($md and $md->{Restricted});
    $restricted = ($md->{restricted} eq 'yes') if ($md and $md->{restricted});
    my $genome_id = $gid;
    $genome_id = $md->{Genome} if ($md and $md->{Genome});
    $genome_id = $md->{genome} if ($md and $md->{genome});
    unless ($name) { #($name && $genome_id) {
        print "log: error: missing required metadata field for file ", basename($file), ", skipping file ...\n";
        return;
    }

    # Run load script
    print "log: ----------------------------------------------------\n";
    print "log: Loading experiment '$name'\n";
    $file = escape($file);
    my $install_dir = $P->{EXPDIR};
    my $cmd = catfile($P->{SCRIPTDIR}, 'load_experiment.pl') . ' ' .
        "-config $config -user_name '".$user->user_name."' -restricted 1 -name '$name' -desc '$description' " .
        "-version '$version' -wid $wid -gid $gid -source_name '$source' " .
        "-staging_dir $path -install_dir $install_dir -data_file '$file' ";
    print "Running: ", $cmd, "\n";
    my $output = qx{ $cmd }; # TODO: use run() here?
    print $output;
    if ( $? != 0 ) {
        print "load_experiment.pl failed with rc=$?\n",
              "log: error: Experiment '$name' was not loaded due to an error\n";
        return;
    }

    # Extract experiment id from output
    my ($id) = $output =~ /experiment id: (\d+)/;
    unless ($id) {
        print "log: error: unable to retrieve experiment id\n";
        exit(-1);
    }
    print "log: Experiment '$name' id$id loaded successfully\n";

    # Add Add user-defined metadata fields to experiment
    my $exp = $coge->resultset('Experiment')->find($id);
    die unless $exp;

    foreach my $column_name (keys %$md) {
        # Skip built-in fields
        next if (grep { /$column_name/i } ('filename', 'name', 'description', 'source', 'version', 'restricted', 'genome'));

        # Add annotation
        my $type = $coge->resultset('AnnotationType')->find_or_create( { name => $column_name } );
        die unless ($type);
        $exp->add_to_experiment_annotations(
            {
                annotation_type_id => $type->id,
                annotation         => $md->{$column_name}
            }
        );
    }

    # TODO add support for "_link" option

    return $id;
}

sub create_notebook { #FIXME use routine CoGe::Core::Notebook
    my %opts    = @_;
    my $name    = $opts{name};
    my $desc    = $opts{desc};
    my $item_list = $opts{item_list};    # optional
    my $type_id = 2; # hardcoded to "mixed" type

    $name = '' unless defined $name;
    $desc = '' unless defined $desc;

    # Create the new list
    my $list = $coge->resultset('List')->create(
        {
            name         => $name,
            description  => $desc,
            #list_type_id => 5, # FIXME hardcoded to 5, # mdb removed 12/14/16 COGE-800
            creator_id   => $user->id,
            restricted   => 1
        }
    );
    return unless $list;

    # Set user as owner
    my $conn = $coge->resultset('UserConnector')->create(
        {
            parent_id   => $user->id,
            parent_type => 5,           #FIXME hardcoded to "user"
            child_id    => $list->id,
            child_type  => 1,           #FIXME hardcoded to "list"
            role_id     => 2            #FIXME hardcoded to "owner"
        }
    );
    return unless $conn;

    # Add selected items to new notebook
    foreach (@$item_list) {
        my $conn = $coge->resultset('ListConnector')->find_or_create(
            {
                parent_id  => $list->id,
                child_id   => $_->{item_id},
                child_type => $_->{item_type}
            }
        );
        return unless $conn;
    }

    # Record in log
    CoGe::Accessory::Web::log_history(
        db          => $coge,
        user_id     => $user->id,
        page        => "LoadBatch",
        description => 'created notebook ' . list->info,
        parent_id   => $list->id,
        parent_type => 1 #FIXME magic number
    );

    return $list;
}

sub run {
    my $cmd = shift;
    print "run:  $cmd\n";
    my $rc = execute($cmd);
    if ( $rc != 0 ) {
        print "log: error: command failed with rc=$rc: $cmd\n";
        exit(-1);
    }
}
