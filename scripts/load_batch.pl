#!/usr/bin/perl -w
use strict;
use CoGeX;
use CoGe::Accessory::Web qw(get_defaults);
use File::Path;
use File::Basename;
use File::Spec::Functions qw(catdir catfile);
use File::Touch;
use Getopt::Long;
use URI::Escape::JavaScript qw(escape unescape);
use Data::Dumper;

use vars qw(
    $staging_dir $data_file $notebook_name $notebook_desc $gid $user_name 
    $config $log_file $user $genome $result_dir
);
  
my $DEBUG = 0;
my $DELIMITER = '\t';

#-------------------------------------------------------------------------------

GetOptions(
    "staging_dir=s" => \$staging_dir,    # temporary staging path
    "result_dir=s"  => \$result_dir,     # results path
    "data_file=s"   => \$data_file,      # input data file (JS escape)
    "name=s"        => \$notebook_name,  # notebook name (JS escaped)
    "desc=s"        => \$notebook_desc,  # notebook description (JS escaped)
    "gid=s"         => \$gid,            # genome id
    "user_name=s"   => \$user_name,      # user name
    "config=s"      => \$config,         # CoGe config file
    "log_file=s"    => \$log_file        # optional log file
);

# Open log file
$| = 1;
$log_file = "$staging_dir/log.txt" unless $log_file;
mkpath($staging_dir) unless -r $staging_dir;
open( my $log, ">>$log_file" ) or die "Error opening log file $log_file";
$log->autoflush(1);
print $log "Starting $0 (pid $$)\n";

# Process and verify parameters
$data_file     = unescape($data_file);
$notebook_name = unescape($notebook_name);
$notebook_desc = unescape($notebook_desc);

if ($user_name eq 'public') {
    print $log "log: error: not logged in\n";
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
    print $log "log: error: couldn't connect to database\n";
    exit(-1);
}

# Retrieve genome
$genome = $coge->resultset('Genome')->find( { genome_id => $gid } );
unless ($genome) {
    print $log "log: error finding genome id$gid\n";
    exit(-1);
}

# Retrieve user
$user = $coge->resultset('User')->find( { user_name => $user_name } );
unless ($user) {
    print $log "log: error finding user '$user_name'\n";
    exit(-1);
}

# Untar data file
my $data_dir = catdir($staging_dir, 'data');
unless (-r $data_file) {
    print $log "log: error: cannot access data file\n";
    exit(-1);
}
unless ($data_file =~ /\.tar\.gz$/) {
    print $log "log: error: data file needs to be a tarball (end with .tar.gz)\n";
    exit(-1);
}
print $log "log: Decompressing/extracting data\n";
mkpath($data_dir);
execute( $P->{TAR}.' -xf '.$data_file.' --directory '.$data_dir );

# Find metadata file
my ($metadata_file) = glob("$data_dir/*.txt");
unless ($metadata_file) { # require a metadata file for now but could be optional in the future
    print $log "log: error: cannot find .txt metadata file\n";
    exit(-1);
}

# Load metadata file
my $metadata;
if ($metadata_file) {
    print $log "log: Loading metadata file '".fileparse($metadata_file)."'\n";
    $metadata = get_metadata($metadata_file);
    unless ($metadata) {
        print $log "log: error: no metadata loaded\n";
        exit(-1);
    }
    print $log Dumper($metadata), "\n";
}

# Load each experiment file
my $exp_count = 0;
my ($notebook, $exp_ids) = process_dir($data_dir, $metadata);

# Save result document
if ($result_dir) {
    mkpath($result_dir);
    CoGe::Accessory::TDS::write(
        catfile($result_dir, '1'),
        {
            genome_id => int($genome->id),
            notebook_id => int($notebook->id),
            experiments => $exp_ids
        }
    );
}

# Yay!
CoGe::Accessory::Web::log_history(
    db          => $coge,
    user_id     => $user->id,
    page        => "LoadBatch",
    description => 'load batch experiments',
    link        => 'NotebookView.pl?nid=' . $notebook->id
);

print $log "log: Loaded $exp_count experiments\n";
close($log);

# Create "log.done" file to indicate completion to JEX
my $logdonefile = "$staging_dir/log.done";
touch($logdonefile);

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
                print $log "log: error: empty header line\n";
                exit(-1);
            }
#            print STDERR Dumper \@header,"\n";
            next;
        }
        
        # Parse data line
        my @tok = map { trim($_) } split($DELIMITER, $line);
        if (@tok < @header) {
            print $log "log: error: missing fields (", scalar(@tok), "<", scalar(@header), ") on line $lineNum\n";
#            print STDERR Dumper \@tok,"\n";
            exit(-1);
        }
        my $i = 0;
        my %fields = map { $_ => $tok[$i++] } @header;
        #print STDERR Dumper \%fields, "\n";
        
        # Make sure required fields are present
        my $filename = $fields{Filename};
        if (!$filename or !$fields{Name}) {
            print $log "log: error: missing required column:\nline $lineNum: $line\n";
            exit(-1);
        }
        if ($data{$filename}) {
            print $log "log: error: duplicate filename '$filename'\n";
            exit(-1);
        }
        
        $data{$filename} = \%fields;
    }
    close($fh);
    return \%data;
}

sub trim {
    my $s = shift;
    $s =~ s/^\s+//;
    $s =~ s/\s+$//;
    $s =~ s/^\"//;
    $s =~ s/\"$//;
    return $s;
}

sub process_dir {
    my $dir = shift;
    my $metadata = shift;
    
    print $log "process_dir: $dir\n";
    
    # Load all experiment files in directory
    my @experiments;
    opendir( my $fh, $dir ) or die;
    my @contents = sort readdir($fh);
    foreach my $item ( @contents ) {
        next if ($item =~ /^\./);
        print $log "item: $dir/$item\n";
        if ( $item =~ /\.csv|\.bam|\.bed/ && -r "$dir/$item" ) { # file
            my $md = ();
            if ($metadata) {
                if ($metadata->{$item}) {
                    $md = $metadata->{$item};
                }
                else {
                    print $log "WARNING: no metadata for $item, skipping ...\n";
                    next;
                }
            }
            my $eid = process_file(
                file     => "$dir/$item",
                metadata => $md
            );
            push @experiments, int($eid) if $eid;
        }
    }
    closedir($fh);
    
    unless (@experiments) {
        print $log "log: error: no experiment files found\n";
        exit(-1);
    }

    # Create notebook of experiments
    my $notebook = create_notebook(name => $notebook_name, desc => $notebook_desc, item_list => \@experiments);
    unless ($notebook) {
        print $log "log: error: failed to create notebook '$notebook_name'\n";
        exit(-1);
    }
    print $log "notebook id: ".$notebook->id."\n"; # !!!! don't change, gets parsed by calling code
    print $log "log: Created notebook '$notebook_name'\n";
    return ($notebook, \@experiments);
}

sub process_file {
    my %opts = @_;
    my $file = $opts{file};
    my $md   = $opts{metadata};
    
    $exp_count++;
    
    my $exp_staging_dir = catdir($staging_dir, $exp_count);
    
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
    die unless ($name and $gid and $source and $user);

    # Run load script
    print $log "log: Loading experiment '$name'\n";
    $file = escape($file);
    my $cmd = catfile($P->{SCRIPTDIR}, 'load_experiment.pl') . ' ' .
        "-config $config -user_name '".$user->user_name."' -restricted 1 -name '$name' -desc '$description' " .
        "-version '$version' -gid $gid -source_name '$source' " .
        "-staging_dir $exp_staging_dir -install_dir ".$P->{EXPDIR}." -data_file '$file' " .
        "-log_file $log_file"; # reuse log so that all experiment loads are present
    print $log "Running: " . $cmd, "\n";
    return if ($DEBUG);
    my $output = qx{ $cmd };
    if ( $? != 0 ) {
        print $log "log: error: load_experiment.pl failed with rc=$?\n";
        exit(-1);
    }
    open( $log, ">>$log_file" ) or die "Error opening log file $log_file"; # Reopen log file
    print $log "log: Experiment '$name' loaded successfully\n";

    # Extract experiment id from output
    my ($eid) = $output =~ /experiment id: (\d+)/;
    if (!$eid) {
        print $log "log: error: unable to retrieve experiment id\n";
        exit(-1);
    }
    print $log "Captured experiment id $eid\n";

    # Add experiment annotations from metadata file
    if ($md) {
        my $exp = $coge->resultset('Experiment')->find($eid);
        die unless $exp;

        foreach my $column_name (keys %$md) {
            # Skip built-in fields
            next if (grep { /$column_name/i } ('filename', 'name', 'description', 'source', 'version', 'restricted'));
            
            # Add annotation
            my $type = $coge->resultset('AnnotationType')->find_or_create( { name => $column_name } );
            $exp->add_to_experiment_annotations(
                {
                    annotation_type_id => $type->id,
                    annotation         => $md->{$column_name}
                }
            );
        }

        # TODO add support for "_link" option
    }

    return $eid;
}

sub create_notebook {
    my %opts    = @_;
    my $name    = $opts{name};
    my $desc    = $opts{desc};
    my $item_list = $opts{item_list};    # optional
    my $type_id = 2; # hardcoded to "experiment" type

    $name = '' unless defined $name;
    $desc = '' unless defined $desc;

    # Create the new list
    my $list = $coge->resultset('List')->create(
        {
            name         => $name,
            description  => $desc,
            list_type_id => $type_id,
            restricted => 1
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
                child_id   => $_,
                child_type => 3 # hardcoded to "experiment" type
            }
        );
        return unless $conn;
    }    

    # Record in log
    CoGe::Accessory::Web::log_history(
        db          => $coge,
        user_id     => $user->id,
        page        => "User",
        description => 'create notebook id' . $list->id
    );

    return $list;
}

sub execute { # FIXME move into Util.pm
    my $cmd = shift;
    print $log "$cmd\n";
    my @cmdOut    = qx{$cmd};
    my $cmdStatus = $?;
    if ( $cmdStatus != 0 ) {
        print $log "log: error: command failed with rc=$cmdStatus: $cmd\n";
        exit(-1);
    }
}