#!/usr/bin/perl -w
use strict;
use CoGeX;
use CoGe::Accessory::Web qw(get_defaults);
use File::Path qw(make_path);
use File::Basename;
use URI::Escape::JavaScript qw(escape);
use Data::Dumper;

#-------------------------------------------------------------------------------
# Parameters

my $DEBUG = 0;
my $COGE_DIR      = '/opt/apache/coge';
my $INSTALL_DIR   = '/storage/coge/data/experiments';

my $gid           = 22716;
my $data_dir      = '/home/mbomhoff/tmp/data/';
my $staging_dir   = '/home/mbomhoff/tmp/staging';
my $metadata_file; # none
my $citation_link; # none
my $user          = 'mbomhoff';
my $source        = 'Matt Bomhoff';
my $description   = '';
my $version       = '1';

#-------------------------------------------------------------------------------

# Load config file
my $config = "$COGE_DIR/web/coge.conf";
my $P      = CoGe::Accessory::Web::get_defaults($config);

# Connect to database
my $connstr = "dbi:mysql:dbname=".$P->{DBNAME}.";host=".$P->{DBHOST}.";port=".$P->{DBPORT}.";";
my $coge = CoGeX->connect( $connstr, $P->{DBUSER}, $P->{DBPASS} );

# Load metadata file
my $metadata;
if ($metadata_file) {
    get_metadata($metadata_file);
    unless ($metadata) {
        print STDERR "load_all_experiments: no metadata\n";
        exit(-1);
    }
    #print STDERR Dumper $metadata;
}

# Load each experiment file
my $exp_count = 0;
process_dir($data_dir);

# Yay!
print STDERR "Loaded $exp_count experiments\n";
print STDERR "All done!\n";
exit;

#-------------------------------------------------------------------------------
sub get_metadata {
    my $file = shift;
    open( IN, $file );
    my %data;

    #Name   Mark    Genome  Protocol        Platform ID     Platform name   Citation        Processing      Raw files       Processed files Series  Samples
    while (<IN>) {
        chomp;
        next unless $_;
        next if /^#/;
        my @line = split(/\t/);
        my $processed_file = $line[9];
        $data{ $processed_file } = \@line;
    }
    close IN;
    return \%data;
}

sub process_dir {
    my $dir = shift;
    my $notebook_name = shift;
    
    my @experiments;
    
    print STDERR "process_dir: $dir\n";
    opendir( my $fh, $dir ) or die;
    my @contents = sort readdir($fh);
    foreach my $item ( @contents ) {
        next if ($item =~ /^\./);
	print STDERR "item: $dir/$item\n";
        if (-d "$dir/$item") { # directory
            my $exp_list = process_dir("$dir/$item", $item); # pass in directory name as notebook name
            
            if (!$DEBUG && @$exp_list > 0 && $notebook_name) {
                # Create notebook of experiments
                if (!create_notebook(name => $notebook_name, item_list => $exp_list)) {
                    print STDERR "Failed to create notebook '$notebook_name'\n";
                    exit(-1);
                }
                print STDERR "Created notebook '$notebook_name'\n";
            }
        }
        elsif ( $item =~ /\.csv|\.bam/ && -r "$dir/$item" ) { # file
            my $md = ();
            if ($metadata) {
                if ($metadata->{$item}) {
                    $md = $metadata->{$item};
                }
                else {
                    print STDERR "WARNING: no metadata for $item, skipping ...\n";
                    next;
                }
            }
            my $eid = process_file(
                file     => "$dir/$item",
                metadata => $md
            );
            push @experiments, $eid if $eid;
        }
    }
    closedir($fh);
    
    return \@experiments;
}

sub process_file {
    my %opts = @_;
    my $file = $opts{file};
    my $md   = $opts{metadata};
    
    my (
        $name,        $mark,            $genome,   $protocol,
        $platform_id, $platform_name,   $citation, $processing,
        $raw_files,   $processed_files, $series,   $samples
    );
    
    if ($md) {
        (
            $name,        $mark,            $genome,   $protocol,
            $platform_id, $platform_name,   $citation, $processing,
            $raw_files,   $processed_files, $series,   $samples
        ) = @$md;
    }
    else {
        ($name) = fileparse($file, qr/\.[^.]*/);
    }

    $exp_count++;
    
    # Create/clear staging directory
    my $staging = "$staging_dir/$exp_count";
    if (-e $staging) {
        print "Clearing staging directory: $staging\n";
        `rm $staging/*`;
    }
    else {
        make_path($staging);
    }

    # Run load script
    $file = escape($file);
    my $cmd = "$COGE_DIR/scripts/load_experiment.pl " .
        "-config $config -user_name $user -restricted 1 -name '$name' -desc '$description' " .
        "-version '$version' -gid $gid -source_name '$source' " .
	    "-allow_negative 1 " .
        "-staging_dir $staging -install_dir $INSTALL_DIR -data_file '$file'";
    print "Running: " . $cmd, "\n";
    return if ($DEBUG);
    my $output = qx{ $cmd };
    if ( $? != 0 ) {
        print STDERR "load_experiment.pl failed with rc=$?\n";
        exit(-1);
    }

    # Extract experiment id from output
    my ($eid) = $output =~ /experiment id: (\d+)/;
    if (!$eid) {
        print STDERR "Unable to retrieve experiment id\n";
        exit(-1);
    }
    print STDERR "Captured experimemt id $eid\n";

    # Add experiment annotations from metadata file
    if ($md) {
        my $exp = $coge->resultset('Experiment')->find($eid);
        my $exp_type = $coge->resultset('ExperimentType')->find_or_create( { name => $protocol } );
        $coge->resultset('ExperimentTypeConnector')->create(
            {
                experiment_type_id => $exp_type->id,
                experiment_id      => $exp->id
            }
        );
    
        $exp_type = $coge->resultset('ExperimentType')->find_or_create( { name => $platform_name } );
        $coge->resultset('ExperimentTypeConnector')->create(
            {
                experiment_type_id => $exp_type->id,
                experiment_id      => $exp->id
            }
        );
    
        my $etype =
          $coge->resultset('AnnotationType')
          ->find_or_create( { name => "Citation" } );
        $exp->add_to_experiment_annotations(
            {
                annotation_type_id => $etype->id,
                annotation         => $citation,
                link               => $citation_link
            }
        );
        $etype =
          $coge->resultset('AnnotationType')
          ->find_or_create( { name => "Mark" } );
        $exp->add_to_experiment_annotations(
            {
                annotation_type_id => $etype->id,
                annotation         => $mark,
            }
        );
        $etype =
          $coge->resultset('AnnotationType')
          ->find_or_create( { name => "Platform" } );
        $exp->add_to_experiment_annotations(
            {
                annotation_type_id => $etype->id,
                annotation         => $platform_name,
            }
        );
        $etype =
          $coge->resultset('AnnotationType')
          ->find_or_create( { name => "Platform ID" } );
        $exp->add_to_experiment_annotations(
            {
                annotation_type_id => $etype->id,
                annotation         => $platform_id,
            }
        );
        $etype =
          $coge->resultset('AnnotationType')
          ->find_or_create( { name => "Processed files" } );
        $exp->add_to_experiment_annotations(
            {
                annotation_type_id => $etype->id,
                annotation         => $processed_files,
            }
        );
        $etype =
          $coge->resultset('AnnotationType')
          ->find_or_create( { name => "Processing" } );
        $exp->add_to_experiment_annotations(
            {
                annotation_type_id => $etype->id,
                annotation         => $processing,
            }
        );
        $etype =
          $coge->resultset('AnnotationType')
          ->find_or_create( { name => "Protocol" } );
        $exp->add_to_experiment_annotations(
            {
                annotation_type_id => $etype->id,
                annotation         => $protocol,
            }
        );
        $etype =
          $coge->resultset('AnnotationType')
          ->find_or_create( { name => "Series" } );
        $exp->add_to_experiment_annotations(
            {
                annotation_type_id => $etype->id,
                annotation         => $series,
            }
        );
        $etype =
          $coge->resultset('AnnotationType')
          ->find_or_create( { name => "Samples" } );
        $exp->add_to_experiment_annotations(
            {
                annotation_type_id => $etype->id,
                annotation         => $samples,
            }
        );
        $etype =
          $coge->resultset('AnnotationType')
          ->find_or_create( { name => "Gene mutant" } );
        $exp->add_to_experiment_annotations(
            {
                annotation_type_id => $etype->id,
                annotation         => $genome,
            }
        );
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

    # Get user
    my ($user) = $coge->resultset('User')->search({user_name => $user});
    return unless $user;

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

    return 1;
}
