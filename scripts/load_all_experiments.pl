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
my $DELIMITER     = '\s*\"?\t\"?\s*';

my $gid           = 16904; # ID for genome to assign experiments to
my $data_dir      = '/home/mbomhoff/tmp/test';  # location of data files
my $staging_dir   = '/home/mbomhoff/tmp/staging';           # temporary staging path
my $metadata_file = '/home/mbomhoff/tmp/jmccurdy/51genos/dmr_genetic_influences_ip_array.txt'; # optional metadata file
my $user          = 'mbomhoff';     # owner
my $source        = 'Matt Bomhoff'; # default source if not specified in metadata file
my $version       = '1';            # default version if not specified in metadata file

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
    $metadata = get_metadata($metadata_file);
    unless ($metadata) {
        print STDERR "load_all_experiments: error: no metadata loaded\n";
        exit(-1);
    }
    #print STDERR Dumper $metadata;
}

# Load each experiment file
my $exp_count = 0;
process_dir($data_dir, $metadata);

# Yay!
print STDERR "Loaded $exp_count experiments\n";
print STDERR "All done!\n";
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
                print STDERR "load_all_experiments: error: empty header line\n";
                exit(-1);
            }
            #print STDERR join(",\n", @header),"\n";
            next;
        }
        
        # Parse data line
        my @tok = map { trim($_) } split($DELIMITER, $line);
        die if (@tok != @header);
        my $i = 0;
        my %fields = map { $_ => $tok[$i++] } @header;
        #print STDERR Dumper \%fields, "\n";
        
        # Make sure required fields are present
        my $filename = $fields{Filename};
        if (!$filename or !$fields{Name}) {
            print STDERR "load_all_experiments: error: missing required column:\nline $lineNum: $line\n";
            exit(-1);
        }
        if ($data{$filename}) {
            print STDERR "load_all_experiments: error: duplicate filename '$filename'\n";
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
    my $notebook_name = shift;
    my @experiments;
    
    print STDERR "process_dir: $dir ", ($notebook_name ? $notebook_name : ''), "\n";
    opendir( my $fh, $dir ) or die;
    my @contents = sort readdir($fh);
    foreach my $item ( @contents ) {
        next if ($item =~ /^\./);
        print STDERR "item: $dir/$item\n";
        if (-d "$dir/$item") { # directory
            my $exp_list = process_dir("$dir/$item", $metadata, $item); # pass in directory name as notebook name
            if (!$DEBUG && @$exp_list > 0) {
                # Create notebook of experiments
                if (!create_notebook(name => $item, item_list => $exp_list)) {
                    print STDERR "Failed to create notebook '$item'\n";
                    exit(-1);
                }
                print STDERR "Created notebook '$item'\n";
            }
        }
        elsif ( $item =~ /\.csv|\.bam|\.bed/ && -r "$dir/$item" ) { # file
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
    
    # Check params and set defaults
    my $name;
    $name = $md->{Name} if ($md and $md->{Name});
    $name = $md->{name} if ($md and $md->{name});
    my $description = '';
    $description = $md->{Description} if ($md and $md->{Description});
    $description = $md->{description} if ($md and $md->{description});
    $source      = $md->{Source} if ($md and $md->{Source});
    $source      = $md->{source} if ($md and $md->{source});
    $version     = $md->{Version} if ($md and $md->{Version});
    $version     = $md->{version} if ($md and $md->{version});
    my $restricted = 1;
    $restricted = ($md->{Restricted} eq 'yes') if ($md and $md->{Restricted});
    $restricted = ($md->{restricted} eq 'yes') if ($md and $md->{restricted});
    die unless ($name and $gid and $source and $user);

    # Run load script
    $file = escape($file);
    my $cmd = "$COGE_DIR/scripts/load_experiment.pl " .
        "-config $config -user_name $user -restricted 1 -name '$name' -desc '$description' " .
        "-version '$version' -gid $gid -source_name '$source' " .
	    #"-allow_negative 0 " .
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
