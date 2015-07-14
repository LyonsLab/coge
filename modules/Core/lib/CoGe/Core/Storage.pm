package CoGe::Core::Storage;

=head1 NAME

CoGe::Core::Storage

=head1 SYNOPSIS

Abstraction layer on top of storage sub-system created for issues
77 and 157.  All accesses to genome FASTA sequences, experiment data files, etc. 
should go through this module.

=head1 DESCRIPTION

=head1 AUTHOR

Matt Bomhoff

=head1 COPYRIGHT

The full text of the license can be found in the
LICENSE file included with this module.

=head1 SEE ALSO

=cut

use strict;
use warnings;

use CoGe::Accessory::Web qw(get_defaults);
use CoGe::Accessory::TDS qw(read);
use CoGe::Accessory::Jex;
use CoGe::Accessory::Workflow;
use CoGe::Accessory::IRODS qw(irods_iget irods_ils);
use File::Basename;
use POSIX qw(floor);
use File::Spec::Functions;
use File::Path qw(mkpath);
use File::Find qw(find);
use File::Slurp;
use List::Util qw[min max];
use File::Slurp;
use JSON::XS qw(encode_json decode_json);
use Data::Dumper;
use POSIX qw(ceil);

BEGIN {
    use vars qw ($VERSION @ISA @EXPORT @EXPORT_OK $DATA_TYPE_QUANT 
                 $DATA_TYPE_POLY $DATA_TYPE_ALIGN $DATA_TYPE_MARKER);
    require Exporter;

    $VERSION = 0.1;
    @ISA     = qw(Exporter);
    @EXPORT = qw(
      get_tiered_path get_workflow_paths get_upload_path get_log
      get_genome_file index_genome_file get_genome_seq get_genome_path
      get_genome_cache_path get_workflow_results add_workflow_result
      get_workflow_results_file get_workflow_log_file get_download_path
      get_experiment_path get_experiment_files get_experiment_data
      create_experiments_from_batch
      create_genome_from_file create_genome_from_NCBI
      create_annotation_dataset reverse_complement
      get_irods_file get_irods_path
      $DATA_TYPE_QUANT $DATA_TYPE_POLY $DATA_TYPE_ALIGN $DATA_TYPE_MARKER
    );
    @EXPORT_OK = qw(data_type);

    # Experiment Data Types -- move to CoGe::Core::Experiment?
    $DATA_TYPE_QUANT  = 1; # Quantitative data
    $DATA_TYPE_POLY   = 2; # Polymorphism data
    $DATA_TYPE_ALIGN  = 3; # Alignments
    $DATA_TYPE_MARKER = 4; # Markers
}

################################################ subroutine header begin ##
=head2 get_tiered_path
 Usage     :
 Purpose   : This method determines the correct directory structure for storing
             the files for a genome/experiment.
 Returns   :
 Argument  :
 Throws    : none
 Comments  : The idea is to build a dir structure that holds large amounts of
             files, and is easy to lookup based on genome ID number.
             The strucuture is three levels of directorys, and each dir holds
             1000 files and/or directorys.
             Thus:
             ./0/0/0/ will hold files 0-999
             ./0/0/1/ will hold files 1000-1999
             ./0/0/2/ will hold files 2000-2999
             ./0/1/0 will hold files 1000000-1000999
             ./level0/level1/level2/
See Also   :
=cut
################################################## subroutine header end ##

sub get_tiered_path {
    my $id = shift;
    return unless $id;

    my $level0 = floor( $id / 1000000000 ) % 1000;
    my $level1 = floor( $id / 1000000 ) % 1000;
    my $level2 = floor( $id / 1000 ) % 1000;
    my $path   = catdir( $level0, $level1, $level2, $id );

    return $path;
}

sub get_genome_cache_path {
    my $gid = shift;
    return unless $gid;
    
    my $cache_dir = CoGe::Accessory::Web::get_defaults()->{'CACHEDIR'};
    unless ($cache_dir) {
        print STDERR "Storage::get_genome_cache_path missing CACHEDIR\n";
        return;
    }
    
    return catdir($cache_dir, $gid, "fasta");
}

sub get_genome_path {
    my $gid = shift;
    return unless $gid;

    my $seqdir = CoGe::Accessory::Web::get_defaults()->{'SEQDIR'};
    unless ($seqdir) {
        print STDERR "Storage::get_genome_path: WARNING, conf file parameter SEQDIR is blank!\n";
    }

    my $path = $seqdir . '/' . get_tiered_path($gid) . '/';
    unless ( -r $path ) {
        print STDERR "Storage::get_genome_path: genome path '$path' doesn't exist!\n";
    }

    return $path;
}

sub get_genome_file {
    my $gid = shift;
    return unless $gid;
    my $base_path = shift;    # optional

    $base_path = get_genome_path($gid) unless ($base_path);

    # First check for legacy filename
    my $file_path = $base_path . "$gid.faa";
    return $file_path if ( -r $file_path );

    # Not there, so check for new filename
    $file_path = $base_path . 'genome.faa';
    if ( not -r $file_path ) {
        print STDERR "Storage::get_genome_file: genome file '$file_path' not found or not readable!\n";
    }

    return $file_path;
}

sub index_genome_file {
    my %opts      = @_;
    my $gid       = $opts{gid};
    my $file_path = $opts{file_path};
    return unless ( $file_path or $gid );
    my $compress = $opts{compress};    # optional flag

    unless ($file_path) {
        $file_path = get_genome_file($gid);
        return unless $file_path;
    }

    # Index fasta file
    my $samtools = CoGe::Accessory::Web::get_defaults()->{'SAMTOOLS'};
    unless ($samtools) {
        print STDERR "Storage::index_genome_file: WARNING, conf file parameter SAMTOOLS is blank!\n";
    }
    my $cmd = "$samtools faidx $file_path";
    qx{ $cmd };
    if ( $? != 0 ) {
        print STDERR "Storage::index_genome_file: command failed with rc=$?: $cmd\n";
        return $?;
    }

    # Optionally generate compressed version of fasta/index files
    if ($compress) {
        my $razip = CoGe::Accessory::Web::get_defaults()->{'RAZIP'};
        unless ($razip) {
            print STDERR "Storage::index_genome_file: WARNING, conf file parameter RAZIP is blank!\n";
        }

        # Compress fasta file into RAZF using razip
        $cmd = "$razip -c $file_path > $file_path.razf";
        qx{ $cmd };
        if ( $? != 0 ) {
            print STDERR "Storage::index_genome_file: command failed with rc=$?: $cmd\n";
            return $?;
        }

        # Index compressed fasta file
        $cmd = "$samtools faidx $file_path.razf";
        qx{ $cmd };
        if ( $? != 0 ) {
            print STDERR "Storage::index_genome_file: command failed with rc=$?: $cmd\n";
            return $?;
        }
    }

    return 0;
}

sub get_genome_seq {
    my %opts = @_;
    my $gid  = $opts{gid};    # required
    unless ($gid) {
        print STDERR "Storage::get_genome_seq: genome id not specified!\n";
        return;
    }
    my $chr   = $opts{chr};
    my $start = $opts{start};
    $start = 1 unless ( $start and $start > 0 );
    my $stop = $opts{stop};
    $stop = $opts{end} if ( not defined $stop );
    $stop = $start unless $stop;
    my $strand            = $opts{strand};
    my $format            = $opts{format};
    my $fasta             = ( defined $format and $format eq 'fasta' );
    #my $FASTA_LINE_LENGTH = 60;
    my $seq;
    #print STDERR "Storage::get_genome_seq gid=$gid chr=" . ($chr ? $chr : '') . " start=" . (defined $start ? $start : '') . " stop=" . (defined $stop ? $stop : '') . "\n";

    # Validate params
    my $len;
    if ( defined $start and defined $stop ) {
        ( $start, $stop ) = ( $stop, $start ) if ( $stop < $start );
        $len = abs( $stop - $start ) + 1;
    }

    # No chromosome specified, return whole genome fasta file
    unless (defined $chr and $chr ne '') {
        # mdb added 4/9/14 issue 359
        my $file_path = get_genome_file($gid);
        my $seq = read_file($file_path);
        return $seq;

# mdb removed 4/9/14 issue 359
#        my $fh;
#        if ( !open( $fh, $file_path ) ) {
#            print STDERR "Storage::get_genome_seq: can't open fasta file '$file_path'!\n";
#            return;
#        }
#        my $readLen = read( $fh, $seq, -s $fh );
#        close($fh);
#        return $seq;
    }

    # Determine file path - first try new indexed method, otherwise
    # revert to old method
    my $file_path = get_genome_file($gid);
    my $file_index_sz = -s "$file_path.fai";
    #print STDERR "file_path=$file_path file_index_sz=$file_index_sz\n";
    if ($file_index_sz) {    # new indexed method
                             # Kludge chr/contig name
        if   ( $chr =~ /\D+/ ) { $chr = 'lcl|' . $chr; }
        else                   { $chr = 'gi|' . $chr; }

        # Extract requested piece of sequence file
        my $region =
          $chr . ( defined $start && defined $stop ? ":$start-$stop" : '' );
        my $samtools = CoGe::Accessory::Web::get_defaults()->{'SAMTOOLS'};
        unless ($samtools) {
            print STDERR "Storage::get_genome_seq: WARNING, conf file parameter SAMTOOLS is blank!\n";
        }
        my $cmd = "$samtools faidx $file_path '$region'";

        #print STDERR "$cmd\n";
        $seq = qx{$cmd};
        unless ($fasta) {
            	# remove header line
            	$seq =~ s/^(?:.*\n)//;
#2/17/14:  Note by EL:  THere is a problem where the following type sof regex sbustitutions fail if the string is longer then about 1G (http://www.perlmonks.org/?node_id=754854).  Need to take these strings and divide them into smaller pieces for processing

        	my @groups;
        	my $seq_length = length($seq);
        	if ($seq_length > 1000000) {
                	my $n = ceil($seq_length/1000000);
                	@groups = unpack "a$n" x (($seq_length/$n)-1) . "a*", $seq;
        	}
        	else {
                	push @groups, $seq;
        	}
        	my $new_seq;
        	foreach my $item (@groups) {
            	# remove end-of-lines
                	$item =~ s/\n//g;
                	$new_seq .= $item;
        	}
        	$seq = $new_seq;
        	}

        #print STDERR "$seq\n";
        my $cmdStatus = $?;
        if ( $cmdStatus != 0 ) {
            print STDERR "Storage::get_genome_seq: command failed with rc=$cmdStatus: $cmd\n";
            return;
        }
    }
    else {
        print STDERR "Storage::get_genome_seq: ERROR, file not found or zero length: $file_path\n";
        return;
    }
#    else {    # old method
#        $file_path =
#            CoGe::Accessory::Web::get_defaults()->{'SEQDIR'} . '/'
#          . get_tiered_path($gid)
#          . "/chr/$chr";
#
#        # Extract requested piece of sequence file
#        open( my $fh, $file_path ) or die "File not found: '$file_path'";
#        seek( $fh, $start - 1, 0 ) if ( $start > 1 );
#
#        $len = -s $file_path unless ( defined $len );
#        if ($fasta) {
#            $seq = ">$chr\n";
#            while ( $len > 0 ) {
#                my $count = read( $fh, my $s, min( $len, $FASTA_LINE_LENGTH ) );
#                $seq .= $s . "\n";
#                $len -= $count;
#            }
#        }
#        else {
#            read( $fh, $seq, $len );
#        }
#
#        close($fh);
#    }

    $seq = reverse_complement($seq)
      if ( defined $strand and $strand =~ /-/ );  #FIXME broken for fasta format

    return $seq;
}

sub get_experiment_path {
    my $eid = shift;
    return unless $eid;

    my $expdir = CoGe::Accessory::Web::get_defaults()->{'EXPDIR'};
    unless ($expdir) {
        print STDERR "Storage::get_experiment_path: WARNING, conf file parameter EXPDIR is blank!\n";
    }

    my $path = $expdir . '/' . get_tiered_path($eid) . '/';
    unless ( -r $path ) {
        print STDERR "Storage::get_experiment_path: experiment path '$path' doesn't exist!\n";
        return;
    }

    return $path;
}

sub get_experiment_files {
    my ($eid, $data_type)  = @_;
    return unless ($eid and $data_type);

    my $file_path = get_experiment_path($eid);
    my @files;

    if ($data_type == $DATA_TYPE_QUANT) {
        @files = glob "$file_path/*.csv";
    }
    elsif ($data_type == $DATA_TYPE_POLY) {
        @files = glob "$file_path/*.vcf";
    }
    elsif ($data_type == $DATA_TYPE_ALIGN) {
        @files = glob "$file_path/*.bam";
    }
    else {
        print STDERR "Storage::get_experiment_files: unknown data type $data_type!\n";
        return;
    }

    return \@files;
}

sub get_fastbit_format {
    my $eid = shift;
    my $data_type = shift;

    my $storage_path = get_experiment_path($eid);
    my $format_file = $storage_path . '/format.json';

    # Backward compatibility, see issue 352
    # FIXME: remove someday by adding format.json files to old experiments
    if (not -r $format_file) {
        if (!$data_type || $data_type == $DATA_TYPE_QUANT) {
            return {
                columns => [
                    { name => 'chr',    type => 'key' },
                    { name => 'start',  type => 'unsigned long' },
                    { name => 'stop',   type => 'unsigned long' },
                    { name => 'strand', type => 'byte' },
                    { name => 'value1', type => 'double' },
                    { name => 'value2', type => 'double' }
                ]
            };
        }
        elsif ($data_type == $DATA_TYPE_POLY) {
            return {
                columns => [
                    { name => 'chr',   type => 'key' },
                    { name => 'start', type => 'unsigned long' },
                    { name => 'stop',  type => 'unsigned long' },
                    { name => 'type',  type => 'key' },
                    { name => 'id',    type => 'text' },
                    { name => 'ref',   type => 'key' },
                    { name => 'alt',   type => 'key' },
                    { name => 'qual',  type => 'double' },
                    { name => 'info',  type => 'text' }
                ]
            };
        }
        elsif ( $data_type == $DATA_TYPE_MARKER ) {
            return {
                columns => [
                    { name => 'chr',    type => 'key' },
                    { name => 'start',  type => 'unsigned long' },
                    { name => 'stop',   type => 'unsigned long' },
                    { name => 'strand', type => 'key' },
                    { name => 'type',   type => 'key' },
                    { name => 'score',  type => 'double' },
                    { name => 'attr',   type => 'text' }
                ]
            };
        }
        return; # should never happen!
    }

    # Otherwise read format json file
    return CoGe::Accessory::TDS::read($format_file);
}

sub parse_fastbit_line {
    my $format = shift;
    my $line = shift;
    my $chr = shift;

    $line =~ s/"//g;
    my @items = split(/,\s*/, $line);
    return if ( $items[0] !~ /^\"?$chr/); # make sure it's a row output line

    my %result;
    foreach (@{$format->{columns}}) {
        my $item = shift @items;
        my $name = $_->{name};
        my $type = $_->{type};
        if ($type =~ /long|byte|double/) { $result{$name} = 0 + $item } # convert to numeric
        else { $result{$name} = '' . $item; }
    }
    return \%result;
}

sub get_experiment_data {
    my %opts = @_;
    my $eid  = $opts{eid};    # required
    my $data_type = $opts{data_type};   # required
    unless ($eid) {
        print STDERR "Storage::get_experiment_data: experiment id not specified!\n";
        return;
    }
    my $chr   = $opts{chr};
    my $start = $opts{start};
    my $stop  = $opts{stop};
    $stop = $opts{end} if ( not defined $stop );
    $start = 0 if ($start < 0);
    $stop = 0 if ($stop < 0);

    my $storage_path = get_experiment_path($eid);
    my $cmd;

    if (!$data_type ||
        $data_type == $DATA_TYPE_QUANT ||
        $data_type == $DATA_TYPE_POLY ||
        $data_type == $DATA_TYPE_MARKER)
    {
        my $pFormat = get_fastbit_format($eid, $data_type);
        my $columns = join(',', map { $_->{name} } @{$pFormat->{columns}});
        my $cmdpath = CoGe::Accessory::Web::get_defaults()->{FASTBIT_QUERY};
        $cmd = "$cmdpath -v 1 -d $storage_path -q \"select $columns where 0.0=0.0 and chr='$chr' and start <= $stop and stop >= $start order by start limit 999999999\" 2>&1";

        #print STDERR "\n$cmd\n";
        my @cmdOut = qx{$cmd};
        #print STDERR @cmdOut;
        my $cmdStatus = $?;
        if ( $? != 0 ) {
            print STDERR "Storage::get_experiment_data: error $? executing command: $cmd\n";
            return;
        }

        # Parse output into a hash result
        my @results;
        foreach (@cmdOut) {
            chomp;
            if (/^\"/) { #if (/^\"$chr\"/) { # potential result line
                my $pResult = parse_fastbit_line($pFormat, $_, $chr);
                push @results, $pResult if ($pResult);
            }
        }

        return \@results;
    }
    elsif ( $data_type == $DATA_TYPE_ALIGN ) { # FIXME move output parsing from Storage.pm to here
        my $cmdpath = CoGe::Accessory::Web::get_defaults()->{SAMTOOLS};
        $cmd = "$cmdpath view $storage_path/alignment.bam $chr:$start-$stop 2>&1";
        #print STDERR "$cmd\n";
        my @cmdOut = qx{$cmd};
        #print STDERR @cmdOut;
        my $cmdStatus = $?;
        if ( $? != 0 ) {
            print STDERR "Storage::get_experiment_data: error $? executing command: $cmd\n";
            return;
        }
        
        # Return if error message detected (starts with '[')
        map { return if (/^\[/) } @cmdOut; # mdb added 5/6/15 COGE-594
        
        return \@cmdOut;
    }
    else {
        print STDERR "Storage::get_experiment_data: unknown data type\n";
        return;
    }
}

sub get_log {
    my %opts = @_;
    my $item_id = $opts{item_id}; # required
    my $item_type  = $opts{item_type}; # required
    my $getEverything = $opts{getEverything}; # optional - not finished
    my $html = $opts{html}; # optional
    #print STDERR Dumper \%opts, "\n";
    return unless ($item_id and $item_type);
    
    my $storage_path;
    if ($item_type eq 'genome') { #TODO use hash of function refs for this instead
        $storage_path = get_genome_path($item_id);
    }
    elsif ($item_type eq 'experiment') {
        $storage_path = get_experiment_path($item_id);
    }
        
    my $old_log_file = catfile($storage_path, 'log.txt');
    
    my $log_file;
    if (-e $old_log_file) { # Old logging method pre-JEX
        $log_file = $old_log_file;
    }
    else { # Current JEX logging method
        # Get workflow id
        my $metadata_file = catfile($storage_path, 'metadata.json');
        if (! -e $metadata_file) {
            print STDERR "Storage::get_log cannot find log file\n";
            return;
        }
        my $md = CoGe::Accessory::TDS::read($metadata_file);
        my $workflow_id = $md->{workflow_id};
        
        # Determine path
        my (undef, $results_path) = get_workflow_paths(undef, $workflow_id);
        $log_file = catfile($results_path, 'debug.log');
    }
    
    # Read-in log file
    open(my $fh, $log_file);
    my $log = read_file($fh); # should be a relatively small file, read it all at once
    #my @log = <$fh>;
    #print STDERR \@log, "\n";
    close($fh);
    
    if ($html) {
        # Convert newlines/tabs to HTML equivalents
        $log =~ s/\t/&nbsp;&nbsp;&nbsp;&nbsp;/g;
        $log =~ s/\\t/&nbsp;&nbsp;&nbsp;&nbsp;/g;
        $log =~ s/\n/<br>/g;
        $log =~ s/\\n/<br>/g;
    }
    
    # Parse out content - NOT WORKING YET, LEFT OFF HERE
#    if (!$getEverything) {
#        #my @sections = $log =~ m{Command Output:(.*)#########################}sg;
#        
#        my $i = 0;
#        my @sections;
#        my $capture = 0;
#        my $line;
#        foreach (@log) {
#            if (/^Command Output:/) {
#                if ($capture) {
#                    $capture = 1;
#                    $line = '';
#                }
#            }
#            elsif (/#########################$/) {
#                $capture = 0;
#                $i++;
#            }
#            elsif ($capture) {
#                $line .= $_;
#                print $line, "\n";
#                push @{$sections[$i]}, $line;
#            }
#        }
#        
#        print STDERR \@sections, "\n";
#        #my @lines = split("/\/\n", ${sections[0]});
#    }
    
    return $log;
}

# mdb removed 5/13/15 - deprecated, instead use CoGe::Builder::Load::Experiment
#sub create_experiment {
#    my %opts = @_;
#    my $genome = $opts{genome}; # genome object or id
#    my $user = $opts{user};
#    my $irods = $opts{irods};
#    my $files = $opts{files};
#    #my $file_type = $opts{file_type};
#    my $metadata = $opts{metadata};
#    my $options = $opts{options};
#    #print STDERR "create_experiment ", Dumper $metadata, "\n";
#
#    my $conf = CoGe::Accessory::Web::get_defaults();
#
#    my $gid = $genome =~ /^\d+$/ ? $genome : $genome->id;
#
#    # Connect to workflow engine and get an id
#    my $jex = CoGe::Accessory::Jex->new( host => $conf->{JOBSERVER}, port => $conf->{JOBPORT} );
#    unless (defined $jex) {
#        return (undef, "Could not connect to JEX");
#    }
#
#    # Create the workflow
#    my $workflow = $jex->create_workflow( name => 'Create Experiment', init => 1 );
#    unless ($workflow and $workflow->id) {
#        return (undef, 'Could not create workflow');
#    }
#
#    # Setup log file, staging, and results paths
#    my ($staging_dir, $result_dir) = get_workflow_paths($user->name, $workflow->id);
#    $workflow->logfile( catfile($result_dir, 'debug.log') );
#
#    # Create list of files to load
#    my @staged_files;
#    push @staged_files, @$files if ($files);
#
#    # Create jobs to retrieve irods files
#    my %load_params;
#    foreach my $item (@$irods) {
#        next unless ($item->{type} eq 'irods');
#        %load_params = _create_iget_job($conf, $item->{path}, $staging_dir);
#        unless ( %load_params ) {
#            return (undef, "Could not create iget task");
#        }
#        $workflow->add_job(\%load_params);
#        push @staged_files, $load_params{outputs}[0];
#    }
#
#    # Create load job
#    my $ignoreMissing = ( $options->{ignoreMissing} ? 1 : 0 );
#    %load_params = _create_load_experiment_job(
#        conf => $conf, 
#        metadata => $metadata, 
#        gid => $gid, 
#        wid => $workflow->id, 
#        user_name => $user->name, 
#        staging_dir => $staging_dir, 
#        files => \@staged_files, 
#        #file_type => $file_type, 
#        ignoreMissing => $ignoreMissing
#    );
#    unless ( %load_params ) {
#        return (undef, "Could not create load task");
#    }
#    $workflow->add_job(\%load_params);
#
#    # Submit the workflow
#    my $result = $jex->submit_workflow($workflow);
#    if ($result->{status} =~ /error/i) {
#        return (undef, "Could not submit workflow");
#    }
#
#    return ($result->{id}, undef);
#}

sub create_experiments_from_batch {
    my %opts = @_;
    my $genome = $opts{genome};     # genome object or id
    my $notebook = $opts{notebook}; # optional notebook object or id
    my $user = $opts{user};         # user running the job
    my $assignee = $opts{assignee}; # user object or id to assign to
    my $irods = $opts{irods};
    my $files = $opts{files};
    my $metadata = $opts{metadata};

    my $conf = CoGe::Accessory::Web::get_defaults();

    my $gid = $genome =~ /^\d+$/ ? $genome : $genome->id;
    my $nid;
    if ($notebook) {
        $nid = $notebook =~ /^\d+$/ ? $notebook : $notebook->id;
    }

    # Connect to workflow engine and get an id
    my $jex = CoGe::Accessory::Jex->new( host => $conf->{JOBSERVER}, port => $conf->{JOBPORT} );
    unless (defined $jex) {
        return (undef, "Could not connect to JEX");
    }

    # Create the workflow
    my $workflow = $jex->create_workflow( name => 'Create Experiments', init => 1 );
    unless ($workflow and $workflow->id) {
        return (undef, 'Could not create workflow');
    }

    # Setup log file, staging, and results paths
    my ($staging_dir, $result_dir) = get_workflow_paths($user->name, $workflow->id);
    $workflow->logfile( catfile($result_dir, 'debug.log') );

    # Create list of files to load
    my @staged_files;
    push @staged_files, @$files;

    # Create jobs to retrieve irods files
    my %load_params;
    foreach my $item (@$irods) {
        next unless ($item->{type} eq 'irods');
        %load_params = _create_iget_job($conf, $item->{path}, $staging_dir);
        unless ( %load_params ) {
            return (undef, "Could not create iget task");
        }
        $workflow->add_job(\%load_params);
        push @staged_files, $load_params{outputs}[0];
    }
    
    # Change user if "assign to user" specified
    my $user_name = $user->name;
    $user_name = $assignee->name if ($assignee);

    # Create load job
    %load_params = _create_load_batch_job($conf, $metadata, $gid, $workflow->id, $user->name, \@staged_files, $staging_dir, $result_dir, $nid);
    unless ( %load_params ) {
        return (undef, "Could not create load task");
    }
    $workflow->add_job(\%load_params);

    # Submit the workflow
    my $result = $jex->submit_workflow($workflow);
    if ($result->{status} =~ /error/i) {
        return (undef, "Could not submit workflow");
    }

    return ($result->{id}, undef);
}

# Note: Using user names in file path instead of user ID.  This is okay because
# user names are guaranteed to only consist of letters, numbers, hyphens,
# and underscores.
sub get_workflow_paths {
    my ( $user_name, $workflow_id ) = remove_self(@_); # required because this routine is called internally and externally, is there a better way?
    unless ($workflow_id) {
        print STDERR "Storage::get_workflow_paths ERROR: missing required param\n";
        return;
    }
    
    my $tmp_path = CoGe::Accessory::Web::get_defaults()->{SECTEMPDIR};
    my ($staging_path, $results_path);
    if (!$user_name) {
        my $staging_dir = catdir($tmp_path, 'staging');
        my @tmp = read_dir($staging_dir);
        my @wdir = grep { -d "$staging_dir/$_/$workflow_id" } read_dir($staging_dir);
        print STDERR Dumper "wdir:\n", \@wdir, "\n";
        if (@wdir != 1) {
            print STDERR "Storage::get_workflow_paths ERROR: ambiguous user directory\n";
            return;
        }
        $user_name = $wdir[0];
    }

    $staging_path = catdir($tmp_path, 'staging', $user_name, $workflow_id);
    $results_path = catdir($tmp_path, 'results', $user_name, $workflow_id);
    
    return ($staging_path, $results_path);
}

sub add_workflow_result {
    my ( $user_name, $workflow_id, $result ) = @_;
    unless ($user_name && $workflow_id && $result) {
        print STDERR "Storage::add_result ERROR: missing required param\n";
        return;
    }
    
    my (undef, $results_path) = get_workflow_paths($user_name, $workflow_id);
    my $results_file = catfile($results_path, '.results');
    
    return CoGe::Accessory::TDS::append($results_file, { results => [ $result ] });
}

sub get_workflow_results {
    my ( $user_name, $workflow_id ) = @_;
    unless ($user_name && $workflow_id) {
        print STDERR "Storage::get_workflow_results ERROR: missing required param\n";
        return;
    }
    
    my (undef, $results_path) = get_workflow_paths($user_name, $workflow_id);
    my $results_file = catfile($results_path, '.results');
    
    my $results = CoGe::Accessory::TDS::read($results_file);
    return unless ($results && $results->{results} && @{$results->{results}});
    return $results->{results};
}

sub get_workflow_results_file {
    my ( $user_name, $workflow_id ) = remove_self(@_); # required because this routine is called internally and externally, is there a better way?
    unless ($user_name && $workflow_id) {
        print STDERR "Storage::get_workflow_results_file ERROR: missing required param\n";
        return;
    }
    
    my (undef, $results_path) = get_workflow_paths($user_name, $workflow_id);
    my $results_file = catfile($results_path, '.results');
    return $results_file;
}

sub get_workflow_log_file {
    my ( $user_name, $workflow_id ) = remove_self(@_); # required because this routine is called internally and externally, is there a better way?
    unless ($user_name && $workflow_id) {
        print STDERR "Storage::get_workflow_log_file ERROR: missing required param\n";
        return;
    }
    
    my (undef, $results_path) = get_workflow_paths($user_name, $workflow_id);
    my $results_file = catfile($results_path, 'debug.log');
    return $results_file;
}

sub get_upload_path {
    my ( $user_name, $load_id ) = remove_self(@_); # required because this routine is called internally and externally, is there a better way?
    unless ($user_name && $load_id) {
         print STDERR "Storage::get_upload_path ERROR: missing required param\n";
        return;
    }
    
    my $conf = CoGe::Accessory::Web::get_defaults();
    return catdir($conf->{SECTEMPDIR}, 'uploads', $user_name, $load_id);
}

sub get_download_path {
    my ($type, $id, $uuid) = @_;
    $uuid = '' unless $uuid; # optional uuid
    my $conf = CoGe::Accessory::Web::get_defaults();
    return catfile($conf->{SECTEMPDIR}, 'downloads', $type, $id, $uuid);
}

sub remove_self { # TODO move to Utils.pm
    return @_ unless ( ref($_[0]) );
    my @a = @_[1..$#_];
    return wantarray ? @a : \@a;
}

sub create_genome_from_file {
    my %opts = @_;
    my $user = $opts{user};
    my $irods = $opts{irods};
    my $files = $opts{files};
    my $metadata = $opts{metadata};

    # Connect to workflow engine and get an id
    my $conf = CoGe::Accessory::Web::get_defaults();
    my $jex = CoGe::Accessory::Jex->new( host => $conf->{JOBSERVER}, port => $conf->{JOBPORT} );
    unless (defined $jex) {
        return (undef, "Could not connect to JEX");
    }

    # Create the workflow
    my $workflow = $jex->create_workflow( name => 'Create Genome', init => 1 );
    unless ($workflow and $workflow->id) {
        return (undef, 'Could not create workflow');
    }

    # Setup log file, staging, and results paths
    my ($staging_dir, $result_dir) = get_workflow_paths($user->name, $workflow->id);
    $workflow->logfile( catfile($result_dir, 'debug.log') );

    # Create list of files to load
    my @staged_files;
    push @staged_files, @$files if ($files);

    # Create jobs to retrieve irods files
    my %load_params;
    foreach my $item (@$irods) {
        next unless ($item->{type} eq 'irods');
        %load_params = _create_iget_job($conf, $item->{path}, $staging_dir);
        unless ( %load_params ) {
            return (undef, "Could not create iget task");
        }
        $workflow->add_job(\%load_params);
        push @staged_files, $load_params{outputs}[0];
    }

    # Create load job
    %load_params = _create_load_genome_job($conf, $metadata, $user->name, $staging_dir, \@staged_files, $result_dir, $irods);
    unless ( %load_params ) {
        return (undef, "Could not create load task");
    }
    $workflow->add_job(\%load_params);

    # Submit the workflow
    my $result = $jex->submit_workflow($workflow);
    if ($result->{status} =~ /error/i) {
        return (undef, "Could not submit workflow");
    }

    return ($result->{id}, undef);
}

sub create_genome_from_NCBI {
    my %opts = @_;
    my $user = $opts{user};
    my $accns = $opts{accns};

    # Connect to workflow engine and get an id
    my $conf = CoGe::Accessory::Web::get_defaults();
    my $jex = CoGe::Accessory::Jex->new( host => $conf->{JOBSERVER}, port => $conf->{JOBPORT} );
    unless (defined $jex) {
        return (undef, "Could not connect to JEX");
    }

    # Create the workflow
    my $workflow = $jex->create_workflow( name => 'Create Genome (NCBI)', init => 1 );
    unless ($workflow and $workflow->id) {
        return (undef, 'Could not create workflow');
    }

    # Setup log file, staging, and results paths
    my ($staging_dir, $result_dir) = get_workflow_paths($user->name, $workflow->id);
    $workflow->logfile( catfile($result_dir, 'debug.log') );

    # Create load job
    my %load_params = _create_load_genome_from_NCBI_job($conf, $user->name, $accns, $staging_dir, $result_dir);
    unless ( %load_params ) {
        return (undef, "Could not create load task");
    }
    $workflow->add_job(\%load_params);

    # Submit the workflow
    my $result = $jex->submit_workflow($workflow);
    if ($result->{status} =~ /error/i) {
        return (undef, "Could not submit workflow");
    }

    return ($result->{id}, undef);
}

sub _create_iget_job {
    my ($conf, $irods_path, $staging_dir) = @_;

    my $dest_file = catdir($staging_dir, 'irods', $irods_path);
    my $dest_path = dirname($dest_file);
    mkpath($dest_path);
    my $cmd = irods_iget( $irods_path, $dest_path, { no_execute => 1 } );

    return (
        cmd => $cmd,
        script => undef,
        args => [],
        inputs => [],
        outputs => [
            $dest_file
        ],
        description => "Fetching $irods_path..."
    );
}

# mdb removed 5/13/15 - deprecated
#sub _create_load_experiment_job {
#    my %opts = @_;
#    my $conf = $opts{conf};
#    my $metadata = $opts{metadata};
#    my $gid = $opts{gid};
#    my $wid = $opts{wid};
#    my $user_name = $opts{user_name};
#    my $staging_dir = $opts{staging_dir};
#    my $files = $opts{files};
#    #my $file_type = $opts{file_type};
#    my $ignoreMissing = $opts{ignoreMissing};
#    #print STDERR "_create_load_experiment_job ", Dumper $metadata, "\n";
#
#    my $cmd = catfile($conf->{SCRIPTDIR}, "load_experiment.pl");
#    return unless $cmd; # SCRIPTDIR undefined
#
#    my $file_str = join(',', map { basename($_) } @$files);
#    #$file_type = 'csv' unless $file_type;
#
#    return (
#        cmd => $cmd,
#        script => undef,
#        args => [
#            ['-user_name', $user_name, 0],
#            ['-name', '"' . $metadata->{name} . '"', 0],
#            ['-desc', '"' . $metadata->{description} . '"', 0],
#            ['-version', '"' . $metadata->{version} . '"', 0],
#            ['-restricted', ( $metadata->{restricted} ? 1 : 0 ), 0],
#            ['-gid', $gid, 0],
#            ['-wid', $wid, 0],
#            ['-source_name', '"' . $metadata->{source_name} . '"', 0],
#            #['-types', qq{"Expression"}, 0], # FIXME
#            #['-annotations', $ANNOTATIONS, 0],
#            ['-staging_dir', "'".$staging_dir."'", 0],
#            #['-file_type', $file_type, 0], # FIXME
#            ['-data_file', "'".$file_str."'", 0],
#            ['-config', $conf->{_CONFIG_PATH}, 1],
#            #['-result_dir', "'".$result_dir."'", 0],
#            ['-ignore-missing-chr', $ignoreMissing, 0]
#        ],
#        inputs => [
#            ($conf->{_CONFIG_PATH}, @$files)
#        ],
#        outputs => [
#            [$staging_dir, 1],
#            catdir($staging_dir, 'log.done')
#        ],
#        description => "Loading experiment data..."
#    );
#}

sub _create_load_batch_job {
    my ($conf, $metadata, $gid, $wid, $user_name, $files, $staging_dir, $result_dir, $nid) = @_;
    my $cmd = catfile($conf->{SCRIPTDIR}, "load_batch.pl");
    return unless $cmd; # SCRIPTDIR undefined

    #my $file_str = join(',', map { basename($_) } @$files);
    my $file_str = join(',', @$files);

    my %params = (
        cmd => $cmd,
        script => undef,
        args => [
            ['-user_name', $user_name, 0],
            ['-name', '"' . $metadata->{name} . '"', 0],
            ['-desc', '"' . $metadata->{description} . '"', 0],
            ['-gid', $gid, 0],
            ['-wid', $wid, 0],
            ['-staging_dir', "'".$staging_dir."'", 0],
            ['-result_dir', "'".$result_dir."'", 0],
            ['-files', "'".$file_str."'", 0],
            ['-config', $conf->{_CONFIG_PATH}, 1]
        ],
        inputs => [
            ($conf->{_CONFIG_PATH}, @$files)
        ],
        outputs => [
            [$staging_dir, 1],
            catdir($staging_dir, 'log.done')
        ],
        description => "Loading batch experiments..."
    );
    
    push $params{args}, ['-nid', $nid, 0] if ($nid);
    
    return %params;
}

sub _create_load_genome_job {
    my ($conf, $metadata, $user_name, $staging_dir, $files, $result_dir, $irods) = @_;
    my $cmd = catfile($conf->{SCRIPTDIR}, "load_genome.pl");
    return unless $cmd; # SCRIPTDIR undefined

    my $file_str = join(',', map { basename($_) } @$files);
    my $irods_str = join(',', map { basename($_) } @$irods);

    return (
        cmd => $cmd,
        script => undef,
        args => [
            ['-user_name', $user_name, 0],
            ['-name', '"' . $metadata->{name} . '"', 0],
            ['-desc', '"' . $metadata->{description} . '"', 0],
            ['-version', '"' . $metadata->{version} . '"', 0],
            ['-restricted', ( $metadata->{restricted} ? 1 : 0 ), 0],
            ['-source_name', '"' . $metadata->{source_name} . '"', 0],
            ['-organism_id', '"' . $metadata->{organism_id} . '"', 0],
            ['-type_id', '"' . $metadata->{type_id} . '"', 0],
            ['-staging_dir', "'".$staging_dir."'", 0],
            ['-fasta_files', "'".$file_str."'", 0],
            ['-irods_files', "'".$irods_str."'", 0],
            ['-config', $conf->{_CONFIG_PATH}, 1],
            ['-result_dir', "'".$result_dir."'", 0]
        ],
        inputs => [
            ($conf->{_CONFIG_PATH}, @$files)
        ],
        outputs => [
            [$staging_dir, 1],
            catdir($staging_dir, 'log.done')
        ],
        description => "Loading genome data ..."
    );
}

sub _create_load_genome_from_NCBI_job {
    my ($conf, $user_name, $accns, $staging_dir, $result_dir, $files) = @_;
    my $cmd = catfile($conf->{SCRIPTDIR}, "genbank_genome_loader.pl");
    return unless $cmd; # SCRIPTDIR undefined

    my %p = (
        cmd => $cmd,
        script => undef,
        args => [
            ['-user_name', $user_name, 0],
            ['-staging_dir', "'".$staging_dir."'", 0],
            ['-result_dir', "'".$result_dir."'", 0],
            ['-config', $conf->{_CONFIG_PATH}, 1],
            ['-GO', 1, 0]
        ],
        inputs => [
            $conf->{_CONFIG_PATH}
        ],
        outputs => [
            [$staging_dir, 1],
            catdir($staging_dir, 'log.done')
        ],
        description => "Loading genome data ..."
    );

    foreach (@$accns) {
        push @{ $p{args} }, ['-accn', "'".$_."'", 0];
    }

    return %p;
}

sub create_annotation_dataset {
    my %opts = @_;
    my $user = $opts{user};
    my $irods = $opts{irods};
    my $files = $opts{files};
    my $metadata = $opts{metadata};

    # Connect to workflow engine and get an id
    my $conf = CoGe::Accessory::Web::get_defaults();
    my $jex = CoGe::Accessory::Jex->new( host => $conf->{JOBSERVER}, port => $conf->{JOBPORT} );
    unless (defined $jex) {
        return (undef, "Could not connect to JEX");
    }

    # Create the workflow
    my $workflow = $jex->create_workflow( name => 'Load Genome Annotation', init => 1 );
    unless ($workflow and $workflow->id) {
        return (undef, 'Could not create workflow');
    }

    # Setup log file, staging, and results paths
    my ($staging_dir, $result_dir) = get_workflow_paths($user->name, $workflow->id);
    $workflow->logfile( catfile($result_dir, 'debug.log') );

    # Create list of files to load
    my @staged_files;
    push @staged_files, @$files;

    # Create jobs to retrieve irods files
    my %load_params;
    foreach my $item (@$irods) {
        next unless ($item->{type} eq 'irods');
        %load_params = _create_iget_job($conf, $item->{path}, $staging_dir);
        unless ( %load_params ) {
            return (undef, "Could not create iget task");
        }
        $workflow->add_job(\%load_params);
        push @staged_files, $load_params{outputs}[0];
    }

    # Create load job
    %load_params = _create_load_annotation_job($conf, $metadata, $user->name, $staging_dir, \@staged_files, $result_dir);
    unless ( %load_params ) {
        return (undef, "Could not create load task");
    }
    $workflow->add_job(\%load_params);

    # Submit the workflow
    my $result = $jex->submit_workflow($workflow);
    if ($result->{status} =~ /error/i) {
        return (undef, "Could not submit workflow");
    }

    return ($result->{id}, undef);
}

sub _create_load_annotation_job {
    my ($conf, $metadata, $user_name, $staging_dir, $files, $result_dir) = @_;
    my $cmd = catfile($conf->{SCRIPTDIR}, "load_annotation.pl");
    return unless $cmd; # SCRIPTDIR undefined

    my $file_str = join(',', map { basename($_) } @$files);

    return (
        cmd => $cmd,
        script => undef,
        args => [
            ['-user_name', $user_name, 0],
            ['-name', '"' . $metadata->{name} . '"', 0],
            ['-desc', '"' . $metadata->{description} . '"', 0],
            ['-link', '"' . $metadata->{link} . '"', 0],
            ['-version', '"' . $metadata->{version} . '"', 0],
            ['-restricted', ( $metadata->{restricted} ? 1 : 0 ), 0],
            ['-source_name', '"' . $metadata->{source_name} . '"', 0],
            ['-gid', '"' . $metadata->{genome_id} . '"', 0],
            ['-staging_dir', "'".$staging_dir."'", 0],
            ['-data_file', "'".$file_str."'", 0],
            ['-config', $conf->{_CONFIG_PATH}, 1],
            ['-result_dir', "'".$result_dir."'", 0]
        ],
        inputs => [
            ($conf->{_CONFIG_PATH}, @$files)
        ],
        outputs => [
            [$staging_dir, 1],
            catdir($staging_dir, 'log.done')
        ],
        description => "Loading annotation data ..."
    );
}

sub get_irods_path {
    my ($path) = @_;
    #print STDERR "irods_get_path: ", $path, "\n";
    
# mdb removed 6/17/15 COGE-313
#    if ( $path !~ /^$basepath/ ) {
#        print STDERR "Attempt to access '$path' denied (basepath='$basepath')\n";
#        return;
#    }

    # mdb added 6/17/15 COGE-313
    if ( $path eq '/iplant/home/' ) {
        #print STDERR "Attempt to access '$path' denied\n";
        return;
    }

    my $result = CoGe::Accessory::IRODS::irods_ils($path, escape_output => 1);
    #print STDERR Dumper $result, "\n";

    return $result;
}

sub get_irods_file {
    my ($src_path, $dest_path) = @_;
    my ($filename)   = $src_path =~ /([^\/]+)\s*$/;
    my ($remotepath) = $src_path =~ /(.*)$filename$/;

    my $localpath     = 'irods/' . $remotepath;
    my $localfullpath = catdir($dest_path . $localpath);
    $localpath = catfile($localpath, $filename);
    my $localfilepath = catfile($localfullpath, $filename);
    #print STDERR "get_file $path $filename $localfilepath\n";

    my $do_get = 1;

    #   if (-e $localfilepath) {
    #       my $remote_chksum = irods_chksum($path);
    #       my $local_chksum = md5sum($localfilepath);
    #       $do_get = 0 if ($remote_chksum eq $local_chksum);
    #       print STDERR "$remote_chksum $local_chksum\n";
    #   }

    if ($do_get) {
        mkpath($localfullpath);
        CoGe::Accessory::IRODS::irods_iget( $src_path, $localfullpath );
    }

    return { path => $localpath, size => -s $localfilepath };
}

sub reverse_complement { #TODO move into Util.pm
    my $seq   = shift;
    my $rcseq = reverse($seq);
    $rcseq =~ tr/ATCGatcg/TAGCtagc/;
    return $rcseq;
}

sub data_type {
    my $data_type = shift;

    # Experiment Data Types
    return "Quantitative" if $data_type == $DATA_TYPE_QUANT;
    return "Polymorphism" if $data_type == $DATA_TYPE_POLY;
    return "Alignment" if $data_type == $DATA_TYPE_ALIGN;
    return "Marker" if $data_type == $DATA_TYPE_MARKER;
    return "Unknown";
}

1;
