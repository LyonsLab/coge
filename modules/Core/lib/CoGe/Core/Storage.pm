package CoGe::Core::Storage;

=head1 NAME

CoGe::Core::Storage

=head1 SYNOPSIS

Abstraction layer on top of genome storage sub-system created for issues
77 and 157.  All accesses to genome FASTA sequences should go through this
module.

=head1 DESCRIPTION

=head1 AUTHOR

Matt Bomhoff

=head1 COPYRIGHT

The full text of the license can be found in the
LICENSE file included with this module.

=head1 SEE ALSO

=cut

#use v5.14;
use strict;
use warnings;

use CoGe::Accessory::Web qw(get_defaults);
use CoGe::Accessory::TDS qw(read);
use CoGe::Accessory::Jex;
use CoGe::Accessory::Workflow;
use CoGe::Accessory::IRODS qw(irods_iget);
use File::Basename;
use POSIX qw(floor);
use File::Spec::Functions;
use File::Path qw(mkpath);
use List::Util qw[min max];
use File::Slurp;
use JSON::XS qw(decode_json);
use Data::Dumper;
use POSIX qw(ceil);

BEGIN {
    use vars qw ($VERSION @ISA @EXPORT $DATA_TYPE_QUANT $DATA_TYPE_POLY $DATA_TYPE_ALIGN $DATA_TYPE_MARKER);
    require Exporter;

    $VERSION = 0.1;
    @ISA     = qw (Exporter);
    @EXPORT = qw(
      get_tiered_path get_workflow_paths
      get_genome_file index_genome_file get_genome_seq get_genome_path
      get_experiment_path get_experiment_files get_experiment_data
      create_experiment create_experiments_from_batch
      create_genome_from_file create_genome_from_NCBI
      create_annotation_dataset reverse_complement
    );

    # Experiment Data Types
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
        return;
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
        print STDERR
          "Storage::index_genome_file: command failed with rc=$?: $cmd\n";
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
            print STDERR
              "Storage::index_genome_file: command failed with rc=$?: $cmd\n";
            return $?;
        }

        # Index compressed fasta file
        $cmd = "$samtools faidx $file_path.razf";
        qx{ $cmd };
        if ( $? != 0 ) {
            print STDERR
              "Storage::index_genome_file: command failed with rc=$?: $cmd\n";
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
        return \@cmdOut;
    }
    else {
        print STDERR "Storage::get_experiment_data: unknown data type\n";
        return;
    }
}

sub create_experiment {
    my %opts = @_;
    my $genome = $opts{genome}; # genome object or id
    my $user = $opts{user};
    my $irods = $opts{irods};
    my $files = $opts{files};
    my $file_type = $opts{file_type};
    my $metadata = $opts{metadata};

    print STDERR (caller(0))[3], "\n";

    my $conf = CoGe::Accessory::Web::get_defaults();

    my $gid = $genome =~ /^\d+$/ ? $genome : $genome->id;

    # Connect to workflow engine and get an id
    my $jex = CoGe::Accessory::Jex->new( host => $conf->{JOBSERVER}, port => $conf->{JOBPORT} );
    unless (defined $jex) {
        return (undef, "Could not connect to JEX");
    }

    # Create the workflow
    my $workflow = $jex->create_workflow( name => 'Create Experiment', init => 1 );
    unless ($workflow and $workflow->id) {
        return (undef, 'Could not create workflow');
    }

    # Setup log file, staging, and results paths
    my ($staging_dir, $result_dir) = get_workflow_paths($user->name, $workflow->id);
    $workflow->logfile( catfile($staging_dir, 'log_main.txt') );

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
        $workflow->add_job(%load_params);
        push @staged_files, $load_params{outputs}[0];
    }

    # Create load job
    %load_params = _create_load_experiment_job($conf, $metadata, $gid, $user->name, $staging_dir, \@staged_files, $file_type, $result_dir);
    unless ( %load_params ) {
        return (undef, "Could not create load task");
    }
    $workflow->add_job(%load_params);

    # Submit the workflow
    my $result = $jex->submit_workflow($workflow);
    if ($result->{status} =~ /error/i) {
        return (undef, "Could not submit workflow");
    }

    return ($result->{id}, undef);
}

sub create_experiments_from_batch {
    my %opts = @_;
    my $genome = $opts{genome}; # genome object or id
    my $user = $opts{user};
    my $irods = $opts{irods};
    my $files = $opts{files};
    my $metadata = $opts{metadata};

    print STDERR (caller(0))[3], "\n";

    my $conf = CoGe::Accessory::Web::get_defaults();

    my $gid = $genome =~ /^\d+$/ ? $genome : $genome->id;

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
    $workflow->logfile( catfile($staging_dir, 'workflow.log') );

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
        $workflow->add_job(%load_params);
        push @staged_files, $load_params{outputs}[0];
    }

    # Create load job
    %load_params = _create_load_batch_job($conf, $metadata, $gid, $user->name, \@staged_files, $staging_dir, $result_dir);
    unless ( %load_params ) {
        return (undef, "Could not create load task");
    }
    $workflow->add_job(%load_params);

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
    unless ($user_name and $workflow_id) {
        print STDERR "Storage::get_workflow_paths ERROR: missing required param\n";
        return;
    }

    my $tmp_path = CoGe::Accessory::Web::get_defaults()->{SECTEMPDIR};
    my $staging_path = catdir($tmp_path, 'staging', $user_name, $workflow_id);
    my $results_path = catdir($tmp_path, 'results', $user_name, $workflow_id);
    return ($staging_path, $results_path);
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

    print STDERR (caller(0))[3], "\n";

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
    $workflow->logfile( catfile($staging_dir, 'log_main.txt') );

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
        $workflow->add_job(%load_params);
        push @staged_files, $load_params{outputs}[0];
    }

    # Create load job
    %load_params = _create_load_genome_job($conf, $metadata, $user->name, $staging_dir, \@staged_files, $result_dir, $irods);
    unless ( %load_params ) {
        return (undef, "Could not create load task");
    }
    $workflow->add_job(%load_params);

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

    print STDERR (caller(0))[3], "\n";

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
    $workflow->logfile( catfile($staging_dir, 'log_main.txt') );

    # Create load job
    my %load_params = _create_load_genome_from_NCBI_job($conf, $user->name, $accns, $staging_dir, $result_dir);
    unless ( %load_params ) {
        return (undef, "Could not create load task");
    }
    $workflow->add_job(%load_params);

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

sub _create_load_experiment_job {
    my ($conf, $metadata, $gid, $user_name, $staging_dir, $files, $file_type, $result_dir) = @_;
    my $cmd = catfile($conf->{SCRIPTDIR}, "load_experiment.pl");
    return unless $cmd; # SCRIPTDIR undefined

    my $file_str = join(',', map { basename($_) } @$files);
    $file_type = 'csv' unless $file_type;

    return (
        cmd => $cmd,
        script => undef,
        args => [
            ['-user_name', $user_name, 0],
            ['-name', '"' . $metadata->{name} . '"', 0],
            ['-desc', '"' . $metadata->{description} . '"', 0],
            ['-version', '"' . $metadata->{version} . '"', 0],
            ['-restricted', ( $metadata->{restricted} ? 1 : 0 ), 0],
            ['-gid', $gid, 0],
            ['-source_name', '"' . $metadata->{source_name} . '"', 0],
            #['-types', qq{"Expression"}, 0], # FIXME
            #['-annotations', $ANNOTATIONS, 0],
            ['-staging_dir', "'".$staging_dir."'", 0],
            ['-file_type', $file_type, 0], # FIXME
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
        description => "Loading experiment data..."
    );
}

sub _create_load_batch_job {
    my ($conf, $metadata, $gid, $user_name, $files, $staging_dir, $result_dir) = @_;
    my $cmd = catfile($conf->{SCRIPTDIR}, "load_batch.pl");
    return unless $cmd; # SCRIPTDIR undefined

    my $file_str = join(',', map { basename($_) } @$files);

#    my $cmd =
#        "$BINDIR/load_batch.pl "
#      . "-user_name $user_name "
#      . '-name "' . escape($name) . '" '
#      . '-desc "' . escape($description) . '" '
#      . "-gid $gid "
#      . "-staging_dir $stagepath "
#      . '-data_file "' . escape( join( ',', @files ) ) . '" '
#      . "-config $CONFIGFILE";

    return (
        cmd => $cmd,
        script => undef,
        args => [
            ['-user_name', $user_name, 0],
            ['-name', '"' . $metadata->{name} . '"', 0],
            ['-desc', '"' . $metadata->{description} . '"', 0],
            ['-gid', $gid, 0],
            ['-staging_dir', "'".$staging_dir."'", 0],
            ['-result_dir', "'".$result_dir."'", 0],
            ['-data_file', "'".$file_str."'", 0],
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
}

sub _create_load_genome_job {
    my ($conf, $metadata, $user_name, $staging_dir, $files, $result_dir, $irods) = @_;
    my $cmd = catfile($conf->{SCRIPTDIR}, "load_genome.pl");
    return unless $cmd; # SCRIPTDIR undefined

    my $file_str = join(',', map { basename($_) } @$files);
    my $irods_str = join(',', map { basename($_) } @$irods);

    $file_str = '' if ($irods_str);

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
    my $cmd = catfile($conf->{SCRIPTDIR}, "load_genomes_n_stuff", "genbank_genome_loader.pl");
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

    print STDERR (caller(0))[3], "\n";

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
    $workflow->logfile( catfile($staging_dir, 'log_main.txt') );

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
        $workflow->add_job(%load_params);
        push @staged_files, $load_params{outputs}[0];
    }

    # Create load job
    %load_params = _create_load_annotation_job($conf, $metadata, $user->name, $staging_dir, \@staged_files, $result_dir);
    unless ( %load_params ) {
        return (undef, "Could not create load task");
    }
    $workflow->add_job(%load_params);

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

sub reverse_complement { #TODO move into Util.pm
    my $seq   = shift;
    my $rcseq = reverse($seq);
    $rcseq =~ tr/ATCGatcg/TAGCtagc/;
    return $rcseq;
}

1;
