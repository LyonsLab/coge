package CoGe::Builder::Buildable;

use Moose;

use Array::Utils qw(array_minus);
use File::Basename qw(basename dirname);
use File::Spec::Functions qw(catdir catfile);
use File::Path qw(make_path);
use JSON qw(encode_json);
use URI::Escape::JavaScript qw(escape);
use String::ShellQuote;
use Data::Dumper;

use CoGe::Accessory::IRODS qw(irods_set_env irods_iput);
use CoGe::Accessory::Web qw(get_command_path get_tiny_link url_for);
use CoGe::Accessory::Utils;
use CoGe::Core::Storage;
use CoGe::Core::Metadata qw(tags_to_string);
use CoGe::Builder::Data::Extractor;
use CoGe::Exception::Generic;

# Public attributes
has 'request'       => ( is => 'ro', isa => 'CoGe::Request::Request', required => 1 );
has 'jex'           => ( is => 'rw', isa => 'CoGe::JEX::Jex' );
has 'workflow'      => ( is => 'rw', isa => 'CoGe::JEX::Workflow' );
has 'site_url'      => ( is => 'rw' );
has 'page'          => ( is => 'rw' );

# Private attributes
has 'tasks'         => ( is => 'rw', isa => 'ArrayRef', default => sub { [] });
has 'staging_dir'   => ( is => 'rw');#, traits => ['Private'] );
has 'result_dir'    => ( is => 'rw');#, traits => ['Private'] );
has 'metadata'      => ( is => 'rw');#, traits => ['Private'] );

# requires 'build'; # must be implemented in sub-class

my $previous_outputs = [];

my $PICARD;
sub BUILD { # called immediately after constructor
    my $self = shift;
    $PICARD = $self->conf->{PICARD};
    unless ($PICARD) {
        CoGe::Exception::Generic->throw(message => 'Missing PICARD in config file');
    }
}

# This allows us to instantiate subclasses with a single arg $self.
# Called before constructor.
around BUILDARGS => sub {
    my $orig = shift;
    my $class = shift;

    if ( @_ == 1 && ref $_[0]) { # $self is only argument
        return $class->$orig(
            request     => $_[0]->request,
            workflow    => $_[0]->workflow,
            staging_dir => $_[0]->staging_dir,
            result_dir  => $_[0]->result_dir
        );
    }
    else { # canonical arguments
        return $class->$orig(@_);
    }
};

sub type    { shift->request->payload->{type}       }
sub params  { shift->request->payload->{parameters} }
sub user    { shift->request->user }
sub db      { shift->request->db   }
sub conf    { shift->request->conf }

sub get_name { # override this
    return '';
}

sub get_site_url { # override this
    return '';
}

sub submit {
    my $self = shift;

    # Add tasks to workflow
    unless ($self->workflow->add_jobs($self->tasks)) {
        return {
            success => JSON::false,
            error   => { Error => "failed to build workflow" }
        };
    }

    # Dump info to file for debugging
    if ($self->result_dir) {
        #my $cmd = 'chmod g+rw ' . $self->result_dir;
        #`$cmd`;

        # Dump raw workflow
        open(my $fh, '>', catfile($self->result_dir, 'workflow.log'));
        print $fh Dumper $self->workflow, "\n";
        close($fh);

        # Dump params
        open($fh, '>', catfile($self->result_dir, 'params.log'));
        print $fh Dumper $self->request->payload, "\n";
        close($fh);
    }

    # Submit workflow to JEX
    my $resp = $self->jex->submit_workflow($self->workflow);
    my $success = $self->jex->is_successful($resp);
    unless ($success) {
        print STDERR 'JEX response: ', Dumper $resp, "\n";
        return {
            success => JSON::false,
            error => { JEX => 'failed to submit workflow' }
            #TODO return $resp error message from JEX
        }
    }

    # Success
    my $response = {
        id => $resp->{id},
        success => JSON::true
    };
    $response->{site_url} = $self->site_url if ($self->site_url);
    return $response;
}

sub pre_build { # Default method, SynMap & SynMap3D override this
    my $self = shift;

    # Connect to JEX
    my $jex = CoGe::JEX::Jex->new( host => $self->conf->{JOBSERVER}, port => $self->conf->{JOBPORT} );
    unless ($jex) {
        CoGe::Exception::Generic->throw(message => "Couldn't connect to JEX");
    }
    $self->jex($jex);

    # Initialize workflow -- NOTE: init => 1 means that a new workflow will be created right away
    $self->workflow( $jex->create_workflow(name => $self->get_name, init => 1 ) );
    unless ($self->workflow && $self->workflow->id) {
        CoGe::Exception::Generic->throw(message => "Couldn't create workflow");
    }

    # Setup workflow paths
    my ($staging_dir, $result_dir) = get_workflow_paths(($self->user ? $self->user->name : 'public'), $self->workflow->id);
    make_path($staging_dir);
    make_path($result_dir);
    $self->staging_dir($staging_dir);
    $self->result_dir($result_dir);
    $self->workflow->logfile( catfile($result_dir, "debug.log") );

    # Set "page" and "site_url" attributes
    my $site_url = $self->get_site_url();
    my $requester = $self->request->requester;
    if ($requester) { # requester is from internal web page - external API requests will not have a 'requester' field
        my $page = $requester->{page}; # page name used for logging
        my $url  = $requester->{url};  # used to set site_url with special query params
        $self->page($page) if $page;
        unless ($site_url) {
            $url = $page unless $url;
            $site_url = url_for($url, wid => $self->workflow->id);
        }
    }
    $self->site_url( get_tiny_link(url => $site_url) ) if $site_url;

    # Add data retrieval tasks
    my $data = $self->params->{source_data};
    if ($data) {
        my $dr = CoGe::Builder::Data::Extractor->new($self);
        $dr->build($data);
        $self->add($dr);
        return $dr;
    }

    return;
}

sub post_build {
    my $self = shift;

    # Add task to add results to notebook
    if ($self->params->{notebook} || $self->params->{notebook_id}) {
        #TODO add_items_to_notebook_job and create_notebook_job and their respective scripts can be consolidated
        if ($self->params->{notebook_id}) { # Use existing notebook
            $self->add_to_all(
                $self->add_items_to_notebook( notebook_id => $self->params->{notebook_id} )
            );
        }
        else { # Create new notebook
            $self->add_to_all(
                $self->create_notebook( metadata => $self->params->{metadata} )
            );
        }
    }

    # Add task to send notification email
    if ( $self->params->{email} ) {
        my $email = $self->params->{email};
        $email = $self->user->email if ((!$email || length($email) < 3) && $self->user && $self->user->email);
        if ($email && length($email) >= 3) {
            $self->add_to_all(
                $self->send_email(
                    to      => $email,
                    subject => 'CoGe Workflow Notification',
                    body    => 'Your workflow finished: '.$self->get_name().
                        ($self->site_url ? "\nLink: ".$self->site_url : '').
                        "\n\nNote: you received this email because you submitted a job on ".
                        "CoGe (http://genomevolution.org) and selected the option to be notified."
                )
            );
        }
    }

    # Add task to send notification to callback url
    if ( $self->params->{callback_url} ) {
        $self->add_to_all(
            $self->curl_get( url => $self->params->{callback_url} )
        );
    }
}

# Add independent task(s)
sub add {
    my ($self, $tasks, $dependencies) = @_;
    unless ($tasks) {
        CoGe::Exception::Generic->throw(message => "Missing tasks");
    }

    if (ref($tasks) eq 'HASH') { # single task hash ref
        $tasks = [ $tasks ];
    }
    elsif (ref($tasks) eq 'ARRAY') { # array of task hash refs
        # nothing to do
    }
    else { # assume Buildable object
        unless ($tasks->can('tasks')) {
            CoGe::Exception::Generic->throw(message => "Invalid tasks", details => Dumper $tasks);
        }
        $tasks = $tasks->tasks;
    }

    # Chain task(s) to given dependencies
    if ($dependencies) {
        $dependencies = [ $dependencies ] unless (ref($dependencies) eq 'ARRAY');
        foreach (@$tasks) {
            my @new_dep = array_minus(@$dependencies, @{$_->{inputs}}); # prevent duplicates
            push @{$_->{inputs}}, @new_dep;
        }
    }

    # Add tasks
    push @{$self->tasks}, @$tasks;

    # Save outputs
    $previous_outputs = [];
    foreach (@$tasks) {
        push @$previous_outputs, @{$_->{outputs}};
    }

    return wantarray ? @$previous_outputs : $previous_outputs;
}

# Chain task(s) to the previous task's outputs
sub add_to_previous {
    my ($self, $tasks) = @_;
    return $self->add($tasks, $previous_outputs);
}

# Chain task(s) to all previous outputs
sub add_to_all {
    my ($self, $tasks) = @_;
    return $self->add($tasks, $self->outputs);
}

sub previous_output {
    my ($self, $index) = @_;
 
    $index = 0 unless defined $index;
    if ($previous_outputs && @$previous_outputs >= $index) {
        return $previous_outputs->[$index];    
    }
    
    return;
}

sub previous_outputs {
    my $self = shift;
    if ($previous_outputs && @$previous_outputs) {
        return wantarray ? @$previous_outputs : $previous_outputs;
    }
    return;
}

sub outputs {
    my $self = shift;
    my @outputs;
    foreach (@{$self->tasks}) {
        push @outputs, @{$_->{outputs}} if $_->{outputs};
    }
    return \@outputs;
}

###############################################################################
# Common Tasks
###############################################################################

# Generate synchronization dependency for chaining pipelines
sub wait {
    my $self = shift;
    my $dependencies = shift // [];

    my $wait_file = catfile($self->staging_dir, 'wait_' . get_unique_id() . '.done');

    return {
        cmd     => "touch $wait_file",
        args    => [],
        inputs  => [ @$dependencies ],
        outputs => [ $wait_file ],
        description => 'Waiting for tasks to complete'
    };
}

# Generate GFF file of genome annotations
sub create_gff {
    my $self = shift;
    my %params = @_;
    
    # Build argument list
    my $args = [
        ['-gid',     $params{gid},                0],
        ['-f',       $params{output_file},        1],
        ['-config',  $self->conf->{_CONFIG_PATH}, 0],
        ['-cds',     $params{cds} // 0,           0],
        ['-annos',   $params{annos} // 0,         0],
        ['-nu',      $params{nu} // 0,            0],
        ['-id_type', $params{id_type} // 0,       0],
        ['-upa',     $params{upa} // 0,           0]
    ];
    push @$args, ['-chr',     $params{chr},     0] if (defined $params{chr});
    push @$args, ['-add_chr', $params{add_chr}, 0] if (defined $params{add_chr});
    
    return {
        cmd         => catfile($self->conf->{SCRIPTDIR}, "coge_gff.pl"),
        args        => $args,
        outputs     => [ $params{output_file} ],
        description => "Generating GFF"
    };
}

# Generate BED file of genome annotations
sub create_bed {
    my ($self, %params) = @_;

    return {
        cmd  => catfile($self->conf->{SCRIPTDIR}, "coge2bed.pl"), #FIXME script name makes no sense
        args => [
            ['-gid',    $params{gid},                0],
            ['-f',      $params{output_file},        1],
            ['-config', $self->conf->{_CONFIG_PATH}, 0]
        ],
        outputs => [ $params{output_file} ]
    };
}

# Generate TBL file of genome annotations
sub create_tbl {
    my ($self, %params) = @_;

    return {
        cmd     => catfile($self->conf->{SCRIPTDIR}, "export_NCBI_TBL.pl"),
        args    => [
            ['-gid',    $params{gid},                0],
            ['-f',      $params{output_file},        1],
            ["-config", $self->conf->{_CONFIG_PATH}, 0]
        ],
        outputs => [ $params{output_file} ]
    };
}

sub export_to_irods {
    my ($self, %params) = @_;

    irods_set_env(catfile($self->conf->{_HOME_PATH}, 'irodsEnv')); # mdb added 2/9/16 -- for hypnotoad, use www-data's irodsEnvFile
    my $cmd = irods_iput($params{src_file}, $params{dest_file}, { no_execute => 1, overwrite => ($params{overwrite} // 0) });

    my $done_file = catfile($self->staging_dir, basename($params{src_file})) . '.iput.done';

    return {
        cmd => qq[$cmd && touch $done_file],
        description => "Exporting file to IRODS " . $params{dest_file},
        args => [],
        inputs => [ $params{src_file} ],
        outputs => [ $done_file ]
    };
}

sub irods_imeta {
    my ($self, %params) = @_;
    
    my $done_file = catfile($self->staging_dir, basename($params{dest_file})) . '.imeta.done';
    
    my $cmd = catdir($self->conf->{SCRIPTDIR}, 'irods.pl') .
        " -cmd metadata" .
        " -dest " . $params{dest_file} . 
        " -metafile " . $params{metadata_file} .
        " -env " . catfile($self->conf->{_HOME_PATH}, 'irodsEnv') .
        " && touch $done_file";
    
    return {
        cmd => $cmd,
        description => "Setting IRODS metadata for " . $params{dest_file},
        args => [],
        inputs => [],
        outputs => [ $done_file ]
    };
}

sub join_files {
    my ($self, %params) = @_;
    my $input_files = $params{input_files};
    my $input_dir   = $params{input_dir};
    my $output_file = $params{output_file};
    
    my @files;
    push @files, @$input_files   if $input_files;
    push @files, $input_dir.'/*' if $input_dir;
    
    my $cmd = "mkdir -p \$(dirname $output_file) && cat " . join(' ', shell_quote(@files)) . ' > ' . $output_file;
    
    return {
        cmd => $cmd,
        args => [],
        inputs => [
            @$input_files,
#            @$done_files
        ],
        outputs => [
            $output_file
        ],
        description => 'Joining files'
    };
}

sub send_email {
    my ($self, %params) = @_;
    my $from    = 'CoGe Support <coge.genome@gmail.com>';
    my $to      = $params{to};
    my $subject = $params{subject};
    my $body    = $params{body};
    
    my $done_file = catfile($self->staging_dir, "send_email.done");
    
    my $args = [
        ['-from',      '"'.escape($from).'"',    0],
        ['-to',        '"'.escape($to).'"',      0],
        ['-subject',   '"'.escape($subject).'"', 0],
        ['-body',      '"'.escape($body).'"',    0],
    ];

    return {
        cmd => catfile($self->conf->{SCRIPTDIR}, "send_email.pl"),
        args => $args,
        inputs => [],
        outputs => [],
        description => "Sending notification email"
    };
}

sub curl_get {
    my ($self, %params) = @_;
    my $url = $params{url};

    my $cmd = get_command_path('CURL');
    my $output_file = catfile($self->staging_dir, 'curl_output.log');

    return {
        cmd => "$cmd -s -o $output_file $url",
        args => [],
        inputs => [],
        outputs => [ $output_file ],
        description => "Sending GET request to $url"
    };
}

sub reheader_fasta {
    my $self = shift;
    my $gid = shift;

    my $fasta     = get_genome_file($gid);
    my $cache_dir = get_genome_cache_path($gid);

    my $output_file = to_filename($fasta) . '.reheader.faa';

    return {
        cmd => catfile($self->conf->{SCRIPTDIR}, "fasta_reheader.pl"),
        args => [
            ["", $fasta, 1],
            ["", $output_file, 0]
        ],
        inputs => [
            $fasta
        ],
        outputs => [
            catfile($cache_dir, $output_file)
        ],
        description => "Reheader fasta file",
    };
}

sub index_fasta {
    my $self = shift;
    my $fasta = shift;

    return {
        cmd => get_command_path('SAMTOOLS'),
        args => [
            ["faidx", $fasta, 1],
        ],
        inputs => [
            $fasta,
        ],
        outputs => [
            $fasta . '.fai',
        ],
        description => "Indexing FASTA file",
    };
}

sub sam_to_bam {
    my $self = shift;
    my $samfile = shift;

    my $filename = to_filename($samfile);
    my $cmd = get_command_path('SAMTOOLS');

    return {
        cmd => $cmd,
        args => [
            ["view", '', 0],
            ["-bS", $samfile, 1],
            [">", $filename . ".bam", 0]
        ],
        inputs => [
            $samfile
        ],
        outputs => [
            catfile($self->staging_dir, $filename . ".bam")
        ],
        description => "Converting SAM to BAM"
    };
}

sub index_bam {
    my $self = shift;
    my $bamfile = shift;

    return {
        cmd => get_command_path('SAMTOOLS'),
        args => [
            ["index", $bamfile, 1],
        ],
        inputs => [
            $bamfile,
        ],
        outputs => [
            $bamfile . '.bai'
        ],
        description => "Indexing BAM file",
    };
}

sub add_result {
    my ($self, %params) = @_;
    my $username = $self->user->name;
    my $wid      = $self->workflow->id;
    my $result   = encode_json($params{result});
    
#    my $result_file = get_workflow_results_file($username, $wid);

    return {
        cmd  => catfile($self->conf->{SCRIPTDIR}, "add_workflow_result.pl"),
        args => [
            ['-user_name', $username,            0],
            ['-wid',       $wid,                 0],
            ['-result',    shell_quote($result), 0] #TODO pass via temp file instead
        ],
        inputs  => [],
        outputs => [
#            $result_file,  # force this to run (for case of multiple results)
        ],
        description => "Adding workflow result"
    };
}

sub add_items_to_notebook {
    my ($self, %params) = @_;
    my $notebook_id = $params{notebook_id};

    my $result_file = get_workflow_results_file($self->user->name, $self->workflow->id);

    my $log_file = catfile($self->staging_dir, 'add_items_to_notebook', 'log.txt');

    return {
        cmd => catfile($self->conf->{SCRIPTDIR}, 'add_items_to_notebook.pl'),
        args => [
            ['-uid', $self->user->id, 0],
            ['-wid', $self->workflow->id, 0],
            ['-notebook_id', $notebook_id, 0],
            ['-config', $self->conf->{_CONFIG_PATH}, 0],
            ['-log', $log_file, 0]
        ],
        inputs => [],
        outputs => [
            $result_file,
            $log_file
        ],
        description => "Adding experiment to notebook"
    };
}

sub create_notebook {
    my ($self, %params) = @_;
    my $metadata = $params{metadata};
    my $annotations = $params{annotations}; # array ref

    my $result_file = get_workflow_results_file($self->user->name, $self->workflow->id);

    my $log_file = catfile($self->staging_dir, 'create_notebook', 'log.txt');

    my $annotations_str = '';
    $annotations_str = join(';', @$annotations) if (defined $annotations && @$annotations);

    return {
        cmd => catfile($self->conf->{SCRIPTDIR}, 'create_notebook.pl'),
        args => [
            ['-uid', $self->user->id, 0],
            ['-wid', $self->workflow->id, 0],
            ['-name', shell_quote($metadata->{name}), 0],
            ['-desc', shell_quote($metadata->{description}), 0],
            ['-type', 2, 0],
            ['-restricted', $metadata->{restricted}, 0],
            ['-annotations', qq{"$annotations_str"}, 0],
            ['-config', $self->conf->{_CONFIG_PATH}, 0],
            ['-log', $log_file, 0]
        ],
        inputs => [],
        outputs => [
            $result_file,
            $log_file
        ],
        description => "Creating notebook of results"
    };
}

sub picard_deduplicate {
    my $self = shift;
    my $bam_file = shift;

    my $cmd = 'java -jar ' . $PICARD;

    my $output_file = $bam_file . '-deduplicated.bam';

    return {
        cmd => "$cmd MarkDuplicates REMOVE_DUPLICATES=true INPUT=$bam_file METRICS_FILE=$bam_file.metrics OUTPUT=$output_file.tmp ; mv $output_file.tmp $output_file",
        args => [
#            ['MarkDuplicates', '', 0],
#            ['REMOVE_DUPLICATES=true', '', 0],
#            ["INPUT=$bam_file", '', 0],
#            ["METRICS_FILE=$bam_file.metrics", '', 0],
#            ["OUTPUT=$output_file", '', 0],
        ],
        inputs => [
            $bam_file
        ],
        outputs => [
            $output_file
        ],
        description => "Deduplicating PCR artifacts using Picard"
    };
}

sub sort_bam {
    my $self = shift;
    my $bam_file = shift;

    my $filename = to_filename($bam_file);
    my $cmd = get_command_path('SAMTOOLS');

    return {
        cmd => $cmd,
        args => [
            ["sort", '', 0],
            ["", $bam_file, 1],
            ["-o", $filename . "-sorted.bam", 1] # mdb changed 1/5/17 -- added -o for SAMtools 1.3.1
        ],
        inputs => [
            $bam_file
        ],
        outputs => [
            catfile($self->staging_dir, $filename . "-sorted.bam")
        ],
        description => "Sorting BAM file"
    };
}

sub load_bam { #TODO combine with load_experiment
    my $self = shift;
    my %opts = @_;
    my $annotations = $opts{annotations};
    my $bam_file = $opts{bam_file};
    my $metadata = $opts{metadata} || $self->params->{metadata};
    my $additional_metadata = $self->params->{additional_metadata};

    my $cmd = 'perl ' . catfile($self->conf->{SCRIPTDIR}, "load_experiment.pl");

    my $output_name = "load_bam_" . to_filename_base($bam_file);
    my $output_path = catdir($self->staging_dir, $output_name);

    my $result_file = get_workflow_results_file($self->user->name, $self->workflow->id);

    # Add tags
    my @tags = ( 'BAM' ); # add BAM tag
    push @tags, @{$metadata->{tags}} if $metadata->{tags};
    my $tags_str = tags_to_string(\@tags);

    my $args = [
        ['-user_name',   $self->user->name, 0],
        ['-name',        ($metadata->{name} ? shell_quote($metadata->{name} . " (BAM alignment)") : '""'), 0],
        ['-desc',        shell_quote($metadata->{description}), 0],
        ['-version',     shell_quote($metadata->{version}), 0],
        ['-link',        shell_quote($metadata->{link}), 0],
        ['-restricted',  shell_quote($metadata->{restricted}), 0],
        ['-gid',         $self->request->genome->id, 0],
        ['-wid',         $self->workflow->id, 0],
        ['-source_name', shell_quote($metadata->{source_name}), 0],
        ['-tags',        shell_quote($tags_str), 0],
        ['-staging_dir', $output_name, 0],
        ['-file_type',   'bam', 0],
        ['-data_file',   $bam_file, 0],
        ['-config',      $self->conf->{_CONFIG_PATH}, 0]
    ];

    # Add additional metadata
    if ($additional_metadata && @$additional_metadata) { # new method using metadata file
        my $metadata_file = catfile($output_path, 'metadata.dump');
        make_path($output_path);
        CoGe::Accessory::TDS::write($metadata_file, $additional_metadata);
        push @$args, ['-metadata_file', $metadata_file, 0];
    }
    if ($annotations && @$annotations) { # legacy method
        my $annotations_str = join(';', @$annotations);
        push @$args, ['-annotations', shell_quote($annotations_str), 0] if ($annotations_str);
    }

    return {
        cmd => $cmd,
        args => $args,
        inputs => [
            $bam_file
        ],
        outputs => [
            [$output_path, '1'],
            catfile($output_path, "log.done"),
            $result_file
        ],
        description => "Loading alignment as new experiment"
    };
}

sub load_experiment {
    my $self = shift;
    my %opts = @_;
    my $metadata = $opts{metadata} || $self->params->{metadata};
    my $additional_metadata = $self->params->{additional_metadata};
    my $annotations = $opts{annotations};
    my $gid = $opts{gid};
    my $input_file = $opts{input_file};
    my $normalize = $opts{normalize} || 0;
    my $name = $opts{name} || to_filename_base($input_file); # optional name for this load

    my $output_name = "load_experiment_$name";
    my $output_path = catdir($self->staging_dir, $output_name);


    my $args = [
        ['-gid',         $gid,                0],
        ['-wid',         $self->workflow->id, 0],
        ['-user_name',   $self->user->name,   0],
        ['-name',        ($metadata->{name}        ? shell_quote($metadata->{name}) : '""'),        0],
        ['-desc',        ($metadata->{description} ? shell_quote($metadata->{description}) : '""'), 0],
        ['-version',     ($metadata->{version}     ? shell_quote($metadata->{version}) : '""'),     0],
        ['-link',        ($metadata->{link}        ? shell_quote($metadata->{link}) : '""'),        0],
        ['-restricted',  shell_quote($metadata->{restricted}), 0],
        ['-source_name', ($metadata->{source_name} ? shell_quote($metadata->{source_name}) : '""'), 0],
        ['-staging_dir', $output_name, 0],
        ['-data_file',   $input_file,  0],
        ['-normalize',   $normalize,   0],
        ['-disable_range_check', '',   0], # mdb added 8/26/16 COGE-270 - allow values outside of [-1, 1]
        ['-config',      $self->conf->{_CONFIG_PATH}, 0]
    ];

    # Add tags
    if ($metadata->{tags}) {
        my $tags_str = tags_to_string($metadata->{tags});
        push @$args, ['-tags', shell_quote($tags_str), 0];
    }

    # Add additional metadata
    if ($additional_metadata && @$additional_metadata) { # new method using metadata file
        my $metadata_file = catfile($output_path, 'metadata.dump');
        make_path($output_path); #TODO maybe this file should be located somewhere else
        open(my $fh, ">$metadata_file");
        print $fh Dumper $additional_metadata;
        close($fh);
        push @$args, ['-metadata_file', $metadata_file, 0];
    }
    if ($annotations && @$annotations) { # legacy method
        my $annotations_str = join(';', @$annotations);
        push @$args, ['-annotations', shell_quote($annotations_str), 0] if ($annotations_str);
    }

    return {
        cmd => catfile($self->conf->{SCRIPTDIR}, "load_experiment.pl"),
        args => $args,
        inputs => [
            $input_file
        ],
        outputs => [
            [$output_path, '1'],
            catfile($output_path, "log.done")
        ],
        description => "Loading" . ($name ? " $name" : '') . " experiment"
    };
}

__PACKAGE__->meta->make_immutable;

1;
