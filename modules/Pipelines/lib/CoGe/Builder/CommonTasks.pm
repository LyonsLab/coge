package CoGe::Builder::CommonTasks;

#######################################################################################
# LEGACY MODULE: continue to move functions out and eventually remove this module
#######################################################################################

use strict;
use warnings;

use File::Spec::Functions qw(catdir catfile);
use File::Basename qw(basename dirname);
use JSON qw(encode_json);
use Data::Dumper;

use CoGe::Accessory::Utils;
use CoGe::Accessory::Web qw(get_defaults);
use CoGe::Accessory::TDS;
use CoGe::Core::Storage;

require Exporter;
our @ISA = qw(Exporter);
our @EXPORT = qw(
    generate_gff generate_features copy_and_mask create_gff_generation_job add_workflow_result
    create_transdecoder_longorfs_job create_transdecoder_predict_job
);

our $CONF = CoGe::Accessory::Web::get_defaults();

sub add_workflow_result {
    my %opts = @_;
    my $username = $opts{username};
    my $wid = $opts{wid};
    my $result = encode_json($opts{result});
    my $dependency = $opts{dependency};
    
    my $result_file = get_workflow_results_file($username, $wid);

    return {
        cmd  => catfile($CONF->{SCRIPTDIR}, "add_workflow_result.pl"),
        args => [
            ['-user_name', $username, 0],
            ['-wid', $wid, 0],
            ['-result', "'".$result."'", 0]
        ],
        inputs  => [$dependency],
        outputs => [$result_file],
        description => "Adding workflow result"
    };
}

sub copy_and_mask {
    my %args = @_;

    my $desc = $args{mask} ? "Copying and masking genome" : "Copying genome";
    $desc .= " (no annotations)" if $args{seq_only};

    my $cmd = "/copy_genome/copy_load_mask_genome.pl";

    return (
        cmd   => catfile($CONF->{SCRIPTDIR}, $cmd),
        args  => [
            ["-conf", $CONF->{_CONFIG_PATH}, 0],
            ["-gid", $args{gid}, 0],
            ["-uid", $args{uid}, 0],
            ["-wid", $args{wid}, 0],
            ["-mask", $args{mask}, 0],
            ["-staging_dir", $args{staging_dir}, 0],
            ["-sequence_only", $args{seq_only}, 0]
        ],
        description => $desc
    );
}

sub generate_features {
    my %args = @_;

    my $filename = $args{basename} . "-gid-" . $args{gid};
    $filename .= "-prot" if $args{protein};
    $filename .= ".fasta";
    my $path = get_genome_cache_path($args{gid});
    my $output_file = catfile($path, $filename);

    return $output_file, (
        cmd    => catfile($CONF->{SCRIPTDIR}, "export_features_by_type.pl"),
        args   => [
            ["-config", $CONF->{_CONFIG_PATH}, 0],
            ["-f", $filename, 0],
            ["-gid", $args{gid}, 0],
            ["-ftid", $args{fid}, 0],
            ["-prot", $args{protein}, 0],
        ],
        outputs => [$output_file]
    );
}

#FIXME dup'ed below -- mdb 11/7/16
sub create_gff_generation_job {
    my %opts = @_;

    # Required params
    my $gid = $opts{gid};
    my $organism_name = $opts{organism_name};
    #my $validated = $opts{validated}; # mdb removed 1/5/15 -- why is this here?

    my $cmd = catfile($CONF->{SCRIPTDIR}, "coge_gff.pl");
    my $name = sanitize_name($organism_name) . "-1-name-0-0-id-" . $gid . "-1.gff";

    return {
        cmd => $cmd,
        script => undef,
        args => [
            ['-f', $name, 1],
            ['-staging_dir', '.', 0],
            ['-gid', $gid, 0],
            ['-upa', 1, 0],
            ['-id_type', "name", 0],
            ['-cds', 0, 0],
            ['-annos', 0, 0],
            ['-nu', 1, 0],
            ['-config', $CONF->{_CONFIG_PATH}, 0],
        ],
        inputs => [],
        outputs => [
            catfile(get_genome_cache_path($gid), $name)
        ],
        description => "Generating genome annotations GFF file"
    };
}

#FIXME dup'ed above -- mdb 11/7/16
sub generate_gff {
    my %inputs = @_;
    my %args = (
        annos   => 0,
        id_type => 0,
        cds     => 0,
        nu      => 0,
        upa     => 0,
        add_chr => 0
    );
    @args{(keys %inputs)} = (values %inputs);

    # Check for a genome or dataset id
    return unless $args{gid};

    # Set the default basename as the id if the basename is not set
    $args{basename} = $args{gid} unless $args{basename};

    # Generate the output filename
    my $param_string = join( "-", map { $_ . $args{$_} } qw(annos cds id_type nu upa add_chr) );
    my $filename = $args{basename} . "_" . $param_string;
    if ($args{chr}) {
    	$filename .= "_" . $args{chr};
    }
    $filename .= ".gid" . $args{gid};
    $filename .= ".gff";
    $filename =~ s/\s+/_/g;
    $filename =~ s/\)|\(/_/g;
    my $path = get_genome_cache_path($args{gid});
    my $output_file = catfile($path, $filename);
    
    # Build argument list
    my $args = [
        ['-gid', $args{gid}, 0],
        ['-f', $filename, 1],
        ['-config', $CONF->{_CONFIG_PATH}, 0],
        ['-cds', $args{cds}, 0],
        ['-annos', $args{annos}, 0],
        ['-nu', $args{nu}, 0],
        ['-id_type', $args{id_type}, 0],
        ['-upa', $args{upa}, 0]
    ];
    push @$args, ['-chr', $args{chr}, 0] if (defined $args{chr});
    push @$args, ['-add_chr', $args{add_chr}, 0] if (defined $args{add_chr});
    
    return $output_file, {
        cmd     => catfile($CONF->{SCRIPTDIR}, "coge_gff.pl"),
        script  => undef,
        args    => $args,
        inputs  => [],
        outputs => [ $output_file ],
        description => "Generating GFF"
    };
}

sub create_transdecoder_longorfs_job {
    my $input_file = shift;

    my $cmd = catfile($CONF->{TRANSDECODER}, 'TransDecoder.LongOrfs');
    
    my $done_file = $input_file . '.longorfs.done';
    my $output_path = $input_file . '.transdecoder_dir';

    return {
        cmd => "$cmd -t $input_file && touch $done_file",
        script => undef,
        args => [],
        inputs => [
            $input_file,
        ],
        outputs => [
            [$output_path, '1'],
            $done_file
        ],
        description => "Running TransDecoder.LongOrfs"
    };
}

sub create_transdecoder_predict_job {
    my $input_file = shift;
    my $input_dir = shift; 
    my $dependency_file = shift;

    my $cmd = catfile($CONF->{TRANSDECODER}, 'TransDecoder.Predict');
    
    my $done_file = $input_file . '.predict.done';
    my $output_file = $input_file . '.transdecoder.gff3';

    return {
        cmd => "cd $input_dir && $cmd -t $input_file && touch $done_file", # transdecoder won't work unless run from the dir that contains the input dir
        script => undef,
        args => [],
        inputs => [
            $input_file,
            $dependency_file
        ],
        outputs => [
            $output_file,
            $done_file
        ],
        description => "Running TransDecoder.Predict"
    };
}


1;
