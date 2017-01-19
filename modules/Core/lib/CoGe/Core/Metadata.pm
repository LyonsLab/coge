package CoGe::Core::Metadata;

use v5.14;
use strict;
use warnings;
use Switch;
use File::Basename;
# use File::Slurp;
use File::Spec::Functions qw(catfile);
use Data::Dumper;
use Text::Unidecode qw(unidecode);

use CoGeX;
use CoGe::Accessory::IRODS qw(irods_get_base_path irods_imeta_ls irods_iput);

BEGIN {
    our (@ISA, $VERSION, @EXPORT);
    require Exporter;

    $VERSION = 0.0.1;
    @ISA = qw(Exporter);
    @EXPORT = qw( create_annotation create_annotations export_annotations tags_to_string to_annotations update_annotation );
}

sub create_annotations {
    my %opts = @_;
    my $db = $opts{db};
    my $target = $opts{target};           # experiment, genome, or list object
    my ($target_id, $target_type) = ($opts{target_id}, $opts{target_type}); # or an experiment/genome/list id and type
    my $annotations = $opts{annotations}; # semicolon-separated list of annotations (image_file|link|group|type|text;[...])
    my $anno_file = $opts{anno_file};     # metadata file
    my $locked = $opts{locked};           # boolean flag to indicate locked (not editable) annotations

    my @result = ();

    if ($annotations) { # legacy method
        foreach (split(/\s*;\s*/, $annotations)) {
            my @tok = split(/\s*\|\s*/, $_);
            my ($image_file, $link, $group_name, $type_name, $anno_text);
            $image_file = shift @tok if (@tok == 5);
            $link = shift @tok if (@tok == 4);
            $group_name = shift @tok if (@tok == 3);
            $type_name = shift @tok if (@tok == 2);
            $anno_text = shift @tok if (@tok == 1);
            unless ($anno_text and $type_name) {
                print STDERR "CoGe::Core::Metadata: missing required annotation type and text fields\n";
                print STDERR Dumper [ $image_file, $link, $group_name, $type_name, $anno_text ], "\n";
                return;
            }

            my $anno = create_annotation(
                db          => $db,
                target      => $target,
                target_id   => $target_id,
                target_type => $target_type,
                group_name  => $group_name,
                type_name   => $type_name,
                text        => $anno_text,
                link        => $link,
                image_file  => $image_file,
                locked      => $locked
            );
            unless ($anno) {
                print STDERR "CoGe::Core::Metadata: error creating annotation\n";
                return;
            }

            push @result, $anno;
        }
    }

    if ($anno_file) { # new method
        my $annos = CoGe::Accessory::TDS::read($anno_file);
        if ($annos) {
            foreach (@$annos) {
                my $anno = create_annotation(
                    db          => $db,
                    target      => $target,
                    target_id   => $target_id,
                    target_type => $target_type,
                    group_name  => $_->{group},
                    type_name   => $_->{type},
                    text        => $_->{text},
                    link        => $_->{link},
                    locked      => $locked
                );
                unless ($anno) {
                    print STDERR "CoGe::Core::Metadata: error creating annotation\n";
                    return;
                }

                push @result, $anno;
            }
        }
    }

    return \@result;
}

sub _init {
    my $opts = shift;
    my $type_name = $opts->{type_name};
    return 'type_name required' unless $type_name;

    my $db = $opts->{db};
    my $group_name = $opts->{group_name};
    my $group_id;
    $group_id = $db->resultset('AnnotationTypeGroup')->find_or_create({ name => $group_name })->id if defined $group_name;

    my $type_id = $db->resultset('AnnotationType')->find_or_create({ name => $type_name, annotation_type_group_id => $group_id })->id;
    return 'error creating annotation type' unless $type_id;

    my $link = $opts->{link};
    if ($link) {
        $link =~ s/^\s+//;
        $link = 'http://' . $link if !($link =~ /^(\w+)\:\/\//);
    }

    my $bisque_id;
    my $image_fh = $opts->{image_fh};
    my $image_file = $opts->{image_file};
    if ($image_file) {
        $bisque_id = create_image(fh => $image_fh, filename => $image_file, tmp_filename => $opts->{image_tmp_file}, db => $db, user => $opts->{user});
        return 'error creating image' unless $image;
    }

    return (undef, $type_id, $link, $bisque_id);
}

sub create_annotation {
    my %opts = @_;

    my ($error, $type_id, $link, $bisque_id) = _init(\%opts);
    if ($error) {
        warn 'create_annotation: ' . $error;
        return;
    }

    my $db          = $opts{db};
    my $target      = $opts{target};    # DBIX experiment/genome/notebook
    my $target_id   = $opts{target_id}; # or id/type
    my $target_type = $opts{target_type};
    if ($target_id && $target_type) {
        switch (lc($target_type)) {
            case 'experiment' { $target = $db->resultset('Experiment')->find($target_id); }
            case 'genome' { $target = $db->resultset('Genome')->find($target_id); }
            case 'notebook' { $target = $db->resultset('List')->find($target_id); }
        }
    }
    unless ($target) {
        warn 'create_annotation: target not found';
        return;
    }

    my $locked = $opts{locked};
    my $text   = $opts{text};
    if (ref($target) =~ /Experiment/) {
        return $db->resultset('ExperimentAnnotation')->find_or_create({
            experiment_id => $target->id,
            annotation_type_id => $type_id,
            annotation => $text,
            link => $link,
            bisque_id => $bisque_id,
            locked => $locked
        });
    }
    if (ref($target) =~ /Genome/) {
        return $db->resultset('GenomeAnnotation')->find_or_create({
            genome_id => $target->id,
            annotation_type_id => $type_id,
            annotation => $text,
            link => $link,
            bisque_id => $bisque_id,
            locked => $locked
        });
    }
    if (ref($target) =~ /List/) {
        return $db->resultset('ListAnnotation')->find_or_create({
            list_id => $target->id,
            annotation_type_id => $type_id,
            annotation => $text,
            link => $link,
            bisque_id => $bisque_id,
            locked => $locked
        });
    }
    warn 'create_annotation: unknown target type';

    return;
}

sub update_annotation {
    my %opts = @_;

    my $annotation_id = $opts{annotation_id};
    unless ($annotation_id) {
        warn 'update_annotation: annotation_id required';
        return;
    }
 
    my ($error, $type_id, $link, $bisque_id) = _init(\%opts);
    if ($error) {
        warn 'update_annotation: ' . $error;
        return;
    }

    my $db = $opts{db};
    my $annotation;
    switch (lc($opts{target_type})) {
        case 'experiment' { $annotation = $db->resultset('ExperimentAnnotation')->find($annotation_id); }
        case 'genome' { $annotation = $db->resultset('GenomeAnnotation')->find($annotation_id); }
        case 'notebook' { $annotation = $db->resultset('ListAnnotation')->find($annotation_id); }
    }
    $annotation->annotation($opts{text});
    $annotation->link($link);
    $annotation->annotation_type_id($type_id);
    $annotation->bisque_id($bisque_id) if ($bisque_id);
    $annotation->update;

    return;
}

sub export_annotations {
    my %opts = @_;
    my $annotations = $opts{annotations}; # array ref of DBIX Annotation objects
    my $export_path = $opts{export_path}; # export directory
    
    return () unless (defined($annotations) and @$annotations);

    my @files = ();
    my $annotation_file = File::Spec->catdir($export_path, "annotations.csv");
    push @files, basename($annotation_file);

    unless (-r $annotation_file) {
        open(my $fh, ">", $annotation_file);

        say $fh "#Type Group, Type, Annotation, Link, Image filename, BisQue URL";
        foreach my $a ( @$annotations ) {
            my $group = (
                defined $a->type->group
                ? '"' . $a->type->group->name . '","' . $a->type->name . '"'
                : '"' . $a->type->name . '",""'
            );

            my $info = $a->info;
            my $url = defined($a->link) ? $a->link : "";

            # Escape quotes
            $info =~ s/"/\"/g;

            my $filename;
            if ($a->image) {
                $filename = $a->image->filename;
                my $img =  File::Spec->catdir($export_path, $filename);

                eval {
                    open(my $imh, ">", $img) or die "image=$filename could not be generated";
                    print $imh $a->image->image;
                    close($imh);

                    push @files, $filename;
                };

                say $fh "log: error: $@" if ($@);
            }
            my $bisque_id = defined($a->bisque_id) ? 'http://bisque.iplantcollaborative.org/data_service/' . $a->bisque_id : '';
            say $fh qq{$group,$info,$url,$filename,$bisque_id};
        }
        close($fh);
    }

    return @files;
}

# Convert metadata hash to a string that can be passed to load scripts (eg, load_experiment.pl).
# The output of this function is eventually parsed by create_annotations() above.
sub to_annotations {
    my $md_array = shift; # array of metadata hash refs
    # [
    #    {
    #       type_group: "string",
    #       type: "string",
    #       text: "string",
    #       link: "string"
    #    }
    # ]
    # will be converted to
    # link|type_group|type|text;link|type_group|type|text;...
    
    my @annotations;
    foreach my $md (@$md_array) {
        next unless ($md->{type} && $md->{text}); # minimum required fields
        
        my @fields = map { unidecode($_) } # mdb added unidecode 12/9/16 COGE-693 - handle UTF-8 unicode chars from NCBI-SRA
                     map { $md->{$_} || '' } ('link', 'type_group', 'type', 'text');

        shift @fields unless $fields[0]; # remove undefined link
        shift @fields unless $fields[0]; # remove undefined type_group
        
        push @annotations, join('|', @fields);
    }
    
    return wantarray ? @annotations : \@annotations;
}

sub tags_to_string {
    my $tags_array = shift;
    
    my $tags_str = '';
    if ($tags_array && @$tags_array) {
        my %tags = map { $_ => 1 } (@$tags_array);
        $tags_str = join(';', map { $_ =~ s/\;//g; $_ } sort keys %tags);
    }
    
    return $tags_str;
}

sub create_image {
    my %opts = @_;
    my $fh = $opts{fh};
    my $filename = $opts{filename};
    my $db = $opts{db};
    return unless (($fh || $filename) && $db);

    my $dest = catfile(dirname(irods_get_base_path($opts{user}->name)), 'bisque_data', basename($filename));
    irods_iput($opts{tmp_filename}, $dest);
    for my $i (0..9) {
        sleep 5;
        my $result = irods_imeta_ls($dest, 'ipc-bisque-id');
        if (@$result == 4 && substr($result->[2], 0, 6) eq) 'value:') {
            return substr($result->[2], 6);
        }
    }
    warn 'unable to get bisque id';
}

1;
