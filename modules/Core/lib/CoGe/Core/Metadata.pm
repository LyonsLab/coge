package CoGe::Core::Metadata;

use v5.14;
use strict;
use warnings;

use Data::Dumper;
use File::Basename qw(basename);
use LWP::UserAgent;
use Mojo::JSON;
use Switch;
use Text::Unidecode qw(unidecode);

use CoGe::Accessory::BisQue qw(create_bisque_image delete_bisque_image get_bisque_data_url);
use CoGe::Accessory::TDS;
use CoGeX;

BEGIN {
    our (@ISA, $VERSION, @EXPORT);
    require Exporter;

    $VERSION = 0.0.1;
    @ISA = qw(Exporter);
    @EXPORT = qw(
        create_annotation create_annotations delete_annotation export_annotations get_annotation get_annotations
        get_type_groups search_annotation_types tags_to_string to_annotations update_annotation
    );
}

sub create_annotation {
    my %opts = @_;
    my ($error, $type_id, $link, $image_id, $bisque_id, $bisque_file) = _init(\%opts);
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

    my $locked = $opts{locked} ? 1 : 0;
    my $text   = $opts{text};
    my $args = {
        annotation_type_id => $type_id,
        annotation => $text,
        link => $link,
        image_id => $image_id,
        locked => $locked
    };
    if ($bisque_file) {
        $args->{bisque_file} = $bisque_file;
        $args->{bisque_id} = $bisque_id;
        $args->{bisque_user} = $opts{user}->id;
    }
    if (ref($target) =~ /Experiment/) {
        $args->{experiment_id} = $target->id;
        return $db->resultset('ExperimentAnnotation')->find_or_create($args);
    }
    if (ref($target) =~ /Genome/) {
        $args->{genome_id} = $target->id;
        return $db->resultset('GenomeAnnotation')->find_or_create($args);
    }
    if (ref($target) =~ /List/) {
        $args->{list_id} = $target->id;
        return $db->resultset('ListAnnotation')->find_or_create($args);
    }
    warn 'create_annotation: unknown target type';
}

sub create_annotations {
    my %opts = @_;
    my $db = $opts{db};
    my $user = $opts{user};
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
                user        => $user,
                db          => $db,
                target      => $target,
                target_id   => $target_id,
                target_type => $target_type,
                group_name  => $group_name,
                type_name   => $type_name,
                text        => $anno_text,
                link        => $link,
                image_file  => $image_file,
                locked      => $locked,
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
                    user        => $user,
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

sub delete_annotation {
    my ($aid, $object_id, $object_type, $db, $user) = @_;
    return 'delete_annotation: missing parameter' unless $aid && $object_id && $object_type && $db;

    my ($error, $object) = _get_object(undef, $object_id, $object_type, 1, $db, $user);
    return $error if $error;

    my $annotation = $db->resultset($object_type . 'Annotation')->find( { lc($object_type) . '_annotation_id' => $aid } );
    $object_type = lc($object_type);
    $object_type = 'notebook' if $object_type eq 'list';
    delete_bisque_image($object_type, $object_id, $annotation->bisque_file, $annotation->bisque_id, $user) if $annotation->bisque_file;
    $annotation->delete();
    return undef;
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
            my $bisque_link = defined($a->bisque_id) ? get_bisque_data_url($a->bisque_id) : '';
            say $fh qq{$group,$info,$url,$filename,$bisque_link};
        }
        close($fh);
    }

    return @files;
}

sub get_annotation {
    my ($aid, $object_type, $db) = @_;
    return unless $aid && $object_type && $db;

    my $annotation = $db->resultset($object_type . 'Annotation')->find($aid);
    return unless $annotation;

    my $type       = '';
    my $type_group = '';
    if ( $annotation->type ) {
        $type = $annotation->type->name;
        $type_group = $annotation->type->group->name if ( $annotation->type->group );
    }
    return {
        annotation  => $annotation->annotation,
        link        => $annotation->link,
        type        => $type,
        type_group  => $type_group,
        bisque_file => $annotation->bisque_file,
        bisque_id   => $annotation->bisque_id
    };
}

sub get_annotations {
    my ($id, $object_type, $db, $for_api) = @_;
    return unless $id && $object_type && $db;

    my $object = $db->resultset($object_type)->find($id);
    return unless $object;

    # Categorize annotations based on type group and type
    my %groups;
    foreach my $a ( $object->annotations ) {
        my $group = ( $a->type->group ? $a->type->group->name : '');
        my $type = $a->type->name;
        push @{ $groups{$group}{$type} }, {
            annotation => $a->info,
            bisque_id => $a->bisque_id,
            id => $a->id,
            image_id => $a->image_id,
            link => $a->link,
            locked => $for_api ? ($a->locked ? Mojo::JSON->true : Mojo::JSON->false) : $a->locked
        };
    }
    return \%groups;
}

sub get_type_groups {
    my $db = shift;
    my %unique;

    my $rs = $db->resultset('AnnotationTypeGroup');
    while (my $atg = $rs->next) {
        $unique{ $atg->name }++;
    }
    return [ sort keys %unique ];
}

sub search_annotation_types {
    my ($search_term, $type_group, $db) = @_;
    return '' unless $search_term && $db;

    $search_term = '%' . $search_term . '%';

    my $group;
    if ($type_group) {
        $group = $db->resultset('AnnotationTypeGroup')->find( { name => $type_group } );
    }

    my @types;
    if ($group) {
        @types = $db->resultset("AnnotationType")->search(
            \[
                'annotation_type_group_id = ? AND (name LIKE ? OR description LIKE ?)',
                [ 'annotation_type_group_id', $group->id ],
                [ 'name',                     $search_term ],
                [ 'description',              $search_term ]
            ]
        );
    }
    else {
        @types = $db->resultset("AnnotationType")->search(
            \[
                'name LIKE ? OR description LIKE ?',
                [ 'name',        $search_term ],
                [ 'description', $search_term ]
            ]
        );
    }

    my %unique;
    map { $unique{ $_->name } = 1 } @types;
    return [ sort keys %unique ];
}

sub update_annotation {
    my %opts = @_;

    my $annotation_id = $opts{annotation_id};
    unless ($annotation_id) {
        warn 'update_annotation: annotation_id required';
        return;
    }
 
    my ($error, $type_id, $link, $image_id, $bisque_id, $bisque_file, $object) = _init(\%opts);
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
    $annotation->image_id($image_id) if ($image_id);
    if ($bisque_file) {
        $annotation->bisque_file($bisque_file);
        $annotation->bisque_id($bisque_id);
        $annotation->bisque_user($opts{user}->id);
    }
    elsif ($opts{'delete_bisque_image'}) {
        $annotation->bisque_file(undef);
        $annotation->bisque_id(undef);
        $annotation->bisque_user(undef);
        delete_bisque_image(lc($opts{target_type}), $object->id, $annotation->bisque_file, $annotation->bisque_id, $opts{user});
    }
    $annotation->update;

    return;
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

sub _create_image {
    my %opts = @_;
    my $fh = $opts{fh};
    my $filename = $opts{filename};
    my $db = $opts{db};
    return unless (($fh || $filename) && $db);

    unless (defined $fh) {
        unless (open($fh, $filename)) {
            print STDERR "Storage::create_image: ERROR, couldn't open file '$filename'\n";
            return;
        }
    }

    my $contents = read_file($fh, binmode => ':raw');
    unless ($contents) {
        print STDERR "Storage::create_image: ERROR, couldn't read file '$filename'\n";
        return;
    }

    my $image = $db->resultset('Image')->create({
        filename => $filename,
        image    => $contents
    });
    return unless $image;
}

sub _get_object {
    my ($object, $id, $object_type, $own_or_edit, $db, $user) = @_;
    $object = $db->resultset($object_type eq 'notebook' ? 'List' : ucfirst($object_type))->find($id) unless $object;
    return 'Resource not found' unless $object;
    if ($own_or_edit) {
        return 'User not logged in' unless $user;
        unless ($object_type) {
            $object_type = lc(ref($object));
            $object_type = substr($object_type, rindex($object_type, ':') + 1);
            $object_type = 'notebook' if $object_type eq 'list';
        }
        return 'Access denied' unless $user->is_owner_editor(lc($object_type) => $object->id);
    }
    return 'Access denied' unless !$object->restricted || ($user && $user->has_access_to($object));
    return undef, $object;
}

sub _init {
    my $opts = shift;
    my $type_name = $opts->{type_name};
    return 'type_name required' unless $type_name;

    my $db = $opts->{db};
    my ($error, $object) = _get_object($opts->{target}, $opts->{target_id}, $opts->{target_type}, 1, $db, $opts->{user});
    return $error if $error;
    my $group_name = $opts->{group_name};
    my $group_id;
    $group_id = $db->resultset('AnnotationTypeGroup')->find_or_create({ name => $group_name })->id if defined $group_name;

    my $type_id = $db->resultset('AnnotationType')->find_or_create({ name => $type_name, annotation_type_group_id => $group_id })->id;
    return 'Error creating annotation type' unless $type_id;

    my $link = $opts->{link};
    if ($link) {
        $link =~ s/^\s+//;
        $link = 'http://' . $link if !($link =~ /^(\w+)\:\/\//);
    }

    my $image_id;
    my $bisque_id;
    my $bisque_file;
    my $filename = $opts->{filename};
    if ($filename) {
        my $upload = $opts->{image};
        if ($upload) {
            ($bisque_id, $bisque_file) = create_bisque_image($object, $upload, $opts->{user});
            return 'Error creating image' unless $bisque_id;
        }
        else {
            my $image = _create_image($opts);
            return 'Error creating image' unless $image;
            $image_id = $image->id;
        }
    }

    return (undef, $type_id, $link, $image_id, $bisque_id, $bisque_file, $object);
}

1;
