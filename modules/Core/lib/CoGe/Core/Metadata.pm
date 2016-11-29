package CoGe::Core::Metadata;

use v5.14;
use strict;
use warnings;

use File::Basename;
use File::Slurp;
use Data::Dumper;

use CoGeX;

BEGIN {
    our (@ISA, $VERSION, @EXPORT);
    require Exporter;

    $VERSION = 0.0.1;
    @ISA = qw(Exporter);
    @EXPORT = qw( create_annotations export_annotations to_annotations 
                  tags_to_string create_image );
}

sub create_annotations {
    my %opts = @_;
    my $db = $opts{db};
    my $target = $opts{target};           # experiment, genome, or list object
    my ($target_id, $target_type) = ($opts{target_id}, $opts{target_type}); # or an experiment/genome/list id and type
    my $annotations = $opts{annotations}; # semicolon-separated list of annotations (image_file|link|group|type|text;[...])
    my $locked = $opts{locked};           # boolean flag to indicate locked (not editable) annotations

    my @result = ();
    foreach ( split(/\s*;\s*/, $annotations) ) {
        my @tok = split(/\s*\|\s*/, $_);
        my ($image_file, $link, $group_name, $type_name, $anno_text);
        $image_file = shift @tok if (@tok == 5);
        $link       = shift @tok if (@tok == 4);
        $group_name = shift @tok if (@tok == 3);
        $type_name  = shift @tok if (@tok == 2);
        $anno_text  = shift @tok if (@tok == 1);
        unless ($anno_text and $type_name) {
            print STDERR "CoGe::Core::Metadata: missing required annotation type and text fields\n";
            print STDERR Dumper [ $image_file, $link, $group_name, $type_name, $anno_text ], "\n";
            return;
        }
        
        my $image_id;
        if ($image_file) {
            my $image = create_image(filename => $image_file, db => $db);
            unless($image) {
                print STDERR "CoGe::Core::Metadata: error creating image\n";
                return;
            }
            $image_id = $image->id;
        }

        # Create type group - first try to find a match by name only
        my ($group, $type, $anno);
        if ($group_name) {
            $group = $db->resultset('AnnotationTypeGroup')->find({ name => $group });
            if (!$group) {
                $group = $db->resultset('AnnotationTypeGroup')->create({ name => $group_name }); # null description
            }
            unless ($group) {
                print STDERR "CoGe::Core::Metadata: error creating annotation type group\n";
                return;
            }
        }

        # Create type - first try to find a match by name and group
        $type = $db->resultset('AnnotationType')->find({
            name => $type_name,
            annotation_type_group_id => ($group ? $group->id : undef) }
        );
        if (!$type) {
            $type = $db->resultset('AnnotationType')->create({
                name => $type_name,
                annotation_type_group_id => ($group ? $group->id : undef)
            }); # null description
        }
        unless ($type) {
            print STDERR "CoGe::Core::Metadata: error creating annotation type\n";
            return;
        }

        # Fetch target DBIX object if id/type given
        if ($target_id && $target_type) {
            if (lc($target_type) eq 'experiment') {
                $target = $db->resultset('Experiment')->find($target_id);
            }
            if (lc($target_type) eq 'genome') {
                $target = $db->resultset('Genome')->find($target_id);
            }
            if (lc($target_type) eq 'notebook') {
                $target = $db->resultset('List')->find($target_id);
            }
        }
        
        # Create annotation
        if ($target) {
            if (ref($target) =~ /Experiment/) {
                $anno = $db->resultset('ExperimentAnnotation')->find_or_create({
                    experiment_id => $target->id,
                    annotation_type_id => ($type ? $type->id : undef),
                    annotation => $anno_text,
                    link => $link,
                    image_id => $image_id,
                    locked => $locked
                }); # null description
            }
            elsif (ref($target) =~ /Genome/) {
                $anno = $db->resultset('GenomeAnnotation')->find_or_create({
                    genome_id => $target->id,
                    annotation_type_id => ($type ? $type->id : undef),
                    annotation => $anno_text,
                    link => $link,
                    image_id => $image_id,
                    locked => $locked
                }); # null description
            }
            elsif (ref($target) =~ /List/) {
                $anno = $db->resultset('ListAnnotation')->find_or_create({
                    list_id => $target->id,
                    annotation_type_id => ($type ? $type->id : undef),
                    annotation => $anno_text,
                    link => $link,
                    image_id => $image_id,
                    locked => $locked
                }); # null description
            }
            else {
                print STDERR "CoGe::Core::Metadata: unknown target type\n";
                return;
            }
        }
        else {
            print STDERR "CoGe::Core::Metadata: target not found\n";
            return;
        }
        
        unless ($anno) {
            print STDERR "CoGe::Core::Metadata: error creating annotation\n";
            return;
        }

        push @result, $anno;
    }

    return \@result;
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

        say $fh "#Type Group, Type, Annotation, Link, Image filename";
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

            if ($a->image) {
                my $filename = $a->image->filename;
                my $img =  File::Spec->catdir($export_path, $filename);

                eval {
                    open(my $imh, ">", $img) or die "image=$filename could not be generated";
                    print $imh $a->image->image;
                    close($imh);

                    push @files, $filename;
                };

                say $fh "log: error: $@" if ($@);
                say $fh qq{$group,"$info","$url","$filename"};
            } else {
                say $fh qq{$group,"$info","$url",""};
            }
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
        
        my @fields = map { $md->{$_} || '' } ('link', 'type_group', 'type', 'text');
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

sub create_image { # mdb: don't have a better place for this atm
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

1;
