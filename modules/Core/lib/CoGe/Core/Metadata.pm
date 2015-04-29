package CoGe::Core::Metadata;

use v5.14;
use strict;
use warnings;

use File::Basename;

use CoGeX;

BEGIN {
    our (@ISA, $VERSION, @EXPORT);
    require Exporter;

    $VERSION = 0.0.1;
    @ISA = qw(Exporter);
    @EXPORT = qw( create_annotations export_annotations );
}

sub create_annotations {
    my %opts = @_;
    my $db = $opts{db};
    my $target = $opts{target};           # experiment, genome, or list object
    my ($target_id, $target_type) = ($opts{target_id}, $opts{target_type}); # or an experiment/genome/list id and type
    my $annotations = $opts{annotations}; # semicolon-separated list of annotations (link:group:type:text;...)
    my $locked = $opts{locked};           # boolean

    my @result = ();
    foreach ( split(/\s*;\s*/, $annotations) ) {
        my @tok = split(/\s*\|\s*/, $_);
        my ($link, $group_name, $type_name, $anno_text);
        $link  = shift @tok if (@tok == 4);
        $group_name = shift @tok if (@tok == 3);
        $type_name  = shift @tok if (@tok == 2);
        $anno_text  = shift @tok if (@tok == 1);
        unless ($anno_text and $type_name) {
            print STDERR "missing required annotation type and text fields\n";
            return;
        }

        # Create type group - first try to find a match by name only
        my ($group, $type, $anno);
        if ($group_name) {
            $group = $db->resultset('AnnotationTypeGroup')->find({ name => $group });
            if (!$group) {
                $group = $db->resultset('AnnotationTypeGroup')->create({ name => $group }); # null description
            }
            unless ($group) {
                print STDERR "error creating annotation type group\n";
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
            print STDERR "error creating annotation type\n";
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
                $anno = $db->resultset('ExperimentAnnotation')->create({
                    experiment_id => $target->id,
                    annotation_type_id => ($type ? $type->id : undef),
                    annotation => $anno_text,
                    link => $link,
                    locked => $locked
                }); # null description
            }
            elsif (ref($target) =~ /Genome/) {
                $anno = $db->resultset('GenomeAnnotation')->create({
                    genome_id => $target->id,
                    annotation_type_id => ($type ? $type->id : undef),
                    annotation => $anno_text,
                    link => $link,
                    locked => $locked
                }); # null description
            }
            elsif (ref($target) =~ /List/) {
                $anno = $db->resultset('ListAnnotation')->create({
                    list_id => $target->id,
                    annotation_type_id => ($type ? $type->id : undef),
                    annotation => $anno_text,
                    link => $link,
                    locked => $locked
                }); # null description
            }
            else {
                print STDERR "unknown target type\n";
                return;
            }
        }
        else {
            print STDERR "target not found\n";
            return;
        }
        
        unless ($anno) {
            print STDERR "error creating annotation\n";
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

1;
