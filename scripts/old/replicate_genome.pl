#!/usr/bin/perl -w

use strict;
use CoGeX;
use Data::Dumper;
use Getopt::Long;

my ($dsgid);

GetOptions( "dsgid=i" => \$dsgid, );

unless ($dsgid) {
    print qq{
Usage:
$0 -dsgid <dataset_group id>

This program generates a copy of a coge dataset_group (genome).  The genomic sequence is NOT copied and that from the original genome is used.

Options:

 -dsgid              coge database dataset_group

};
    exit;
}

my $connstr = 'dbi:mysql:dbname=DB;host=HOST;port=PORT';
my $coge = CoGeX->connect($connstr, 'USER', 'PASSWORD' );

#$coge->storage->debugobj(new DBIxProfiler());
#$coge->storage->debug(1);

my $dsg = $coge->resultset('DatasetGroup')->find($dsgid);

unless ($dsg) {
    print
"Problem with retrieving the coge database object for the genometo be copied.\nExiting.\n";
    exit;
}

my $name;
$name = $dsg->name if $dsg->name;
$name = $dsg->organism->name unless $name;

print "Replicating dataset_group_object $name\n";

$name .= " (copy)";

my $desc;
$desc = $dsg->description if $dsg->description;
$desc .= " " if $desc;
$desc .= "Copied from dataset_group id: " . $dsg->id;
my $new_dsg = $coge->resultset('DatasetGroup')->create(
    {
        name                     => $name,
        description              => $desc,
        version                  => $dsg->version,
        organism_id              => $dsg->organism->id,
        genomic_sequence_type_id => $dsg->genomic_sequence_type_id,
        file_path                => $dsg->file_path
    }
);

foreach my $gs ( $dsg->genomic_sequences ) {
    $new_dsg->add_to_genomic_sequences(
        {
            sequence_length => $gs->sequence_length,
            chromosome      => $gs->chromosome,
        }
    );
}
my @dss;
print "Replicating datasets\n";
foreach my $ds ( $dsg->datasets ) {
    my $new_ds = $coge->resultset('Dataset')->create(
        {
            name           => $ds->name,
            description    => $ds->description,
            version        => $ds->version,
            link           => $ds->link,
            data_source_id => $ds->data_source_id
        }
    );

    #create linker
    $coge->resultset('DatasetConnector')->create(
        {
            dataset_id       => $new_ds->id,
            dataset_group_id => $new_dsg->id,
        }
    );
    push @dss, [ $ds, $new_ds ];
}
print "Replicating features\n";
foreach my $item (@dss) {
    my ( $ds1, $ds2 ) = @$item;

    foreach my $feat1 ( $ds1->features ) {

        #    sleep (.3);
        my $chr   = $feat1->chromosome;
        my $feat2 = $ds2->add_to_features(

            #    print Dumper
            {
                feature_type_id => $feat1->type->id,
                start           => $feat1->start,
                stop            => $feat1->stop,
                strand          => $feat1->strand,
                chromosome      => $chr,
            }
        );
        foreach my $name ( $feat1->feature_names ) {
            $feat2->add_to_feature_names(

                #	print Dumper
                {
                    name         => $name->name,
                    description  => $name->description,
                    primary_name => $name->primary_name,

                    #				     };
                }
            );
        }
        foreach my $loc ( $feat1->locs ) {
            $feat2->add_to_locations(

                #	print Dumper
                {
                    start  => $loc->start,
                    stop   => $loc->stop,
                    strand => $loc->strand,

                    #				  chromosome=>$loc->chromosome,
                    chromosome => $chr,

                    #				 };
                }
            );
        }
        foreach my $anno ( $feat1->annotations ) {
            $feat2->add_to_annotations(

                #	print Dumper
                {
                    annotation         => $anno->annotation,
                    annotation_type_id => $anno->annotation_type_id,
                    link               => $anno->link,

                    #				   };
                }
            );
        }
        foreach my $seq ( $feat1->sequences ) {
            $feat2->add_to_sequences(

                #	print Dumper
                {
                    sequence_type_id => $seq->sequence_type_id,
                    sequence_data    => $seq->sequence_data,

                    #				  };
                }
            );
        }
    }
}

print "Finished\n";
print "New DatasetGroup id: " . $new_dsg->id, "\n";
print "Link:\n";
print "\t"
  . "http://genomevolution.org/CoGe/GenomeView.pl?dsgid="
  . $new_dsg->id . "\n";
