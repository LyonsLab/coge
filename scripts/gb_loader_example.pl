#!/usr/bin/perl -w

use strict;

use GB;
use Genome;

#initialize genome database data source information
my $gb_file = "NC_003074.gbk";    #Arabidopsis chromosome 3;
my $gb;                           # = GB->new(file=>$gb_file);

##############
#HEADER OBJECT
##############
my $head_obj = $gb->header();

my $genome;                       # = Genome->new();
my $organism = $genome->get_organism_obj->find_or_create(
    {
        name        => $head_obj->organism,    #this adds "Arabidopsis thaliana"
        description => $head_obj
          ->organism_long,    #this adds the complete kingdom to species list
    }
);
my $data_source = $genome->get_data_source_obj->find_or_create(
    {
        name        => "NCBI",
        description => "National Center for Biotechnology Information",
        link        => "www.ncbi.nlm.nih.gov",
    }
);

my $data_information_desc = "LOCUS: " . $head_obj->locus . "\n";
$data_information_desc .= "ACCESSION: " . $head_obj->accession . "\n";
$data_information_desc .= "VERSION: " . $head_obj->version . "\n";

my $data_information = $genome->get_data_information_obj->find_or_create(
    {
        name => "$gb_file",
        desc => $data_information_desc,
        link =>
          "ftp://ftp.ncbi.nih.gov/genomes/Arabidopsis_thaliana/NC_003074.gbk",
        data_source_id => $data_source->id,
    }
);

##############
#FEATURES
##############
my $source;    #need to keep the source object for processing features
while ( my $gb_feature = $gb->get_next_feature ) {
    my $feat_type = $genome->get_feature_type_obj->find_or_create(
        { name => $gb_feature->type } );

    my $db_feature = $genome->get_feature_obj->create(
        {
            feature_type_id     => $feat_type->id,
            data_information_id => $data_information->id,
            organism_id         => $organism->id,
        }
    );

    #expect first feature to be the source feature!
    if ( $gb_feature->type =~ /source/i ) {
        $source = $gb_feature;    #set source;
             #generate name based on organism name and chromosome
        my $feat_name = $genome->get_feature_name_obj->create(
            {
                name        => $organism->name,
                description => "Chromosome " . $source->chromosome,
                feature_id  => $db_feature->id
            }
        );

        #generate name for accession
        $feat_name = $genome->get_feature_name_obj->create(
            {
                name       => $head_obj->accession,
                feature_id => $db_feature->id
            }
        );

        #generate name for version
        $feat_name = $genome->get_feature_name_obj->create(
            {
                name       => $head_obj->version,
                feature_id => $db_feature->id
            }
        );
    }

    #must have a source object to get chromosome info
    unless ($source) {
        die "Error -> can't proceed without a source feature";
    }

    ######################
    #LOCATION
    ######################
    foreach my $loc ( $gb_feature->locations ) {
        my $location = $genome->get_location_obj->create(
            {
                feature_id => $db_feature->id,
                start      => $loc->start,
                stop       => $loc->stop,
                strand     => $loc->strand,       #expect this to be "1" or "-1"
                chromosome => $source->chromosome
                , #make sure that chromosomes are in real numbers (not roman numerals [IX])
            }
        );
    }

    ######################
    #ANNOTATIONS
    ######################
    while ( my $anno = $gb_feature->next_Annotation() ) {

        #deal with db_xref
        #/db_xref=taxon:3702
        #/db_xref=GeneID:821318
        #/db_xref=GI:18379324
        if ( $anno->type =~ /xref/i ) {
            my $anno_type_group =
              $genome->get_annotation_type_group_obj->find_or_create(
                { name => $anno->type } );
            while ( my $xref_anno = $anno->next_Annotation() ) {
                my $anno_type =
                  $genome->get_annotation_type_obj->find_or_create(
                    {
                        name                     => $xref_anno->type,
                        annotation_type_group_id => $anno_type_group->id(),
                    }
                  );
                while ( my $sub_anno = $xref_anno->next_Annotation() ) {
                    $genome->get_annotation_obj->create(
                        {
                            name               => $sub_anno->to_string(),
                            feature_id         => $db_feature->id,
                            annotation_type_id => $anno_type->id
                        }
                    );
                }
            }
        }
        elsif (
               $anno->type =~ /locus_tag/i
            || $anno->type =~ /transcript_id/i
            || $anno->type =~ /protein_id/i
            || $anno->gene =~ /gene/i
            || $anno->type =~
            /synonym/i    ##synonyms are embedded in the /note= tag!
          )               #these are names
        {
            my $feat_name = $genome->get_feature_name_obj->create(
                {
                    name       => $anno->to_string,
                    feature_id => $db_feature->id
                }
            );
        }
        elsif ( $anno->type =~
            /translation/i )   #this needs to be entered into the sequence table
        {
            my $seq_type = $genome->get_sequence_type_obj->find_or_create(
                {
                    name        => "protein",
                    description => "translation",
                }
            );
            my $sequence = $genome->get_sequence_obj->create(
                {
                    sequence_type_id => $seq_type->id(),
                    sequence_data    => $anno->to_string(),
                    feature_id       => $db_feature->id(),
                }
            );
        }
        elsif ( $anno->type =~ /go/i
          ) #go will need special treatment in the parser as the GO annotation is embedded in the /note= tag.  Remove it from the /note tag and leave everything else in the /note tag.
            #example:
            #go_function: nucleic acid binding [goid 0003676]
        {
            my $anno_type_group =
              $genome->get_annotation_type_group_obj->find_or_create(
                { name => $anno->type } );    #this should be "go_function"
            while ( my $go_cat_anno = $anno->next_Annotation() ) {
                my $anno_type =
                  $genome->get_annotation_type_obj->find_or_create(
                    {
                        name => $go_cat_anno->type,    #this should be "0003676"
                        annotation_type_group_id => $anno_type_group->id(),
                    }
                  );
                while ( my $sub_anno = $go_cat_anno->next_Annotation() ) {
                    $genome->get_annotation_obj->create(
                        {
                            name => $sub_anno->to_string()
                            ,    #this should be "nucleic acid binding"
                            feature_id         => $db_feature->id,
                            annotation_type_id => $anno_type->id
                        }
                    );
                }
            }
        }
        else                     ##everything else
        {
            my $anno_type = $genome->get_annotation_type_obj->find_or_create(
                { name => $anno->type, } );
            my $name = $anno->to_string();
            warn
"Name is longer than 255 characters.  Will be truncated in the database:\n$name\n"
              if length $name > 255;
            $genome->get_annotation_obj->create(
                {
                    name               => $name,
                    feature_id         => $db_feature->id(),
                    annotation_type_id => $anno_type->id(),
                }
            );
        }
    }
}
