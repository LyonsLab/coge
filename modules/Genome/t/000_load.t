# -*- perl -*-

# t/001_load.t - check module loading and create testing directory
use Data::Dumper;
use lib "../lib/";
use Test::More qw(no_plan);

BEGIN { use_ok( 'CoGe::Genome' ); }

my $o = CoGe::Genome->new ();
isa_ok ($o, 'CoGe::Genome');

my $ao = $o->get_annotation_obj();
isa_ok ($ao, 'CoGe::Genome::DB::Annotation');
$ao = $o->get_anno_obj();
isa_ok ($ao, 'CoGe::Genome::DB::Annotation');
my $ato = $o->get_annotation_type_obj();
isa_ok ($ato, 'CoGe::Genome::DB::Annotation_type');
$ato = $o->get_anno_type_obj();
isa_ok ($ato, 'CoGe::Genome::DB::Annotation_type');
my $atgo = $o->get_annotation_type_group_obj();
isa_ok ($atgo, 'CoGe::Genome::DB::Annotation_type_group');
$atgo = $o->get_anno_type_group_obj();
isa_ok ($atgo, 'CoGe::Genome::DB::Annotation_type_group');
my $dso = $o->get_data_source_obj();
isa_ok ($dso, 'CoGe::Genome::DB::Data_source');
$dso = $o->get_source_obj();
isa_ok ($dso, 'CoGe::Genome::DB::Data_source');
my $dio = $o->get_data_information_obj();
isa_ok ($dio, 'CoGe::Genome::DB::Data_information');
$dio = $o->get_information_obj();
isa_ok ($dio, 'CoGe::Genome::DB::Data_information');
my $fo = $o->get_feature_obj();
isa_ok ($fo, 'CoGe::Genome::DB::Feature');
$fo = $o->get_feat_obj();
isa_ok ($fo, 'CoGe::Genome::DB::Feature');
my $fno = $o->get_feature_name_obj();
isa_ok ($fno, 'CoGe::Genome::DB::Feature_name');
$fno = $o->get_feat_name_obj();
isa_ok ($fno, 'CoGe::Genome::DB::Feature_name');
my $fto = $o->get_feature_type_obj();
isa_ok ($fto, 'CoGe::Genome::DB::Feature_type');
$fto = $o->get_feature_type_obj();
isa_ok ($fto, 'CoGe::Genome::DB::Feature_type');
my $gso = $o->get_genomic_sequence_obj();
isa_ok ($gso, 'CoGe::Genome::DB::Genomic_sequence');
$gso = $o->get_genomic_seq_obj();
isa_ok ($gso, 'CoGe::Genome::DB::Genomic_sequence');
my $lo = $o->get_location_obj();
isa_ok ($lo, 'CoGe::Genome::DB::Location');
$lo = $o->get_loc_obj();
isa_ok ($lo, 'CoGe::Genome::DB::Location');
my $so = $o->get_organism_obj();
isa_ok ($so, 'CoGe::Genome::DB::Organism');
$so = $o->get_org_obj();
isa_ok ($so, 'CoGe::Genome::DB::Organism');
$so = $o->get_species_obj();
isa_ok ($so, 'CoGe::Genome::DB::Organism');
my $seqo = $o->get_sequence_obj();
isa_ok ($seqo, 'CoGe::Genome::DB::Sequence');
$seqo = $o->get_seq_obj();
isa_ok ($seqo, 'CoGe::Genome::DB::Sequence');
my $seqto = $o->get_sequence_type_obj();
isa_ok ($seqto, 'CoGe::Genome::DB::Sequence_type');
$seqto = $o->get_seq_type_obj();
isa_ok ($seqto, 'CoGe::Genome::DB::Sequence_type');

my $ds = $dso->find_or_create({name=>"test",
			       description=>"created by testing script",
			       link=>"www.test.com"
			      });
isa_ok($ds, 'CoGe::Genome::DB::Data_source');
my $org = $so->find_or_create({name=>"test org",
			       description=>"created by testing script",});
isa_ok($org, 'CoGe::Genome::DB::Organism');
my $di = $dio->find_or_create({name=>"test",
			       description=>"created by testing script",
			       link=>"www.test.com",
			       data_source_id=>$ds->id,
			       organism_id=>$org->id,
			      });
isa_ok($di, 'CoGe::Genome::DB::Data_information');
my $ft = $fto->find_or_create({name=>"test feature type",
			       description=>"created by testing script",});
isa_ok($ft, 'CoGe::Genome::DB::Feature_type');
my $atg = $atgo->find_or_create({name=>"test annotation type group",
				 description=>"created by testing script",});
isa_ok($atg, 'CoGe::Genome::DB::Annotation_type_group');
my $at = $ato->find_or_create({name=>"test annotation type",
			       description=>"created by testing script",
			       annotation_type_group_id=>$atg->id});
isa_ok($at, 'CoGe::Genome::DB::Annotation_type');
my $seqt = $seqto->find_or_create({name=>"test sequence type",
				   description=>"created by testing script",
				  });
isa_ok($seqt, 'CoGe::Genome::DB::Sequence_type');


#insert genomic sequence
for (my $c=1; $c<=5; $c++)
  {
    my $start = 1;
    my $stop = 0;
    for (my $s=1; $s<=3;$s++)
      {
	$stop += 100000; 
	my $seq = generate_random_seq();
	$gso->create({start=>$start,
		      stop=>$stop,
		      chromosome=>$c,
		      sequence_data=>$seq,
		      data_information_id=>$di->id,
		  });
	$start = $stop+1;
      }
    my $seq = $gso->get_sequence(start=>10,
				 stop =>20,
				 chr  => $c,
				 info_id=> $di->id,
				 strand => 1,
				 );
    my $seq2 = $gso->get_sequence(start=>10,
				 stop =>20,
				 chr  => $c,
				  info_id=> $di->id,
				 strand => -1,
				 );
}

my $c = 1;
my $f = $fo->find_or_create({
                             feature_type_id=>$ft->id,
                             data_information_id=>$di->id,
                            });
my $fn = $fno->find_or_create({feature_id=>$f->id,
			       name=>"test feature",
			       description=>"created by testing script",
			      });


for (my $i=1; $i<5; $i++)
  {
    $lo->find_or_create({feature_id=>$f->id,
			 start=>$i*100,
			 stop=>$i*100+50,
			 strand=>-1,
			 chromosome=>$c
			     
			});
    my $seq = generate_random_seq(1000);
    $seqo->find_or_create({sequence_type_id=>$seqt->id,
			 feature_id=>$f->id,
			 sequence_data=>$seq,
			});
    $fno->find_or_create({name=>"test feature name",
			  description=>"created by testing script",
			  feature_id=>$f->id,
			 });
    $ao->find_or_create({annotation=>"test annotation",
			 feature_id=>$f->id,
			 annotation_type_id=>$at->id,
			});
  }
$f->annotation_pretty_print;
$f->genbank_location_string;
foreach my $feat ($o->get_feature_by_name("test feature name"))
{
  isa_ok($feat, 'CoGe::Genome::DB::Feature');
}
		  
#print Dumper $f->organism, $f->feature_type, $f->data_information, $f->feature_names, $f->locations, $f->sequences, $f->annotations;
#print Dumper $org->id, $org->name, $org->desc, $org->genomic_sequences;
#print Dumper $org->features;
#$gso->delete_data_information($di->id);
#$ao->delete_data_information($di->id);
#$lo->delete_data_information($di->id);
#$fno->delete_data_information($di->id);
#$seqo->delete_data_information($di->id);
#$fo->delete_data_information($di->id);
can_ok ($ds, "delete");
$ds->delete;
can_ok ($ft, "delete");
$ft->delete;
can_ok ($atg, "delete");
$atg->delete;
can_ok ($seqt, "delete");
$seqt->delete;
can_ok ($org, "delete");
$org->delete;

#insert 3 billion bases
# for (my $c=1; $c<=23; $c++)
# {
#     my $start = 1;
#     for (my $s=1; $s<=1300;$s++)
#     {


sub generate_random_seq
  {
    my $length = shift || 100000;
    my $GC_bias = shift || 60;
    my $seq;
    for (1..$length)
      {
	if (int(rand(100))+1 <= $GC_bias)
	  {
	    if (int(rand(2)))
	      {
		$seq .= "G";
	      }
	    else
	      {
		$seq .= "C";
	      }
	  }
	else
	  {
	    if (int(rand(2)))
	      {
		$seq .= "A";
	      }
	    else
	      {
		$seq .= "T";
	      }
	    
	  }
      }
    return $seq;
  }
