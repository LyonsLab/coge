# -*- perl -*-

# t/001_load.t - check module loading and create testing directory

use Test::More "no_plan";

BEGIN { use_ok( 'CoGe::ECNCS' ); }

my $object = CoGe::ECNCS->new ();
isa_ok ($object, 'CoGe::ECNCS');

my $algo_obj = $object->get_algorithm_obj();
isa_ok($algo_obj, 'CoGe::ECNCS::DB::Algorithm');
$algo_obj = $object->get_alg_obj();
isa_ok($algo_obj, 'CoGe::ECNCS::DB::Algorithm');
my $algdata_obj = $object->get_algorithm_data_obj();
isa_ok($algdata_obj, 'CoGe::ECNCS::DB::Algorithm_data');
$algdata_obj = $object->get_alg_data_obj();
isa_ok($algdata_obj, 'CoGe::ECNCS::DB::Algorithm_data');
my $algrun_obj = $object->get_algorithm_run_obj();
isa_ok($algrun_obj, 'CoGe::ECNCS::DB::Algorithm_run');
$algrun_obj = $object->get_alg_run_obj();
isa_ok($algrun_obj, 'CoGe::ECNCS::DB::Algorithm_run');
my $auth_obj = $object->get_author_obj();
isa_ok($auth_obj, 'CoGe::ECNCS::DB::Author');
$auth_obj = $object->get_auth_obj();
isa_ok($auth_obj, 'CoGe::ECNCS::DB::Author');
my $class_obj = $object->get_classification_obj();
isa_ok($class_obj, 'CoGe::ECNCS::DB::Classification');
$class_obj = $object->get_class_obj();
isa_ok($class_obj, 'CoGe::ECNCS::DB::Classification');
my $dm_obj = $object->get_data_mask_obj();
isa_ok($dm_obj, 'CoGe::ECNCS::DB::Data_mask');
$dm_obj = $object->get_mask_obj();
isa_ok($dm_obj, 'CoGe::ECNCS::DB::Data_mask');
my $ecncs_obj = $object->get_ecncs_obj();
isa_ok($ecncs_obj, 'CoGe::ECNCS::DB::Ecncs');
$ecncs_obj = $object->get_cns_obj();
isa_ok($ecncs_obj, 'CoGe::ECNCS::DB::Ecncs');
my $loc_obj = $object->get_location_obj();
isa_ok($loc_obj, 'CoGe::ECNCS::DB::Location');
$loc_obj = $object->get_loc_obj();
isa_ok($loc_obj, 'CoGe::ECNCS::DB::Location');
my $spike_obj = $object->get_spike_obj();
isa_ok($spike_obj, 'CoGe::ECNCS::DB::Spike');
my $status_obj = $object->get_status_obj();
isa_ok($status_obj, 'CoGe::ECNCS::DB::Status');
$status_obj = $object->get_stat_obj();
isa_ok($status_obj, 'CoGe::ECNCS::DB::Status');

my $alg = $algo_obj->find_or_create(name=>'test algorithm',
				   description=>'auto created by test script',
				   link=>'www.fake.net',
				   );
my $spike = $spike_obj->find_or_create(
				       length=>9,
				       sequence_data=>'test_spike',
				       );
my $class = $class_obj->find_or_create(name=>'test classification',
				       description=>'auto created by test script',
				       );
my $stat = $status_obj->find_or_create(name=>'test stat',
				       description=>'auto created by test script',
				       );
my $auth = $auth_obj->find_or_create(name=>'test author',
				     description=>'auto created by test script',
				     );
#create 100 algorithm runs
for (my $i = 1; $i <=100; $i++)
{
    my $algrun = $algrun_obj->create({algorithm_id     =>$alg->id,
				      algorithm_options=>"test options",
				  status_id        =>$stat->id,
				  classification_id=>$class->id,
				  author_id        =>$auth->id,
				  note             =>"auto_created by test script",
#				  date
				  });
    #make algo_data
    for (my $j=1; $j<=2; $j++)
    {
	my $loc = $loc_obj->create({start=>int(rand(10000)),
				    stop =>int(rand(10000))+1000,
				    start=>1,
				    data_information_id => "100",
				});
	my $algdata = $algdata_obj->create({
	    algorithm_run_id =>$algrun->id,
	    location_id      =>$loc->id,
	    feature_id       =>"100",
	    spike_id         =>$spike->id,
	});
    }
    #make data_mask
    for (my $j=1; $j<=5; $j++)
    {
	my $loc = $loc_obj->create({start=>int(rand(10000)),
				    stop =>int(rand(10000))+100,
				    start=>1,
				    data_information_id => "100",
				});
	my $datamask = $dm_obj->create({
	    algorithm_run_id =>$algrun->id,
	    location_id      =>$loc->id,
	});
    }

    #make ecncs
    for (my $j=1; $j<=5; $j++)
    {
	my $loc = $loc_obj->create({start=>int(rand(10000)),
				    stop =>int(rand(10000))+100,
				    start=>1,
				    data_information_id => "100",
				});
	my $ecncs = $ecncs_obj->create({
	    algorithm_run_id =>$algrun->id,
	    location_id      =>$loc->id,
	    name             => $j,
	    'pval'           =>"1.4e-5",
	    alignment        =>"ATTA-TGCTAGATCGGCTCTGATAGC-ACCGTGTGC",
	    similarity            =>"120/150",
	    identity         =>"89",
	    note             =>"auto_created by test script",
	    valid            => "Y",
	});
    }

}

$alg->delete;
$spike->delete;
$class->delete;
$stat->delete;
$auth->delete;
