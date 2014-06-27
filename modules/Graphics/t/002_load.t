# -*- perl -*-

# t/002_load.t - check module loading and create testing directory

use Test::More qw(no_plan);

BEGIN {
    use_ok( 'CoGe::Graphics::Chromosome' );
    use_ok( 'CoGe::Graphics::Feature' );
    use_ok( 'CoGe::Graphics::Feature::Gene' );
    use_ok( 'CoGe::Graphics::Feature::NucTide' );
}

my $o = CoGe::Graphics::Chromosome->new ();
isa_ok ($o, 'CoGe::Graphics::Chromosome');

# my $a = CoGe::Graphics::Feature->new({start=>50000});
# isa_ok ($a, 'CoGe::Graphics::Feature');
# $a->ih(5);
# $a->iw(5);
# $a->gd->fill(0,0,$a->get_color(200,200, 255));
# $a->label("A");
# $a->strand("1");
# $a->fill(1);

# my $t = CoGe::Graphics::Feature->new({start=>50000});
# isa_ok ($t, 'CoGe::Graphics::Feature');
# $t->ih(5);
# $t->iw(5);
# $t->gd->fill(0,0,$t->get_color(200,200, 255));
# $t->label("T");
# $t->strand("-1");
# $t->fill(1);

# my $c = CoGe::Graphics::Feature->new({start=>50001});
# isa_ok ($c, 'CoGe::Graphics::Feature');
# $c->ih(5);
# $c->iw(5);
# $c->gd->fill(0,0,$c->get_color(200,255, 200));
# $c->label("C");
# $c->strand("1");
# $c->fill(1);

# my $g = CoGe::Graphics::Feature->new({start=>50001});
# isa_ok ($g, 'CoGe::Graphics::Feature');
# $g->ih(5);
# $g->iw(5);
# $g->gd->fill(0,0,$g->get_color(200,255, 200));
# $g->label("G");
# $g->strand("-1");
# $g->fill(1);

my $ge = CoGe::Graphics::Feature->new({start=>49990, stop=>50010});
isa_ok ($ge, 'CoGe::Graphics::Feature');
$ge->ih(10);
$ge->iw(10);
$ge->gd->fill(0,0,$ge->get_color(250,20, 200));
$ge->label("Gene");
$ge->strand("1");

my $ge2 = CoGe::Graphics::Feature->new({start=>49980, stop=>50000});
isa_ok ($ge2, 'CoGe::Graphics::Feature');
$ge2->ih(10);
$ge2->iw(10);
$ge2->gd->fill(0,0,$ge2->get_color(255,20, 25));
$ge2->label("Gene");
$ge2->strand("-1");

my $ge3 = CoGe::Graphics::Feature::Gene->new();
isa_ok ($ge3, 'CoGe::Graphics::Feature::Gene');
$ge3->strand("1");
$ge3->color([0,255,0]);
$ge3->add_segment(start=>49900, end=>50000);
$ge3->add_segment(start=>50100, end=>50200);
$ge3->add_segment(start=>50300, end=>50400);
$ge3->label("multi-part_gene");

my $ge4 = CoGe::Graphics::Feature::Gene->new();
isa_ok ($ge4, 'CoGe::Graphics::Feature::Gene');
$ge4->strand("-1");
$ge4->color([255,0,255]);
$ge4->add_segment(start=>49900, end=>50000);
$ge4->add_segment(start=>50100, end=>50200);
$ge4->add_segment(start=>50300, end=>50400);
$ge4->label("multi-part_gene");

my $ge5 = CoGe::Graphics::Feature::Gene->new();
isa_ok ($ge5, 'CoGe::Graphics::Feature::Gene');
$ge5->strand("-1");
$ge5->color([255,0,255]);
$ge5->add_segment(start=>48900, end=>50000);
$ge5->add_segment(start=>50500, end=>51200);
$ge5->add_segment(start=>52300, end=>54400);
$ge5->label("multi-part_gene");

my $ge6 = CoGe::Graphics::Feature::Gene->new();
isa_ok ($ge6, 'CoGe::Graphics::Feature::Gene');
$ge6->strand("1");
$ge6->color([0,0,255]);
$ge6->add_segment(start=>48900, end=>50000);
$ge6->add_segment(start=>50500, end=>51200);
$ge6->add_segment(start=>52300, end=>54400);
$ge6->label("multi-part_gene");

my $ge7 = CoGe::Graphics::Feature::Gene->new();
isa_ok ($ge7, 'CoGe::Graphics::Feature::Gene');
$ge7->strand("1");
$ge7->color([255,0,0]);
$ge7->add_segment(start=>48900, end=>50000);
$ge7->add_segment(start=>50500, end=>51200);
$ge7->add_segment(start=>52300, end=>54400);
$ge7->label("multi-part_gene");
$ge7->order(10);

my $ge8 = CoGe::Graphics::Feature::Gene->new();
isa_ok ($ge8, 'CoGe::Graphics::Feature::Gene');
$ge8->strand("1");
$ge8->color([0,100,100]);
$ge8->add_segment(start=>55000, end=>56000);
$ge8->add_segment(start=>56500, end=>57200);
$ge8->add_segment(start=>58300, end=>59400);
$ge8->label("multi-part_gene");
$ge8->order(1);

#$o->add_feature($a);
#$o->add_feature($t);
#$o->add_feature($c);
#$o->add_feature($g);
$o->add_feature($ge);
$o->add_feature($ge2);
$o->add_feature($ge3);
$o->add_feature($ge4);
$o->add_feature($ge5);
$o->add_feature($ge6);
$o->add_feature($ge7);
$o->add_feature($ge8);

my %nt = (
	  0=>"A",
	  1=>"T",
	  2=>"C",
	  3=>"G",
	 );

my %comp = (
	    "A"=>"T",
	    "T"=>"A",
	    "C"=>"G",
	    "G"=>"C",
	   );

 for my $i (54300..54500)
   {
     my $nt = $nt{int rand(4)};
     my $cnt = $comp{$nt};
     my $f1 = CoGe::Graphics::Feature::NucTide->new({nt=>$nt, strand=>1, start=>$i});
     my $f2 = CoGe::Graphics::Feature::NucTide->new({nt=>$cnt, strand=>-1, start=>$i});
     $o->add_feature($f1);
     $o->add_feature($f2);
   }

$o->DEBUG(1);
$o->chr_length(250000000);
#$o->set_point (125000000);
$o->set_point (54400);

#$o->set_region(start=>5000, stop=>6000);
$o->mag(14);
$o->iw(800);
$o->draw_chromosome(1);

$o->generate_png();
