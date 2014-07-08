#!/usr/bin/perl -w

use strict;
use CoGeX;
use CoGe::Genome::Accessory::Annotation;
use DBIxProfiler;
use Benchmark qw/timethese/;

$| = 1;

my $connstr = 'dbi:mysql:dbname=DB;host=HOST;port=PORT';
my $coge = CoGeX->connect($connstr, 'USER', 'PASSWORD' );
$coge->storage->debugobj(new DBIxProfiler());
$coge->storage->debug(1);
#timethese(1000, { 'anno' => \&annotation_pretty_print,});

print annotation_pretty_print();

sub annotation_pretty_print
  {
    my ($self) = $coge->resultset('Feature')->esearch( {'me.feature_id' => '108234'})->next;
    my $anno_obj = new CoGe::Genome::Accessory::Annotation(Type=>"anno");
    $anno_obj->Val_delimit("\n");
    $anno_obj->Val_delimit("\n");
    $anno_obj->Add_type(0);
    $anno_obj->String_end("\n");
    my $start = $self->start;
    my $stop = $self->stop;
    my $chr = $self->chr;
    my $strand = $self->strand;
    #look into changing this to set_id
    my $info_id = $self->dataset->id;
    my $location = "Chr ".$chr." ";
    $location .= join (", ", map {$_->start."-".$_->stop} $self->locs);
    $location .="(".$strand.")";
    #my $location = "Chr ".$chr. "".$start."-".$stop.""."(".$strand.")";
    $anno_obj->add_Annot(new CoGe::Genome::Accessory::Annotation(Type=>"Location", Values=>[$location], Type_delimit=>": ", Val_delimit=>" "));
    my $anno_type = new CoGe::Genome::Accessory::Annotation(Type=>"Name(s)");
    $anno_type->Type_delimit(": ");
    $anno_type->Val_delimit(", ");
    foreach my $name ($self->names)
      {
	$anno_type->add_Annot($name);
      }

    $anno_obj->add_Annot($anno_type);
    foreach my $anno (sort {$b->type->name cmp $a->type->name} $self->annos)
      {
	my $type = $anno->type();
	my $group = $type->group();
	my $anno_type = new CoGe::Genome::Accessory::Annotation(Type=>$type->name);
	$anno_type->Val_delimit("\n");

	$anno_type->add_Annot($anno->annotation);
	if (ref ($group) =~ /group/i)
	  {
	    my $anno_g = new CoGe::Genome::Accessory::Annotation(Type=>$group->name);
	    $anno_g->add_Annot($anno_type);
	    $anno_g->Type_delimit(": ");
	    $anno_g->Val_delimit(", ");
	    $anno_obj->add_Annot($anno_g);
	  }
	else
	  {
	    $anno_type->Type_delimit(": ");
	    $anno_obj->add_Annot($anno_type);
	  }
      }
    return $anno_obj->to_String;
  }
