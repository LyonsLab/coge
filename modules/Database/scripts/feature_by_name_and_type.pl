# -*- perl -*-


use CoGeX;
my $connstr = 'dbi:mysql:genomes:biocon:3306';
my $s = CoGeX->connect($connstr, 'cnssys', 'CnS' );

#$s->storage->debug(1);


my $rs = $s->resultset('Feature')->search(
                    { 
                    'feature_names.name' => { 'like' => 'At1%' },
                    'feature_type.name' =>   'CDS' 
                    },
                    {
                        join => [ 'feature_names','feature_type'] 
                    },

                );


my %seen;
while (my $feat =$rs->next()){
    my $fn = $feat->feature_names;
    my $type = $feat->feature_type->name;
    my @names = grep { ! $seen{$_->name} } $fn->next();
    map { $seen{$_->name}=1 }  $feat->feature_names->next(); 
    map { print $_->name  . "\n" } @names;
}
   

