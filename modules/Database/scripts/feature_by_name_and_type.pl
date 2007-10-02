# -*- perl -*-


use CoGeX;
my $connstr = 'dbi:mysql:genomes:biocon:3306';
my $s = CoGeX->connect($connstr, 'cnssys', 'CnS' );

#$s->storage->debug(1);

# real    3m35.392s
# after
#real    2m41.003s
# real 2m33


my $rs = $s->resultset('Feature')->search(
                    { 'feature_names.name' => { 'like' => 'At1g16%' },
                    'feature_type.name' =>   'CDS' }, {
                        join => [ 'feature_names','feature_type'] 
                    });
open(OUT,">","/tmp/seq_new.txt");
while (my $feat =$rs->next()){
    my $s =  $feat->genomic_sequence();
    print OUT $s . "\n";
}
   

