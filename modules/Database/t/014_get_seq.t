use CoGeX::Fasta;
use Test::More tests => 1;

my $db = "t/asdf.fasta";
my $seq = CoGeX::Fasta->get_seq({'seqid'=>'2', 'start'=>1, 'stop'=>11, blastdb=>$db});

like($seq, qr/[actgnx]/i, "sequence");
