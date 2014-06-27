#!/usr/bin/env perl

use strict;
use CGI;
use DBI;
use Data::Dumper;
use POSIX;
use JSON::Syck;
use LWP::Simple;

my $q = new CGI;
print "Content-Type: text/html\n\n";

my $tmpdir = "/opt/apache/CoGe/tmp/GEvo/";

if($ENV{SERVER_NAME} !~ /(toxic|synteny)/){
    $tmpdir = "tmp/";
}

my $db = "$tmpdir/" . $q->param('db') . ".sqlite";
unless (-r $db) {
    print STDERR $q->url(-query=>1),"\n";
    warn "database file $db does not exist or cannot be read!\n";
    exit;
}

our $dbh = DBI->connect("dbi:SQLite:dbname=$db") || die "cant connect to db";
my $sth;

sub get_cns {
    # get all the stuff that had been saved in the mysql db as CNS
    # (and here as well). use those with image_id == 1 to find pair in
    # image_id == 2, if we cant find the pair, then it must be
    # off-screen.
    $sth = $dbh->prepare(qq! SELECT xmin, xmax, ymin, ymax, id, pair_id, image_id
                             FROM image_data WHERE bpmin || "," || bpmax IN
                            (SELECT bpmin || "," || bpmax
                                FROM image_data WHERE type = "CNS")
                            AND type = "HSP"; !);
    $sth->execute();
    my @cnss;
    my %seen_ids;
    while (my $cns  = $sth->fetchrow_hashref()) {
        push(@cnss, $cns);
        $seen_ids{$cns->{id}}=1;
    }

    #print STDERR "CNS: " .  scalar(@cnss) . "\n";
    # remove any that dont have their partner in the visible window.
    @cnss = grep { $seen_ids{$_->{pair_id}} } @cnss;

    # make it easier to look up the pair.
    my %cnss; map { $cnss{$_->{id}} = $_ } @cnss;
    my @coordslist;

    foreach my $cns (values %cnss){
        # we'll get anything iwth image_id == 2 as a pair.
        if ($cns->{image_id} != 1) { next;}
        my $pair = $cnss{$cns->{pair_id}};
        push(@coordslist, {
             'img1' => [$cns->{xmin},  $cns->{ymin},  $cns->{xmax},  $cns->{ymax},  $cns->{id}]
            ,'img2' => [$pair->{xmin}, $pair->{ymin}, $pair->{xmax}, $pair->{ymax}, $pair->{id}]

        });
    }
    #print STDERR "CNS: " .  scalar(@cnss) . "\n";
    #print STDERR "coords: " .  Dumper @coordslist;
    return \@coordslist;
}

if ($q->param('get_info')){
    my %result;
    my %data;
    $sth = $dbh->prepare("SELECT * FROM image_info order by display_id");
    my $sth2 = $dbh->prepare("select min(xmin), max(xmax), image_id from image_data where type='anchor' group by image_id order by image_id;");
#    my $sth3 = $dbh->prepare("select * from image_info order by iname;");
    $sth->execute();
    while( my $title = $sth->fetchrow_hashref() ) {
        $data{$title->{display_id}}{img}{image_name}=$title->{iname};
        $data{$title->{display_id}}{img}{title}=$title->{title};
        $data{$title->{display_id}}{img}{width}=$title->{px_width};
        $data{$title->{display_id}}{img}{bpmin}=$title->{bpmin};
        $data{$title->{display_id}}{img}{bpmax}=$title->{bpmax};
        $data{$title->{display_id}}{img}{id}=$title->{id};
    }

    $sth2->execute();
    foreach my $anchor (@{$sth2->fetchall_arrayref()}) {
        $data{$anchor->[2]}{anchor}{max}=$anchor->[0];
        $data{$anchor->[2]}{anchor}{min}=$anchor->[1];
    }

    my $i = 0;
    foreach my $item (sort {$a<=>$b} keys %data) {

        my $img = $data{$item}{img};
        my $name = $img->{image_name};
        $result{$name}{title} = $img->{title};
        $result{$name}{i}=$i;
        $result{$name}{extents} = {img_width=>$img->{width},
                        bpmin=>$img->{bpmin},
                        bpmax=>$img->{bpmax}
                       };
        my $anc = $data{$item}{anchor};
        $result{$name}{anchors} = {
                        idx=>$data{$item}{img}{id},
                        xmax=>$anc->{min},
                        xmin=>$anc->{max},
                       };
        $i++;
    }

#    print STDERR Dumper %result;
    my ($coordslist, $rand)  = &get_cns();
    $result{'CNS'} = $coordslist;
    print JSON::Syck::Dump(\%result);
    undef $coordslist ;
    exit();
}

if ($q->param('predict')){
    my $log = $db;
    $log =~ s/sqlite$/log/g;
    open(FH, "<", $log);
    my ($bl2seq, $eval1, $eval2);
    my $seen = 0;
    while (my $line = <FH>){
        last if $seen == 2;
        chomp $line;
        if ($line =~ /bl2seq/){ $bl2seq = $line; }
        elsif ($line =~ /cutoff/i){
            print STDERR $line . "\n";
            if ($seen == 0 ){ $eval1 = $line; $seen++; }
            elsif ($seen == 1){ $eval2 = $line; $seen++; }
        }
    }

    $bl2seq =~ s/.+(\/usr\/bin\/bl2seq.+)/$1/;
    $eval1  =~ s/.+\s([^\s]+)$/$1/;
    $eval2  =~ s/.+\s([^\s]+)$/$1/;
    chomp $bl2seq;

    # if necessary, fix for dev machine...
    if( $tmpdir =~ /^\/var\/www\//){ # TODO check -e (exists)
        $bl2seq =~ s/\/opt\/apache\/CoGe\/tmp\//$tmpdir\/tmpdir\//g;
    }
    $bl2seq =~ s/\-o\s[^\s]+//;
    # use tab-delimited output, and only the top strand
    $bl2seq .= " -D 1 -S 1 ";
    my $outfile = $log;
    $outfile =~ s/log$/blast/;
    print STDERR $outfile;
    `$bl2seq | grep -v '#' > $outfile`;
    # use tab-delimited output, and only the top strand
    my @res = `/usr/bin/python predict_cns.py "$outfile" $eval1 $eval2`;
    print STDERR "FROM PYTHON:" . scalar(@res)  . " pairs\n";
    close(FH);
    exit();
}

my $x    = $q->param('x');
my $y    = $q->param('y');
my $all  = $q->param('all') || 0;
my $img_id = $q->param('img');

my $statement;
if ($q->param('follow')){
}
elsif ($q->param('bbox')) {
    my @bbox = split(/,/, $q->param('bbox'));
    my $query ="SELECT * FROM image_data WHERE ? + 1 > xmin AND ? - 1 < xmax AND ? - 1 > ymin AND ? + 1 < ymax AND image_id = ? and pair_id != -99 and type = 'HSP'";
    $sth = $dbh->prepare($query);
    $sth->execute($bbox[2], $bbox[0], $bbox[3], $bbox[1], $img_id);
}
elsif($all){
    $sth = $dbh->prepare("SELECT distinct(image_track) FROM image_data WHERE ? BETWEEN ymin and ymax and image_id = ? order by abs(image_track) DESC");
    $sth->execute($y, $img_id);
    my ($track) = $sth->fetchrow_array();

    $statement = qq{ SELECT id, xmin, xmax, ymin, ymax, image_id, image_track, pair_id, color FROM image_data
                    WHERE ( (image_track = ?) or (image_track = (? * -1) ) )
                    and image_id = ? and pair_id != -99 and type = "HSP" };
    $sth = $dbh->prepare($statement);

    $sth->execute($track, $track, $img_id);
}
else{
  my $query ="SELECT * FROM image_data WHERE ? + 3 > xmin AND ? - 3 < xmax AND ? BETWEEN ymin and ymax and image_id = ?";
    $sth = $dbh->prepare($query);
  my @params = ($x, $x, $y, $img_id);
#  print STDERR $query,"\n";
#  print STDERR Dumper \@params;
    $sth->execute(@params);
}

my @results;
while( my $result = $sth->fetchrow_hashref() ){
    my $sth2 = $dbh->prepare("SELECT id, xmin, xmax, ymin, ymax, image_id, image_track, pair_id, color FROM image_data where id = ?");
    $sth2->execute( $result->{pair_id} );
    my $pair = $sth2->fetchrow_hashref();

    my $annotation = $result->{annotation};
    $annotation = $annotation =~ /http/ ? get($annotation) : $annotation;
#    $annotation =~ s/"/'/g;
#    $annotation =~ s/<br\/?>/\n/ig;
#    $annotation =~ s/<.*?>//g;
    my $f1name = $result->{image_id}; # GEvo_rIKDAf4x_1.png -> 1
    my $f2name = $pair->{image_id};

    # TODO: clean this up. we should know if there's a pair or not.
    my @f1pts = map {floor  $result->{$_} + 0.5 } qw/xmin ymin xmax ymax/; push(@f1pts, $result->{'id'}); push(@f1pts, $result->{'image_track'});
    my @f2pts = map { floor $pair->{$_} + 0.5 } qw/xmin ymin xmax ymax/;   push(@f2pts, $pair->{'id'}); push(@f2pts, $pair->{'image_track'});
    my $has_pair = 0;
    map { $has_pair += $_ } @f2pts;

    my $link = $result->{link};
    my $color = ($result->{color} ne 'NULL' && $result->{color} || $pair->{color}) ;
    $color =~ s/#/0x/;

    push(@results, {  link       => "/CoGe/$link"
                    , annotation => $annotation
                    , has_pair   => $has_pair
                    , color      => $color
                    , features   => {'key' . $f1name => \@f1pts,'key'. $f2name => \@f2pts}
                 });
    #print STDERR Dumper @results;
}
print JSON::Syck::Dump({resultset => \@results});
