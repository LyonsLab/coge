package CoGeX::Result::Dataset;

use strict;
use warnings;
use Data::Dumper;
use POSIX;
use Carp qw (cluck);
use CoGe::Core::Storage qw(reverse_complement);
use CoGeDBI qw(get_features get_feature_names get_feature_annotations get_locations);

use base 'DBIx::Class::Core';

use Text::Wrap;
use Carp;

=head1 NAME

CoGeX::Dataset

=head1 SYNOPSIS

This object uses the DBIx::Class to define an interface to the C<dataset> table in the CoGe database.

=head1 DESCRIPTION

=head1 AUTHORS

 Eric Lyons
 Brent Pedersen

=head1 COPYRIGHT

This program is free software; you can redistribute
it and/or modify it under the same terms as Perl itself.

The full text of the license can be found in the
LICENSE file included with this module.

=cut

__PACKAGE__->table("dataset");
__PACKAGE__->resultset_class("CoGeX::ResultSet::Dataset");
__PACKAGE__->add_columns(
    "dataset_id",
    {
        data_type     => "INT",
        default_value => undef,
        is_nullable   => 0,
        size          => 11
    },
    "data_source_id",
    { data_type => "INT", default_value => 0, is_nullable => 0, size => 11 },
    "name",
    {
        data_type     => "VARCHAR",
        default_value => "",
        is_nullable   => 0,
        size          => 100
    },
    "description",
    {
        data_type     => "VARCHAR",
        default_value => undef,
        is_nullable   => 1,
        size          => 255,
    },
    "version",
    {
        data_type     => "VARCHAR",
        default_value => undef,
        is_nullable   => 1,
        size          => 50,
    },
    "link",
    {
        data_type     => "TEXT",
        default_value => undef,
        is_nullable   => 1,
        size          => 65535,
    },
    "date",
    {
        data_type     => "DATETIME",
        default_value => "",
        is_nullable   => 0,
        size          => 19
    },
    "restricted",
    { data_type => "int", default_value => "0", is_nullable => 0, size => 1 },
    "deleted",
    { data_type => "int", default_value => "0", is_nullable => 0, size => 1 },
    "creator_id",
    { data_type => "INT", default_value => 0, is_nullable => 0, size => 11 },
);

__PACKAGE__->set_primary_key("dataset_id");
__PACKAGE__->has_many( "features" => "CoGeX::Result::Feature", 'dataset_id' );
__PACKAGE__->has_many(
    "dataset_connectors" => "CoGeX::Result::DatasetConnector",
    'dataset_id'
);
__PACKAGE__->belongs_to(
    "data_source" => "CoGeX::Result::DataSource",
    'data_source_id'
);

################################################ subroutine header begin ##

=head2 genomes

 Usage     :
 Purpose   :
 Returns   :
 Argument  :
 Throws    :
 Comments  :

See Also   :

=cut

################################################## subroutine header end ##

sub genomes {
    my $self = shift;
    my %opts = @_;
    my $chr  = $opts{chr};

    my @genomes;
    foreach my $dsc ( $self->dataset_connectors() ) {
        if ( defined $chr ) {
            my %chrs = map { $_, 1 } $dsc->genome->chromosomes;
            next unless $chrs{$chr};
        }
        push @genomes, $dsc->genome;
    }
    return wantarray ? @genomes : \@genomes;
}

sub first_genome {
    my $self = shift;
    my %opts = @_;
    my $skip_deleted = $opts{skip_deleted};
    
    foreach my $dsc ( $self->dataset_connectors() ) {
        next if ($skip_deleted && $dsc->genome->deleted);
        return $dsc->genome;
    }
}

sub dataset_groups {
    cluck "Dataset::dataset_groups is obselete, please use ->genomes\n";
    shift->genomes(@_);
}

################################################ subroutine header begin ##

=head2 organism

 Usage     :
 Purpose   :
 Returns   :
 Argument  :
 Throws    :
 Comments  :

See Also   :

=cut

################################################## subroutine header end ##

sub organism {
    my $self = shift;
    my %opts = @_;
    my %orgs = map { $_->id, $_ } map { $_->organism } $self->genomes;
    if ( keys %orgs > 1 ) {
        warn "sub organism in Dataset.pm fetched more than one organism!  Very odd:\n";
        warn join( "\n", map { $_->name } values %orgs ), "\n";
        warn "Only one will be returned\n";
    }
    my ($org) = values %orgs;
    return $org;
}

sub datasource {
    shift->data_source(@_);
}

sub source {
    shift->data_source(@_);
}

sub desc {
    shift->description(@_);
}

################################################ subroutine header begin ##

=head2 info

 Usage     : $self->info
 Purpose   : returns a string of information about the data set.

 Returns   : returns a string
 Argument  : none
 Throws    :
 Comments  : To be used to quickly generate a string about the data set

See Also   :

=cut

################################################## subroutine header end ##

sub info {
    my $self = shift;
    my $info;
    $info .= $self->name if $self->name;
    $info .= ": " . $self->description if $self->description;
    $info .=
      " (v" . $self->version . ", " . $self->date . ", id" . $self->id . ")";
    return $info;
}

################################################ subroutine header begin ##

=head2 get_genomic_sequence

 Usage     :
 Purpose   :
 Returns   :
 Argument  :
 Throws    :
 Comments  :

See Also   :

=cut

################################################## subroutine header end ##

# mdb replaced 7/31/13, issue 77
#sub get_genomic_sequence
#{
#	my $self = shift;
#	my %opts = @_;
#
#	my $start = $opts{start} || $opts{begin};
#	my $stop  = $opts{stop}  || $opts{end};
#	my $chr   = $opts{chr};
#	$chr = $opts{chromosome} unless defined $chr;
#	my $strand   = $opts{strand};
#	my $seq_type = $opts{seq_type} || $opts{gstid};
#	my $debug    = $opts{debug};
#	my $dsgid    = $opts{dsgid};
#	my $server   = $opts{server};                     #server from which to retrieve genomic sequence if not stored on local machine.  Web retrieval from CoGe/GetSequence.pl
#	my $dsg;
#	$dsg = $dsgid if $dsgid && ref($dsgid) =~ /Genome/;
#	return $dsg->genomic_sequence( start => $start, stop => $stop, chr => $chr, strand => $strand, debug => $debug, server => $server ) if $dsg;
#	my $seq_type_id = ref($seq_type) =~ /GenomicSequenceType/i ? $seq_type->id : $seq_type;
#	$seq_type_id = 1 unless $seq_type_id && $seq_type_id =~ /^\d+$/;
#
#	foreach my $tmp_dsg ( $self->groups )
#	{
#		if ( ( $dsgid && $tmp_dsg->id == $dsgid ) || ( $seq_type_id && $tmp_dsg->genomic_sequence_type->id == $seq_type_id ) )
#		{
#			return $tmp_dsg->genomic_sequence( start => $start, stop => $stop, chr => $chr, strand => $strand, debug => $debug, server => $server );
#		}
#	}
#
#	#hmm didn't return -- perhaps the seq_type_id was off.  Go ahead and see if anything can be returned
#	#    carp "In Dataset.pm, sub get_genomic_sequence.  Did not return sequence from a genome with a matching sequence_type_id.  Going to try to return some sequence from any genome.\n";
#	($dsg) = $self->groups;
#	return $dsg->genomic_sequence( start => $start, stop => $stop, chr => $chr, strand => $strand, debug => $debug, server => $server );
#}
sub get_genomic_sequence {
    my $self = shift;
    my %opts = @_;

    my $start  = $opts{start};
    my $stop   = $opts{stop};
    my $chr    = $opts{chr};
    $chr = $opts{chromosome} unless defined $opts{chr};
    my $strand = $opts{strand};
    my $gstid  = $opts{gstid};
    #print STDERR "Dataset->get_genomic_sequence:  start $start stop $stop chr $chr strand $strand\n";
    $gstid = 1 unless $gstid;    #FIXME hardcoded type value
    my $gid = $opts{gid};
    my $genome = $opts{genome}; #genome object
    $genome = $gid if ( $gid && ref($gid) =~ /Genome/ );
    unless ($genome && ref($genome) =~ /Genome/) {
        foreach ( $self->genomes ) {
            if ( $_->genomic_sequence_type_id == $gstid ) {
                $genome = $_;
                last;
            }
        }

# Hmmm didn't return -- perhaps the seq_type_id was off.  Go ahead and see if anything can be returned.
        unless ($genome) {
            ($genome) = $self->genomes;
        }
    }
    return $genome->get_genomic_sequence(
        chr    => $chr,
        start  => $start,
        stop   => $stop,
        strand => $strand
    );
}

# mdb removed 7/31/13, issue 77
#sub get_genome_sequence
#{
#	return shift->get_genomic_sequence(@_);
#}
#sub genomic_sequence
#{
#	return shift->get_genomic_sequence(@_);
#}

################################################ subroutine header begin ##

=head2 trim_sequence

 Usage     :
 Purpose   :
 Returns   :
 Argument  :
 Throws    :
 Comments  :

See Also   :

=cut

################################################## subroutine header end ##

sub trim_sequence {
    my $self = shift;
    my ( $seq, $seqstart, $seqend, $newstart, $newend ) = @_;

    my $start = $newstart - $seqstart;
    my $stop = length($seq) - ( $seqend - $newend ) - 1;
    $seq = substr( $seq, $start, $stop - $start + 1 );

    return ($seq);
}

################################################## subroutine header start ##

=head2 last_chromsome_position

 Usage     : my $last = $genome_seq_obj->last_chromosome_position($chr);
 Purpose   : gets the last genomic sequence position for a dataset given a chromosome
 Returns   : an integer that refers to the last position in the genomic sequence refered
             to by a dataset given a chromosome
 Argument  : string => chromsome for which the last position is sought
 Throws    :
 Comments  :

See Also   :

=cut

################################################## subroutine header end ##

sub last_chromosome_position {
    my $self = shift;
    my $chr  = shift;
    return 0 unless defined $chr;

    my $dsg = $self->first_genome; #my ($dsg) = $self->genomes; # mdb changed 4/23/14 issue 364
#    my ($item) = $dsg->get_chromosome($chr);
#    unless ($item) {
#        warn "Dataset::last_chromosome_position: unable to find genomic_sequence object for '$chr'";
#        return 0;
#    }
#    my $stop = $item->sequence_length();
#    unless ($stop) {
#        warn "No genomic sequence for ", $self->name, " for chr $chr\n";
#        return 0;
#    }
#    return $stop;
	return $dsg->get_chromosome_length($chr);
}

################################################ subroutine header begin ##

=head2 last_chromosome_position_old

 Usage     :
 Purpose   :
 Returns   :
 Argument  :
 Throws    :
 Comments  :

See Also   :

=cut

################################################## subroutine header end ##

#sub last_chromosome_position_old {
#    my $self = shift;
#    my $chr  = shift;
#    my $stop = $self->genomic_sequences( { chromosome => "$chr", }, )->get_column('stop')->max;
#    unless ($stop) {
#        warn "No genomic sequence for ", $self->name, " for chr $chr\n";
#        return;
#    }
#    return $stop;
#}

################################################ subroutine header begin ##

=head2 sequence_type

 Usage     :
 Purpose   :
 Returns   :
 Argument  :
 Throws    :
 Comments  :

See Also   :

=cut

################################################## subroutine header end ##

sub total_length {
    my $self = shift;
    my %opts = @_;
    my $ftid = $opts{ftid};
    my $search;
    my $join = {
        select => [ { sum => 'stop' } ],
        as     => ['total_length']
        ,    # remember this 'as' is for DBIx::Class::ResultSet not SQL
    };
    if ($ftid) {
        $search = { 'feature_type_id' => $ftid };
    }
    else {
        $search = { 'name' => 'chromosome' };
        $join->{join} = 'feature_type';
    }

    my $rs = $self->features( $search, $join );
    my $total_length = $rs->first->get_column('total_length');
    return ( defined $total_length ? $total_length : 0);
}

############################################### subroutine header begin ##

=head2 chromosome_count

 Usage     : $self->chromosome_count
 Purpose   : get count of chromosomes in the dataset group
 Returns   : number
 Argument  :
 Throws    :
 Comments  :

See Also   :

=cut

################################################## subroutine header end ##

sub chromosome_count {
    my $self = shift;
    my %opts = @_;
    my $ftid = $opts{ftid};
    my $search;
    my $join;
    if ($ftid) {
        $search = { 'feature_type_id' => $ftid };
    }
    else {
        $search = { 'name' => 'chromosome' };
        $join->{join} = 'feature_type';
    }

    my $count = $self->features( $search, $join )->count;
    return $count;
}

sub has_gene_annotation {
    return shift->features( { 'feature_type_id' => { -in => [ 1, 2, 3 ] } } )->count;
}

sub sequence_type {
    my $self   = shift;
    my (@dsgs) = $self->genomes;
    my %types  = map { $_->id, $_ } map { $_->genomic_sequence_type } @dsgs;
    my @types  = values %types;

    #    my ($type) = $self->genomic_sequences->slice(0,0);
    #    return $type ? $type->genomic_sequence_type : undef;
    if ( @types == 1 ) {
        return shift @types;
    }
    elsif ( @types > 1 ) {
        return wantarray ? @types : \@types;
    }
    else {
        return undef;
    }
}

sub genomic_sequence_type {
    my $self = shift;
    return $self->sequence_type(@_);
}

################################################ subroutine header begin ##

=head2 get_chromosomes

 Usage     :
 Purpose   :
 Returns   :
 Argument  :
 Throws    :
 Comments  :

See Also   :

=cut

################################################## subroutine header end ##

sub get_chromosomes {
    my $self   = shift;
    my %opts   = @_;
    my $ftid   = $opts{ftid};   #feature_type_id for feature_type of name "chromosome";
    my $length = $opts{length}; #option to return length of chromosomes as well
    my $limit  = $opts{limit};  #optional number of chromosomes to return, sorted by size
    my $max    = $opts{max};    #optional number of chromosomes to avoid slow query, sorted by size
    my @data;

    #this query is faster if the feature_type_id of feature_type "chromosome" is known.
    #features of this type refer to the entire stored sequence which may be a fully
    #assembled chromosome, or a contig, supercontig, bac, etc.
    my $search = {};
    my $search_type = { order_by => { -desc => 'stop' } };
    if ($ftid) {
        $search->{feature_type_id} = $ftid;
    }
    else {
        $search->{name}      = "chromosome";
        $search_type->{join} = "feature_type";
    }
    
    if ($max && $self->features( $search, $search_type )->count() > $max ) {
        return;
    }

    if ($limit) {
        $search_type->{rows} = $limit;
    }
    
    if ($length) {
        @data = $self->features( $search, $search_type );
    }
    else {
        @data = map { $_->chromosome } $self->features( $search, $search_type );
    }
    unless (@data) {
        my %seen;
        foreach my $feat ( $self->features({}, { distinct => 'chromosome' }) ) {
	       $seen{$feat->chromosome} = 1 if $feat->chromosome;
        }
        @data = sort keys %seen;
    }
    return wantarray ? @data : \@data;
}

sub chromosomes {
    my $self = shift;
    $self->get_chromosomes(@_);
}

################################################ subroutine header begin ##

=head2 has_chromosome

 Usage     : $ds->has_chromosome(chr=>"12")
 Purpose   : test to see if a dataset has a particular chromsome
 Returns   : 1 if yes, 0 if no
 Argument  :
 Throws    :
 Comments  :

See Also   :

=cut

################################################## subroutine header end ##

sub has_chromosome {
    my $self = shift;
    my %opts = @_;
    my $chr  = $opts{chr};
    my ($res) =
      $self->features->count( { "feature_type.name" => "chromosome", },
        { join => ["feature_type"] } );
    if ($res) {
        my ($res) = $self->features(
            {
                "feature_type.name" => "chromosome",
                "chromosome"        => "$chr",
            },
            { join => ["feature_type"] }
        );
        return 1 if $res;
        return 0;
    }
    else {
        my ($res) = $self->features->count( { "chromosome" => "$chr", }, );
        return 1 if $res;
        return 0;
    }
}

################################################ subroutine header begin ##

=head2 percent_gc

 Usage     :
 Purpose   :
 Returns   :
 Argument  :
 Throws    :
 Comments  :

See Also   :

=cut

################################################## subroutine header end ##

sub percent_gc {
    my $self  = shift;
    my %opts  = @_;
    my $count = $opts{count};

    #    my $chr = $opts{chr};
    my $seq    = $self->get_genomic_sequence(%opts);
    my $length = length $seq;
    return unless $length;
    my ($gc) = $seq =~ tr/GCgc/GCgc/;
    my ($at) = $seq =~ tr/ATat/ATat/;
    my ($n)  = $seq =~ tr/nN/nN/;
    my ($x)  = $seq =~ tr/xX/xX/;
    return ( $gc, $at, $n, $x ) if $count;
    return sprintf( "%.4f", $gc / $length ), sprintf( "%.4f", $at / $length ),
      sprintf( "%.4f", $n / $length ), sprintf( "%.4f", $x / $length );
}

sub gc_content {
    shift->percent_gc(@_);
}

################################################ subroutine header begin ##

=head2 fasta

 Usage     :
 Purpose   :
 Returns   :
 Argument  :
 Throws    :
 Comments  :

See Also   :

=cut

################################################## subroutine header end ##

sub fasta {
    my $self = shift;
    my %opts = @_;
    my $col  = $opts{col};

    #$col can be set to zero so we want to test for defined variable
    $col = $opts{column} unless defined $col;
    $col = $opts{wrap}   unless defined $col;
    $col = 100           unless defined $col;
    my $chr = $opts{chr};
    ($chr) = $self->get_chromosomes unless defined $chr;
    my $strand = $opts{strand} || 1;
    my $start  = $opts{start}  || 1;
    $start = 1 if $start < 1;
    my $stop  = $opts{stop} || $self->last_chromosome_position($chr);
    my $prot  = $opts{prot};
    my $rc    = $opts{rc};
    my $gstid = $opts{gstid};
    $strand = -1 if $rc;
    my $seq = $self->get_genomic_sequence(
        chr   => $chr,
        start => $start,
        stop  => $stop,
        gstid => $gstid
    );
    $stop = $start + length($seq) - 1 if $stop > $start + length($seq) - 1;
    my $head = ">" . $self->organism->name . " (" . $self->name;
    $head .= ", " . $self->description if $self->description;
    $head .= ", v"
      . $self->version . ")"
      . ", Location: "
      . $start . "-"
      . $stop
      . " (length: "
      . ( $stop - $start + 1 )
      . "), Chromosome: "
      . $chr
      . ", Strand: "
      . $strand;

    $Text::Wrap::columns = $col;
    my $fasta;

    $seq = reverse_complement($seq) if $rc;
    if ($prot) {
        my $trans_type = $self->trans_type;
        my $feat       = new CoGeX::Result::Feature;
        my ( $seqs, $type ) = $feat->frame6_trans(
            seq        => $seq,
            trans_type => $trans_type,
            gstid      => $gstid
        );
        foreach my $frame ( sort { length($a) <=> length($b) || $a cmp $b }
            keys %$seqs )
        {
            $seq = $seqs->{$frame};
            $seq = reverse_complement($seq) if $rc;
            $seq = join( "\n", wrap( "", "", $seq ) ) if $col;
            $fasta .= $head . " $type frame $frame\n" . $seq . "\n";
        }
    }
    else {
        $seq = join( "\n", wrap( "", "", $seq ) ) if $col;
        $fasta = $head . "\n" . $seq . "\n";
    }
    return $fasta;
}

################################################ subroutine header begin ##
=head2 gff
 Usage     : $ds->gff(print=>1)
 Purpose   : generating a gff file for a dataset from all the features it contains
 Returns   : a string
 Argument  : name_re     =>    regular expression for only displaying features containing a name that matches
             print       =>    print the gff file as the lines are retrieved
             annos       =>    print annotations as well (takes longer)
             debug       =>    prints some debugging stuff
             no_gff_head =>    won't print "gff-version 3".  Used when this function is called by Genome->gff(); (default 0)
             ds          =>    Dataset object.  Uses $self if none is specified
             id          =>    Starting number to be used for ID tag.  This in incremented by one with each entry.
             cds         =>    Only print CDS gene features (skip all ncRNA and other features).  Will print genes, mRNA, and CDS entries
             id_type     =>    Specify if the GFF entry IDs are going to be unique numbers or unique names.
             unique_parent_annotations => Flag to NOT print redundant annotations in children entries.  E.g. if parent has an annotation, a child will not have that annotation
             name_unique =>   Flag for specifying that the name tag of an entry will be unique
See Also   : genome->gff
=cut
################################################## subroutine header end ##
our ($FEATURES, $FEATURE_NAMES, $FEATURE_ANNOS, $FEATURE_LOCS, %FEATURES_BY_NAME);
sub feature_sort { $a->{start} <=> $b->{start} || $a->{type_id} <=> $b->{type_id} || $a->{id} <=> $b->{id} };
sub gff {
    my $self = shift;
    my %opts = @_;
    my $name_re = $opts{name_re};  #regular expression to search for a specific name
    my $debug   = $opts{debug};    #debug flag
    my $print   = $opts{print};    #flag to print gff as it is being retrieved
    my $annos   = $opts{annos};    #flag to retrieve and add annotations
    my $no_gff_head = $opts{no_gff_head}; #flag to NOT print gff headers (used in conjunction with genome->gff which generates its own headers
    my $ds    = $opts{ds};  #ds object, uses self if one is not passed in
    my $count = $opts{id};  #number to be used for unique identification of each id.  starts at 0 unless one is passed in
    my $cds   = $opts{cds}; #flag to only print protein coding genes
    my $name_unique = $opts{name_unique}; #flag for making Name tag of output unique by appending type and occurrence to feature name
    my $unique_parent_annotations = $opts{unique_parent_annotations}; #flag so that annotations are not propogated to children if they are contained by their parent
    my $id_type  = $opts{id_type};  #type of ID (name, num):  unique number; unique name
    my $cds_exon = $opts{cds_exon}; #option so that CDSs are used for determining an exon instead of the mRNA.  This keeps UTRs from being called an exon
    my $chromosome = $opts{chr};    #optional, set to only include features on a particular chromosome
    my $add_chr = $opts{add_chr};   #option to add "chr" before chromosome name
    $id_type = "name" unless defined $id_type;
    $count = 0 unless ($count && $count =~ /^\d+$/);
    $ds = $self unless $ds;

    # mdb added 7/28/15 - performance improvement
    my $dbh = $self->{_result_source}{schema}{storage}->dbh;
    print STDERR "Pre-caching features\n" if $debug;
    $FEATURES      = get_features($dbh, undef, $self->id);
    print STDERR "Pre-caching feature names\n" if $debug;
    $FEATURE_NAMES = get_feature_names($dbh, undef, $self->id);
    print STDERR "Pre-caching feature annotations\n" if $debug;
    $FEATURE_ANNOS = get_feature_annotations($dbh, undef, $self->id);
    print STDERR "Pre-caching feature locations\n" if $debug;
    $FEATURE_LOCS  = get_locations($dbh, undef, $self->id);
    print STDERR "Indexing features by name\n" if $debug;
    foreach (values %$FEATURE_NAMES) {
        foreach my $fname (values %$_) {
            my $name = $fname->{name};
            my $fid  = $fname->{fid};
            my $chr  = $fname->{chr};
            $FEATURES_BY_NAME{$chr}{$name}{$fid} = 1;
        }
    }

    # Generate GFF header
    print STDERR "Generating GFF header\n" if $debug;
# mdb removed 8/26/15 -- performance improvement
#    my %chrs;
#    foreach my $chr ( $ds->get_chromosomes() ) {
#        $chrs{$chr} = $ds->last_chromosome_position($chr);
#    }
    my $chrs = $ds->first_genome->chromosome_lengths_by_name(); # mdb added 8/26/15 -- performance improvement
    my @chrs = keys %FEATURES_BY_NAME;
    @chrs = ($chromosome) if ($chromosome);
    
    my $tmp;
    $tmp = "##gff-version\t3\n" unless $no_gff_head;
    $tmp .= "##CoGe Dataset ID: ".$ds->id."\n";
    $tmp .= "##CoGe Dataset Name: ". $ds->name."\n" if $ds->name;
    my $desc = $ds->desc;
    $desc =~ s/\n/ /g if $desc;
    $tmp .= "##CoGe Dataset Desc: ". $desc."\n" if $desc;
    my $output;
    $output .= $tmp if $tmp;
    print $tmp if $print && !$no_gff_head;
    foreach my $chr (@chrs) {
        next unless $chrs->{$chr};
        $tmp = "##sequence-region $chr 1 " . $chrs->{$chr} . "\n";
        $output .= $tmp;
        print $tmp if $print;
    }

    # Generate GFF entries
    print STDERR "Generating GFF entries\n" if $debug;
    my %fids = ();  #skip fids that we have processed
    my %types;      #track the number of different feature types encountered
    my %ids2names;  #lookup table for unique id numbers to unique id names (determined by $id_type)
    my %unique_ids; #place to make sure that each ID used is unique;
    my %prev_annos; #hash to store previously used annotations by parents.  Used in conjunction with the $unique_parent_annotations flag
    my %prior_genes;  #place to store previous gene models for alternatively spliced transcript lookup of parents.  keys is the primary name of the transcript.  value is coge feature object
    my %orphaned_RNA; #place to store RNAs without a parent feature
    #my $prior_gene;  #place to store the prior gene, if needed.  Some alternatively spliced transcripts have one gene per transcript, some have one gene for all transcripts.
    #my $prior_gene_id; #place to store the prior gene's ID for use the GFF file.  Use of this is tied to having and using a prior_gene

    #my $prefetch = [ 'feature_type', 'feature_names' ];
    #push @$prefetch, {'annotations' => 'annotation_type'} if $annos;
    foreach my $chr (sort @chrs) {
        my %seen = (); #for storage of seen names organized by $feat_name{$name}
#        my $rs_feat = $ds->features(
#            { chromosome => $chr }, # mdb added 4/23/14 issue 364
#            {   'prefetch' => $prefetch,
#                'order_by' => [ 'me.chromosome', 'me.start', 'me.feature_type_id'
#                  ] #go by order in genome, then make sure that genes (feature_type_id == 1) is first
#            });
        my $rs_feat = [ sort feature_sort values %{$FEATURES->{$chr}} ];

        #main: while ( my $feat = $rs_feat->next ) {
        main: foreach my $feat (@$rs_feat) {
	        print STDERR "Processing feature: ", $feat->{type_name}, " id=", $feat->{id}, " chr=", $feat->{chromosome}, " start=", $feat->{start}, "\n" if $debug;
            next if ( $fids{ $feat->{id} } );
            next unless ($feat->{start} && defined $feat->{chromosome});
            
            #my $ft = $feat->feature_type;
            #print STDERR join ("\t", $ft->name, $feat->chromosome, $feat->start),"\n" if $debug;
            $types{ $feat->{type_name} }++;
            $count++;
            
            my @out;   #store hashref of items for output
            my %notes; #additional annotations to add to gff line
            my @feat_names = _feat_names($feat->{id});
            if ($name_re) {
                @feat_names = grep { $_ =~ /$name_re/i } @feat_names;#$feat->names;
                next unless @feat_names;
            }
             
            if ($feat->{type_name} =~ /gene/i && !$prior_genes{ $feat_names[0] }) {
                $prior_genes{ $feat_names[0] } = $feat;
		        $prior_genes{$count} = $feat;
		    }
            if ( $feat->{type_name} =~ /RNA/ ) { #perhaps an alternatively spliced transcript
                #check for name congruence
                my $match = 0;
                my $prior_gene;
                my $prior_gene_id;
                name_search: foreach my $name ( _feat_names($feat->{id}) ) { #$feat->names ) {
                    if ( $prior_genes{$name} ) {
                        $match         = 1;
                        $prior_gene    = $prior_genes{$name};
                        $prior_gene_id = $name;
                        last;
                    }
                }
                if ($match) {
                    if ($self->_search_rna(
                            name_search => [ _feat_names($prior_gene->{id}) ], #[ $prior_gene->names ],
                            notes       => \%notes,
                            fids        => \%fids,
                            types       => \%types,
                            count       => \$count,
                            out         => \@out,
                            name_re     => $name_re,
                            parent_feat => $prior_gene,
                            parent_id   => $prior_gene_id,
                            chr         => $chr,
                            ds          => $ds,
                            cds_exon    => $cds_exon,
			                prior_genes => \%prior_genes))
                    {
                        my $tmp = $self->_format_gff_line(
                            out         => \@out,
                            notes       => \%notes,
                            cds         => $cds,
                            seen        => \%seen,
                            print       => $print,
                            annos       => $annos,
                            name_unique => $name_unique,
                            ids2names   => \%ids2names,
                            id_type     => $id_type,
                            unique_ids  => \%unique_ids,
                            unique_parent_annotations => $unique_parent_annotations,
                            prev_annos  => \%prev_annos,
			                add_chr     => $add_chr
                        );
                        $output .= $tmp if $tmp;
                        next main;
                    }
                }
                else {
                    $orphaned_RNA{ $feat->{id} } = $feat;
                    next main;
                }
            }

            foreach my $loc ( _feat_locs($feat) ) { #$feat->locs() ) {
                push @out,
                  {
                    f           => $feat,
                    start       => $loc->{start},
                    stop        => $loc->{stop},
                    name_re     => $name_re,
                    id          => $count,
                    parent_id   => 0,
                    type        => $feat->{type_name},
		            prior_genes => \%prior_genes
                  };
                $count++;
            }
            $fids{ $feat->{id} } = 1;    #feat_id has been used

#			unless ( $ft->id == 1 && @feat_names )    #if not a gene, don't do the next set of searches.
            unless ( $feat->{type_name} =~ /gene/i && @feat_names ) {
                my $tmp = $self->_format_gff_line(
                    out                       => \@out,
                    notes                     => \%notes,
                    cds                       => $cds,
                    seen                      => \%seen,
                    print                     => $print,
                    annos                     => $annos,
                    name_unique               => $name_unique,
                    ids2names                 => \%ids2names,
                    id_type                   => $id_type,
                    unique_ids                => \%unique_ids,
                    unique_parent_annotations => $unique_parent_annotations,
                    prev_annos                => \%prev_annos,
		            add_chr                   => $add_chr
                );
                $output .= $tmp if $tmp;
                next;
            }

            #does this gene have an RNA?
            if ($self->_search_rna(
                    name_search => \@feat_names,
                    notes       => \%notes,
                    fids        => \%fids,
                    types       => \%types,
                    count       => \$count,
                    out         => \@out,
                    name_re     => $name_re,
                    parent_feat => $feat,
                    chr         => $chr,
                    ds          => $ds,
                    cds_exon    => $cds_exon,
		            prior_genes => \%prior_genes
                ))
            {
                my $tmp = $self->_format_gff_line(
                    out                       => \@out,
                    notes                     => \%notes,
                    cds                       => $cds,
                    seen                      => \%seen,
                    print                     => $print,
                    annos                     => $annos,
                    name_unique               => $name_unique,
                    ids2names                 => \%ids2names,
                    id_type                   => $id_type,
                    unique_ids                => \%unique_ids,
                    unique_parent_annotations => $unique_parent_annotations,
                    prev_annos                => \%prev_annos,
		            add_chr 		          => $add_chr
                );
                $output .= $tmp if $tmp;
                next main;
            }
            #dump other stuff for gene that does not have mRNAs.
            my $sub_rs = $self->_feat_search(
                name_search => \@feat_names,
                skip_ftids  => [1, 2], #1 == gene; 2 == mRNA
                #ds          => $ds,
                chr         => $chr,
            );
            
            my $parent_id = $count - 1;
            #while ( my $f = $sub_rs->next() ) {
            foreach my $f (@$sub_rs) {
                next if ( $fids{ $f->{id} } );
		        
                my $ftn = $self->process_feature_type_name( $f->{type_name} );
                push @{ $notes{gene}{"Encoded_feature"} }, $self->escape_gff($ftn);
                foreach my $loc ( _feat_locs($f) ) { #$f->locs() ) {
                    # Need to do some special stuff for such cases where a CDS retrieved in the absense of an mRNA

                    if ( $f->{type_name} eq "CDS" ) {
                        #let's add the mRNA, change the parent and count (id)
                        push @out, {
                            f         => $f,
                            start     => $loc->{start},
                            stop      => $loc->{stop},
                            name_re   => $name_re,
                            id        => $count,
                            parent_id => $parent_id,
                            type      => "mRNA",
		    	            prior_genes => \%prior_genes
                        };
                        $parent_id = $count;
                        $count++;
                    }
                    push @out, {
                        f         => $f,
                        start     => $loc->{start},
                        stop      => $loc->{stop},
                        name_re   => $name_re,
                        id        => $count,
                        parent_id => $parent_id,
                        type      => "exon",
		                prior_genes => \%prior_genes
                    } unless $f->{type_name} =~ /pseudogene/;
                    $count++;
                    push @out, {
                        f         => $f,
                        start     => $loc->{start},
                        stop      => $loc->{stop},
                        name_re   => $name_re,
                        id        => $count,
                        parent_id => $parent_id,
                        type      => $f->{type_name},
		                prior_genes => \%prior_genes
                    };
                    $fids{ $f->{id} } = 1;    #feat_id has been used;
                    $types{ $f->{type_name} }++;
                }
                my $tmp = $self->_format_gff_line(
                    out                       => \@out,
                    notes                     => \%notes,
                    cds                       => $cds,
                    seen                      => \%seen,
                    print                     => $print,
                    annos                     => $annos,
                    name_unique               => $name_unique,
                    ids2names                 => \%ids2names,
                    id_type                   => $id_type,
                    unique_ids                => \%unique_ids,
                    unique_parent_annotations => $unique_parent_annotations,
                    prev_annos                => \%prev_annos,
		            add_chr                   => $add_chr
                );
                $output .= $tmp if $tmp;
                last;
            }
        }
    }
    if ($print) {
        print "#Orphaned RNAs\n";
        foreach my $fid ( sort keys %orphaned_RNA ) {
            next if $fids{$fid};
            print "#", join( "\t", $fid, $orphaned_RNA{$fid}->names ), "\n";
        }
    }
    return $output, $count;
}

sub _feat_names {
    my $fid = shift;
    my (@primary_names, @secondary_names);
    foreach my $feature_name ( values %{$FEATURE_NAMES->{$fid}} ) {
        if ( $feature_name->{primary_name} ) {
            push @primary_names, $feature_name->{name};
        }
        else {
            push @secondary_names, $feature_name->{name};
        }
    }
    my @names = ( sort @primary_names, sort @secondary_names );
    return wantarray ? @names : \@names;
}

sub _feat_locs {
    my $feature = shift;
    my @locs;
    foreach my $loc (values %{$FEATURE_LOCS->{$feature->{id}}}) {
        next if ($loc->{strand} ne $feature->{strand});
        next if ($loc->{chr} ne $feature->{chromosome});
        next if ($loc->{start} < $feature->{start} || $loc->{start} > $feature->{stop});
        next if ($loc->{stop} < $feature->{start} || $loc->{stop} > $feature->{stop});
        push @locs, $loc;
    }
    my @sorted = sort { $a->{start} <=> $b->{start} || $a->{lid} <=> $b->{lid} } @locs;
    return wantarray ? @sorted : \@sorted;
}

sub _feat_annos {
    my $fid = shift;
    #print STDERR "_feat_annos\n";
    my @annos = sort { $a->{faid} <=> $b->{faid} } values %{$FEATURE_ANNOS->{$fid}};
    return wantarray ? @annos : \@annos;
}

sub _feat_search {
    my $self        = shift;
    my %opts        = @_;
    my $name_search = $opts{name_search};
    my $skip_ftids  = $opts{skip_ftids};
    #my $ds          = $opts{ds};
    my $chr         = $opts{chr};
    #print STDERR "_feat_search\n";
    
    my %feats;
    foreach my $name (@$name_search) {
        foreach my $fid (keys %{$FEATURES_BY_NAME{$chr}{$name}}) {
            my $f = $FEATURES->{$chr}{$fid};
            next if (grep { $f->{type_id} eq $_ } @$skip_ftids);
            $feats{$fid} = $f;
        }
    }
    
    return [ sort feature_sort values %feats ];
    
#    return $ds->features(
#        {
#            'me.chromosome'      => $chr,
#            'feature_names.name' => { 'IN' => $name_search },
#            'me.feature_type_id' => { 'NOT IN' => $skip_ftids },
#        },
#        {
#            'join'     => 'feature_names',
#            'prefetch' => [ 'feature_type', 'locations' ],
#            'order_by' => [ 'me.start', 'locations.start', 'me.feature_type_id' ]
#        }
#    );
}

sub _search_rna {
    my $self        = shift;
    my %opts        = @_;
    my $name_search = $opts{name_search};
    my $notes       = $opts{notes};
    my $fids        = $opts{fids};
    my $types       = $opts{types};
    my $count       = $opts{count};
    my $out         = $opts{out};
    my $parent_feat = $opts{parent_feat};
    my $parent_id   = $opts{parent_id};
    my $name_re     = $opts{name_re};
    my $chr         = $opts{chr};
    my $ds          = $opts{ds};
    my $cds_exon    = $opts{cds_exon}; #CDSs are used for determining an exon instead of the mRNA.  This keeps UTRs from being called an exon
    my $prior_genes = $opts{prior_genes};
    #print STDERR "_search_rna\n";

    my $rna_rs = $self->_feat_search(
        name_search => $name_search,
        skip_ftids  => [1], # gene type
        #ds          => $ds,
        chr         => $chr
    );

    #assemble RNA info
    #while ( my $f = $rna_rs->next() ) {
    foreach my $f (@$rna_rs) {
        next if ( $fids->{ $f->{id} } );
        next unless $f->{type_name} =~ /RNA/i; #searching for feat_types of RNA
        
        #process the RNAs
        $parent_id = $self->_process_rna(
            notes       => $notes,
            fids        => $fids,
            types       => $types,
            count       => $count,
            out         => $out,
            f           => $f,
            name_re     => $name_re,
            parent_feat => $parent_feat,
            parent_id   => $parent_id,
            cds_exon    => $cds_exon,
	        prior_genes => $prior_genes
        );

        #get CDSs (mostly)
        my $sub_rs = $self->_feat_search(
            name_search => [ _feat_names($f->{id}) ], #[ $f->names ],
            skip_ftids  => [ 1, $f->{type_id} ],
            #ds          => $ds,
            chr         => $chr,
        );

        my %tmp_types;    #only want to process one of each type.
        #while ( my $f = $sub_rs->next() ) {
        foreach my $f (@$sub_rs) {
            next if ( $fids->{ $f->{id} } );
            next if $tmp_types{ $f->{type_name} };
            
            $tmp_types{ $f->{type_name} }++;
            my $ftn = $self->process_feature_type_name( $f->{type_name} );
            $fids->{ $f->{id} } = 1;    #feat_id has been used;
            $types->{ $f->{type_name} }++;
            
	        foreach my $loc ( _feat_locs($f) ) { #$f->locs() ) {
                push @$out, {
                    f         => $f,
                    start     => $loc->{start},
                    stop      => $loc->{stop},
                    name_re   => $name_re,
                    id        => $$count,
                    parent_id => $parent_id,
                    type      => $f->{type_name},
		            prior_genes => $prior_genes
                };
                $$count++ if $cds_exon;
                push @$out, {
                    f         => $f,
                    start     => $loc->{start},
                    stop      => $loc->{stop},
                    name_re   => $name_re,
                    id        => $$count,
                    parent_id => $parent_id,
                    type      => "exon",
		            prior_genes => $prior_genes,
                }
                if $cds_exon;
                $$count++;
            }
        }
        return 1;    #return after doing this once
    }
    return 0;
}

sub _process_rna {
    my $self        = shift;
    my %opts        = @_;
    my $notes       = $opts{notes};
    my $fids        = $opts{fids};
    my $types       = $opts{types};
    my $count       = $opts{count};
    my $out         = $opts{out};
    my $f           = $opts{f};
    my $parent_feat = $opts{parent_feat};
    my $parent_id   = $opts{parent_id};
    my $name_re     = $opts{name_re};
    my $prior_genes = $opts{prior_genes};
    my $cds_exon = $opts{cds_exon}; #option so that CDSs are used for determining an exon instead of the mRNA.  This keeps UTRs from being called an exon
    #print STDERR "_process_rna\n";
    
    my $ftn = $self->process_feature_type_name( $f->{type_name} );
    push @{ $notes->{gene}{"encoded_feature"} }, $self->escape_gff($ftn);
    $fids->{ $f->{id} } = 1;    #feat_id has been used;
    $types->{ $f->{type_name} }++;

    #have mRNA.  mRNA in CoGe translates to what most people have settled on calling exons.  the output mRNA therefore needs to be a replicate of the gene
    $parent_id = $$count - 1 unless $parent_id;
    push @$out, {
        f         => $f,
        start     => $f->{start},
        stop      => $f->{stop},
        name_re   => $name_re,
        id        => $$count,
        parent_id => $parent_id,
        type      => $f->{type_name},
	    prior_genes => $prior_genes
    };    #need to add a mRNA from its start to stop, no exons/introns
    
    #dump exons for mRNA
    $parent_id = $$count;
    $$count++;
    foreach my $loc ( _feat_locs($f) ) { #$f->locs() ) {
        unless ($cds_exon) {
            push @$out, {
                f         => $f,
                start     => $loc->{start},
                stop      => $loc->{stop},
                name_re   => $name_re,
                id        => $$count,
                parent_id => $parent_id,
                type      => "exon"
            };
            $$count++;
        }
    }
    
    return $parent_id;
}

sub _format_gff_line {
    my $self  = shift;
    my %opts  = @_;
    my $out   = $opts{out};   #array of hashref of output items
    my $notes = $opts{notes}; #hashref of notes; keyed by feature type
    my $cds   = $opts{cds};   #only print CDS genes
    my $print = $opts{print}; #are lines printed here?
    my $annos = $opts{annos}; #are annotations retrieved?
    my $unique_parent_annos = $opts{unique_parent_annotations}; #parent annotations are NOT propogated to children
    my $prev_annos  = $opts{prev_annos};  #hash for storing previously seen annotations by parents and children -- used in conjuction with $unique_parent_annos flag
    my $seen        = $opts{seen};        #general var for checking if a simlar feature has been seen before (looked up by type and name string)
    my $ids2names   = $opts{ids2names};   #hash to looking up names for a particular id.  May be a number (same as the id) or a unique name;
    my $id_type     = $opts{id_type};     #type of ID (name, num):  unique number; unique name
    my $name_unique = $opts{name_unique}; #flag for making Name tag of output unique by appending type and occurrence to feature name
    my $unique_ids  = $opts{unique_ids};  #hash for making sure that each used ID happens once for each ID
    my $add_chr     = $opts{add_chr};     #flag to add "chr" before the chromosome
    #print STDERR "_format_gff_line\n";
    
    my $output;
    foreach my $item (@$out) {
        my $f         = $item->{f};         #feature object
        my $type      = $item->{type};      #feature type.  may be retrieved from feature object, but these may differ
        my $name_re   = $item->{name_re};   #regex for searching for a specific name
        my $id        = $item->{id};        #unique id for the gff feature;
        my $parent_id = $item->{parent_id}; #unique id for the parent of the gff feature
        my $start     = $item->{start};     #start of entry.  May be the feat start, may be a loc start.  Need to declare in logic outside of this routine
        my $stop      = $item->{stop};      #stop of entry.  May be the feat stop, may be a loc stop.  Need to declare in logic outside of this routine
	    my $prior_genes = $item->{prior_genes}; #ref to genes list.  Used to look up a previous gene if needed.
        my $parsed_type = $self->process_feature_type_name($type);
        
        #check to see if we are only printing CDS genes
        return
          if $cds
              && (   $type ne "gene"
                  && $type ne "mRNA"
                  && $type ne "CDS"
                  && $type ne "exon" );

        my @feat_names = _feat_names($f->{id});#$f->names;
        if ($name_re) {
            @feat_names = grep { $_ =~ /$name_re/i } @feat_names;
            next unless @feat_names;
        }
        my ($alias) = join( ",", map { $self->escape_gff($_) } @feat_names );
        my ($name) = @feat_names;
        $name = $parsed_type unless $name;
        $seen->{$parsed_type}{$name}++;

        #create a unique name for the feature type given the name of the feature
        my $unique_name = $name;
        $unique_name .= "." . $parsed_type . $seen->{$parsed_type}{$name}
          if $unique_ids->{$unique_name};

        #store the unqiue name and associate it with the unique ID number
        warn "ERROR!  $id is already in use in \$ids2names lookup table: " . $ids2names->{$id} . " fid: ". $f->{id}."\n"
            if ( $ids2names->{$id} );
        $ids2names->{$id} = $unique_name;

        #if unique names are requested for the Name tag, use it
        $name = $unique_name if ($name_unique);

        #for the feature and parent ids, are we using names or numbers?
        if ( $id_type eq "name" ) {
            $id        = $ids2names->{$id}        if $ids2names->{$id};
            $parent_id = $ids2names->{$parent_id} if $ids2names->{$parent_id};
        }
        warn "ERROR:  ID $id has been previously used!" if ( $unique_ids->{$id} );
        $unique_ids->{$id}++;

    	#need to check and correct for an edge case when item is something like a miRNA, and parent is defined, and parent is not a gene.  EHL: 6/2/15
    	my $is_RNA = $parsed_type =~ /RNA/i ? 1 : 0;
    	my $is_mRNA = $parsed_type =~ /mRNA/i ? 1 : 0;
    	my $parent_is_gene = $prior_genes->{$parent_id}{type_name} =~ /gene/ ? 1 : 0 if $prior_genes->{$parent_id};
    	$parent_is_gene = 0 unless defined $parent_is_gene;
    	if ($parent_id && $is_RNA && !$is_mRNA &&! $prior_genes->{$parent_id}  && !$parent_is_gene) {
    	    my $temp_id = $parent_id;
    	    while (!$parent_is_gene && $temp_id !~ /[A-Za-z]/ && $temp_id > 0) { # mdb added $temp_id > 0 8/25/15 to prevent infinite loop
    		    $temp_id--;
    		    $parent_is_gene = $prior_genes->{$temp_id}{type_name} =~ /gene/ ? 1 : 0 if $prior_genes->{$temp_id};
    	    }
    	    $parent_id = $temp_id;
    	}
    	  
        my $attrs;
        my $test_var = $type =~ /gene/ ? 1 : 0;
        $attrs .= "Parent=$parent_id;" if $parent_id && $test_var == 0;#!$type=~/gene/i;
        $attrs .= "ID=$id";
        $attrs .= ";Name=$name"        if $name;
        $attrs .= ";Alias=$alias"      if $alias;
        $attrs .= ";$parsed_type=$name" if $name;
        $attrs .= ";coge_fid=" . $f->{id} . "";
        foreach my $key ( sort keys %{ $notes->{$type} } ) {
            $attrs .= ";$key=" . join( ",", sort @{ $notes->{$type}{$key} } );
        }
    
        my $anno_stuff;
        if ($annos) {
            my %annos;
            foreach my $anno ( _feat_annos($f->{id}) ) { #$f->annotations ) {
                next unless defined $anno->{annotation};

                if ( $unique_parent_annos && $prev_annos->{$parent_id}{$anno->{annotation} } ) {
                    #we have used this annotation in a parent annotation
                    $prev_annos->{$id}{ $anno->{annotation} } = 1;
                    next;
                }
                if ($unique_parent_annos) {
                    $prev_annos->{$id}{ $anno->{annotation} } = 1;
                    $prev_annos->{$parent_id}{ $anno->{annotation} } = 1;
                }
                
                my $atn;      #attribute name
                my $value;    #annotation;
    
                #if there is a group it should be the attr name otherwise it is the anno type
                if ($anno->{atgname}) { #$anno_type->annotation_type_group ) {
                    $atn = $anno->{atgname}; #$anno_type->annotation_type_group->name;
                    $value .= $anno->{atname} . ', '; #$anno_type->name . ", ";
                }
                else {
                    $atn = $anno->{atname}; #$anno_type->name;
                }
                $atn = $self->escape_gff($atn);
                $atn =~ s/\s+/_/g;
                $atn = 'Note' unless $atn;
                $value .= $anno->{annotation};
                $value = $self->escape_gff($value);
                $annos{$atn}{$value} = 1;
            }
    
            foreach my $key ( sort keys %annos ) {
                $anno_stuff .= $key . '=' . join( ',', sort keys %{ $annos{$key} } ) . ';';
            }
        }

        #assemble GFF line for printing
	    my $chr = $f->{chromosome};
	    $chr = 'chr' . $chr if $add_chr;
        my $anno_str = join(
            "\t",
            (   $chr, 'CoGe', $parsed_type, $start, $stop, '.',
                ($f->{strand} == 1 ? '+' : '-'), 
                '.', $attrs
            )
        );
        $anno_str .= ";$anno_stuff" if $anno_stuff;
        $anno_str =~ s/;$//g;
        $output .= $anno_str . "\n";
    }
    
    print $output if ($print and $output);
    return $output;
}

sub process_feature_type_name {
    my $self = shift;
    my $ftn  = shift;
    $ftn = lc($ftn);
    $ftn =~ s/\s+/_/;
    $ftn =~ s/rna/RNA/i;
    $ftn =~ s/cds/CDS/i;
    return $ftn;
}

sub escape_gff {
    my $self = shift;
    my $tmp  = shift;

    #escape for gff
    $tmp =~ s/\t/\%09/g;
    $tmp =~ s/\n/\%0A/g;
    $tmp =~ s/\r/\%0D/g;
    $tmp =~ s/;/\%3B/g;
    $tmp =~ s/=/\%3D/g;
    $tmp =~ s/\%/\%25/g;
    $tmp =~ s/&/\%26/g;
    $tmp =~ s/,/\%2C/g;
    return $tmp;
}

sub trans_type {
    my $self = shift;
    my $trans_type;
    foreach
      my $feat ( $self->features( { feature_type_id => 3 }, { rows => 10 } ) )
    {

        #	next unless $feat->type->name =~ /cds/i;
        my ( $code, $type ) = $feat->genetic_code;
        ($type) = $type =~ /transl_table=(\d+)/ if $type =~ /transl_table/;
        return $type if $type;
    }
    return 1;    #universal genetic code type;
}

sub distinct_feature_type_ids {
    my $self = shift;
    my %ids = map { $_->feature_type_id => 1 } $self->features(
        {},
        {
            columns  => ['feature_type_id'],
            distinct => 1
        }
    );
    return wantarray ? keys %ids : [ keys %ids ];
}

sub distinct_feature_type_names {
    my $self = shift;
    my %names = map { $_->feature_type->name => 1 } $self->features(
        {},
        {
            join     => 'feature_type',
            columns  => ['feature_type.name'],
            distinct => 1
        }
    );
    return wantarray ? keys %names : [ keys %names ];
}

sub translation_type {
    my $self = shift;
    foreach my $feat (
        $self->features(
            { 'annotation_type_id' => 10973 },
            { join                 => 'feature_annotations' }
        )
      )
    {
        foreach
          my $anno ( $feat->annotations( { annotation_type_id => 10973 } ) )
        {
            return $anno->annotation;
        }
    }
}

sub user_groups() {
    my $self = shift;

    my @groups = ();
    foreach ( $self->user_group_data_connectors() ) {
        push( @groups, $_->user_group() );
    }

    return @groups;
}

1;
