package CoGeX::Result::Genome;

use strict;
use warnings;
use base 'DBIx::Class::Core';

use CoGe::Core::Chromosomes;
use CoGe::Core::Storage qw( get_genome_seq get_genome_file );
use CoGe::Accessory::Utils qw( commify );
use JSON qw(encode_json);
use Data::Dumper;
use Text::Wrap;
use base 'Class::Accessor';
use Carp;

use constant ERROR   => -1;
use constant LOADING => 1;
use constant LOADED  => 2;

BEGIN {
	use Exporter 'import';
	our @EXPORT = qw( 
	   ERROR LOADING LOADED
	);
}

=head1 NAME

CoGeX::Genome

=head1 SYNOPSIS

This object uses the DBIx::Class to define an interface to the C<genome> table in the CoGe database.

=head1 DESCRIPTION

=head1 USAGE

 use CoGeX;

=head1 METHODS

=head1 AUTHORS

 Eric Lyons
 Brent Pedersen
 Daniel Hembry
 Matt Bomhoff

=head1 COPYRIGHT

This program is free software; you can redistribute
it and/or modify it under the same terms as Perl itself.

The full text of the license can be found in the
LICENSE file included with this module.

=head1 SEE ALSO

=cut

my $node_types = CoGeX::node_types();

__PACKAGE__->table("genome");
__PACKAGE__->add_columns(
    "genome_id",
    {
        data_type     => "INT",
        default_value => undef,
        is_nullable   => 0,
        size          => 11
    },
    "name",
    {
        data_type     => "VARCHAR",
        default_value => "",
        is_nullable   => 0,
        size          => 200
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
        is_nullable   => 0,
        size          => 50,
    },
    "organism_id",
    { data_type => "INT", default_value => 0, is_nullable => 0, size => 11 },
    "genomic_sequence_type_id",
    { data_type => "INT", default_value => 0, is_nullable => 0, size => 11 },
    "restricted",
    { data_type => "INT", default_value => "0", is_nullable => 0, size => 1 },
    "message",
    {
        data_type     => "text",
        default_value => undef,
        is_nullable   => 1,
    },
    "link",
    {
        data_type     => "text",
        default_value => undef,
        is_nullable   => 1,
    },
    "deleted",
    { data_type => "INT", default_value => "0", is_nullable => 0, size => 1 },
    "certified", # mdb added 7/27/16 COGE-536
    { data_type => "INT", default_value => "0", is_nullable => 0, size => 1 },
    "creator_id",
    { data_type => "INT", default_value => 0, is_nullable => 0, size => 11 },
    "date",
    { data_type => "TIMESTAMP", default_value => undef, is_nullable => 0 },
    "status",
    { data_type => "INT", default_value => undef, is_nullable => 1, size => 1 },
);

__PACKAGE__->set_primary_key("genome_id");
__PACKAGE__->has_many(
    "dataset_connectors" => "CoGeX::Result::DatasetConnector",
    'genome_id'
);
__PACKAGE__->belongs_to(
    "organism" => "CoGeX::Result::Organism",
    'organism_id'
);
__PACKAGE__->belongs_to(
    "genomic_sequence_type" => "CoGeX::Result::GenomicSequenceType",
    'genomic_sequence_type_id'
);
__PACKAGE__->belongs_to(
    "creator" => "CoGeX::Result::User", 
    { 'foreign.user_id' => 'self.creator_id' }
);
__PACKAGE__->has_many(
    "genome_annotations" => "CoGeX::Result::GenomeAnnotation",
    'genome_id'
);
__PACKAGE__->has_many(    # parent lists
    'list_connectors' => 'CoGeX::Result::ListConnector',
    { 'foreign.child_id' => 'self.genome_id' },
    { where => [ -and => [ child_type => $node_types->{genome} ] ] }
);
__PACKAGE__->has_many(    # parent users
    'user_connectors' => 'CoGeX::Result::UserConnector',
    { "foreign.child_id" => "self.genome_id" },
    {
        where => [
            -and => [
                parent_type => $node_types->{user},
                child_type  => $node_types->{genome}
            ]
        ]
    }
);
__PACKAGE__->has_many(    # parent groups
    'group_connectors' => 'CoGeX::Result::UserConnector',
    { "foreign.child_id" => "self.genome_id" },
    {
        where => [
            -and => [
                parent_type => $node_types->{group},
                child_type  => $node_types->{genome}
            ]
        ]
    }
);
__PACKAGE__->has_many(
    "experiments" => "CoGeX::Result::Experiment",
    "genome_id"
);
__PACKAGE__->has_many(
    "favorite_connectors" => "CoGeX::Result::FavoriteConnector",
    { "foreign.child_id" => "self.genome_id" },
    { where => [ -and => [ child_type  => $node_types->{genome} ] ] }
);

__PACKAGE__->mk_accessors('_chromosomes');

sub page {
    return 'GenomeInfo.pl';
}

sub object_type {
    return 'genome';
}

sub item_type {
    return $node_types->{genome};   
}

sub desc {
    return shift->description(@_);
}

sub annotations {
    shift->genome_annotations(@_);
}

################################################

sub to_hash {
	my $self = shift;
	return {
		id            => $self->id,
		type          => 'genome',
		name          => $self->name,
		description   => $self->description,
		version       => $self->version,
		organism      => $self->organism->name,
		sequence_type => $self->genomic_sequence_type->name,
		restricted    => $self->restricted,
		message       => $self->message,
		link          => $self->link,
		deleted       => $self->deleted,
		certified     => $self->certified,
        status        => $self->status
	}
}

sub to_json {
	return encode_json( shift->to_hash );
}

################################################ subroutine header begin ##

=head2 lists

 Usage     : $self->lists
 Purpose   : Get lists that contain this genome
 Returns   : array or array ref of parent lists
 Argument  :
 Throws    :
 Comments  :

See Also   :

=cut

################################################## subroutine header end ##

sub lists {
    my $self = shift;

    my %lists;
    foreach my $conn ( $self->list_connectors ) {
        my $list = $conn->parent_list;
        $lists{ $list->id } = $list if ($list);
    }

    return wantarray ? values %lists : [ values %lists ];
}

sub notebooks {
    shift->lists(@_);
}

sub notebooks_desc {
    my $self = shift;
    return join(',', map { $_->name } $self->notebooks) || '';
}

##################################################

sub datasets {
    my $self = shift;
    my %opts = @_;
    my $chr  = $opts{chr};
    $chr = $opts{chromosome} unless defined $chr;
    my $restricted = $opts{restricted}; #FIXME why is this here? restricted datasets are not used

    my %datasets;
    foreach my $dsc ( $self->dataset_connectors() ) {
        my $ds = $dsc->dataset;
        next unless $ds;
        next if ( $restricted and not $ds->restricted );

        if ( defined $chr ) {
            if ( $ds->has_chromosome(chr => $chr) ||
                 $ds->has_chromosome(chr => 'lcl|'.$chr) || # mdb added 8/18/16 for COGE-735
                 $ds->has_chromosome(chr => 'gi|'.$chr) )   # mdb added 8/18/16 for COGE-735
            {
                $datasets{ $ds->id } = $ds;
            }
        }
        else {
            $datasets{ $ds->id } = $ds;
        }
    }
    return wantarray ? values %datasets : [ values %datasets ];
}

sub first_dataset {
    my $self = shift;

    foreach my $dsc ( $self->dataset_connectors() ) {
        return $dsc->dataset;
    }
}

################################################ subroutine header begin ##
#
#=head2 owner
#
#Usage     :
#Purpose   : Returns user object
#Returns   : user object
#Argument  : None
#Throws    : None
#Comments  :
#
#=cut
#
################################################### subroutine header end ##

sub owner {
    my $self = shift;

    foreach ( $self->user_connectors( { role_id => 2 } ) ) {    #FIXME hardcoded
        return $_->parent;
    }
}

################################################ subroutine header begin ##

=head2 get_chromosome

 Usage     :
 Purpose   : replace DBIX one-to-many relation with Chromosomes class
 Returns   :
 Argument  : name
 Throws    :
 Comments  :

See Also   :

=cut

################################################## subroutine header end ##

sub get_chromosome {
    my $self = shift;
    my $name = shift;
	my $c = CoGe::Core::Chromosomes->new($self->id);
	if ($c->find($name)) {
		return (chromosome=>$c->name, sequence_length=>$c->length);
	}
	print STDERR "CoGeX::Result::Genome::get_chromosome ERROR, chromosome '$name' not found\n";
	return 0;
}

################################################ subroutine header begin ##

=head2 get_chromosome_length

 Usage     :
 Purpose   : replace DBIX one-to-many relation with Chromosomes class
 Returns   :
 Argument  : name
 Throws    :
 Comments  :

See Also   :

=cut

################################################## subroutine header end ##

sub get_chromosome_length {
    my $self = shift;
    my $name = shift;
	my $c = CoGe::Core::Chromosomes->new($self->id);
	return $c->length if $c->find($name);
	warn "CoGeX::Result::Genome::get_chromosome_length ERROR, chromosome '$name' not found for genome ", $self->id;
	return 0;
}

################################################

sub get_genomic_sequence {
    my $self  = shift;
    my %opts  = @_;
    my $start = $opts{start} || $opts{begin};
    my $stop  = $opts{stop} || $opts{end};
    my $chr   = $opts{chr};
    $chr = $opts{chromosome} unless defined $chr;
    $chr = "1" unless defined $chr;
    my $strand = $opts{strand};
    return
      if ( defined $start && defined $stop && $start < 1 && $stop < 1 )
      ;    #asking for sequence beyond the start
    $start = 1 unless $start;
#    my $last = $self->sequence_length($chr);
    my $last = $self->get_chromosome_length($chr);
    return 0 unless $last; # chromosome not found
    $stop = $last unless $stop;
    $stop  = $last if $stop > $last;
    $start = $last if $start > $last;
    return
      if ( $start > $last && $stop > $last )
      ;    #asking for sequence beyond the end;
    return if $start > $last && $stop > $last;    #outside of chromosome
    return if $start < 1 && $stop < 1;

    if ( defined $start && defined $stop ) {
        $start = 1 if $start < 1;
        $stop = $start unless $stop;
        return undef unless ( $start =~ /^\d+$/ and $stop =~ /^\d+$/ );
        ( $start, $stop ) = ( $stop, $start ) if $stop < $start;
    }
    return get_genome_seq(    #$self->get_seq( # mdb changed 7/31/13 issue 77
        gid    => $self->id,
        chr    => $chr,
        start  => $start,
        stop   => $stop,
        strand => $strand
    );
}

# mdb added 7/29/13, issue #77
sub file_path {
    my $self = shift;
    return CoGe::Core::Storage::get_genome_file( $self->id );
}

sub storage_path {
    return shift->file_path;
}

# mdb added 8/6/13, issue #157
sub create_index {
    my $self = shift;
    return CoGe::Core::Storage::index_genome_file( gid => $self->id );
}

sub is_indexed {
    my $self = shift;
    return ( -e $self->file_path . '.fai' );
}

##################################################

sub is_loading {
    my $self = shift;
    return ($self->status && $self->status == LOADING);
}

sub is_loaded {
    my $self = shift;
    return (!$self->status || $self->status == LOADED); # legacy rows will have NULL status
}

sub is_error {
    my $self = shift;
    return ($self->status && $self->status == ERROR);
}

################################################ subroutine header begin ##

=head2 is_editable

 Usage     : is this genome editable by the specified user?
 Purpose   :
 Returns   : 0 or 1
 Argument  :
 Throws    : None
 Comments  :

=cut

################################################## subroutine header end ##

sub is_editable {
    my $self = shift;
    my $user = shift;

    return ( $user->is_admin || $user->is_owner_editor( dsg => $self->id ) ) ? 1 : 0;
}

sub is_deletable {
    my $self = shift;
    my $user = shift;

    return ( $user->is_admin || $user->is_owner( dsg => $self->id ) ) ? 1 : 0;
}

##################################################

sub sequence_type {
    shift->genomic_sequence_type(@_);
}

################################################

sub type {
    shift->genomic_sequence_type(@_);
}

################################################

sub chromosomes {
    my $self = shift;
	return CoGe::Core::Chromosomes->new($self->id)->names;
}

################################################

sub chromosomes_all {
    my $self = shift;
    return CoGe::Core::Chromosomes->new($self->id)->all;
}

################################################

sub chromosome_lengths {
    my $self = shift;
	return CoGe::Core::Chromosomes->new($self->id)->lengths;
}

sub chromosome_lengths_by_name {
    my $self = shift;
    return CoGe::Core::Chromosomes->new($self->id)->lengths_by_name;
}

################################################

sub percent_gc {
    my $self     = shift;
    my %opts     = @_;
    my $count    = $opts{count};
    my $sent_chr = $opts{chr};

    my @chr;
    push @chr, $sent_chr if $sent_chr;

    my $gc     = 0;
    my $at     = 0;
    my $n      = 0;
    my $x      = 0;
    my $length = 0;

    unless ($sent_chr) {
        foreach my $chr ( $self->chromosomes ) {
            push @chr, $chr;
        }
    }
    foreach my $chr (@chr) {
        my $seq = $self->get_genomic_sequence( chr => $chr );
        $length += length $seq;
        $gc     += $seq =~ tr/GCgc/GCgc/;
        $at     += $seq =~ tr/ATat/ATat/;
        $n      += $seq =~ tr/nN/nN/;
        $x      += $seq =~ tr/xX/xX/;
    }
    return unless $length;
    return ( $gc, $at, $n, $x ) if $count;
    return sprintf( "%.4f", $gc / $length ), sprintf( "%.4f", $at / $length ),
      sprintf( "%.4f", $n / $length ), sprintf( "%.4f", $x / $length );
}

################################################ subroutine header begin ##

=head2 fasta

 Usage     :
 Purpose   :
 Returns   :
 Argument  :
            col      =>   number of sequence characters per line (default 100)
            chr_name =>   fasta header contains only the chromosome name (default 0)
            start    =>  start position (default 1)
            stop     =>  stop position  (default $self->sequence_legnth($chr)
            chr      =>  chromosome for which to get sequence (default:  whatever $self->chromosomes gets first)
            rc       =>  generate the reverse complement (default: 0)
            prot     =>  translate to protein, will do 6 frame automatically if it is not in a proper reading frame (default: 0)

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
    my $chr_name = $opts{chr_name};   #makes header contain only chromosome name
    my $chr      = $opts{chr};
    ($chr) = $self->chromosomes unless defined $chr;
    my $strand = $opts{strand} || 1;
    my $start  = $opts{start}  || 1;
    $start = 1 if $start < 1;
#    my $stop = $opts{stop} || $self->sequence_length($chr);
    my $stop = $opts{stop} || $self->get_chromosome_length($chr);
    my $prot = $opts{prot};
    my $rc   = $opts{rc};
    $strand = -1 if $rc;

    my $seq = $self->get_genomic_sequence(
        start => $start,
        stop  => $stop,
        chr   => $chr
    );
    $stop = $start + length($seq) - 1 if $stop > $start + length($seq) - 1;
    my $head;

    if ($chr_name) {
        $head .= ">" . $chr;
    }
    else {
        $head .= ">" . $self->organism->name;
        $head .= $self->name if $self->name;
        $head .= ", " . $self->description if $self->description;
        $head .= " (v"
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
    }

    $Text::Wrap::columns = $col;
    my $fasta;

    $seq = $self->reverse_complement($seq) if $rc;
    if ($prot) {
        my $trans_type = $self->trans_type;
        my $feat       = new CoGeX::Result::Feature;
        my ( $seqs, $type ) =
          $feat->frame6_trans( seq => $seq, trans_type => $trans_type );
        foreach my $frame (
            sort { length($a) <=> length($b) || $a cmp $b }
            keys %$seqs
          )
        {
            $seq = $seqs->{$frame};
            $seq = $self->reverse_complement($seq) if $rc;
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

 Usage     : $dsg->gff(print=>1)
 Purpose   : generating a gff file for a genome from all the datasets it contains
 Returns   : a string
 Argument  : name_re     =>    regular expression for only displaying features containing a name that matches
             print       =>    print the gff file as the lines are retrieved
             annos       =>    print annotations as well (takes longer)
             cds         =>    Only print CDS gene features (skip all ncRNA and other features).  Will print genes, mRNA, and CDS entries
             id_type     =>    Specify if the GFF entry IDs are going to be unique numbers or unique names.
             unique_parent_annotations => Flag to NOT print redundant annotations in children entries.  E.g. if parent has an annotation, a child will not have that annotation
             name_unique =>   Flag for specifying that the name tag of an entry will be unique
 Throws    :
 Comments  :

See Also   : dataset->gff

=cut

################################################## subroutine header end ##

sub gff {
    my $self    = shift;
    my %opts    = @_;
    my $name_re = $opts{name_re};
    my $debug   = $opts{debug};
    my $print   = $opts{print};
    my $annos   = $opts{annos};
    my $cds     = $opts{cds}; #only print CDS gene features
    my $chr		= $opts{chr}; #optional, set to only include features on a particular chromosome
    my $add_chr = $opts{add_chr};
    my $unique_parent_annotations = $opts{unique_parent_annotations}; #parent annotations are NOT propogated to children
    my $id_type = $opts{id_type}; #type of ID (name, num):  unique number; unique name
    $id_type = "name" unless defined $id_type;

    my $name_unique = $opts{name_unique}; #flag for making Name tag of output unique by appending type and occurrence to feature name
    my $base_url = $opts{base_url}; #server for coge to link back to GenomeInfo

    my $output; #store the goodies
    $output .= "##gff-version\t3\n";
    $output .= "##Generated by CoGe\n";
    $output .= "##Organism name: " . $self->organism->name . "\n";
    $output .= "##Organism desc: " . $self->organism->description . "\n"
      if $self->organism->description;
    $output .= "##Version: " . $self->version . "\n";
    $output .= "##CoGe Genome ID (gid): " . $self->id . "\n";
    $output .= "##CoGe GenomeInfo Link: $base_url/GenomeInfo.pl?gid=" . $self->id . "\n" if $base_url;
    $output .= "##\n";
    $output .= "##\n";
    print $output if $print;
    my $id = 0;
    foreach my $ds ( $self->datasets ) {
        my $tmp;
        ( $tmp, $id ) = $ds->gff(
            name_re                   => $name_re,
            debug                     => $debug,
            print                     => $print,
            annos                     => $annos,
            no_gff_head               => 1,
            id                        => $id,
            cds                       => $cds,
            name_unique               => $name_unique,
            id_type                   => $id_type,
            unique_parent_annotations => $unique_parent_annotations,
            chr						  => $chr,
	        add_chr                   => $add_chr
        );
        $output .= $tmp if $tmp;
    }
    return $output;
}

################################################

sub trans_type {
    my $self = shift;
    foreach my $ds ( $self->datasets() ) {
        my $type = $ds->trans_type();
        return $type if $type;
    }
    return 1;    #universal genetic code type;
}

# mdb removed 7/30/13 issue 77, moved into Accessory::Storage
#sub get_path {
#    my $self      = shift;
#    my $gid = $self->id;
#    my $level0    = floor( $gid / 1000000000 ) % 1000;
#    my $level1    = floor( $gid / 1000000 ) % 1000;
#    my $level2    = floor( $gid / 1000 ) % 1000;
#    my $path      = catdir( $level0, $level1, $level2, $gid )
#      ; #adding $gid for final directory.  blast's formatdb will be run on the faa file and this will help keep that stuff organized
#    return $path;
#}

sub chr_info {
    my $self    = shift;
    my %opts    = @_;
    my $summary = $opts{summary};

    my $html;
    my $total_length;
    my $count = 0;
    my $chr_list;
#    my @gs = sort {
#             $a->chromosome =~ /(\d+)/ <=> $b->chromosome =~ /(\d+)/
#          || $a->chromosome cmp $b->chromosome
#    } $self->genomic_sequences;
#
#    foreach my $gs (@gs) {
#        my $chr    = $gs->chromosome;
#        my $length = $gs->sequence_length;
#        $total_length += $length;
#        $length = commify($length);
#        $chr_list .= qq{$chr:  $length bp<br>};
#        $count++;
#    }
    my $c = CoGe::Core::Chromosomes->new($self->id);
    while ($c->next) {
        $total_length += $c->length;
        my $chr    = $c->name;
        my $length = commify($c->length);
        $chr_list .= qq{$chr:  $length bp<br>};
        $count++;
    }
    $html .=
        qq{Chromosome count: $count<br>}
      . qq{Total length: }
      . commify($total_length) . " bp";
    $html .= "<br>" . qq{-----------<br>Chr:   (length)<br>} . $chr_list
      unless $summary;
    return $html;
}

############################################### subroutine header begin ##

=head2 length

 Usage     : $self->length
 Purpose   : get total length of sequence in dataset group
 Returns   : number
 Argument  :
 Throws    :
 Comments  :

See Also   :

=cut

################################################## subroutine header end ##

sub length {
    my $self = shift;
    return CoGe::Core::Chromosomes->new($self->id)->total_length;
}

############################################### subroutine header begin ##

=head2 chromosome_count

 Usage     :
 Purpose   :
 Returns   : number of chromosomes in the genome
 Argument  :
 Throws    :
 Comments  :

See Also   :

=cut

################################################## subroutine header end ##

sub chromosome_count {
	my $self = shift;
    return CoGe::Core::Chromosomes->new($self->id)->count;
}

################################################ subroutine header begin ##

=head2 features

 Usage     : $self->features
 Purpose   : run through associated datasets and get their features
 Returns   : array of feature objects
 Argument  :
 Throws    :
 Comments  :

See Also   :

=cut

################################################## subroutine header end ##

sub features {
    my $self = shift;
    my @feats;
    foreach my $ds ( $self->datasets ) {
        my @tmp = $ds->features(@_);
        push @feats, @tmp;
    }
    return wantarray ? @feats : \@feats;
}

sub feature_count {
	my $self = shift;
	my %opts = @_;
    my $chr = $opts{chr}; # optional
    my $feature_type_id = $opts{feature_type_id}; # optional

    my $total = 0;
    foreach my $ds ( $self->datasets ) {
    	my $n = $ds->features( chromosome => $chr, feature_type_id => $feature_type_id )->count;
    	$total += $n if $n;
    }

    return $total;
}

sub translation_type {
    my $self = shift;
    foreach my $ds ( $self->datasets ) {
        my $trans_type = $ds->translation_type;
        return $trans_type if $trans_type;
    }
}

sub distinct_feature_type_ids {
    my $self = shift;
    my %opts = @_;
    my %ids =
      map { $_ => 1 } map { $_->distinct_feature_type_ids } $self->datasets;
    return wantarray ? keys %ids : [ keys %ids ];
}

sub distinct_feature_type_names {
    my $self = shift;
    my %opts = @_;
    my %names =
      map { $_ => 1 } map { $_->distinct_feature_type_names } $self->datasets;
    return wantarray ? keys %names : [ keys %names ];
}

sub has_gene_features {
    my $self = shift;
    
    foreach ($self->distinct_feature_type_ids) { # FIXME use grep instead
        return 1 if ($_ == 1 || $_ == 2 || $_ == 3); # FIXME hardcoded feature types
    }
    
    return 0;
}

sub source {
    my $self = shift;
    my %sources;
    foreach my $ds ( $self->datasets ) {
        $sources{ $ds->data_source->id } = $ds->data_source;
    }
    return wantarray ? values %sources : [ values %sources ];
}

################################################ subroutine header begin ##

=head2 info

 Usage     : $self->info
 Purpose   : returns a string of information about the genome.

 Returns   : returns a string
 Argument  : none
 Throws    :
 Comments  : To be used to quickly generate a string about the genome

See Also   :

=cut

################################################## subroutine header end ##

sub info {
    my $self = shift;
    my %opts = @_;
    
    my $info;
    $info .= "&#x1f512; "              if ($self->restricted && !$opts{hideRestrictedSymbol}); #TODO move this into view code
    $info .= $self->organism->name     if $self->organism;
    $info .= " (" . $self->name . ")"  if $self->name;
    $info .= ": " . $self->description if $self->description;
    $info .= " (v" . $self->version . ", id" . $self->id . ")";
    $info .= ': ' . $self->genomic_sequence_type->name
      if $self->genomic_sequence_type;
    return $info;
}

############################################### subroutine header begin ##

=head2 info_html

 Usage     :
 Purpose   : provides quick information about the genome wrapped with a link to LIstView
 Returns   : a string
 Argument  :
 Throws    :
 Comments  : name, description, restricted, type
           :

See Also   :

=cut

################################################## subroutine header end ##

sub info_html {
    my $self = shift;
    my $info = $self->info;
    return
        qq{<span class=link onclick='window.open("GenomeInfo.pl?gid=}
      . $self->id
      . qq{")'>}
      . $info
      . "</span>";
}

sub info_file {
    my $self = shift;
    
    my $restricted = ($self->restricted) ? "yes" : "no";
    my $genome_name = $self->genome->info(hideRestrictedSymbol=>1);

    my @lines = (
        qq{"Name","} . $self->name . '"',
        qq{"Description","} . $self->description . '"',
        qq{"Source","} . $self->source->info . '"',
        qq{"Version","} . $self->version . '"',
        qq{"Organism","} . $self->organism->name . '"',
        qq{"Sequence Type", "} . $self->genomic_sequence_type->name . '"',
        qq{"Notebooks","} . $self->notebooks_desc . '"',
        qq{"Restricted","$restricted"},
    );
    push @lines, qq{"Link","} . $self->link . '"' if ($self->link);

    return join("\n", @lines);
}

############################################### subroutine header begin ##

=head2 date

 Usage     :
 Purpose   : returns load date from first dataset entry
 Returns   : a string
 Argument  :
 Throws    :
 Comments  :

See Also   :

=cut

################################################## subroutine header end ##

sub get_date {
    my $self = shift;
    return $self->date if ($self->date && $self->date ne '0000-00-00 00:00:00');
    my ($ds) = $self->datasets;
    return $ds ? $ds->date : '';
}

##################################################

sub get_codon_usage {
    my $self = shift;
    my $chr = shift;

    my %seqs; # prefetch the sequences with one call to genomic_sequence (slow for many seqs)
    my @chrs = (defined $chr and $chr) ? ($chr) : $self->chromosomes;
    for my $chr (@chrs) {
        $seqs{$chr} = $self->get_genomic_sequence(chr => $chr);
    }

    my %codons;
    my $codon_total = 0;
    my ( $code, $code_type );
    my $search = { "feature_type.name" => "CDS" };
    $search->{"me.chromosome"} = $chr if defined $chr;

    foreach my $ds ($self->datasets) {
        foreach my $feat (
            $ds->features(
                $search,
                {
                    join => [
                        "feature_type", 'locations',
                        { 'dataset' => { 'dataset_connectors' => 'genome' } }
                    ],
                    prefetch => [
                        'locations',
                        { 'dataset' => { 'dataset_connectors' => 'genome' } }
                    ]
                }
            )
          )
        {
            my $seq = substr(
                $seqs{ $feat->chromosome },
                $feat->start - 1,
                $feat->stop - $feat->start + 1
            );
            $feat->genomic_sequence( seq => $seq );
            ( $code, $code_type ) = $feat->genetic_code() unless $code;
            my ($codon) = $feat->codon_frequency( counts => 1 );
            grep { $codon_total += $_ } values %$codon;
            grep { $codons{$_}  += $codon->{$_} } keys %$codon;
        }
    }
    %codons = map { $_, $codons{$_} / $codon_total } keys %codons;

    return "Codon Usage: $code_type" .
        CoGe::Accessory::genetic_code->html_code_table(
            data => \%codons,
            code => $code
        );
}
1;
