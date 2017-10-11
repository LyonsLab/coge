package CoGe::Core::Chromosomes;

=head1 NAME

CoGe::Core::Chromosomes

=head1 SYNOPSIS

provides class for accessing chromosome data from FASTA index files (.fai)

=head1 DESCRIPTION

=head1 AUTHOR

Sean Davey

=head1 COPYRIGHT

The full text of the license can be found in the
LICENSE file included with this module.

=head1 SEE ALSO

=cut

use strict;
use warnings;

use CoGe::Core::Storage qw(get_genome_file);
use Data::Dumper;

################################################ subroutine header begin ##

=head2 new

 Usage     :
 Purpose   :
 Returns   : newly instantiated object for this class
 Argument  : genome id (required)
 Throws    :
 Comments  :

See Also   :

=cut

################################################## subroutine header end ##

sub new {
	my ($class, $gid) = @_;
	my $self = {};
	my $genome_file = get_genome_file($gid);
	if ($genome_file) {
	    my $fh;
	    if (!open($fh, $genome_file . '.fai')) { # sd added 12/8/2015 COGE-687
	        #warn 'error opening index file in Chromosomes::new()';
	        sleep 1;
	        open($fh, $genome_file . '.fai');
	    }
        $self->{fh} = $fh if open($fh, $genome_file . '.fai');
	}
	$self->{lines} = 0;
	return bless $self, $class;
}

################################################ subroutine header begin ##

=head2 DESTROY

 Usage     :
 Purpose   : closes file if still open when object goes out of scope
 Returns   :
 Argument  :
 Throws    :
 Comments  :

See Also   :

=cut

################################################## subroutine header end ##

sub DESTROY {
	my $self = shift;
	if ($self->{fh}) {
		close($self->{fh});
	}	
}

################################################ subroutine header begin ##

=head2 count

 Usage     :
 Purpose   :
 Returns   : number of chromosomes for the genome
 Argument  :
 Throws    :
 Comments  :

See Also   :

=cut

################################################## subroutine header end ##

sub count {
	my $self = shift;
	while ($self->next) {
	}
	return $self->{lines};
}

################################################ subroutine header begin ##

=head2 find

 Usage     : 
 Purpose   : iterates through list and stops when matching chromosome is found
 Returns   : 1 if found, 0 otherwise
 Argument  : name of chromosome to find (required)
 Throws    :
 Comments  :

See Also   :

=cut

################################################## subroutine header end ##

sub find {
	my $self = shift;
	my $name = shift;
	while ($self->next) {
		if ($name eq $self->name) {
			return 1;
		}
	}
	return 0;
}

################################################ subroutine header begin ##

=head2 length

 Usage     :
 Purpose   : 
 Returns   : length of current chromosome
 Argument  :
 Throws    :
 Comments  :

See Also   :

=cut

################################################## subroutine header end ##

sub length {
	my $self = shift;
	return $self->{tokens}[1];
}

################################################ subroutine header begin ##

=head2 lengths

 Usage     :
 Purpose   : 
 Returns   : array of lengths of all chromosomes for the genome
 Argument  :
 Throws    :
 Comments  :

See Also   :

=cut

################################################## subroutine header end ##

sub lengths {
	my $self = shift;
	my @a;
	while ($self->next) {
		push @a, $self->length;
	}
    return wantarray ? @a : \@a;
}

sub lengths_by_name {
    my $self = shift;
    my %lengths;
    while ($self->next) {
        $lengths{$self->name} = $self->length;
    }
    return \%lengths;
}

################################################ subroutine header begin ##

=head2 name

 Usage     :
 Purpose   : 
 Returns   : name of current chromosome
 Argument  :
 Throws    :
 Comments  :

See Also   :

=cut

################################################## subroutine header end ##

sub name {
	my $self = shift;
	my $name = $self->{tokens}[0];
	my $index = index($name, '|');
	$name = substr($name, $index + 1) if $index != -1;
	return $name;
}

################################################ subroutine header begin ##

=head2 names

 Usage     :
 Purpose   : 
 Returns   : array of names of all chromosomes for the genome
 Argument  :
 Throws    :
 Comments  :

See Also   :

=cut

################################################## subroutine header end ##

sub names {
	my $self = shift;
	my @a;
	while ($self->next) {
		push @a, $self->name;
	}
    return wantarray ? @a : \@a;
}

################################################ subroutine header begin ##

=head2 all

 Usage     :
 Purpose   : 
 Returns   : array of hashes of all names/lengths of all chromosomes for the genome
 Argument  :
 Throws    :
 Comments  :

See Also   :

=cut

################################################## subroutine header end ##

sub all {
    my $self = shift;
    my @a;
    while ($self->next) {
        push @a, { name => $self->name, length => int($self->length) };
    }
    return wantarray ? @a : \@a;
}

################################################ subroutine header begin ##

=head2 hash

 Usage     :
 Purpose   : 
 Returns   : hash by name lengths and offsets of all chromosomes for the genome
 Argument  :
 Throws    :
 Comments  :

See Also   :

=cut

################################################## subroutine header end ##

sub hash {
    my $self = shift;
    my $hash = {};
    while ($self->next) {
        $hash->{$self->name} = { length => int($self->length), offset => int($self->offset) };
    }
    return $hash;
}

################################################ subroutine header begin ##

=head2 next

 Usage     :
 Purpose   : set the current chromosome to be the next one in the list
 Returns   : 1 if new chromosome is current, 0 if no more chromosomes available
 Argument  :
 Throws    :
 Comments  :

See Also   :

=cut

################################################## subroutine header end ##

sub next {
	my $self = shift;
	if (!$self->{fh}) {
		#warn caller;
		return 0;
	}
	my $line = readline($self->{fh});
	if ($line) {
		my @tokens = split('\t', $line);
		@{$self->{tokens}} = @tokens;
		$self->{lines}++;
		return 1;
	}
	close($self->{fh});
	$self->{fh} = 0;
	return 0;
}

################################################ subroutine header begin ##

=head2 offset

 Usage     :
 Purpose   : 
 Returns   : offset of current chromosome
 Argument  :
 Throws    :
 Comments  :

See Also   :

=cut

################################################## subroutine header end ##

sub offset {
	my $self = shift;
	return $self->{tokens}[2];
}

################################################ subroutine header begin ##

=head2 total_length

 Usage     :
 Purpose   : 
 Returns   : length of all chromosomes for the genome
 Argument  :
 Throws    :
 Comments  :

See Also   :

=cut

################################################## subroutine header end ##

sub total_length {
	my $self = shift;
	my $length = 0;
	while ($self->next) {
		$length += $self->length;
	}
    return $length;
}

1;
