package Data::Stag::HashDB;

=head1 NAME

  Data::Stag::HashDB

=head1 SYNOPSIS

  # parsing a file into a hash
  my $hdb = Data::Stag::HashDB->new;
  $hdb->unique_key("ss_details/social_security_no");
  $hdb->record_type("person");
  my $obj = {};
  $hdb->index_hash($obj);
  Data::Stag->parse(-file=>$fn, -handler=>$hdb);
  my $person = $obj->{'999-9999-9999'};
  print $person->xml;

  # indexing an existing stag tree into a hash
  my $personset = Data::Stag->parse($fn);
  my $hdb = Data::Stag::HashDB->new;
  $hdb->unique_key("ss_details/social_security_no");
  $hdb->record_type("person");
  my $obj = {};
  $hdb->index_hash($obj);
  $personset->sax($hdb);
  my $person = $obj->{'999-9999-9999'};
  print $person->xml;


=cut

=head1 DESCRIPTION

Used for building indexes over Stag files or objects

You need to provide a B<record_type> - this is the type of element
that will be indexed

You need to provide a N<unique_key> - this is a single value used to
index the B<record_type>s

For example, if we have data in the stag structure below, and if ss_no
is unique (we assume it is) then we can index all the people in the
database using the code above

  publicinfo:
    persondata:
      person:
        ss_details:
          social_security_no:
        name:
        address: 

There is a subclass of this method callsed Data::Stag::StagDB, which
makes the hash persistent

=head1 PUBLIC METHODS -

=cut

use strict;
use base qw(Data::Stag::BaseHandler);
use Data::Stag qw(:all);

use vars qw($VERSION);
$VERSION="0.11";

sub init {
    my $self = shift;
    $self->SUPER::init(@_);
    $self->nextid(0);
    return $self;
}


=head2 record_type

  Usage   -
  Returns -
  Args    -

=cut

sub record_type {
    my $self = shift;
    $self->{_record_type} = shift if @_;
    return $self->{_record_type} || '';
}

=head2 unique_key

  Usage   -
  Returns -
  Args    -

=cut

sub unique_key {
    my $self = shift;
    $self->{_unique_key} = shift if @_;
    return $self->{_unique_key};
}


=head2 index_hash

  Usage   -
  Returns -
  Args    -

=cut

sub index_hash {
    my $self = shift;
    $self->{_index_hash} = shift if @_;
    if (!$self->{_index_hash}) {
	$self->{_index_hash} = {};
    }
    return $self->{_index_hash};
}


sub nextid {
    my $self = shift;
    if (@_) {
	$self->{_nextid} = shift;
    }
    else {
	$self->{_nextid} = 0 unless $self->{_nextid};
	$self->{_nextid}++;
    }
    return $self->{_nextid};
}

sub end_event {
    my $self = shift;
    my $ev = shift;
    if ($ev eq $self->record_type) {
	my $topnode = $self->popnode;
	$self->add_record(stag_stagify($topnode));
#	my $name_elt = $self->unique_key;
#	my $name;
#	if ($name_elt) {
#	    $name = stag_get($topnode, $name_elt);
#	}
#	if (!$name) {
#	    $name = $ev."_".$self->nextid;
#	}
#	$self->index_hash->{$name} = stag_stagify($topnode);
	return [];
    }
    else {
	return $self->SUPER::end_event($ev, @_);
    }
}

sub add_record {
    my $self = shift;
    my $record = shift;
    
    my $idx = $self->index_hash;
    my $ukey = $self->unique_key;
    my $keyval;
    if ($ukey) {
	$keyval = stag_get($record, $ukey);
    }
    if (!$keyval) {
	$keyval = $record->name."_".$self->nextid;
    }
    $idx->{$keyval} = [] unless $idx->{$keyval};
    my $vals = $idx->{$keyval};
    push(@$vals, $record);
    $idx->{$keyval} = $vals;
    return;
}

sub get_record {
    my $self = shift;
    my $keyval = shift;
    my $records = $self->index_hash->{$keyval} || [];
    if (wantarray) {
	return @$records;
    }
    else {
	return $records->[0];
    }
}

sub reset {
    my $self = shift;
    %{$self->index_hash} = ();
    return;
}

1;

