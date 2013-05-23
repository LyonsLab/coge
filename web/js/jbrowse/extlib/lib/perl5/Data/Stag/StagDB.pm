package Data::Stag::StagDB;

=head1 NAME

  Data::Stag::StagDB - persistent storage and retrieval of stag nodes

=head1 SYNOPSIS

  # parsing a file into a file based index
  my $sdb = Data::Stag::StagDB->new;
  $sdb->unique_key("ss_details/social_security_no");
  $sdb->record_type("person");
  $sdb->indexfile("./person_by_ss-idx");
  Data::Stag->parse(-file=>$fn, -handler=>$sdb);
  my $obj = $sdb->index_hash;
  my $person = $obj->{'999-9999-9999'};
  print $person->xml;

  # indexing an existing stag tree into a file based index
  my $personset = Data::Stag->parse($fn);
  my $sdb = Data::Stag::StagDB->new;
  $sdb->unique_key("ss_details/social_security_no");
  $sdb->record_type("person");
  $sdb->indexfile("./person_by_ss-idx");
  $personset->sax($sdb);
  my $obj = $sdb->index_hash;
  my $person = $obj->{'999-9999-9999'};
  print $person->xml;


=cut

=head1 DESCRIPTION

This module is an extension of L<Data::Stag::HashDB> - you can use it
in the same way.

It creates a simple file based database of stag nodes

This is useful if you want your data to persist; or if you want to use L<Data::Stag::HashDB> but your data will not fit into memory

=head1 PUBLIC METHODS -

=cut

use strict;
use base qw(Data::Stag::HashDB);
use Data::Stag qw(:all);
use MLDBM qw(DB_File Storable);
use Fcntl;

use vars qw($VERSION);
$VERSION="0.11";

sub init {
    my $self = shift;
    $self->SUPER::init(@_);
    return $self;
}


=head2 indexfile

  Usage   -
  Returns -
  Args    -

=cut

sub indexfile {
    my $self = shift;
    $self->{_indexfile} = shift if @_;
    return $self->{_indexfile};
}

sub index_hash {
    my $self = shift;
    $self->{_index_hash} = shift if @_;
    if (!$self->{_index_hash}) {
	my %obj = ();
	my $f = $self->indexfile || 'tmp-stagdb';
	my $dbm = tie %obj, 'MLDBM', $f, O_CREAT|O_RDWR, 0640 or die $!;
	$self->dbm($dbm);
	$self->{_index_hash} = \%obj;
    }
    return $self->{_index_hash};
}

sub dbm {
    my $self = shift;
    $self->{_dbm} = shift if @_;
    return $self->{_dbm};
}


1;

