package CoGe::Accessory::IRODS;

=head1 NAME

CoGe::Accessory::IRODS

=head1 SYNOPSIS

Interface to iCommands

=head1 DESCRIPTION

=head1 AUTHOR

Matt Bomhoff

=head1 COPYRIGHT

The full text of the license can be found in the
LICENSE file included with this module.

=head1 SEE ALSO

=cut

use strict;
use warnings;

use CoGe::Accessory::Web;
#use Data::Dumper;
use IPC::System::Simple qw(capture system $EXITVAL EXIT_ANY);

BEGIN {
    use vars qw ($VERSION @ISA @EXPORT $IRODS_METADATA_PREFIX $IRODS_ENV);
    require Exporter;

    $VERSION = 0.1;
    @ISA     = qw (Exporter);
    @EXPORT = qw( irods_ils irods_imeta irods_iget irods_chksum irods_iput $IRODS_METADATA_PREFIX );

    $IRODS_METADATA_PREFIX = 'ipc-coge-';
}

sub irods_ils {
    my $path = shift;
    $path = '' unless $path;

    #print STDERR "irods_ils: path=$path\n";
    my $env_file = _irods_get_env_file();
    return { error => "Error: iRODS env file missing" } unless $env_file;

    my $cmd = "export irodsEnvFile='$env_file'; ils -l $path 2>&1";

#	print STDERR "cmd: $cmd\n";
#	my @ils = `$cmd`; # old way of executing command, replaced by better error checking below
    my @ils = capture( EXIT_ANY, $cmd );
    if ($EXITVAL) {
        return { error => "Error: ils rc=$EXITVAL" };
    }

    $path = shift @ils;

    #	if ($path =~ /^ERROR/) { # iRODS error message
    #		my $result = { type => 'error', name => $path };
    #		return wantarray ? ($result) : [$result];
    #	}
    chomp($path);
    chop($path);

    my @result;
    foreach my $line (@ils) {
        my ( $type, $backup, $size, $timestamp, $name );

        chomp $line;
        if ( $line =~ /^\s*C\-/ ) {    # directory
            $type = 'directory';
            ($name) = $line =~ /([^\/\s]+)\s*$/;
            if ($name) { $name .= '/'; }
            else       { $name = 'error' }
            ( $size, $timestamp ) = ( '', '' );
        }
        else {                         # file
            $type = 'file';
            ( undef, undef, $backup, undef, $size, $timestamp, undef, $name ) =
              split( /\s+/, $line );
        }
		next if $backup;

        push @result,
          {
            type      => $type,
            size      => $size,
            timestamp => $timestamp,
            name      => $name,
            path      => $path . '/' . $name
          };
    }
    @result =
      sort { $a->{type} cmp $b->{type} } @result;    # directories before files

    return { items => \@result };
}

sub irods_chksum {
    my $path = shift;
    return 0 unless ($path);

    my $env_file = _irods_get_env_file();
    return unless $env_file;

    my $cmd = "export irodsEnvFile='$env_file'; ichksum $path";
    #print STDERR "cmd: $cmd\n";
    my @output = `$cmd`;
    my ($chksum) = $output[0] =~ /\s*\S+\s+(\S+)/;

    #print STDERR "chksum: $chksum\n";
    return $chksum;
}

sub irods_iget {
    my ( $src, $dest, $opts ) = @_;
    my $no_execute = ( $opts and $opts->{no_execute} ); # mdb added 4/10/14 for REST API, get command name but don't execute

    #print STDERR "irods_iget $src $dest\n";
    my $env_file = _irods_get_env_file();
    return unless $env_file;

    my $cmd = "export irodsEnvFile='$env_file'; iget -fT " . ($src || '') . ' ' . ($dest || '');
    return $cmd if $no_execute;
    #print STDERR "cmd: $cmd\n";
    my @result = `$cmd`;
    #print STDERR "@result";

    return;
}

sub irods_iput {
    my ( $src, $dest ) = @_;

    #print STDERR "irods_iput $src $dest\n";
    my $env_file = _irods_get_env_file();
    return unless $env_file;

    my $cmd = "export irodsEnvFile='$env_file'; iput -fT $src $dest";
    #print STDERR "cmd: $cmd\n";
    my @result = `$cmd`;
    #print STDERR "@result";

    return;
}

sub irods_imeta {
    my ( $dest, $pHash ) = @_;
    #print STDERR "irods_imeta $dest\n";
    return unless ( defined $pHash and keys %$pHash );

    my $env_file = _irods_get_env_file();
    return unless $env_file;

	while( my($k, $v) = each %$pHash ) {
		if (not $k or not $v) {
			print STDERR "CoGe::Accessory::IRODS::imeta: improper key/value pair\n";
			next;
		}

	    my $cmd = "export irodsEnvFile='$env_file'; imeta add -d '" . $dest . "' '" . $k . "' '" . $v . "'";
	    print STDERR "cmd: $cmd\n";
	    my @result = `$cmd`;
	    #print STDERR "@result";
	}

    return;
}

sub irods_set_env {
    $IRODS_ENV = shift;
}

sub _irods_get_env_file {
    my $env_file = $IRODS_ENV;
    $env_file //= CoGe::Accessory::Web::get_defaults()->{IRODSENV};

    if ( not defined $env_file or not -e $env_file ) {
        print STDERR "CoGe::Accessory::IRODS: fatal error: iRODS env file missing!\n";
        return;
    }
    return $env_file;
}

1;
