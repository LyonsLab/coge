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
#use Data::Dumper;
use IPC::System::Simple qw(capture system $EXITVAL EXIT_ANY);

BEGIN {
    use vars qw ($VERSION @ISA @EXPORT);
    require Exporter;

    $VERSION = 0.1;
    @ISA     = qw (Exporter);
    @EXPORT = qw( irods_ils );
}

sub irods_ils {
    my $path = shift;
    $path = '' unless $path;

    #	print STDERR "irods_ils: path=$path\n";
    my $env_file = CoGe::Accessory::Web::get_defaults()->{IRODSENV};
    if ( not defined $env_file or not -e $env_file ) {
        print STDERR "fatal error: iRODS env file missing!\n";
        return { error => "Error: iRODS env file missing" };
    }

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

    my $env_file = CoGe::Accessory::Web::get_defaults()->{IRODSENV};
    if ( not defined $env_file or not -e $env_file ) {
        print STDERR "fatal error: iRODS env file missing!\n";
        return;
    }

    my $cmd = "export irodsEnvFile='$env_file'; ichksum $path";
    #print STDERR "cmd: $cmd\n";
    my @output = `$cmd`;
    my ($chksum) = $output[0] =~ /\s*\S+\s+(\S+)/;

    #print STDERR "chksum: $chksum\n";
    return $chksum;
}

sub irods_iget {
    my ( $src, $dest ) = @_;

    #print STDERR "irods_iget $src $dest\n";
    my $env_file = CoGe::Accessory::Web::get_defaults()->{IRODSENV};
    if ( not defined $env_file or not -e $env_file ) {
        print STDERR "fatal error: iRODS env file missing!\n";
        return;
    }

    my $cmd = "export irodsEnvFile='$env_file'; iget -fT $src $dest";

    #print STDERR "cmd: $cmd\n";
    my @result = `$cmd`;
    #print STDERR "@result";

    return;
}

sub irods_iput {
    my ( $src, $dest ) = @_;

    #print STDERR "irods_iput $src $dest\n";
    my $env_file = CoGe::Accessory::Web::get_defaults()->{IRODSENV};
    if ( not defined $env_file or not -e $env_file ) {
        print STDERR "fatal error: iRODS env file missing!\n";
        return;
    }

    my $cmd = "export irodsEnvFile='$env_file'; iput -fT $src $dest";

    print STDERR "cmd: $cmd\n";
    my @result = `$cmd`;
    print STDERR "@result";

    return;
}

1;