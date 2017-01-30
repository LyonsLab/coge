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
use URI::Escape qw(uri_unescape);
use URI::Escape::JavaScript qw(escape unescape);
use File::Basename qw( basename );
use File::Spec::Functions qw( catdir catfile );
use Data::Dumper;
use IPC::System::Simple qw(capture system $EXITVAL EXIT_ANY);

BEGIN {
    use vars qw ($VERSION @ISA @EXPORT @EXPORT_OK $IRODS_METADATA_PREFIX $IRODS_ENV);
    require Exporter;

    $VERSION = 0.1;
    @ISA     = qw (Exporter);
    @EXPORT = qw( irods_ils irods_imeta_add irods_imeta_ls irods_iget irods_chksum irods_imkdir irods_iput irods_irm $IRODS_METADATA_PREFIX );
    @EXPORT_OK = qw( irods_get_base_path irods_set_env );

    $IRODS_METADATA_PREFIX = 'coge-'; #'ipc-coge-'; # mdb changed 9/27/16 -- getting error from IRODS that 'ipc' is protected AVU
}

sub irods_ils {
    my $path = shift;
    $path = '' unless $path;
    my %opts = @_;
    my $escape_output = $opts{escape_output} || 0;
    #print STDERR "irods_ils: path=$path\n";
    my $env_file = _irods_get_env_file();
    return { error => "Error: iRODS env file missing" } unless $env_file;

    $path = uri_unescape($path); # mdb added 8/15/14 issue 441
    $ENV{irodsEnvFile} = $env_file;  # mdb added 2/17/16 for hypnotoad
    my $cmd = "ils -l '$path' 2>&1"; #"export irodsEnvFile='$env_file'; ils -l '$path' 2>&1"; # mdb changed 2/17/16 for hypnotoad

    my @ils = capture( EXIT_ANY, $cmd );
    if ($EXITVAL) {
        return { error => "Error: ils rc=$EXITVAL" };
    }

    #$path = shift @ils; # mdb removed 8/15/14 issue 441
    shift @ils; # skip first line showing path # mdb added 8/15/14 issue 441

#	if ($path =~ /^ERROR/) { # iRODS error message
#		my $result = { type => 'error', name => $path };
#		return wantarray ? ($result) : [$result];
#	}
    
    # mdb removed 8/15/14 issue 441
    #chomp($path);
    #chop($path);

    my @result;
    foreach my $line (@ils) {
        #print STDERR $line, "\n";
        my ( $type, $backup, $size, $timestamp, $name, $order );
        chomp $line;
        
        if ( $line =~ /^\s*C\-/ ) { # directory or link
            if ( $line =~ /(\S+)\s+linkPoint/ ) { # link # mdb added 9/15/15
                $type = 'link';
                $order = 1;
                ($name) = basename($1);
            }
            else { # directory
                $type = 'directory';
                $order = 1;
                ($name) = basename($line);# $line =~ /([^\/]+)\s*$/; # mdb modified 8/14/14 issue 441
                if ($name) { $name =~ s/\s*$//; $name .= '/'; }
                else       { $name = 'error' }
            }
        }
        else { # file
            $type = 'file';
            $order = 2;
            #( undef, undef, $backup, undef, $size, $timestamp, undef, $name ) = split( /\s+/, $line ); # mdb removed 8/14/14 issue 441
            ($backup, $size, $timestamp, $name) = $line =~ /\s+\S+\s+(\d+)\s+\S+\s+(\S+)\s+(\S+)\s+\S+\s+(.+)/; # mdb added 8/14/14 issue 441
        }
		next if $backup;
		
		my $path2 = catfile($path, $name);
		if ($escape_output) {
		    $type      = escape($type);
            $size      = escape($size);
            $timestamp = escape($timestamp);
            $path2     = escape($path2);
            $name      = escape($name);
		}
		
		my %item = (
            type      => $type,
            name      => $name,
            path      => $path2,
            order     => $order
		);
		$item{size} = $size if $size;
		$item{timestamp} = $timestamp if $timestamp;
        push @result, \%item;
    }
    @result = sort { $a->{order} cmp $b->{order} } @result; # directories/links before files

    #print STDERR Dumper \@result, "\n";
    return { path => $path, items => \@result };
}

sub irods_chksum {
    my $path = shift;
    return 0 unless ($path);

    my $env_file = _irods_get_env_file();
    return unless $env_file;

    my $cmd = "export irodsEnvFile='$env_file' && ichksum $path";
    #print STDERR "cmd: $cmd\n";
    my @output = `$cmd`;
    my ($chksum) = $output[0] =~ /\s*\S+\s+(\S+)/;

    #print STDERR "chksum: $chksum\n";
    return $chksum;
}

sub irods_iget {
    my ( $src, $dest, $opts ) = @_;
    my $no_execute = ( $opts and $opts->{no_execute} ); # mdb added 4/10/14 for REST API, get command name but don't execute
    $src = '' unless (defined $src);
    $dest = '' unless (defined $dest);
    $src = unescape($src); # mdb added 8/15/14 issue 441
    $dest = unescape($dest); # mdb added 8/15/14 issue 441

    #$src =~ s/ /\\ /;
    #$dest =~ s/ /\\ /;

    #print STDERR "irods_iget $src $dest\n";
    my $env_file = _irods_get_env_file();
    return unless $env_file;

    my $cmd = "export irodsEnvFile='$env_file' && iget -fT '$src' '$dest'";
    return $cmd if $no_execute;
    #print STDERR "cmd: $cmd\n";
    my @result = `$cmd`;
    #print STDERR "@result";

    return;
}

sub irods_iput {
    my ( $src, $dest, $opts) = @_;
    my $no_execute = ( $opts and $opts->{no_execute} ); # erb added 9/03/14 for REST API, get command name but don't execute
    my $overwrite = ( $opts and $opts->{overwrite} );
    $overwrite = 1 unless defined $overwrite;

    #print STDERR "irods_iput $src $dest\n";
    my $env_file = _irods_get_env_file();
    return unless $env_file;

    my $cmd = "export irodsEnvFile='$env_file' && iput -T";
       $cmd .= " -f " if $overwrite;
       $cmd .= " $src $dest";

    return $cmd if $no_execute;
    #print STDERR "cmd: $cmd\n";
    my @result = `$cmd`;
    warn Dumper \@result;

    return;
}

sub irods_imeta_add {
    my ( $dest, $pHash ) = @_;
    #print STDERR "irods_imeta_add $dest\n";
    return unless ( defined $pHash and keys %$pHash );

    my $env_file = _irods_get_env_file();
    return unless $env_file;

	while( my($k, $v) = each %$pHash ) {
		if (not $k or not $v) {
			print STDERR "CoGe::Accessory::IRODS::imeta_add: improper key/value pair\n";
			next;
		}

	    my $cmd = "export irodsEnvFile='$env_file' && imeta add -d '" . $dest . "' '" . $k . "' '" . $v . "'";
	    #print STDERR "cmd: $cmd\n";
	    my @result = `$cmd`;
	    #print STDERR "@result";
	}

    return;
}

sub irods_imeta_ls {
    my ( $dest, $attribute ) = @_;

    my $env_file = _irods_get_env_file();
    return unless $env_file;

    my $cmd = "export irodsEnvFile='$env_file' && imeta ls -d '" . $dest . "' '" . $attribute . "'";
    my @result = `$cmd`;

    return \@result;
}

sub irods_imkdir {
    my $path = shift;
    return 'path not specified' unless $path;

    my $env_file = _irods_get_env_file();
    return 'irods env file missing' unless $env_file;

    my $cmd = "export irodsEnvFile='$env_file' && imkdir '" . $path . "'";
    my @result = `$cmd`;
    return $result[0] if scalar @result;
}

# sub irods_irm {
#     my $path = shift;
#     return 'path not specified' unless $path;

#     my $env_file = _irods_get_env_file();
#     return 'irods env file missing' unless $env_file;

#     $ENV{irodsEnvFile} = $env_file;
#     warn $path;
#     my $cmd = "irm -rf '" . $path . "'";
#     my @result = `$cmd`;
#     return $result[0] if scalar @result;
# }

sub irods_get_base_path {
    my $username = shift;
    my $basepath = CoGe::Accessory::Web::get_defaults()->{IRODSDIR};
    $basepath =~ s/\<USER\>/$username/;
    return $basepath;
}

sub irods_set_env {
    $IRODS_ENV = shift;
}

sub _irods_get_env_file {
    my $conf = CoGe::Accessory::Web::get_defaults();
    
    # There are multiple ways to specify the irodsEnvFile.  For main
    # installations runf rom Upstart the default can be used.  For sandbox 
    # installations run locally there must be two irodsEnv files: one for 
    # the local user and one for the www-data user.  The user's version
    # is "irodsEnv_local" as set in api.sh.  The www-data version is "irodsEnv".
    my $env_file = $IRODS_ENV;       # 1: set with irods_set_env()
    $env_file //= $conf->{IRODSENV}; # 2: config file
    $env_file //= $ENV{IRODSENV};    # 3: environment variable
    $env_file //= catfile($conf->{_HOME_PATH}, 'irodsEnv'); # 4: default
    #print STDERR "env_file=$env_file\n";
    
    if ( not defined $env_file or not -e $env_file ) {
        print STDERR "CoGe::Accessory::IRODS: fatal error: iRODS env file missing! env_file=", ($env_file ? $env_file : ''), "\n";
        return;
    }
    return $env_file;
}

1;
