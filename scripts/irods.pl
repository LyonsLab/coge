use strict;
use warnings;

use v5.14;

use Carp qw(croak);
use Getopt::Long qw(GetOptions);
use Data::Dumper;

use CoGe::Accessory::IRODS;
use CoGe::Accessory::TDS:

our ($envFile, $command, $metaString, $metaFile, $dest, $src);

GetOptions(
    "cmd|c=s"       => \$command,
    "dest|d=s"      => \$dest,
    "src|s=s"       => \$src,
    "meta|m=s"      => \$metaString,
    "metafile|mf=s" => \$metaFile,   # json file
    "env|e=s"       => \$envFile     # IRODS env file
);

croak "The dest file was not specified" unless $dest;
croak "The command was not specified" unless $command;
croak "Environment file does not exist" unless $envFile and -r $envFile;

CoGe::Accessory::IRODS::irods_set_env($envFile);

# Parse metadata arguments
my %metadata;
if ($metaString) {
    %metadata = (split /[:,]/, $metaString);
}
elsif ($metaFile) {
    croak "The specified metadata file does not exist" unless -r $metaFile;
    my $md = TDS::read($metaFile);
    %metadata = %$md;
}

given($command) {
    when(/submit/i) { # iput and imeta
        croak "The src file was not specified" unless $src;
        croak "The specified src file does not exist" unless -r $src;
        CoGe::Accessory::IRODS::irods_iput($src, $dest);
        CoGe::Accessory::IRODS::irods_imeta($dest, \%metadata) if (keys %metadata);
    }

    when(/fetch/i) { # iget
        croak "The src file was not specified" unless $src;
        CoGe::Accessory::IRODS::irods_iget($src, $dest);
        croak "The file $src could not be fetched" unless -r $dest;
    }
    
    when(/metadata/i) { # imeta
        croak "Medata was not specified" unless (keys %metadata);
        CoGe::Accessory::IRODS::irods_imeta($dest, \%metadata);
    }

    default {
        croak "Unknown command: $command";
    }
}

exit;
