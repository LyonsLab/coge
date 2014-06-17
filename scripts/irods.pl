use strict;
use warnings;

use v5.14;

use Carp qw(croak);
use Getopt::Long qw(GetOptions);
use Data::Dumper;

use CoGe::Accessory::IRODS;

our ($envFile, $command, $metaString, $dest, $src);

GetOptions(
    "cmd|c=s"  => \$command,
    "dest|d=s" => \$dest,
    "src|s=s"  => \$src,
    "meta|m=s" => \$metaString,
    "env|e=s"  => \$envFile
);

croak "The src file was not specified" unless $src;
croak "The dest file was not specified" unless $dest;
croak "The command was not specified" unless $command;
croak "Environment file does not exist" unless $envFile and -r $envFile;

CoGe::Accessory::IRODS::irods_set_env($envFile);

given($command) {
    when(/submit/i) {
        #Exit if the file does not exist
        croak "The input file specified does not exist" unless -r $src;

        my %meta;

        if ($metaString) {
            %meta = (split /[:,]/, $metaString);
        }
        CoGe::Accessory::IRODS::irods_iput($src, $dest);
        CoGe::Accessory::IRODS::irods_imeta($dest, \%meta) if %meta;
    }

    when(/fetch/i) {
        CoGe::Accessory::IRODS::irods_iget($src, $dest);
        croak "The file $src could not be fetched" unless -r $dest;
    }

    default {
        croak "Unknown command $command was given";
    }
}

exit;
