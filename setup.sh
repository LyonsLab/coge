#!/bin/bash

coge_root=$PWD
server_name="coge-web"
server_address="localhost"
DB_HOST="localhost"
DB_NAME="coge"
DB_USER="coge_web"
DB_PORT="3306"

function get_pass {
    read -s -p "Enter password: " pass1; echo
    read -s -p "Retype Password: " pass2; echo

    while [[ "$pass1" != "$pass2" ]]; do
        echo -e "\nPasswords did not match."
        read -s -p "Enter password: " pass1; echo
        read -s -p "Retype Password: " pass2; echo
    done

    eval "$1"=$pass1
}

function check_requirements {
    echo -e "Stage 1: Checking for CoGe requirements...";
    for prog in "$@"
    do
        echo -n "Checking for $prog"
        type $prog > /dev/null 2>&1 ||
            { echo -e >&2 " [Not Found].\nPlease install $prog." && exit; }
        echo " [Found]"
    done
    echo "All requirements were found."
}

function build_perl_environment {
    echo -e "\nStage 2: Configuring a local copy of perl..."
    # Install a local copy of perlbrew
    export PERLBREW_ROOT="$PWD/perl"
    export PERLBREW_HOME="$PERLBREW_ROOT/.perl"
    perl_version="perl-5.12.5"

    curl -kL http://install.perlbrew.pl | bash
    source "$PERLBREW_ROOT/etc/bashrc"
    type perlbrew > /dev/null 2>&1 ||
        { echo >&2 "perlbrew was not successfully installed." && exit; }

    perlbrew init
    perlbrew install $perl_version
    perlbrew switch $perl_version
    perlbrew lib create $perl_version
    perlbrew use "$perl_version@coge"

    curl -L http://cpanmin.us | perl - App::cpanminus

    echo "$PERLBREW_PATH/cpanm"

# FIXME: Find best way to check if module is locally installed
#    type "$PERLBREW_PATH/cpanm" > /dev/null 2>&1 ||
#        { echo >&2 "cpanm was not successfully installed." && exit; }

    # Install third-party modules
    while read line;
    do
        perl -$line -e || cpanm $line
    done < "modules.txt"

    # Install local modules
    pushd "$coge_root/coge-db"
        perl Makefile.PL
        make test && make install
    popd

    echo "Perl has been properly configured."
}

function configure_apache2 {
    echo -e "\nStage 3: Configuring apache..."
    dir="/etc/apache2/sites-available"

    #FIXME: Make these options configurable
    #read -p "Specify the root directory: " coge_root
    #read -p "Specify the server name: " server_name

    echo "<VirtualHost *:80>
    ServerName $server
    DocumentRoot $coge_root/$server_name
    ScriptAlias $server/ \"$coge_root/$server_name/\"
    Alias /$server \"$coge_root/$server\"
    <Directory \"$coge_root/$server_name\">
        # PerlInitHandler Apache2::Reload
        Options Includes Indexes ExecCGI FollowSymLinks
        AllowOverride None
        SetEnv HOME \"$coge_root/$server_name/\"
        Order allow,deny
        Allow from all
    </Directory>
    </VirtualHost>" > $coge_root/coge.site

    echo "Adding coge.site to sites-available."
    echo "Requesting administrative privileges."

    sudo mkdir -p $dir
    sudo mv "coge.site" $dir
    echo "Added coge.site to sites-available"
}

function setup_workspace {
    echo -e "\nStage 4: Creating directory structure..."

    tools="CoGeAlign,FeatList,GenomeList,MotifView,SynFind,TreeView,
                  CoGeBlast,FeatMap,GEvo,OrgView,SynMap,tRNA"

    data="pblast/db,experiments,genomic_sequence,perl,
                  bed,distrib,fasta,image_cache,private"

    tool_dirs=$(echo $tools | tr -d " |\t|\r|\n" | tr -s "," "\n")
    data_dirs=$(echo $data | tr -d " |\t|\r|\n" | tr -s "," "\n")

    for folder in ${tool_dirs[@]}
    do
        mkdir -p -m 777 $coge_root/tmp/$folder
    done

    for folder in ${data_dirs[@]}
    do
        mkdir -p -m 777 $coge_root/data/$folder
    done

    mkdir -p -m 777 $coge_root/diags

    echo "Finished setting up directory structure"
}

function import_database {
    echo -e "\nStage 5: Attempting to import the CoGe database.."
    read -p "MySQL Username: " DB_USER
    get_pass "DB_PASS"

    while [[ -n `mysql -u $DB_USER -p $DB_PASS` ]]; do
        echo "Can't connect, please retry."

        read -p "MySQL Username: " DB_USER;
        get_pass "DB_PASS"
    done

    echo "Importing Database..."
    if [[ -n `mysql -u $DB_USER -p $DB_PASS  coge < schema.sql` ]]; then
        echo "Import failed."
        exit
    fi

    echo "CoGe database was imported successfully."
}

function setup_config {
    echo -e "\nStage 6: Generating configuration file..."

    config="##This is a configuration file for CoGe.
    ##Key Value pairs:
    ##<NAME>    <PATH>

    #database configuration
    DBNAME $DB_NAME
    DBHOST $DB_HOST
    DBPORT $DB_PORT
    DBUSER $DB_USER
    DBPASS $DB_PASS

    #basedir for coge
    COGEDIR $coge_root/coge_web/

    #bin dir for coge's programs
    BINDIR $coge_root/$server_name/bin/

    #data dir for coge's programs
    DATADIR $coge_root/$server_name/data/

    #dir for pair-wise whole genome comparisons (e.g. SynMap)
    DIAGSDIR $coge_root/$server_name/diags/

    #fasta dir
    FASTADIR $coge_root/$server_name/data/fasta/

    #TMPL dir for coge's web page templates
    TMPLDIR $coge_root/$server_name/tmpl/

    #temp dir for coge
    TEMPDIR $coge_root/$server_name/tmp/

    #Base URL for web-site
    URL $server_address/

    #URL for temp directory
    TEMPURL $server_name/tmp/

    #blast style scoring matrix dirs
    BLASTMATRIX $coge_root/$server_name/data/blast/matrix/

    #blastable DB
    BLASTDB $coge_root/$server_name/data/blast/db/

    #directory for bed files
    BEDDIR $coge_root/$server_name/data/bed/

    #servername for links
    SERVER $server_address/$server_name/

    #directory for caching genome browser images
    IMAGE_CACHE $coge_root/$server_name/data/image_cache/"

    echo $config > $coge_root/app.conf
    echo "Successfully generated the configuration file."
}

current=$PWD
cd $(dirname "${BASH_SOURCE[0]}")
script=$PWD
cd $current
current="$PWD/scripts"

if [ "$current" != "$script" ]; then
    echo "Must be executed from the projects root directory."
    exit
fi

# Checks requirements
check_requirements curl git apachectl mysql sudo

#running=$(ps aux | grep mysql | grep -c -v "grep mysql")
#if [ -n $running ]; then
#    echo "MySQL is not running, please start mysql."
#    exit
#fi

# Installs Perlbrew, Perl and required modules
build_perl_environment

# Adds configuration to apache
configure_apache2

# Build directory structure
setup_workspace

# Setup Database
import_database

# Build configuration
setup_config

echo "Finished"
