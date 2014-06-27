if [ ! -f /opt/apache/CoGe/tmp/GEvo/${1}.sqlite ]
then
    echo /opt/apache/CoGe/tmp/GEvo/${1}.sqlite does not exist
    exit;
fi

echo "adding new files to the directory"
svn up /opt/apache/CoGe/gobe/
svn rm /opt/apache/CoGe/gobe/tmp/*
cp /opt/apache/CoGe/tmp/GEvo/${1}*.png /opt/apache/CoGe/gobe/tmp/
cp /opt/apache/CoGe/tmp/GEvo/${1}*.faa /opt/apache/CoGe/gobe/tmp/
cp /opt/apache/CoGe/tmp/GEvo/${1}*.log /opt/apache/CoGe/gobe/tmp/
cp /opt/apache/CoGe/tmp/GEvo/${1}.sqlite /opt/apache/CoGe/gobe/tmp/
svn add /opt/apache/CoGe/gobe/tmp/*

perl -p -i -e "s/GEv[^']+/$1/" /opt/apache/CoGe/gobe/index.html
N=`ls /opt/apache/CoGe/tmp/GEvo/${1}*.png | wc -l`
echo $N
perl -p -i -e "s/\('n'\)\s*\|\|\s*[^\)]+/('n') || ${N}/" /opt/apache/CoGe/gobe/index.html
