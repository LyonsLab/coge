This package was generated automatically by CoGe (genomevolution.org).

The contents vary depending on the type of data:  quantitative, polymorphism, or alignment.

QUANTITATIVE
============
This package contains the raw and formatted data for a quantitative experiment.

The raw data filename ends with a .csv extension.

The other files form the FastBit (https://sdm.lbl.gov/fastbit/) indexed database of the raw data file.

To query the FastBit database, install FastBit and use the “ibis” command:
    ibis -d <directory> -q <query>

Example:
    ibis -d . -q \"select chr,start,stop,strand,value1,value2 where chr=1 limit 10”


POLYMORPHISM
============
This package contains the raw and formatted data for a polymorphism experiment.

The raw data filename ends with a .csv extension.

The other files form the FastBit (https://sdm.lbl.gov/fastbit/) indexed database of the raw data file.

To query the FastBit database, install FastBit and use the “ibis” command:
    ibis -d <directory> -q <query>

Example:
    ibis -d . -q \"select chr,start,stop,type,id,ref,alt,qual,info where chr=1 limit 10”


ALIGNMENT
=========
This package contains the raw data file and index file for an alignment experiment.

The raw data filename ends with a .bam extension.

The index filename ends with a .bam.bai extension.

Indexed BAM files can be queried using SAMTools (http://samtools.sourceforge.net/).
    samtools view <filename> <region>

Example:
    samtools view test.bam “1:1-10000”


