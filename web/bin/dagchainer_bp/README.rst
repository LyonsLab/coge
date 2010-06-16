this is a copy of dagchainer:
http://dagchainer.sourceforge.net/


Changes
-------
1. the calling script has been converted to python. dag_chainer.py
   the options should be identical.

2. it removes the need for temporary files (per chr-chr pair) by 
   reading/writing to and from stdout.

3. does less filtering than the perl eqivalent.

4. the Makefile was update to include compile flags that make the resulting
   executable up to 2 times as fast.

5. added a find_nearby.py script. dagchainer finds the *optimal* path through a list
   of hits, but we may also want to include nearby hits that weren't on the optimial
   path. this uses scipy's cKDTree to fast queries to find the close-by points.
   As such scipy>=0.8 is required to use this script.


Example
-------

1. build the dagchainer executable and the python libs ::

  make

2. an example shell script using the atat.dag from arabidopsis ::

    IN=atat.dag
    OUT=a.out
    python dagtools/evalue_adjust.py -c 4 $IN | python dag_chainer.py -g 40000 -D 160000 -A 4 --merge $OUT -
    python dagtools/find_nearby.py -d 20000 --diags ${OUT}.merge --all ${INFILE} > ${OUT}.merge.all
    python dagtools/plot_dag.py -q athaliana_3 -s athaliana_5 -p dotplot.png -d ${OUT}.merge.all


3. and view the output in dotplot.png. the ${OUT}.merge file has the diagonals
   merged, the ${OUT}.merge.all includes points that are near the diagonals, 
   but not on what dagchainer considers the optimal path.


Library
-------
the dag_chainer.py script runs the dagchainer c++ executable. It also uses tools
from the library in the dagtools module in this directory. The most useful function
is reading a dagline into a python object ::

    >>> from dagtools import DagLine
    >>> for line in open('somefile.dag'):
    ...    if line[0] != "#":
    ...    d = DagLine(line)
    ...    print d.a_accn, d.a_start, b.b_accn, d.b_start

which will loop across the file and access the attributes (available 
attributes are (prefixed by either "a_" or "b_") start, stop, accn, seqid.
and "evalue" (with no prefix).
This DagLine object is essentially a C struct which is filled in C, and
as such it is faster to create and use, and uses less memory than a python
equivalent.
Other scripts used in the Example_ section above can also be used as libraries.
For example plotting a chromosome pair of dag hits is done in the plot_dag.py
library.
The module can be installed via ::

    $ cython dagtools/cdagline.pyx
    $ sudo python setup.py install


