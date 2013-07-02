import sys
import sqlite3

f = sys.argv[1]
org = sys.argv[2]
fdb = sys.argv[3]
db = sqlite3.connect(fdb)
c = db.cursor()
c.execute("""\
CREATE TABLE simple (
   id INTEGER PRIMARY KEY,
   organism TEXT,
   seqid TEXT,
   name TEXT,
   start INTEGER,
   end INTEGER,
   strand TEXT,
   type TEXT,
   UNIQUE(organism, seqid, start, end, type)
)""")

def geninserts(f, org):
    for line in open(f):
        seqid, ftype, start, end, strand, name = line.rstrip().split(",")
        start, end = map(int, (start, end))
        yield org, seqid, name, start, end, strand, ftype

c.executemany("INSERT INTO simple(organism, seqid, name, start, end, strand, type) VALUES"
              "(?, ?, ?, ?, ?, ?, ?)", geninserts(f, org))
db.commit()
db.close()
