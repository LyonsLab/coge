import sys
import sqlite3

f = sys.argv[1]
db = sqlite3.connect(f + ".sqlite")
c = db.cursor()
c.execute("""\
CREATE TABLE simple (
   id INTEGER PRIMARY KEY,
   seqid TEXT,
   name TEXT,
   start INTEGER,
   end INTEGER,
   strand TEXT,
   type TEXT,
   UNIQUE(seqid, start, end, type)
)""")

def geninserts(f):
    for line in open(f):
        seqid, ftype, start, end, strand, name = line.rstrip().split(",")
        start, end = map(int, (start, end))
        yield seqid, name, start, end, strand, ftype

c.executemany("INSERT INTO simple(seqid, name, start, end, strand, type) VALUES"
              "(?, ?, ?, ?, ?, ?)", geninserts(f))
db.commit()
db.close()
