cdef extern from *:
    ctypedef char* const_char_star "const char*"

cimport stdlib

cdef const_char_star dag_format_line = "%s\t%s\t%i\t%i\t%s\t%s\t%i\t%i\t%lf"

cdef extern from "stdio.h":
    int sscanf(char* astr, const_char_star format, ...)

DEF MAXSIZE=256

cdef class DagLine:
    r"""
    given a string of tab-delimited dag output, parse it and create
    an object with the usual attrs.
    """
    cdef public int a_start, a_end, b_start, b_end
    cdef public double evalue
    cdef public char a_accn[MAXSIZE]
    cdef public char b_accn[MAXSIZE]
    cdef public char a_seqid[64]
    cdef public char b_seqid[64]

    def __init__(self, char *sline=NULL):
        if sline == NULL: return
        self.evalue = 1.0
        sscanf(sline, dag_format_line, 
               self.a_seqid, self.a_accn, &self.a_start, &self.a_end,
               self.b_seqid, self.b_accn, &self.b_start, &self.b_end,
               &self.evalue)
        if self.evalue < 1e-250: self.evalue = 1e-250

    def __repr__(self):
        return ("DagLine('%s', '%s')" % (self.a_accn, self.b_accn))


    def __str__(self):
        attrs = ('a_seqid', 'a_accn', 'a_start', 'a_end', 
                 'b_seqid', 'b_accn', 'b_start', 'b_end', 
                 'evalue')
        return "%s\t%s\t%i\t%i\t%s\t%s\t%i\t%i\t%g" \
                % tuple([getattr(self, a) for a in attrs])

    @classmethod
    def from_dict(cls, dict d):
        return _factory(d)

    @classmethod
    def from_pair_dict(cls, dict d):
        return _pair_factory(d)

cdef DagLine _pair_factory(dict d):
    cdef DagLine instance = DagLine.__new__(DagLine)
    cdef dict da = d['A']
    cdef dict db = d['B']
    stdlib.strcpy(instance.a_seqid, da['seqid'])
    stdlib.strcpy(instance.b_seqid, db['seqid'])
    stdlib.strcpy(instance.a_accn, da['accn'])
    stdlib.strcpy(instance.b_accn, db['accn'])

    instance.a_start = da['start']
    instance.b_start = db['start']
    instance.a_end = da['end']
    instance.b_end = db['end']
    instance.evalue = d['evalue']
    return instance


cdef DagLine _factory(dict d):
    cdef DagLine instance = DagLine.__new__(DagLine)
    stdlib.strcpy(instance.a_seqid, d['a_seqid'])
    stdlib.strcpy(instance.b_seqid, d['b_seqid'])
    stdlib.strcpy(instance.a_accn, d['a_accn'])
    stdlib.strcpy(instance.b_accn, d['b_accn'])

    instance.a_start = d['a_start']
    instance.b_start = d['b_start']
    instance.a_end = d['a_end']
    instance.b_end = d['b_end']
    instance.evalue = d['evalue']
    return instance

