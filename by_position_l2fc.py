'''
Per base localization activity for sufficency MPRA
- avg log2 fc for each included base by window size
- avg log2 fc for each excluded base by window size?
* starting from pos-shuffle.py code
'''
import csv

class oligo:
    def __init__(self, name):
        self.name = name
        self.gene = ''
        self.start = []
        self.end = []
        self.span = []
        self.log2fc = []
        self.padj = []
        self.pos_inc = []

    def __repr__(self):
        return self.name

def makeoligo(oligo, line):
    '''
    Parameters
    ----------
    line : a line of a csv file describing an oligo (gfp_fc_exp.csv)
    oligo : an initialized oligo

    Returns
    -------
    o: a completed oligo, ready for downstream analysis
    '''
    o.gene = line[1]
    o.start = int(line[2])
    o.end = int(line[3])
    o.span = line[9]
    o.log2fc = line[5]
    o.padj = line[6]
    py_start_idx = o.start - 1 #reindexinf
    py_end_idx = o.end + 1 # adding 1 to compensate for range works in python vs R
    o.pos_inc = list(range(py_start_idx, py_end_idx))

def l2fc_of_pos(an_oligo, output_list, oligo_len = 260):
    s = int(an_oligo.start) - 1 #reindex -> 0 to 259
    e = int(an_oligo.end) - 1
    

#===================== read in data ==============================#
with open("gfp_fc_toy.csv") as input:
    oligos = []
    gfp_l2fc = csv.reader(input, delimiter = ",")
    next(gfp_l2fc) # skip header

    for line in gfp_l2fc:
        o = oligo(name = line[4])
        o1 = makeoligo(oligo=o, line=line)
        oligos.append(o)

l2fc_by_base = [ [] for i in range(260)]
for i in range(0, len(oligos)):
    o = oligos[i]
