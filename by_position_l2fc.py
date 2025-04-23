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
    py_start_idx = o.start - 1 #reindexing 0 to 259
    py_end_idx = o.end  
    o.pos_inc = list(range(py_start_idx, py_end_idx))

def l2fc_of_pos(an_oligo, output_list, oligo_len = 260):
    '''
    Parameters
    ----------
    an_oligo : an object of type oligo    
    output_list :  [[] for i in range(oligo_len)] -> list to put gene, l2fc, and
    span

    Returns
    -------
    list filled out with all the log2FCs for each oligo a position is
    included in; list of lists:
        gene, log2fc, span 
    tidyr style data
    '''
    out = [an_oligo.gene, an_oligo.span, an_oligo.log2fc, an_oligo.padj] # to be added to list
    for i in an_oligo.pos_inc:
        output_list[i].append(out)

    

#===================== read in data ==============================#
with open("gfp_fc_exp.csv") as input:
    oligos = []
    gfp_l2fc = csv.reader(input, delimiter = ",")
    next(gfp_l2fc) # skip header

    for line in gfp_l2fc:
        o = oligo(name = line[4])
        o1 = makeoligo(oligo=o, line=line)
        oligos.append(o)
#===================== convert lines of csv file to list of lists=====#
l2fc_by_base = [[] for i in range(260)]
for i in range(0, len(oligos)):
    o = oligos[i]
    l2fc_of_pos(an_oligo = o, output_list = l2fc_by_base)
#===================== export ========================#
header = ["position", "gene", "span", "log2fc", "padj"]
with open("l2fc_by_pos_span_suff.csv", "w") as output:
    for i in header:
        output.write(i + ",")
    output.write("\n")

    for i in range(260):
        pos = i + 1 # base position, index-1/ R style indexing
        o = l2fc_by_base[i]
        for j in o:    #  position         gene             span
            output.write(str(pos) + "," + str(j[0]) + "," + str(j[1]) + "," 
                         + str(j[2]) + "," + str(j[3]) + '\n')
#                           l2fc                padj