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

    def __repr__(self):
        return self.name

def makeoligo(line):
    '''
    Parameters
    ----------
    line : a line of a csv file describing an oligo (gfp_fc_exp.csv)

    Returns
    -------
    o: a completed oligo, ready for downstream analysis
    '''
    o = oligo(line[4]) #makes object, line[4] = oligo name
    o.gene = line[1]
    o.start = line[2]
    o.end = line[3]
    o.span = line[9]
    o.log2fc = line[5]
    o.padj = line[6]


