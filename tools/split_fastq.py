#!/usr/bin/env python3
import pyfastx
import sys
infile=sys.argv[1] #'02_test_workflow/fastq/gast-T-ChIC-CTCF-i1_R1.fastq.gz'
prefix=sys.argv[2] #'02_test_workflow/outfile'
nlaBC=sys.argv[3] #'/hpc/hub_oudenaarden/vbhardwaj/annotations/cell_barcodes_inhouse/maya_384NLA.bc'
celseqBC=sys.argv[4] #'/hpc/hub_oudenaarden/vbhardwaj/annotations/cell_barcodes_inhouse/celseq2_barcodes.txt'

## check GC content for CS2 barcode filter
def checkGCcontent(seq):
    total_bases = len(seq)
    gc_bases = len([x for x in seq if x == 'C' or x == 'G'])
    gc_frac = float(gc_bases)/total_bases
    return gc_frac
## search matching string and return the position, hamming dist, string
def ham_dist(s1, s2):
    if len(s1) != len(s2):
        raise ValueError("Undefined")
    return sum(ch1 != ch2 for ch1, ch2 in zip(s1, s2))

def search_min_dist(source,search):
    l = len(search)
    index = 0
    min_dist = l
    min_substring = source[:l]
    for i in range(len(source)-l+1):
        d = ham_dist(search, source[i:i+l])
        if d<min_dist:
            min_dist = d
            index = i
            min_substring = source[i:i+l]
    #return (index,min_dist,min_substring)
    return min_dist

nla=[]
with open(nlaBC, 'r') as f:
    nla.extend(f.read().splitlines())
cs2 = [x.split('\t')[0] for x in open(celseqBC, 'r').readlines()]

nl=open(prefix+'.NLA.txt', 'w')
cs=open(prefix+'.CS2.txt', 'w')
both=open(prefix+'.BOTH.txt', 'w')
none=open(prefix+'.NONE.txt', 'w')

for read in pyfastx.Fastq(infile):
    try:
        seq_one=read.seq[0:12]
        seq_two=read.seq[0:14]
        seq_four=read.seq[15:50]
    except UnicodeDecodeError:
         continue
    ## search for nla/cs2 barcode upto hamming dist of 1
    nlaStat = any([search_min_dist(seq_one, b)<=1 for b in nla])
    csStat = any([search_min_dist(seq_two, b)<=1 for b in cs2])
    #nlaStat = any([b in seq_one for b in nla])
    #csStat = any([b in seq_two for b in cs2])

    if nlaStat is True and checkGCcontent(seq_four) >=0.3:
        nlaStat_final = True
    else:
        nlaStat_final = False

    if csStat is True and checkGCcontent(seq_four) <=0.3:
        csStat_final = True
    else:
        csStat_final = False

    if nlaStat_final is True and csStat_final is True:
        both.write(read.name)
        both.write("\n")
        #both.append(read.seq)
    elif nlaStat_final is True and csStat_final is False:
        nl.write(read.name)
        nl.write("\n")
        #nl.append(read.seq)
    elif nlaStat_final is False and csStat_final is True:
        cs.write(read.name)
        cs.write("\n")
        #cs.append(read.seq)
    elif nlaStat_final is False and csStat_final is False:
        none.write(read.name)
        none.write("\n")
        #none.append(read.seq)

for k in [nl, cs, both, none]:
    k.close()
