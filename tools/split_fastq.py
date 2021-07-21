#!/usr/bin/env python3
import sys
import multiprocessing as mp
import click
import gzip
from Bio.SeqIO.QualityIO import FastqGeneralIterator
from Bio.Seq import Seq
from splitfastq import *
from functools import partial

# Note: about 1 in 10k reads have NLA+T7promoter and should be classified as CS2
# instead of NLA (current situation), but I don't want to bother with the extra search time
def process_reads(handle, fh_prefix, nla, cs2):
    read, seq = handle
    cs=[]
    nl=[]
    bo=[]
    no=[]
    try:
        seq_one=seq[0:12]
        seq_two=seq[0:14]
        #seq_three=seq[0:50]
        seq_four=seq[15:40]
    except UnicodeDecodeError:
        pass
    # search for nla/cs2 barcode upto hamming dist of 1
    nlaStat = any([search_min_dist(seq_one, b)<=1 for b in nla])
    csStat = any([search_min_dist(seq_two, b)<=1 for b in cs2])
    # Loook for T7 promoter in the beginning of the fastq if they are CS2+NLA
    #prom='GCCGGTAATACGACTCACTATAGGGAGTTCTACAGTCCGACGATC'
    #csStat_extended = search_min_dist(seq_three, prom) <=2

    if nlaStat is True and checkGCcontent(seq_four) >=0.3:
        nlaStat_final = True
    else:
        nlaStat_final = False# this can still be an NLA+CS2+polyA situation

    if csStat is True and checkGCcontent(seq_four) <0.3:
        csStat_final = True
    #elif csStat_extended is True:
    #    csStat_final = True#NLA+CS2+polyA situation
    else:
        csStat_final = False
    # if you find NLA+CS2 barcodes one after another, check for GC content and
    # re-classify them as cs2

    if nlaStat_final is True and csStat_final is True:
        bo.append(read)
    elif nlaStat_final is True and csStat_final is False:
        nl.append(read)
    elif nlaStat_final is False and csStat_final is True:
        cs.append(read)
    elif nlaStat_final is False and csStat_final is False:
        no.append(read)

    return {'CS2':cs, 'NLA': nl, 'BOTH': bo, 'NONE': no}

@click.command()
@click.option('--infile')#, type=click.File(mode='r'))
@click.option('--nla_bc')#, type=click.File(mode='r'))#'/hpc/hub_oudenaarden/vbhardwaj/annotations/cell_barcodes_inhouse/maya_384NLA.bc'
@click.option('--celseq_bc')#, type=click.File(mode='r'))#'/hpc/hub_oudenaarden/vbhardwaj/annotations/cell_barcodes_inhouse/celseq2_barcodes.txt'
@click.option('--prefix')#, type=click.STRING)
@click.option('--ncpus')
def run(infile, nla_bc, celseq_bc, prefix, ncpus):
    # get barcodes
    nla=[]
    with open(nla_bc, 'r') as f:
        nla.extend(f.read().splitlines())
    cs2 = [x.split('\t')[0] for x in open(celseq_bc, 'r').readlines()]

    def get_data(in_handle):
        for rec_id, seq, _ in FastqGeneralIterator(in_handle):
            yield rec_id, str(Seq(seq))

    func=partial(process_reads, fh_prefix=prefix, nla=nla, cs2=cs2)
    with gzip.open(infile, 'rt') as input_handle:
        p = mp.Pool(processes=int(ncpus))
        g = p.map(func, get_data(input_handle))
    ## combine and print results
    for n in ['CS2', 'NLA', 'BOTH', 'NONE']:
        full_list=[x[n] for x in g if x[n]!=[]]
        print('{}_{}: {}'.format(prefix, n, len(full_list)))
        with open(prefix+'.'+n+".txt", "w",encoding="utf-8",errors='ignore') as f:
            for rec in full_list:
                f.write(rec[0]+'\n')
        f.close()

if __name__ == "__main__":
    run()
