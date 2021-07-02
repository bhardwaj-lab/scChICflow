#!/usr/bin/env python3
import pyfastx
import sys
import multiprocessing as mp
import click
from splitfastq import *

def process_reads(infile, nla, cs2):
    d = {'NLA':[],
         'CS2':[],
         'BOTH':[],
         'NONE':[]
        }
    n=0
    for read in pyfastx.Fastq(infile):
        #if n>99:
        #    break
        try:
            seq=read.seq[0:25]
            seq_one=read.seq[0:12]
            seq_two=read.seq[0:14]
            seq_four=read.seq[15:50]
        except UnicodeDecodeError:
            pass
        ## search for nla/cs2 barcode upto hamming dist of 1
        nlaStat = any([search_min_dist(seq_one, b)<=1 for b in nla])
        csStat = any([search_min_dist(seq_two, b)<=1 for b in cs2])

        if nlaStat is True and checkGCcontent(seq_four) >=0.3:
            nlaStat_final = True
        else:
            nlaStat_final = False

        if csStat is True and checkGCcontent(seq_four) <=0.3:
            csStat_final = True
        else:
            csStat_final = False

        if nlaStat_final is True and csStat_final is True:
            d['BOTH'].append(read.name)
        elif nlaStat_final is True and csStat_final is False:
            d['NLA'].append(read.name)
        elif nlaStat_final is False and csStat_final is True:
            d['CS2'].append(read.name)
        elif nlaStat_final is False and csStat_final is False:
            d['NONE'].append(read.name)
        #n+=1
    return d

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

    #ensure index file has been created in main process
    fa = pyfastx.Fastq(infile)
    #get sequence counts
    c = len(fa)
    #start the process pool
    pool = mp.Pool(int(ncpus))
    #add five task workers to run
    out=pool.apply_async(process_reads, args=(infile, nla, cs2)).get()
    #wait for tasks to finish
    pool.close()
    pool.join()

    ## write results
    for n in ['CS2', 'NLA', 'BOTH', 'NONE']:
        li = out[n]
        with open(prefix+'_'+n+'.txt', 'w') as f:
            for l in li:
                f.write(l+'\n')
        f.close()

if __name__ == "__main__":
    run()
