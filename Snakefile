import os
import glob
import yaml
### snakemake_workflows initialization ########################################
maindir = os.path.dirname(os.path.dirname(workflow.basedir))
#workflow_rscripts=os.path.join(maindir, "shared", "rscripts")

## some internal functions  ###################################################
def load_configfile(configfile):
   with open(configfile, "r") as f:
       config = yaml.load(f, Loader=yaml.FullLoader)
   return(config)

def get_sample_names(infiles, ext, reads):
    """
    Get sample names without file extensions
    """
    s = set()
    lext = len(ext)
    l0 = len(reads[0])
    l1 = len(reads[1])
    for x in infiles:
        x = os.path.basename(x)[:-lext]
        if x.endswith(reads[0]):
            x = x[:-l0]
        elif x.endswith(reads[1]):
            x = x[:-l1]
        else:
            continue
        s.add(x)
    return sorted(list(s))

# load config file
globals().update(load_configfile(workflow.overwrite_configfile))

## load samples
infiles = sorted(glob.glob(os.path.join(indir, '*'+ext)))
samples = get_sample_names(infiles,ext,reads)

### include modules of other snakefiles ########################################
################################################################################
include: os.path.join(workflow.basedir, "rules", "fastq.snakefile")
include: os.path.join(workflow.basedir, "rules", "map_and_methCall.snakefile")

### conditional/optional rules #################################################
################################################################################
def run_Trimming(trim):
    if trim:
        file_list = [
        expand("FASTQ_trimmed/{sample}{read}.fastq.gz", sample = samples, read = reads),
        expand("FASTQ_trimmed/FastQC/{sample}{read}_fastqc.html", sample = samples, read = reads)
        ]
        return(file_list)
    else:
        return([])

### main rule ##################################################################
################################################################################
localrules: FASTQ1, FASTQ2
rule all:
    input:
        expand("FASTQ/umiTrimmed_{sample}{read}.fastq.gz", sample = samples, read = reads),
        run_Trimming(trim),
        expand("FastQC/{sample}{read}_fastqc.html", sample = samples, read=reads),
        expand("bwa_mapped/{sample}.bam", sample = samples),
        expand("tagged_bam/{sample}.bam", sample = samples),
        expand("meth_calls/{sample}_allC.bed.gz", sample = samples),
        expand("meth_calls/{sample}.methCpG.bw", sample = samples)

### execute after workflow finished ############################################
################################################################################
onsuccess:
    if "verbose" in config and config["verbose"]:
        print("\n--- openTAPS workflow finished successfully! --------------------------------\n")
onerror:
    print("\n !!! ERROR in openTAPS workflow! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n")
