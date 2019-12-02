import os
import glob
import yaml
### snakemake_workflows initialization ########################################
maindir = os.path.dirname(os.path.dirname(workflow.basedir))
#workflow_rscripts=os.path.join(maindir, "shared", "rscripts")
#shell.executable("/bin/bash")
## some internal functions  ###################################################
def load_configfile(configfile):
   with open(configfile, "r") as f:
       config = yaml.load(f, Loader=yaml.FullLoader)
   return(config)

def set_condaEnv():
    return{'CONDA_SHARED_ENV': 'env.yaml'}

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

# update envs
globals().update(set_condaEnv())
# load config file
globals().update(load_configfile(workflow.overwrite_configfile))

## load samples
infiles = sorted(glob.glob(os.path.join(indir, '*'+ext)))
samples = get_sample_names(infiles,ext,reads)

### include modules of other snakefiles ########################################
################################################################################
include: os.path.join(workflow.basedir, "rules", "fastq_map.snakefile")

if method in ['chic-taps', 'nla-taps']:
    include: os.path.join(workflow.basedir, "rules", "methCall.snakefile")
    if lambda_phage:
        include: os.path.join(workflow.basedir, "rules", "methCall_lambda_phage.snakefile")

if method in ['chic-taps', 'chic']:
    include: os.path.join(workflow.basedir, "rules", "dedup_and_qc.snakefile")
    include: os.path.join(workflow.basedir, "rules", "peakcall.snakefile")

include: os.path.join(workflow.basedir, "rules", "QC.snakefile")

### conditional/optional rules #################################################
################################################################################
def run_Trimming(trim):
    if trim:
        file_list = [
        expand("FASTQ_trimmed/{sample}{read}.fastq.gz", sample = samples, read = reads),
        expand("QC/FastQC_trimmed/{sample}{read}_fastqc.html", sample = samples, read = reads)
        ]
        return(file_list)
    else:
        return([])

def meth_check(type=method):
    file_list = []
    if type in ['chic-taps', 'nla-taps']:
        file_list.extend([
        expand("tagged_bam/{sample}.bam", sample = samples),
        #expand("meth_calls/{sample}_methylation.bed.gz", sample = samples),
        expand("coverage/{sample}.methRatio.bw", sample = samples),
        expand("meth_counts/{sample}_CpG_binCounts.meth.csv", sample = samples),
        expand("meth_counts/{sample}_CpG_binCounts.unmeth.csv", sample = samples),
        expand("QC/scMultiOmics/{sample}_QCplots/ConversionMatrix.conversions.png", sample = samples)
        ])
        if lambda_phage:
            file_list.extend([
            expand("tagged_bam/lambda_phage/{sample}.bam", sample = samples),
            expand("tagged_bam/lambda_phage/{sample}.bam.bai", sample = samples),
            expand("meth_calls/lambda_phage/{sample}_stats.txt", sample = samples),
            expand("meth_calls/lambda_phage/{sample}.methRatio.bw", sample = samples),
                            "QC/lambda_stats.txt"])
        if len(samples) > 1:
            file_list.extend(["QC/bwSummary_methRatio_10kBins.npz",
                            "QC/cor-spearman_methRatio_10kBins.png"])

    if type in ['chic', 'chic-taps']:
        file_list.extend([
        expand("dedup_bam/{sample}.bam", sample = samples),
        expand("dedup_bam/{sample}.bam.bai", sample = samples),
        expand("coverage/{sample}_dedup.cpm.bw", sample = samples),
        "QC/featureEnrichment.png",
        "QC/featureEnrichment_biotype.png",
        expand("counts/{sample}.per_barcode.tsv", sample = samples),
        expand("counts/{sample}.windows_total.tsv", sample = samples),
        expand("counts/{sample}.windows_per_barcode.tsv", sample = samples),
        expand("macs2_peaks/{sample}_peaks.narrowPeak", sample = samples),
        expand("macs2_peaks/{sample}_peaks.bed", sample = samples),
        ])
        if len(samples) > 1:
            file_list.extend(["QC/bwSummary_10kBins.npz",
                            "QC/cor-spearman_10kBins.png"])

    if type == 'chic-taps':
        file_list.extend([expand("counts/{sample}_peaks.csv", sample = samples)])
    return(file_list)

### main rule ##################################################################
################################################################################
localrules: FASTQ1, FASTQ2
rule all:
    input:
        expand("FASTQ/umiTrimmed_{sample}{read}.fastq.gz", sample = samples, read = reads),
        run_Trimming(trim),
        expand("QC/FastQC/{sample}{read}_fastqc.html", sample = samples, read=reads),
        expand("bwa_mapped/{sample}.bam", sample = samples),
        expand("bwa_mapped/{sample}.bam.bai", sample = samples),
        meth_check(),
        "QC/multiqc_report.html"


### execute after workflow finished ############################################
################################################################################
onsuccess:
    if "verbose" in config and config["verbose"]:
        print("\n--- openTAPS workflow finished successfully! ------------------\n")
onerror:
    print("\n !!!! ERROR in openTAPS workflow! !!!!\n")
