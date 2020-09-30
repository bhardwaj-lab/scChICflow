rule taps_tagger:
    input:
        bam = "bwa_mapped/{sample}.bam",
        idx = "bwa_mapped/{sample}.bam.bai",
        genome = genome_fasta
    output:
        bam = "tagged_bam/{sample}.bam",
        bai = "tagged_bam/{sample}.bam.bai",
        bed = "meth_calls/{sample}_methylation.bed"
    params:
        method = 'chic' if method == 'chic-taps' else 'nla',
        min_mq = min_mapq,
        #cluster = '--cluster' if cluster else '',
        context = '-context '+bedContext if bedContext else ''
    log:
        out = "logs/taps_tagger_{sample}.out",
        err = "logs/taps_tagger_{sample}.err"
    threads: 1
    conda: CONDA_SHARED_ENV
    shell:#{params.cluster}
        "tapsTagger.py {params.context} \
        -ref {input.genome} -method {params.method} -bed {output.bed} \
        -o {output.bam} -min_mq {params.min_mq} {input.bam} \
        > {log.out} 2> {log.err}"
#--cluster 
## methylation counts per bin per cell
rule meth_bincounts:
    input:
        bam = "tagged_bam/{sample}.bam",
        bai = "tagged_bam/{sample}.bam.bai"
    output:
        csv = "meth_counts/{sample}_CpG_binCounts.meth.csv"
    params:
        binsize = binSize
    log: "logs/meth_bincounts_{sample}.log"
    threads: 1
    conda: CONDA_SHARED_ENV
    shell:
        "bamToCountTable.py -sampleTags SM -featureTags sZ -byValue sZ \
        -joinedFeatureTags reference_name -bin {params.binsize} -o {output.csv} {input.bam} > {log} 2>&1"

## all CpG counts per bin per Cell
rule unmeth_bincounts:
    input:
        bam = "tagged_bam/{sample}.bam",
        bai = "tagged_bam/{sample}.bam.bai"
    output:
        csv = "meth_counts/{sample}_CpG_binCounts.unmeth.csv"
    params:
        binsize = binSize
    log: "logs/all_bincounts_{sample}.log"
    threads: 1
    conda: CONDA_SHARED_ENV
    shell:
        "bamToCountTable.py -sampleTags SM -featureTags sz -byValue sz \
        -joinedFeatureTags reference_name -bin {params.binsize} -o {output.csv} {input.bam} > {log} 2>&1"

## ----------------------------- Bulk-QC --------------------------- ##
rule plotCorrelation_meth:
    input: "QC/bwSummary_methRatio_10kBins.npz"
    output: "QC/cor-spearman_methRatio_10kBins.png"
    log: "logs/plotCorrelation_CpG.log"
    threads: 1
    conda: CONDA_SHARED_ENV
    shell:
        "plotCorrelation -p heatmap -c spearman \
        --skipZeros --removeOutliers \
        --plotNumbers -in {input} -o {output} > {log} 2>&1"

rule scMultiOmics_qc:
    input:
        bam = "tagged_bam/{sample}.bam",
        bai = "tagged_bam/{sample}.bam.bai"
    output: "QC/scMultiOmics/{sample}_QCplots/ConversionMatrix.conversions.png"
    params:
        outdir = "QC/scMultiOmics/{sample}_QCplots",
        stat = 'meth-stats'
    log: "logs/scMultiOmics_qc_{sample}.log"
    threads: 1
    conda: CONDA_SHARED_ENV
    shell:
        "libraryStatistics.py --nort -t {params.stat} \
        --plotsOnly -o {params.outdir} {input.bam} > {log} 2>&1"

## ----------------------------- QC data per cell for all samples --------------------------- ##

rule fragsPerCell:
    input:
        bam = "tagged_bam/{sample}.bam",
        bai = "tagged_bam/{sample}.bam.bai"
    output:
        csv = "scQC/fragsPerCell/{sample}.csv"
    params:
        binsize = binSize
    log: "logs/fragsPerCell_{sample}.log"
    threads: 1
    conda: CONDA_SHARED_ENV
    shell:
        "bamToCountTable.py --dedup -sampleTags SM -featureTags RF \
        -o {output.csv} {input.bam} > {log} 2>&1"

rule methCountsPerCell:
    input:
        bam = "tagged_bam/{sample}.bam",
        bai = "tagged_bam/{sample}.bam.bai"
    output:
        csv = "scQC/methCountsPerCell/{sample}_allMeth.csv"
    log: "logs/methCountsPerCell_{sample}.log"
    threads: 1
    conda: CONDA_SHARED_ENV
    shell:
        "bamToCountTable.py --dedup -sampleTags SM -featureTags MC \
        -o {output.csv} {input.bam} > {log} 2>&1"

rule CpGmethCountsPerCell:
    input:
        bam = "tagged_bam/{sample}.bam",
        bai = "tagged_bam/{sample}.bam.bai"
    output:
        csv = "scQC/methCountsPerCell/{sample}_CpGmeth.csv"
    log: "logs/CpGmethCountsPerCell_{sample}.log"
    threads: 1
    conda: CONDA_SHARED_ENV
    shell:
        "bamToCountTable.py --dedup -sampleTags SM -featureTags sZ \
        -o {output.csv} {input.bam} > {log} 2>&1"

rule CpGNOmethCountsPerCell:
    input:
        bam = "tagged_bam/{sample}.bam",
        bai = "tagged_bam/{sample}.bam.bai"
    output:
        csv = "scQC/methCountsPerCell/{sample}_CpG_nometh.csv"
    log: "logs/CpGNOmethCountsPerCell_{sample}.log"
    threads: 1
    conda: CONDA_SHARED_ENV
    shell:
        "bamToCountTable.py --dedup -sampleTags SM -featureTags sz \
        -o {output.csv} {input.bam} > {log} 2>&1"
