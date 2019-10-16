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
        cluster = '--cluster' if cluster else '',
        context = '-context '+bedContext if bedContext else ''
    log:
        out = "logs/taps_tagger_{sample}.out",
        err = "logs/taps_tagger_{sample}.err"
    threads: 1
    conda: CONDA_SHARED_ENV
    shell:
        "tapsTagger.py {params.cluster} {params.context} \
        -ref {input.genome} -method {params.method} -bed {output.bed} \
        -o {output.bam} -min_mq {params.min_mq} {input.bam} \
        > {log.out} 2> {log.err}"

### Cmd to make bigwigs ###
methBigwig_cmd = """
    grep {params.context} {input.bed} | sort -k1,1 -k2,2n |\
    bedtools merge -s -c 4,5,6 -o distinct,sum,distinct -i - |\
    cut -f1-3,5 {params.grep} > {params.sample}.tmp.bg &&
    bedGraphToBigWig {params.sample}.tmp.bg {input.sizes} {output.bw} 2> {log} &&
    rm {params.sample}.tmp.bg
    """

## methylated C
rule meth_bigwig:
    input:
        bed = "meth_calls/{sample}_methylation.bed",
        bam = "tagged_bam/{sample}.bam",
        bai = "tagged_bam/{sample}.bam.bai",
        sizes = "chrom_sizes.txt"
    output:
        bw = "meth_calls/{sample}.methCpG.bw"
    params:
        sample = '{sample}',
        genome = genome_fasta,
        context = "'Z'",
        grep = "| grep -v 'phage\|G\|K\|J'"
    log: "logs/meth_bigwig_{sample}.err"
    threads: 1
    conda: CONDA_SHARED_ENV
    shell: methBigwig_cmd

## all C
rule all_bigwig:
    input:
        bed = "meth_calls/{sample}_methylation.bed",
        bam = "tagged_bam/{sample}.bam",
        bai = "tagged_bam/{sample}.bam.bai",
        sizes = "chrom_sizes.txt"
    output:
        bw = temp("meth_calls/{sample}.allCpG.bw")
    params:
        sample = '{sample}',
        genome = genome_fasta,
        context = "'z\|Z'",
        grep = "| grep -v 'phage\|G\|K\|J'"
    log: "logs/all_bigwig_{sample}.err"
    threads: 1
    conda: CONDA_SHARED_ENV
    shell: methBigwig_cmd

## ratio meth/all
rule bigwigRatio:
    input:
        methbw = "meth_calls/{sample}.methCpG.bw",
        allbw = "meth_calls/{sample}.allCpG.bw"
    output: "meth_calls/{sample}.methRatio.bw"
    params:
        blklist = blacklist_bed
    log: "logs/bigwigRatio_{sample}.err"
    threads: 10
    conda: CONDA_SHARED_ENV
    shell:
        "bigwigCompare -p {threads} -bl {params.blklist} -bs 1 --skipZeroOverZero \
        --operation ratio -o {output} --skipNAs -b1 {input.methbw} -b2 {input.allbw}"

rule meth_gzip:
    input:
        bed = "meth_calls/{sample}_methylation.bed",
        bw = "meth_calls/{sample}.methRatio.bw"
    output: "meth_calls/{sample}_methylation.bed.gz"
    threads: 1
    shell: 'gzip {input}'


## methylation counts per bin per cell
rule meth_bincounts:
    input:
        bam = "tagged_bam/{sample}.bam",
        bai = "tagged_bam/{sample}.bam.bai"
    output:
        csv = "meth_counts/{sample}_CpG_binCounts.csv"
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

rule bwSummary_meth:
    input:
        bw = expand("meth_calls/{sample}.methRatio.bw", sample = samples),
        bl = blacklist_bed
    output: "QC/bwSummary_methRatio_10kBins.npz"
    log: "logs/bwSummary_methCpG.log"
    threads: 8
    conda: CONDA_SHARED_ENV
    shell:
        "multiBigwigSummary bins -p {threads} -bl {input.bl} \
        --smartLabels -b {input.bw} -o {output}  > {log} 2>&1"

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
        --plotsOnly -o {params.outdir} {input} > {log} 2>&1"

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
