rule umi_dedup:
    input:
        bam = "bwa_mapped/{sample}.bam",
        idx = "bwa_mapped/{sample}.bam.bai"
    output:
        bam = "dedup_bam/{sample}.bam"
    params:
        mapq = min_mapq
    log:
        out = "logs/umi_dedup_{sample}.out",
        err = "logs/umi_dedup_{sample}.err"
    threads: 1
    conda: CONDA_SHARED_ENV
    shell:
        "umi_tools dedup --mapping-quality {params.mapq} \
        --per-cell --umi-tag=RX --cell-tag=BC --extract-umi-method=tag \
        -I {input.bam} -L {log.out} > {output.bam} 2> {log.err}"
# --paired --unmapped-reads use

rule index_dedup:
    input: "dedup_bam/{sample}.bam"
    output: "dedup_bam/{sample}.bam.bai"
    threads: 1
    conda: CONDA_SHARED_ENV
    shell: "samtools index {input}"

rule flagstat_dedup:
    input: "dedup_bam/{sample}.bam"
    output: "QC/flagstat_dedup_{sample}.txt"
    threads: 1
    conda: CONDA_SHARED_ENV
    shell: "samtools flagstat {input} > {output}"

rule readfiltering_dedup:
    input:
        bam = "dedup_bam/{sample}.bam",
        bai = "dedup_bam/{sample}.bam.bai",
        blacklist = blacklist_bed
    output: "QC/readfiltering_dedup_{sample}.txt"
    params:
        mapq = min_mapq
    log: "logs/readfiltering_dedup_{sample}.out",
    threads: 10
    conda: CONDA_SHARED_ENV
    shell:
        "estimateReadFiltering -p {threads} --minMappingQuality {params.mapq} \
        -bl {input.blacklist} -b {input.bam} > {output} 2> {log}"

rule bamCoverage_dedup:
    input:
        bam = "dedup_bam/{sample}.bam",
        bai = "dedup_bam/{sample}.bam.bai",
        blacklist = blacklist_bed
    output: "coverage/{sample}_dedup.cpm.bw"
    params:
        ignore = "chrX chrY chrM"
    log: "logs/bamCoverage_dedup_{sample}.out"
    threads: 10
    conda: CONDA_SHARED_ENV
    shell:
        "bamCoverage --normalizeUsing CPM -p {threads} \
        -ignore {params.ignore} -bl {input.blacklist} \
        -b {input.bam} -o {output} > {log} 2>&1"

rule plotEnrichment:
    input:
        bam = expand("dedup_bam/{sample}.bam", sample = samples),
        bai = expand("dedup_bam/{sample}.bam.bai", sample = samples),
        bl = blacklist_bed,
        gtf = gtf_file
    output: "QC/featureEnrichment.png"
    log: "logs/plotEnrichment.log"
    threads: 8
    conda: CONDA_SHARED_ENV
    shell:
        "plotEnrichment -p {threads} --BED {input.gtf} -bl {input.bl} \
        --smartLabels -b {input.bam} -o {output}  > {log} 2>&1"

rule plotEnrichment_biotype:
    input:
        bam = expand("dedup_bam/{sample}.bam", sample = samples),
        bai = expand("dedup_bam/{sample}.bam.bai", sample = samples),
        bl = blacklist_bed,
        gtf = gtf_file
    output: "QC/featureEnrichment_biotype.png"
    log: "logs/plotEnrichment_biotype.log"
    threads: 8
    conda: CONDA_SHARED_ENV
    shell:
        "plotEnrichment -p {threads} --BED {input.gtf} -bl {input.bl} \
        --attributeKey 'gene_biotype' --smartLabels \
        -b {input.bam} -o {output}  > {log} 2>&1"

rule bwSummary:
    input:
        bw = expand("coverage/{sample}_dedup.cpm.bw", sample = samples),
        bl = blacklist_bed
    output: "QC/bwSummary_10kBins.npz"
    log: "logs/bwSummary.log"
    threads: 8
    conda: CONDA_SHARED_ENV
    shell:
        "multiBigwigSummary bins -p {threads} -bl {input.bl} \
        --smartLabels -b {input.bw} -o {output}  > {log} 2>&1"

rule plotCorrelation:
    input: "QC/bwSummary_10kBins.npz"
    output: "QC/cor-spearman_10kBins.png"
    log: "logs/plotCorrelation.log"
    threads: 1
    conda: CONDA_SHARED_ENV
    shell:
        "plotCorrelation -p heatmap -c spearman \
        --skipZeros --removeOutliers \
        --plotNumbers -in {input} -o {output} > {log} 2>&1"
