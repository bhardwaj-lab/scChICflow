## explanation for the parameters
# umis are de-dupped per cell using the CB and RX tags, and only R1 start pos
# due to 3-bp UMIs, I use no edit distance threshold (method = unique)
# I consider reads with 3'- soft-clipping of 2 or more bases as distinct from non-clipped reads

rule umi_dedup:
    input:
        bam = "mapped_bam/{sample}.bam",
        idx = "mapped_bam/{sample}.bam.bai"
    output:
        bam = "dedup_bam/{sample}.bam",
        stats = "QC/umi_dedup/{sample}_per_umi.tsv"
    params:
        mapq = min_mapq,
        sample = "{sample}",
        tmp = tempDir,
#        paired = "--paired --unmapped-reads use" if protocol == "tchic" else ""
    log:
        out = "logs/umi_dedup_{sample}.out",
        err = "logs/umi_dedup_{sample}.err"
    threads: 1
    conda: CONDA_SHARED_ENV
    shell:
        "umi_tools dedup --mapping-quality {params.mapq} \
        --per-cell --umi-tag=RX --cell-tag=BC --extract-umi-method=tag \
        --method unique --spliced-is-unique --soft-clip-threshold 2 \
        --output-stats=QC/umi_dedup/{params.sample} \
        --temp-dir={params.tmp} \
        -I {input.bam} -L {log.out} > {output.bam} 2> {log.err}"


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
        bai = "dedup_bam/{sample}.bam.bai"
    output: "QC/readfiltering_dedup_{sample}.txt"
    params:
        mapq = min_mapq,
        blacklist = "-bl " + blacklist_bed if blacklist_bed else ""
    log: "logs/readfiltering_dedup_{sample}.out",
    threads: 10
    conda: CONDA_SHARED_ENV
    shell:
        "estimateReadFiltering -p {threads} --minMappingQuality {params.mapq} \
        {params.blacklist} -b {input.bam} > {output} 2> {log}"

rule bamCoverage_dedup:
    input:
        bam = "dedup_bam/{sample}.bam",
        bai = "dedup_bam/{sample}.bam.bai"
    output: "coverage/{sample}_dedup.cpm.bw"
    params:
        ignore = "chrX chrY chrM",
        extendReads = "-e --maxFragmentLength 1000" if protocol=="tchic" else "",
        blacklist = "-bl " + blacklist_bed if blacklist_bed else ""
    log: "logs/bamCoverage_dedup_{sample}.out"
    threads: 10
    conda: CONDA_SHARED_ENV
    shell:
        "bamCoverage --normalizeUsing CPM -bs 100 -p {threads} \
        -ignore {params.ignore} {params.blacklist} \
        -b {input.bam} -o {output} > {log} 2>&1"

rule plotEnrichment:
    input:
        bam = expand("dedup_bam/{sample}.bam", sample = samples),
        bai = expand("dedup_bam/{sample}.bam.bai", sample = samples),
        gtf = gtf_file
    output: "QC/featureEnrichment.png"
    log: "logs/plotEnrichment.log"
    params:
        blacklist = "-bl " + blacklist_bed if blacklist_bed else ""
    threads: 8
    conda: CONDA_SHARED_ENV
    shell:
        "plotEnrichment -p {threads} --BED {input.gtf} {params.blacklist} \
        --smartLabels --variableScales --perSample -b {input.bam} -o {output}  > {log} 2>&1"

rule plotEnrichment_biotype:
    input:
        bam = expand("dedup_bam/{sample}.bam", sample = samples),
        bai = expand("dedup_bam/{sample}.bam.bai", sample = samples),
        gtf = gtf_file
    output: "QC/featureEnrichment_biotype.png"
    log: "logs/plotEnrichment_biotype.log"
    params:
        blacklist = "-bl " + blacklist_bed if blacklist_bed else ""
    threads: 8
    conda: CONDA_SHARED_ENV
    shell:
        "plotEnrichment -p {threads} --BED {input.gtf} {params.blacklist} \
        --attributeKey 'gene_biotype' --smartLabels --variableScales --perSample \
        -b {input.bam} -o {output}  > {log} 2>&1"

rule bwSummary:
    input:
        bw = expand("coverage/{sample}_dedup.cpm.bw", sample = samples)
    output: "QC/bwSummary_10kBins.npz"
    log: "logs/bwSummary.log"
    params:
        blacklist = "-bl " + blacklist_bed if blacklist_bed else ""
    threads: 8
    conda: CONDA_SHARED_ENV
    shell:
        "multiBigwigSummary bins -p {threads} {params.blacklist} \
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
