rule macs2:
    input:
        bam = "dedup_bam/{sample}.bam",
        bai = "dedup_bam/{sample}.bam.bai"
    output:
        peak = "macs2_peaks/{sample}_peaks.narrowPeak",
        summits = "macs2_peaks/{sample}_summits.bed"
    params:
        outdir = "macs2_peaks",
        sample = "{sample}",
        size = genomeSize
    log: "logs/macs2_{sample}.out"
    threads: 1
    conda: CONDA_SHARED_ENV
    shell:
        "macs2 callpeak -q 0.1 -t {input.bam} -f BAM -g {params.size} --keep-dup all \
        --outdir {params.outdir} --name {params.sample} 2> {log}"

rule macs2_bed:
    input: "macs2_peaks/{sample}_peaks.narrowPeak",
    output: "macs2_peaks/{sample}_peaks.bed"
    params:
        qvalue = "1"# -log10(0.1)
    threads: 1
    conda: CONDA_SHARED_ENV
    shell:
        """
        awk 'OFS="\\t" {{ if ($9 >= {params.qvalue}) {{print $1,$2,$3,$4,$5,$6}} }}' {input} > {output}
        """

rule macs2_bed_union:
    input: expand("macs2_peaks/{sample}_peaks.bed", sample = samples)
    output: "macs2_peaks/peaks_union.bed"
    threads: 1
    conda: CONDA_SHARED_ENV
    shell:
        """
        cat {input} | sort -k1,1 -k2n,2 | \
        bedtools merge -c 4,5,6 -o collapse,mean,distinct -i - > {output}
        """

#if method == 'chic-taps':
#    rule count_inpeaks:
#        input:
#            peak = "macs2_peaks/{sample}_peaks.bed",
#            bam = "tagged_bam/{sample}.bam"
#        output: "counts/{sample}_peaks.csv"
#        params:
#            min_mq = min_mapq
#        log: "logs/count_inpeaks_{sample}.out"
#        threads: 1
#        conda: CONDA_SHARED_ENV
#        shell:
#            "bamToCountTable.py -o {output} -bedfile {input.peak} \
#            -minMQ {params.min_mq} --dedup --filterXA -sampleTags SM \
#            -joinedFeatureTags reference_name {input.bam}"

## count deduplicated reads per cell
rule countFrags_perCell:
    input:
        bam = "dedup_bam/{sample}.bam",
        bai = "dedup_bam/{sample}.bam.bai"
    output: "counts/{sample}.per_barcode.tsv"
    log: "logs/counts_{sample}.out"
    threads: 1
    conda: CONDA_SHARED_ENV
    shell:
        """
        samtools view {input.bam} | grep -o "[[:space:]]BC:Z:[ATGC]*" | \
        sed 's/[[:space:]]BC:Z://' | sort | uniq -c | \
        awk 'OFS="\\t" {{ print $2, $1 }}' > {output} 2> {log}
        """
