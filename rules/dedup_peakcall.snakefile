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
    log: "logs/bamCoverage_dedup_{sample}.out"
    threads: 10
    conda: CONDA_SHARED_ENV
    shell:
        "bamCoverage --normalizeUsing CPM -p {threads} -bl {input.blacklist} \
        -b {input.bam} -o {output} > {log} 2>&1"

rule maketags_dedup:
    input:
        bam = "dedup_bam/{sample}.bam",
        bai = "dedup_bam/{sample}.bam.bai"
    output: "homer_peaks/{sample}/tagInfo.txt"
    params:
        dir = "homer_peaks/{sample}"
    log: "logs/maketags_dedup_{sample}.out"
    threads: 1
    conda: CONDA_SHARED_ENV
    shell:
        "makeTagDirectory {params.dir} {input.bam} > {log} 2>&1"

rule homer_findpeaks:
    input: "homer_peaks/{sample}/tagInfo.txt"
    output:
        txt = "homer_peaks/{sample}/peaks.txt",
        bed = "homer_peaks/{sample}_peaks.bed"
    params:
        dir = "homer_peaks/{sample}",
        sample = "{sample}"
    log: "logs/homer_findpeaks_{sample}.out"
    threads: 1
    conda: CONDA_SHARED_ENV
    shell:
        """
        findPeaks {params.dir} -style factor -center -fdr 0.05 -tagThreshold 2 \
        -L 2 -F 0 -C 0 > {output.txt} 2> {log} && \
        awk 'OFS="\\t" {{if ($0 !~ "#") {{print $2, $3, $4, $1, $10, $5 }} }}' \
        {output.txt} > {output.bed}
        """

rule homer_cleanup:
    input: "homer_peaks/{sample}_peaks.bed"
    output: directory("QC/homer_peaks_{sample}")
    params:
        dir = "homer_peaks/{sample}"
#        sample = "{sample}"
    threads: 1
    conda: CONDA_SHARED_ENV
    shell:
        """
        rm {params.dir}/*.tags.tsv && \
        ln -sf -r {params.dir} {output}
        """

if method == 'chic-taps':
    rule count_inpeaks:
        input:
            peak = "homer_peaks/{sample}_peaks.bed",
            bam = "tagged_bam/{sample}.bam"
        output: "counts/{sample}_peaks.csv"
        params:
            min_mq = min_mapq
        log: "logs/count_inpeaks_{sample}.out"
        threads: 1
        conda: CONDA_SHARED_ENV
        shell:
            "bamToCountTable.py -o {output} -bedfile {input.peak} \
            -minMQ {params.min_mq} --dedup --filterXA -sampleTags SM \
            -joinedFeatureTags reference_name {input.bam}"
