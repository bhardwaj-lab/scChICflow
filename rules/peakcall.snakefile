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
    threads: 1
    conda: CONDA_SHARED_ENV
    shell:
        "cut -f1-6 {input} > {output}"

if method == 'chic-taps':
    rule count_inpeaks:
        input:
            peak = "macs2_peaks/{sample}_peaks.bed",
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

## split the dedup bam into cells (make temp)
rule split_cells:
    input:
        bam = "dedup_bam/{sample}.bam",
        bai = "dedup_bam/{sample}.bam.bai"
    output: "dedup_bam/{sample}/split_finished.txt"
    params:
        out = "dedup_bam/{sample}/{sample}.bam"
    log: "logs/split_cells_{sample}.out"
    threads: 1
    conda: CONDA_SHARED_ENV
    shell:
        "ln -r -s {input.bam} {params.out} && \
        bamtools split -tag BC -tagPrefix '' -in {params.out} > {log} 2>&1 && \
        touch {output}"

rule idx_cells:
    input: "dedup_bam/{sample}/split_finished.txt"
    output: "dedup_bam/{sample}/idx_finished.txt"
    params:
        dir = "dedup_bam/{sample}"
    log: "logs/idx_cells_{sample}.out"
    threads: 1
    conda: CONDA_SHARED_ENV
    shell:
        "for file in {params.dir}/*.bam; do samtools index $file; done && \
        touch {output}"

## count total fragments per cell
rule countFrags_perCell:
    input: "dedup_bam/{sample}/split_finished.txt"
    output:
        names = temp("dedup_bam/{sample}.n.txt"),
        counts = temp("dedup_bam/{sample}.c.txt")
    params:
        folder = "dedup_bam/{sample}"
    threads: 1
    conda: CONDA_SHARED_ENV
    shell:
        """
        for dat in {params.folder}/*.bam;
        do echo $(basename $dat .bam) | sed 's/.BC_/\t/g' >> {output.names}
        bamtools count -in $dat  >> {output.counts}
        done
        """

rule countFrags_cleanup:
    input:
        names = "dedup_bam/{sample}.n.txt",
        counts = "dedup_bam/{sample}.c.txt"
    output: "counts_perCell/{sample}_total.txt"
    threads: 1
    conda: CONDA_SHARED_ENV
    shell:
        "paste {input.names} {input.counts} > {output}"

## make windows (of given size) in the genome
rule makeWindows:
    input: "chrom_sizes.txt"
    output: temp("dedup_bam/counting_windows.bed")
    params:
        size = binSize
    log: "logs/makeWindows.err"
    threads: 1
    conda: CONDA_SHARED_ENV
    shell:
        "bedtools makewindows -w {params.size} -g {input} > {output} 2> {log}"

## count reads in windows per cell per sample
rule count_windows:
    input:
        windows = "dedup_bam/counting_windows.bed",
        split = "dedup_bam/{sample}/split_finished.txt",
        idx = "dedup_bam/{sample}/idx_finished.txt"
    output: 'counts_perCell/{sample}_windows.txt'
    params:
        outdir = "dedup_bam/{sample}"
    log: "logs/count_windows_{sample}.err"
    threads: 1
    conda: CONDA_SHARED_ENV
    shell:
        "bedtools multicov -D -bed {input.windows} -bams {params.outdir}/*.bam > {output} 2> {log}"

## remove splitBAM
rule cleanup_splitbam:
    input:
        bedcounts = "counts_perCell/{sample}_windows.txt",
        fragcounts = "counts_perCell/{sample}_total.txt"
    output: temp("rm_finished_{sample}.txt")
    params:
        folder = "dedup_bam/{sample}"
    threads: 1
    conda: CONDA_SHARED_ENV
    shell:
        "rm -rf {params.folder} && touch {output}"
