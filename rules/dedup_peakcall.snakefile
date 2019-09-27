rule umi_dedup:
    input:
        bam = "bwa_mapped/{sample}.bam",
        idx = "bwa_mapped/{sample}.bam.bai"
    output:
        bam = "dedup_bam/{sample}.bam"
    log:
        out = "logs/umi_dedup_{sample}.out",
        err = "logs/umi_dedup_{sample}.err"
    threads: 1
    shell:
        "umi_tools dedup --mapping-quality 10 \
        --per-cell --umi-tag=RX --cell-tag=BC --extract-umi-method=tag \
        -I {input.bam} -L {log.out} > {output.bam} 2> {log.err}"
# --paired --unmapped-reads use

rule index_dedup:
    input: "dedup_bam/{sample}.bam"
    output: "dedup_bam/{sample}.bam.bai"
    threads: 1
    shell: "samtools index {input}"

rule flagstat_dedup:
    input: "dedup_bam/{sample}.bam"
    output: "QC/flagstat_dedup_{sample}.txt"
    threads: 1
    shell: "samtools flagstat {input} > {output}"

rule readfiltering_dedup:
    input:
        bam = "dedup_bam/{sample}.bam",
        bai = "dedup_bam/{sample}.bam.bai",
        blacklist = blacklist_bed
    output: "QC/readfiltering_dedup_{sample}.txt"
    log: "logs/readfiltering_dedup_{sample}.out",
    threads: 10
    shell:
        "estimateReadFiltering -p {threads} --minMappingQuality 10 \
        -bl {input.blacklist} -b {input.bam} > {output} 2> {log}"

rule bamCoverage_dedup:
    input:
        bam = "dedup_bam/{sample}.bam",
        bai = "dedup_bam/{sample}.bam.bai",
        blacklist = blacklist_bed
    output: "coverage/{sample}_dedup.cpm.bw"
    log: "logs/bamCoverage_dedup_{sample}.out"
    threads: 10
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
    shell:
        "makeTagDirectory {params.dir} {input.bam} > {log} 2>&1"

rule homer_findpeaks:
    input: "homer_peaks/{sample}/tagInfo.txt"
    output:
        txt = temp("homer_peaks/{sample}/peaks.txt"),
        bed = "homer_peaks/{sample}_peaks.bed"
    params:
        dir = "homer_peaks/{sample}",
        sample = "{sample}"
    log: "logs/homer_findpeaks_{sample}.out"
    threads: 1
    shell:
        """
        findPeaks homer_test -style factor -center -fdr 0.05 -tagThreshold 2 \
        -L 2 -F 0 -C 0 > {output.txt} > {log} 2>&1 && \
        awk 'OFS="\\t" {{if ($0 !~ "#") {{print $2, $3, $4, $1, $10, $5 }} }}' \
        {output.txt} > {output.bed} && \
        ln -s $PWD/{params.dir} QC/homer_peaks_{params.sample}
        """
