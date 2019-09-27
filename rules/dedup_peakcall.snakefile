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
        blacklist = blacklist_bed
    output: "QC/readfiltering_dedup_{sample}.txt"
    threads: 10
    shell:
        "estimateReadFiltering -p {threads} --minMappingQuality 10 \
        -bl {input.blaklist} -b {input.bam} > {output}"
