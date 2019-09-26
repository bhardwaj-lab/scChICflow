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
    "umi_tools dedup --mapping-quality 10 --paired --unmapped-reads use \
    --per-cell --umi-tag=RX --cell-tag=BC --extract-umi-method=tag \
    -I {input.bam} -L {log.out} > {output.bam} 2> {log.err}"

rule index_dedup:
    input: "dedup_bam/{sample}.bam"
    output: "dedup_bam/{sample}.bam.bai"
    threads: 1
    shell: "samtools index {input}"



"""
cutadapt -j {threads} -e 0.1 -q 16 -O 3 --trim-n --minimum-length 20 \
-a AGATCGGAAGAGC -A AGATCGGAAGAGC -u 12 -u -4 -U -4 --nextseq-trim=16 \
--match-read-wildcards {params.opts} \
-o "{output.r1}" -p "{output.r2}" "{input.r1}" "{input.r2}" > {log.out} 2> {log.err}
"""
