rule bwa_map:
    input:
        r1 = "FASTQ_trimmed/{sample}"+reads[0]+".fastq.gz" if trim else "FASTQ/umiTrimmed_{sample}"+reads[0]+".fastq.gz",
        r2 = "FASTQ_trimmed/{sample}"+reads[1]+".fastq.gz" if trim else "FASTQ/umiTrimmed_{sample}"+reads[1]+".fastq.gz",
        idx = bwa_index
    output: "BWA_mapped/{sample}.bam"
    params:
        sample = '{sample}'
    log:
        out = "logs/bwa_map_{sample}.out",
        err = "logs/bwa_map_{sample}.err"
    threads: 10
    #conda: CONDA_scRIA_ENV
    shell:
        """bwa mem -t {threads} {input.idx} {input.r1} {input.r2} |\
        samtools view -F 4 -h | awk -v sample={params.sample} \
        'OFS="\\t" {{ if($0 ~ "^@") {{print $0}} else \
        {{ split($1,a,"_"); print a[1]";BC:Z:"a[2]";RX:Z:"a[3], $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12, $13, $14, "SM:Z:"sample"_"a[2], "BC:Z:"a[2], "RX:Z:"a[3], "MI:Z:"a[2]a[3] }} }}' |\
        samtools sort -@ {threads} -T {params.sample} -o {output} > {log.out} 2> {log.err}"""
