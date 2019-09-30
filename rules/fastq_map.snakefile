rule FASTQ1:
      input:
          indir+"/{sample}"+reads[0]+ext
      output:
          "FASTQ/{sample}"+reads[0]+".fastq.gz"
      shell:
          "( [ -f {output} ] || ln -s -r {input} {output} )"

rule FASTQ2:
      input:
          indir+"/{sample}"+reads[1]+ext
      output:
          "FASTQ/{sample}"+reads[1]+".fastq.gz"
      shell:
          "( [ -f {output} ] || ln -s -r {input} {output} )"

rule FastQC:
    input:
        "FASTQ/{sample}{read}.fastq.gz"
    output:
        "FASTQ/FastQC/{sample}{read}_fastqc.html"
    params:
        outdir = "FASTQ/FastQC"
    log:
        out = "logs/FastQC.{sample}{read}.out",
        err = "logs/FastQC.{sample}{read}.err"
    threads: 2
    #conda: CONDA_SHARED_ENV
    shell:
        "fastqc -o {params.outdir} {input} > {log.out} 2> {log.err} && \
         mkdir -p QC && \
         ln -s -r {params.outdir} QC/FastQC"

if downsample:
    rule FASTQdownsample:
        input:
            r1 = "FASTQ/{sample}"+reads[0]+".fastq.gz",
            r2 = "FASTQ/{sample}"+reads[1]+".fastq.gz"
        output:
            r1 = "FASTQ/downsample_{sample}"+reads[0]+".fastq.gz",
            r2 = "FASTQ/downsample_{sample}"+reads[1]+".fastq.gz"
        params:
            num_reads = downsample
        threads: 10
        #conda: CONDA_SHARED_ENV
        shell:
            """
            seqtk sample -s 100 {input.r1} {params.num_reads} | pigz -p {threads} -9 > {output.r1}
            seqtk sample -s 100 {input.r2} {params.num_reads} | pigz -p {threads} -9 > {output.r2}
            """

rule umi_trimming:
    input:
        r1 = "FASTQ/downsample_{sample}"+reads[0]+".fastq.gz" if downsample else "FASTQ/{sample}"+reads[0]+".fastq.gz",
        r2 = "FASTQ/downsample_{sample}"+reads[1]+".fastq.gz" if downsample else "FASTQ/{sample}"+reads[1]+".fastq.gz"
    output:
        r1 = "FASTQ/umiTrimmed_{sample}"+reads[0]+".fastq.gz",
        r2 = "FASTQ/umiTrimmed_{sample}"+reads[1]+".fastq.gz"
    params:
        barcodes = barcode_list
    log:
        out = "logs/umi_trimming_{sample}.out",
        err = "logs/umi_trimming_{sample}.err"
    threads: 2
    #conda: CONDA_SHARED_ENV
    shell:
        "umi_tools extract --extract-method=string --bc-pattern=NNNCCCCCCCC \
        --whitelist {params.barcodes} --filter-cell-barcode -I {input.r1} --read2-in {input.r2} \
        -L {log.out} --read2-out {output.r2} | gzip - > {output.r1}"


if trim:
    rule cutadapt:
        input:
            r1 = "FASTQ/umiTrimmed_{sample}"+reads[0]+".fastq.gz",
            r2 = "FASTQ/umiTrimmed_{sample}"+reads[1]+".fastq.gz"
        output:
            r1 = "FASTQ_trimmed/{sample}"+reads[0]+".fastq.gz",
            r2 = "FASTQ_trimmed/{sample}"+reads[1]+".fastq.gz"
        params:
            opts = str(trimmerOptions or '')
        log:
            out = "logs/cutadapt.{sample}.out",
            err = "logs/cutadapt.{sample}.err"
        threads: 8
        shell:
            "cutadapt -j {threads} -e 0.1 -q 16 -O 3 --trim-n --minimum-length 20 \
            -b TGGAATTCTCGGGTGCCAAGG -B TGGAATTCTCGGGTGCCAAGG \
            -b ATCTCGTATGCCGTCTTCTGCTTG -B ATCTCGTATGCCGTCTTCTGCTTG \
            -b GTTCAGAGTTCTACAGTCCGACGATC -B GTTCAGAGTTCTACAGTCCGACGATC \
            --nextseq-trim=16 {params.opts} \
            -o {output.r1} -p {output.r2} {input.r1} {input.r2} > {log.out} 2> {log.err}"

    rule FastQC_trimmed:
        input:
            "FASTQ_trimmed/{sample}{read}.fastq.gz"
        output:
            "FASTQ_trimmed/FastQC/{sample}{read}_fastqc.html"
        params:
            outdir = "FASTQ_trimmed/FastQC"
        log:
            out = "logs/FastQC_trimmed.{sample}{read}.out",
            err = "logs/FastQC_trimmed.{sample}{read}.err"
        threads: 2
        #conda: CONDA_SHARED_ENV
        shell:
            "fastqc -o {params.outdir} {input} > {log.out} 2> {log.err} && \
             ln -s -r FASTQ_trimmed/FastQC QC/FastQC_trimmed "



## MAP

rule bwa_map:
    input:
        r1 = "FASTQ_trimmed/{sample}"+reads[0]+".fastq.gz" if trim else "FASTQ/umiTrimmed_{sample}"+reads[0]+".fastq.gz",
        r2 = "FASTQ_trimmed/{sample}"+reads[1]+".fastq.gz" if trim else "FASTQ/umiTrimmed_{sample}"+reads[1]+".fastq.gz",
        idx = bwa_index
    output: "bwa_mapped/{sample}.bam"
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

rule bwa_index:
    input: "bwa_mapped/{sample}.bam"
    output: "bwa_mapped/{sample}.bam.bai"
    shell: "samtools index {input}"

rule flagstat_bwa:
    input: "bwa_mapped/{sample}.bam"
    output: "QC/flagstat_bwa_{sample}.txt"
    threads: 1
    shell: "samtools flagstat {input} > {output}"

rule readfiltering_bwa:
    input:
        bam = "bwa_mapped/{sample}.bam",
        bai = "bwa_mapped/{sample}.bam.bai",
        blacklist = blacklist_bed
    output: "QC/readfiltering_bwa_{sample}.txt"
    params:
        mapq = min_mapq
    threads: 10
    shell:
        "estimateReadFiltering -p {threads} --minMappingQuality {params.mapq} --samFlagInclude 64 \
        -bl {input.blacklist} -b {input.bam} > {output}"
