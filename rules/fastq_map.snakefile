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
        "FastQC/{sample}{read}_fastqc.html"
    log:
        out = "logs/FastQC.{sample}{read}.out",
        err = "logs/FastQC.{sample}{read}.err"
    threads: 2
    #conda: CONDA_SHARED_ENV
    shell: "fastqc -o FastQC {input} > {log.out} 2> {log.err}"

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
        barcodes = barcode_tsv
    log:
        out = "logs/umi_trimming_{sample}.out",
        err = "logs/umi_trimming_{sample}.err"
    threads: 2
    #conda: CONDA_SHARED_ENV
    shell: "umi_tools extract --extract-method=string --bc-pattern=NNNCCCCCCCC --whitelist {params.barcodes} --filter-cell-barcode -I {input.r1} --read2-in {input.r2} -L {log.out} --read2-out {output.r2} | gzip - > {output.r1}"


if trim:
    rule fastp:
        input:
            r1 = "FASTQ/umiTrimmed_{sample}"+reads[0]+".fastq.gz",
            r2 = "FASTQ/umiTrimmed_{sample}"+reads[1]+".fastq.gz"
        output:
            r1 = "FASTQ_trimmed/{sample}"+reads[0]+".fastq.gz",
            r2 = "FASTQ_trimmed/{sample}"+reads[1]+".fastq.gz",
            json = "FASTQ_trimmed/{sample}.fastp.json",
            html = "FASTQ_trimmed/{sample}.fastp.html"
        params:
            opts = str(trimmerOptions or '')
        log:
            out = "logs/fastp.{sample}.out",
            err = "logs/fastp.{sample}.err"
        threads: 8
        shell:
            """
            fastp -q 16 -l 20 -w {threads} -i {input.r1} -I {input.r2} -o {output.r1} \
            -O {output.r2} -j {output.json} -h {output.html} {params.opts} > {log.out} 2> {log.err}
            """

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
        shell: "fastqc -o {params.outdir} {input} > {log.out} 2> {log.err}"



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