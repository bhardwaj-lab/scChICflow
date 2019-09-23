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
            out = "logs/Cutadapt.{sample}.out",
            err = "logs/Cutadapt.{sample}.err"
        threads: 8
        #conda: CONDA_SHARED_ENV
        shell:
            """
            cutadapt -j {threads} -e 0.1 -q 16 -O 3 --trim-n --minimum-length 20 \
            -a AGATCGGAAGAGC -A AGATCGGAAGAGC -u 12 -u -4 -U -4 --nextseq-trim=16 \
            --match-read-wildcards {params.opts} \
            -o "{output.r1}" -p "{output.r2}" "{input.r1}" "{input.r2}" > {log.out} 2> {log.err}
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
