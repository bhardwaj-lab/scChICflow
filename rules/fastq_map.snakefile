
rule FastQC:
    input:
        "FASTQ/{sample}{read}.fastq.gz"
    output:
        "QC/FastQC/{sample}{read}_fastqc.html"
    params:
        outdir = "QC/FastQC"
    log: "logs/FastQC.{sample}{read}.out",
    threads: 2
    conda: CONDA_SHARED_ENV
    shell:
        "fastqc -o {params.outdir} {input} > {log} 2>&1"

if downsample:
    rule FASTQdownsample:
        input:
            r1 = "FASTQ/{sample}"+reads[0]+".fastq.gz",
            r2 = "FASTQ/{sample}"+reads[1]+".fastq.gz"
        output:
            r1 = temp("FASTQ/downsampled/downsample_{sample}"+reads[0]+".fastq.gz"),
            r2 = temp("FASTQ/downsampled/downsample_{sample}"+reads[1]+".fastq.gz")
        params:
            num_reads = downsample
        threads: 10
        conda: CONDA_SHARED_ENV
        shell:
            """
            seqtk sample -s 100 {input.r1} {params.num_reads} | pigz -p {threads} -9 > {output.r1}
            seqtk sample -s 100 {input.r2} {params.num_reads} | pigz -p {threads} -9 > {output.r2}
            """

rule umi_trimming:
    input:
        r1 = lambda wildcards: "FASTQ/downsampled/downsample_{sample}"+reads[0]+".fastq.gz" if downsample else "FASTQ/{sample}"+reads[0]+".fastq.gz",
        r2 = lambda wildcards: "FASTQ/downsampled/downsample_{sample}"+reads[1]+".fastq.gz" if downsample else "FASTQ/{sample}"+reads[1]+".fastq.gz"
    output:
        r1 = temp("FASTQ/umi_trimmed/umiTrimmed_{sample}"+reads[0]+".fastq.gz"),
        r2 = temp("FASTQ/umi_trimmed/umiTrimmed_{sample}"+reads[1]+".fastq.gz")
    params:
        barcodes = barcode_list
    log:
        out = "logs/umi_trimming_{sample}.out",
        err = "logs/umi_trimming_{sample}.err"
    threads: 2
    conda: CONDA_SHARED_ENV
    shell:
        "umi_tools extract --extract-method=string --bc-pattern=NNNCCCCCCCC \
        --whitelist {params.barcodes} --filter-cell-barcode -I {input.r1} --read2-in {input.r2} \
        -L {log.out} --read2-out {output.r2} | gzip - > {output.r1}"


if trim:
    rule cutadapt:
        input:
            r1 = "FASTQ/umi_trimmed/umiTrimmed_{sample}"+reads[0]+".fastq.gz",
            r2 = "FASTQ/umi_trimmed/umiTrimmed_{sample}"+reads[1]+".fastq.gz"
        output:
            r1 = "FASTQ_trimmed/{sample}"+reads[0]+".fastq.gz",
            r2 = "FASTQ_trimmed/{sample}"+reads[1]+".fastq.gz"
        params:
            opts = str(trimmerOptions or '')
        log:
            out = "logs/cutadapt.{sample}.out",
            err = "logs/cutadapt.{sample}.err"
        threads: 8
        conda: CONDA_SHARED_ENV
        shell:
            "cutadapt -j {threads} -e 0.1 -O 5 \
            -b TGGAATTCTCGGGTGCCAAGG -B TGGAATTCTCGGGTGCCAAGG \
            -b ATCTCGTATGCCGTCTTCTGCTTG -B ATCTCGTATGCCGTCTTCTGCTTG \
            -b GTTCAGAGTTCTACAGTCCGACGATC -B GTTCAGAGTTCTACAGTCCGACGATC \
            --nextseq-trim=30 {params.opts} \
            -o {output.r1} -p {output.r2} {input.r1} {input.r2} > {log.out} 2> {log.err}"

    rule FastQC_trimmed:
        input:
            "FASTQ_trimmed/{sample}{read}.fastq.gz"
        output:
            "QC/FastQC_trimmed/{sample}{read}_fastqc.html"
        params:
            outdir = "QC/FastQC_trimmed"
        log: "logs/FastQC_trimmed.{sample}{read}.out",
        conda: CONDA_SHARED_ENV
        threads: 2
        #conda: CONDA_SHARED_ENV
        shell:
            "fastqc -o {params.outdir} {input} > {log} 2>&1"

## prep chromosome sizes for bigwigs and stuff
rule chrSizes:
    input: genome_fasta
    output: temp("chrom_sizes.txt")
    params:
        genome = genome_fasta+".fai"
    threads: 1
    conda: CONDA_SHARED_ENV
    shell:
        "cut -f1-2 {params.genome} > {output}"

## MAP Fastq
#"""
#bwa mem -v 1 -T {params.mapq} -t {threads} {input.idx} {input.r1} {input.r2} 2> {log.err} |\
#samtools view -h {params.samfilter} | awk -v sample={params.sample} \
#'OFS="\\t" {{ if($0 ~ "^@") {{print $0}} else \
#{{ split($1,a,"_"); $1=""; print a[1]";BC:Z:"a[2]";RX:Z:"a[3], $0, "SM:Z:"sample"_"a[2], "BC:Z:"a[2], "RX:Z:"a[3], "MI:Z:"a[2]a[3] }} }}' |\
#samtools sort -@ {threads} -T {params.sample} -o {output} > {log.out} 2>> {log.err}
#"""

filter_cmd = """
    samtools view -h {params.samfilter} | awk -v sample={params.sample} \
    'OFS="\\t" {{ if($0 ~ "^@") {{print $0}} else \
    {{ split($1,a,"_"); $1=""; print a[1]";BC:Z:"a[2]";RX:Z:"a[3], $0, "SM:Z:"sample"_"a[2], "BC:Z:"a[2], "RX:Z:"a[3], "MI:Z:"a[2]a[3] }} }}' |\
    samtools sort -@ {threads} -T {params.sample} -o {output} > {log.out} 2>> {log.err}
    """

if protocol == "chic":
    mapping_cmd = "bwa mem -v 1 -T {params.mapq} -t {threads} {input.idx} {input.r1} {input.r2} 2> {log.err} | " + filter_cmd
else:
    mapping_cmd = "hisat2 --sensitive --no-spliced-alignment --no-mixed --no-discordant --no-softclip -X 1000 -x {params.idx} -1 {input.r1} -2 {input.r2} 2> {log.err} | " + filter_cmd

rule bam_map:
    input:
        r1 = lambda wildcards: "FASTQ_trimmed/{sample}"+reads[0]+".fastq.gz" if trim else "FASTQ/umi_trimmed/umiTrimmed_{sample}"+reads[0]+".fastq.gz",
        r2 = lambda wildcards: "FASTQ_trimmed/{sample}"+reads[1]+".fastq.gz" if trim else "FASTQ/umi_trimmed/umiTrimmed_{sample}"+reads[1]+".fastq.gz",
        idx = bwa_index if protocol == "chic" else hisat2_index+".1.ht2"
    output: "mapped_bam/{sample}.bam"
    params:
        sample = '{sample}',
        mapq = min_mapq,
        samfilter='-F 4 -F 256',
        idx = "" if protocol == "chic" else hisat2_index
    log:
        out = "logs/bam_map_{sample}.out",
        err = "logs/bam_map_{sample}.err"
    threads: 10
    conda: CONDA_SHARED_ENV
    shell: mapping_cmd


rule bam_index:
    input: "mapped_bam/{sample}.bam"
    output: "mapped_bam/{sample}.bam.bai"
    conda: CONDA_SHARED_ENV
    shell: "samtools index {input}"

rule flagstat_bam:
    input: "mapped_bam/{sample}.bam"
    output: "QC/flagstat_bwa_{sample}.txt"
    threads: 1
    conda: CONDA_SHARED_ENV
    shell: "samtools flagstat {input} > {output}"

## get some stats
rule readfiltering_bam:
    input:
        bam = "mapped_bam/{sample}.bam",
        bai = "mapped_bam/{sample}.bam.bai"
    output: "QC/readfiltering_bwa_{sample}.txt"
    params:
        mapq = min_mapq,
        blacklist = "-bl " + blacklist_bed if blacklist_bed else ""
    threads: 10
    conda: CONDA_SHARED_ENV
    shell:
        "estimateReadFiltering -p {threads} --minMappingQuality {params.mapq} --samFlagInclude 64 \
        {params.blacklist} -b {input.bam} > {output}"
