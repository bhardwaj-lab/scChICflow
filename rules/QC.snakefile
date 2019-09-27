def get_fastqc_input():
    file = [expand("FASTQ/FastQC/{sample}{read}_fastqc.html", sample = samples, read=reads)]
    if trim:
        file.append(expand("FASTQ_trimmed/FastQC/{sample}{read}_fastqc.html", sample = samples, read = reads))
    return(file)

def get_samflags_input():
    file = [expand("QC/flagstat_bwa_{sample}.txt", sample = samples)]
    if method == 'chic':
        file.append(expand("QC/flagstat_dedup_{sample}.txt", sample = samples))
    return(file)

def get_readfilter_input():
    file = [expand("QC/readfiltering_bwa_{sample}.txt", sample = samples)]
    if method == 'chic':
        file.append(expand("QC/readfiltering_dedup_{sample}.txt", sample = samples))
    return(file)

rule multiQC:
    input:
        fastqc = get_fastqc_input(),
        samflags = get_samflags_input(),
        readfilter = get_readfilter_input()
    output: "QC/multiqc_report.html"
    params:
        outdir = "QC"
    log:
        out = "logs/multiqc.out",
        err = "logs/multiqc.err"
    threads: 1
    shell:
        "multiqc -o {params.outdir} > {log.out} 2> {log.err}"
