def get_multiqc_input():
    file = [
        expand("QC/FastQC/{sample}{read}_fastqc.html", sample = samples, read=reads),
        expand("QC/flagstat_bwa_{sample}.txt", sample = samples),
        expand("QC/readfiltering_bwa_{sample}.txt", sample = samples)
        ]
    if trim:
        file.append(expand("QC/FastQC_trimmed/{sample}{read}_fastqc.html", sample = samples, read = reads))
    if method in ['chic', 'chic-taps']:
        file.append([
            expand("QC/flagstat_dedup_{sample}.txt", sample = samples),
            expand("QC/readfiltering_dedup_{sample}.txt", sample = samples),
            expand("homer_peaks/{sample}_peaks.bed", sample = samples)
            ])
    return(file)

rule multiQC:
    input: get_multiqc_input()
    output: "QC/multiqc_report.html"
    params:
        outdir = "QC"
    log:
        out = "logs/multiqc.out",
        err = "logs/multiqc.err"
    threads: 1
    shell:
        "multiqc -o {params.outdir} {params.outdir} > {log.out} 2> {log.err}"
