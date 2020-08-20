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
            expand("QC/readfiltering_dedup_{sample}.txt", sample = samples)
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
    conda: CONDA_SHARED_ENV
    shell:
        "multiqc -f -o {params.outdir} {params.outdir} > {log.out} 2> {log.err}"

## count deduplicated reads per cell
rule countFrags_perCell:
    input:
        bam = "dedup_bam/{sample}.bam",
        bai = "dedup_bam/{sample}.bam.bai"
    output: "counts/{sample}.per_barcode.tsv"
    log: "logs/counts_{sample}.out"
    threads: 1
    conda: CONDA_SHARED_ENV
    shell:
        """
        samtools view {input.bam} | grep -o "[[:space:]]BC:Z:[ATGC]*" | \
        sed 's/[[:space:]]BC:Z://' | sort | uniq -c | \
        awk 'OFS="\\t" {{ print $2, $1 }}' > {output} 2> {log}
        """

rule plate_plot:
    input:
        counts = expand("counts/{sample}.per_barcode.tsv", sample = samples),
        barcodes = barcode_list
    output: "QC/plate_plots.pdf"
    params:
        countdir = "counts",
        rscript = os.path.join(workflow.basedir, "tools", "make_plate_plots.R")
    log: "logs/plate_plots.out"
    threads: 1
    conda: CONDA_SHARED_ENV
    shell:
        """
        Rscript {params.rscript} {input.barcodes} {params.countdir} {output}
        """
