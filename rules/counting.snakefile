
if countRegions == "windows":
    ## make windows (of given size) in the genome
    rule make_regions:
        input: "chrom_sizes.txt"
        output: temp("counts/counting_regions.saf")
        params:
            size = binSize
        log: "logs/makeWindows.err"
        threads: 1
        conda: CONDA_SHARED_ENV
        shell:
            """
            bedtools makewindows -w {params.size} -g {input} |\
            awk 'OFS="\\t" {{ print $1"_"$2"_"$3, $1, $2, $3, "." }}' > {output} 2> {log}
            """

elif countRegions == "bed" or countRegions == "peaks":
    rule make_regions:
        input: lambda wildcards: bedFile if countRegions == "bed" else "macs2_peaks/{sample}_peaks.bed"
        output: temp("counts/counting_regions.saf")
        threads: 1
        conda: CONDA_SHARED_ENV
        shell:
            """
            awk 'OFS="\\t" {{ print $4, $1, $2, $3, $6 }}' {input} > {output}
            """

elif countRegions == "genes":
    print("Counting reads in genes per cell")

## count reads in Genes/regions/windows using featurecounts (bulk)
rule count_regions:
    input:
        regions = lambda wildcards: gtf_file if countRegions == "genes" else "counts/counting_regions.saf",
        bam = "dedup_bam/{sample}.bam",
        bai = "dedup_bam/{sample}.bam.bai"
    output:
        counts = "counts/{sample}."+countRegions+"_total.tsv",
        bam = temp("counts/{sample}.bam.featureCounts.bam")
    params:
        filetype = lambda wildcards: "-t gene" if countRegions == "genes" else "-F SAF"
    log: "logs/featurecounts_{sample}_bulk.err"
    threads: 10
    conda: CONDA_SHARED_ENV
    shell:
        "featureCounts -T {threads} -a {input.regions} {params.filetype} \
        -R BAM -s 0 -o {output.counts} {input.bam} > {log} 2>&1"

rule fcount_sort:
    input: "counts/{sample}.bam.featureCounts.bam"
    output: temp("counts/{sample}.featureCounts.bam")
    threads: 10
    conda: CONDA_SHARED_ENV
    shell:
        "samtools sort -@ {threads} -O BAM -o {output} {input}"

rule fcount_index:
    input: "counts/{sample}.featureCounts.bam"
    output: temp("counts/{sample}.featureCounts.bam.bai")
    threads: 1
    conda: CONDA_SHARED_ENV
    shell:
        "samtools index {input}"

rule count_regions_cells:
    input:
        bam = "counts/{sample}.featureCounts.bam",
        bai = "counts/{sample}.featureCounts.bam.bai"
    output: "counts/{sample}."+countRegions+"_per_barcode.tsv"
    log: "logs/umi_counts_{sample}."+countRegions+".err"
    threads: 1
    conda: CONDA_SHARED_ENV
    shell:
        "umi_tools count --per-gene --gene-tag=XT --assigned-status-tag=XS \
        --per-cell --cell-tag=BC --umi-tag=RX --extract-umi-method=tag \
        --method=percentile -I {input.bam} -S {output} -v 4 --log2stderr --log={log}"
