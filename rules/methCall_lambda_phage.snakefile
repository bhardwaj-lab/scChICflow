rule taps_tagger_phage:
    input:
        bam = "bwa_mapped/{sample}.bam",
        idx = "bwa_mapped/{sample}.bam.bai",
        genome = genome_fasta
    output:
        bam = "tagged_bam/lambda_phage/{sample}.bam",
        bai = "tagged_bam/lambda_phage/{sample}.bam.bai",
        bed = temp("meth_calls/lambda_phage/{sample}_mCpG.bed"),
        stats = "meth_calls/lambda_phage/{sample}_stats.txt"
    params:
        method = 'chic' if method == 'chic-taps' else 'nla',
        min_mq = min_mapq,
        cluster = '--cluster' if cluster else '',
        contig = 'Enterobacteria_phage_lambda_NC_001416.1'
    log: "logs/taps_tagger_phage_{sample}.err"
    threads: 1
    conda: CONDA_SHARED_ENV
    shell:
        "tapsTagger.py -contig {params.contig} {params.cluster} -context Z \
         -ref {input.genome} -method {params.method} -bed {output.bed} \
         -o {output.bam} -min_mq {params.min_mq} {input.bam} \
         > {output.stats} 2> {log}"

rule lamda_stats:
    input: expand("meth_calls/lambda_phage/{sample}_stats.txt", sample = samples)
    output: "QC/lambda_stats.txt"
    threads: 1
    conda: CONDA_SHARED_ENV
    shell:
        """
        for file in {input}; \
        do echo $(basename $file _stats.txt) >> files.txt; \
        grep "z" $file | cut -f2 >> z.txt; \
        grep "Z" $file | cut -f2 >> Z.txt; \
        done && paste files.txt z.txt Z.txt > {output} && \
        rm files.txt z.txt Z.txt
        """

rule meth_bigwig_phage:
    input:
        bam = "tagged_bam/lambda_phage/{sample}.bam",
        bai = "tagged_bam/lambda_phage/{sample}.bam.bai",
        bed = "meth_calls/lambda_phage/{sample}_mCpG.bed",
        sizes = "chrom_sizes.txt"
    output:
        bw = temp("meth_calls/lambda_phage/{sample}.methCpG.bw")
    params:
        sample = 'phage_{sample}',
        genome = genome_fasta,
        context = "'Z'",
        grep = ""
    log: "logs/meth_bigwig_phage_{sample}.err"
    threads: 1
    conda: CONDA_SHARED_ENV
    shell: methBigwig_cmd

rule all_bigwig_phage:
    input:
        bam = "tagged_bam/lambda_phage/{sample}.bam",
        bai = "tagged_bam/lambda_phage/{sample}.bam.bai",
        bed = "meth_calls/lambda_phage/{sample}_mCpG.bed",
        sizes = "chrom_sizes.txt"
    output:
        bw = temp("meth_calls/lambda_phage/{sample}.all.bw")
    params:
        sample = 'phage_{sample}',
        genome = genome_fasta,
        context = "'z\|Z'",
        grep = ""
    log: "logs/meth_bigwig_phage_{sample}.err"
    threads: 1
    conda: CONDA_SHARED_ENV
    shell: methBigwig_cmd

## ratio meth/all
rule bigwigRatio_lambda:
    input:
        methbw = "meth_calls/lambda_phage/{sample}.methCpG.bw",
        allbw = "meth_calls/lambda_phage/{sample}.all.bw"
    output: "meth_calls/lambda_phage/{sample}.methRatio.bw"
    params:
        blklist = blacklist_bed
    log: "logs/bigwigRatio_phage_{sample}.err"
    threads: 10
    conda: CONDA_SHARED_ENV
    shell:
        "bigwigCompare -p {threads} -bl {params.blklist} -bs 1 --skipZeroOverZero \
        --operation ratio -o {output} --skipNAs -b1 {input.methbw} -b2 {input.allbw}"
