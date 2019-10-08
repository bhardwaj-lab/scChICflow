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
        cluster = '--cluster' if cluster else ''
    log: "logs/taps_tagger_phage_{sample}.err"
    threads: 1
    conda: CONDA_SHARED_ENV
    shell:
        "tapsTagger.py {params.cluster} -context Z \
         -ref {input.genome} -method {params.method} -bed {output.bed} \
         -o {output.bam} -min_mq {params.min_mq} {input.bam} \
         > {output.stats} 2> {log.err}"

rule meth_bigwig_phage:
    input:
        bam = "tagged_bam/lambda_phage/{sample}.bam",
        bai = "tagged_bam/lambda_phage/{sample}.bam.bai",
        bed = "meth_calls/lambda_phage/{sample}_mCpG.bed"
    output:
        bw = "meth_calls/lambda_phage/{sample}.methCpG.bw"
    params:
        sample = '{sample}',
        genome = genome_fasta
    log: "logs/meth_bigwig_phage_{sample}.err"
    threads: 1
    conda: CONDA_SHARED_ENV
    shell:
        """
        sort -k1,1 -k2,2n | {input.bed} \
        bedtools merge -s -c 4,5,6 -o distinct,sum,distinct -i - >\
        {params.sample}.tmp.bed &&
        bedGraphToBigWig {params.sample}.tmp.bed chrom_sizes.txt {output.bw} 2> {log} &&
        rm {params.sample}.tmp.bed
        """

rule lamda_stats:
    input: expand("meth_calls/lambda_phage/{sample}_stats.txt", sample = samples)
    output: "QC/lambda_stats.txt"
    threads: 1
    conda: CONDA_SHARED_ENV
    shell:
        """
        for file in {input} \
        do echo $(basename $file _stats.txt) >> files.txt; \
        grep "z" $file | cut -f2 >> z.txt; \
        grep "Z" $file | cut -f2 >> Z.txt; \
        done && paste files.txt z.txt Z.txt > {output} && \
        rm files.txt z.txt Z.txt
        """
