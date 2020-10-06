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
        method = 'nla',# phage reads always have nla adapter
        min_mq = min_mapq,
#        cluster = '--cluster' if cluster else '',
        contig = 'Enterobacteria_phage_lambda_NC_001416.1'
    log: "logs/taps_tagger_phage_{sample}.err"
    threads: 1
    conda: CONDA_SHARED_ENV
    shell:
        "tapsTagger.py -contig {params.contig} -context Z \
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
