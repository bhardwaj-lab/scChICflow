rule taps_tagger:
    input:
        bam = "bwa_mapped/{sample}.bam",
        idx = "bwa_mapped/{sample}.bam.bai",
        genome = genome_fasta
    output:
        bam = "tagged_bam/{sample}.bam",
        bai = "tagged_bam/{sample}.bam.bai",
        bed = "meth_calls/{sample}_methylation.bed"
    params:
        method = 'chic' if method == 'chic-taps' else 'nla',
        min_mq = min_mapq,
        cluster = '--cluster' if cluster else '',
        context = '-context '+bedContext if bedContext else ''
    log:
        out = "logs/taps_tagger_{sample}.out",
        err = "logs/taps_tagger_{sample}.err"
    threads: 1
    conda: CONDA_SHARED_ENV
    shell:
        "tapsTagger.py {params.cluster} {params.context} \
        -ref {input.genome} -method {params.method} -bed {output.bed} \
        -o {output.bam} -min_mq {params.min_mq} {input.bam} \
        > {log.out} 2> {log.err}"

rule meth_bigwig:
    input:
        bed = "meth_calls/{sample}_methylation.bed",
        bam = "tagged_bam/{sample}.bam",
        bai = "tagged_bam/{sample}.bam.bai"
    output:
        bw = "meth_calls/{sample}.methCpG.bw",
        gz = "meth_calls/{sample}_methylation.bed.gz"
    params:
        sample = '{sample}',
        genome = genome_fasta
    log:
        out = "logs/meth_bigwig_{sample}.out",
        err = "logs/meth_bigwig_{sample}.err"
    threads: 1
    conda: CONDA_SHARED_ENV
    shell:
        """
        grep 'Z' {input.bed} | sort -k1,1 -k2,2n |\
        bedtools merge -s -c 4,5,6 -o distinct,sum,distinct -i - |\
        cut -f1-3,5 | grep -v 'phage\|G\|K\|J' > {params.sample}.tmp.bed &&
        cut -f1-2 {params.genome}.fai > chrom_sizes.txt 2> {log.err} &&
        bedGraphToBigWig {params.sample}.tmp.bed chrom_sizes.txt {output.bw} 2> {log.err} &&
        gzip {input.bed} &&
        rm {params.sample}.tmp.bed chrom_sizes.txt
        """

rule meth_bincounts:
    input:
        bam = "tagged_bam/{sample}.bam",
        bai = "tagged_bam/{sample}.bam.bai"
    output:
        csv = "meth_calls/{sample}_CpG_binCounts.csv"
    params:
        binsize = methCountsBinSize
    log: "logs/meth_bincounts_{sample}.log"
    threads: 1
    conda: CONDA_SHARED_ENV
    shell:
        "bamToCountTable.py -sampleTags SM -featureTags sZ -byValue sZ \
        -joinedFeatureTags reference_name -bin {params.binsize} -o {output.csv} {input.bam} > {log} 2>&1"
