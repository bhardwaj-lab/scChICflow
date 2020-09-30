
### Cmd to make bigwigs ###
methBigwig_cmd1 = """
    grep {params.context} {input.bed} | sort -k1,1 -k2,2n |\
    bedtools merge -s -c 4,5,6 -o distinct,sum,distinct -i - |\
    cut -f1-3,5 {params.grep} > {output}
    """

methBigwig_cmd2 = "bedGraphToBigWig {input.bg} {input.sizes} {output} 2> {log}"

#methBigwig_cmd = """
#    grep {params.context} {input.bed} | sort -k1,1 -k2,2n |\
#    bedtools merge -s -c 4,5,6 -o distinct,sum,distinct -i - |\
#    cut -f1-3,5 {params.grep} > {params.sample}.tmp.bg &&
#    bedGraphToBigWig {params.sample}.tmp.bg {input.sizes} {output.bw} 2> {log} &&
#    rm {params.sample}.tmp.bg
#    """

## methylated C
rule meth_bigwig1:
    input:
        bed = "meth_calls/{sample}_methylation.bed",
        bam = "tagged_bam/{sample}.bam",
        bai = "tagged_bam/{sample}.bam.bai"
    output: temp("meth_calls/{sample}.meth.bg")
    params:
        sample = '{sample}',
        genome = genome_fasta,
        context = "'Z'",
        grep = "| grep -v 'phage\|G\|K\|J'"
    log: "logs/meth_bigwig1_{sample}.err"
    threads: 1
    conda: CONDA_SHARED_ENV
    shell: methBigwig_cmd1

rule meth_bigwig2:
    input:
        bg = "meth_calls/{sample}.meth.bg",
        sizes = "chrom_sizes.txt"
    output: "meth_calls/{sample}.methCpG.bw"
    log: "logs/meth_bigwig2_{sample}.err"
    threads: 1
    conda: CONDA_SHARED_ENV
    shell: methBigwig_cmd2

## all C
rule all_bigwig1:
    input:
        bed = "meth_calls/{sample}_methylation.bed",
        bam = "tagged_bam/{sample}.bam",
        bai = "tagged_bam/{sample}.bam.bai",
        sizes = "chrom_sizes.txt"
    output:
        bw = temp("meth_calls/{sample}.all.bg")
    params:
        sample = '{sample}',
        genome = genome_fasta,
        context = "'z\|Z'",
        grep = "| grep -v 'phage\|G\|K\|J'"
    log: "logs/all_bigwig1_{sample}.err"
    threads: 1
    conda: CONDA_SHARED_ENV
    shell: methBigwig_cmd1

rule all_bigwig2:
    input:
        bg = "meth_calls/{sample}.all.bg",
        sizes = "chrom_sizes.txt"
    output: "meth_calls/{sample}.allCpG.bw"
    log: "logs/all_bigwig2_{sample}.err"
    threads: 1
    conda: CONDA_SHARED_ENV
    shell: methBigwig_cmd2

## ratio meth/all
rule bigwigRatio:
    input:
        methbw = "meth_calls/{sample}.methCpG.bw",
        allbw = "meth_calls/{sample}.allCpG.bw"
    output: "coverage/{sample}.methRatio.bw"
    params:
        blklist = blacklist_bed
    log: "logs/bigwigRatio_{sample}.err"
    threads: 30
    conda: CONDA_SHARED_ENV
    shell:
        "bigwigCompare -p {threads} -bl {params.blklist} -bs 1 --skipZeroOverZero \
        --operation ratio -o {output} --skipNAs -b1 {input.methbw} -b2 {input.allbw}"

#rule meth_gzip:
#    input:
#        bed = "meth_calls/{sample}_methylation.bed",
#        methbw = "meth_calls/{sample}.methCpG.bw",
#        allbw = "meth_calls/{sample}.allCpG.bw",
#        bw = "coverage/{sample}.methRatio.bw"
#    output: "meth_calls/{sample}_methylation.bed.gz"
#    threads: 1
#    shell: 'gzip {input.bed}'


rule bwSummary_meth:
    input:
        bw = expand("coverage/{sample}.methRatio.bw", sample = samples),
        bl = blacklist_bed
    output: "QC/bwSummary_methRatio_10kBins.npz"
    log: "logs/bwSummary_methCpG.log"
    threads: 8
    conda: CONDA_SHARED_ENV
    shell:
        "multiBigwigSummary bins -p {threads} -bl {input.bl} \
        --smartLabels -b {input.bw} -o {output}  > {log} 2>&1"

###-------------------------- PHAGE -----------
rule meth_bigwig_phage1:
    input:
        bam = "tagged_bam/lambda_phage/{sample}.bam",
        bai = "tagged_bam/lambda_phage/{sample}.bam.bai",
        bed = "meth_calls/lambda_phage/{sample}_mCpG.bed"
    output:
        bw = temp("meth_calls/lambda_phage/{sample}.lmeth.bg")
    params:
        sample = 'phage_{sample}',
        genome = genome_fasta,
        context = "'Z'",
        grep = ""
    log: "logs/meth_bigwig_phage1_{sample}.err"
    threads: 1
    conda: CONDA_SHARED_ENV
    shell: methBigwig_cmd1


rule meth_bigwig_phage2:
    input:
        bg = "meth_calls/lambda_phage/{sample}.lmeth.bg",
        sizes = "chrom_sizes.txt"
    output: temp("meth_calls/lambda_phage/{sample}.methCpG.bw")
    log: "logs/meth_bigwig_phage2_{sample}.err"
    threads: 1
    conda: CONDA_SHARED_ENV
    shell: methBigwig_cmd2


rule all_bigwig_phage1:
    input:
        bam = "tagged_bam/lambda_phage/{sample}.bam",
        bai = "tagged_bam/lambda_phage/{sample}.bam.bai",
        bed = "meth_calls/lambda_phage/{sample}_mCpG.bed"
    output:
        bw = temp("meth_calls/lambda_phage/{sample}.lall.bg")
    params:
        sample = 'phage_{sample}',
        genome = genome_fasta,
        context = "'z\|Z'",
        grep = ""
    log: "logs/meth_bigwig_phage1_{sample}.err"
    threads: 1
    conda: CONDA_SHARED_ENV
    shell: methBigwig_cmd1

rule all_bigwig_phage2:
    input:
        bg = "meth_calls/lambda_phage/{sample}.lall.bg",
        sizes = "chrom_sizes.txt"
    output: temp("meth_calls/lambda_phage/{sample}.all.bw")
    log: "logs/meth_bigwig_phage2_{sample}.err"
    threads: 1
    conda: CONDA_SHARED_ENV
    shell: methBigwig_cmd2


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
