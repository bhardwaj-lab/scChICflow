# @author: Vivek Bhardwaj (@vivekbhr)
## TODO: switch the star_bugfix path with star version 2.7.10b,
import os
import glob

STARsoloCoords = [1, 6, 7, 8] ## UMI start, UMI length, Cell BC start, cell BC length (1-based)

# optional args
indir='FASTQ_RNA'
outdir='RNA_MAPPING_STAR'

if trim:
    trim_dir = os.path.join(outdir, "FASTQ_trimmed")
else:
    trim_dir = indir

if trim:
    rule cutadapt_rna:
        input:
            R1 = os.path.join(indir,"{sample}"+reads[0]+".fastq.gz"),
            R2 = os.path.join(indir,"{sample}"+reads[1]+".fastq.gz")
        output:
            R1 = temp(os.path.join(outdir, "FASTQ_trimmed", "{sample}"+reads[0]+".fastq.gz")),
            R2 = temp(os.path.join(outdir, "FASTQ_trimmed", "{sample}"+reads[1]+".fastq.gz")),
            QC = os.path.join(outdir, "QC", "cutadapt", "{sample}.out")
        params:
            opts = "-A W{'10'}"
        threads: 10
        resources:
            mem_mb=50000
        shell:
            """
            cutadapt -j {threads} {params.opts} -e 0.1 -q 16 -O 3 --trim-n --minimum-length 10 \
            -a AGATCGGAAGAGC -A AGATCGGAAGAGC --nextseq-trim=16 \
            -b TGGAATTCTCGGGTGCCAAGG -B TGGAATTCTCGGGTGCCAAGG \
            -b ATCTCGTATGCCGTCTTCTGCTTG -B ATCTCGTATGCCGTCTTCTGCTTG \
            -b GTTCAGAGTTCTACAGTCCGACGATC -B GTTCAGAGTTCTACAGTCCGACGATC \
            -o "{output.R1}" -p "{output.R2}" "{input.R1}" "{input.R2}" > {output.QC}
            """

rule FastQC_rna:
    input:
        os.path.join(trim_dir, "{sample}{read}.fastq.gz")
    output:
        os.path.join(trim_dir, "FastQC/{sample}{read}_fastqc.html")
    params:
        outdir = os.path.join(trim_dir, "FastQC")
    log:
        out = os.path.join(outdir, "logs/FastQC.{sample}{read}.out"),
        err = os.path.join(outdir, "logs/FastQC.{sample}{read}.err")
    threads: 1
    shell: "fastqc -o {params.outdir} {input} > {log.out} 2> {log.err}"

# mapped using STARsolo so possible to split SAM files per cell type later on
rule mapReads:
    input:
        read1 = os.path.join(trim_dir, "{sample}"+reads[0]+".fastq.gz"),
        read2 = os.path.join(trim_dir, "{sample}"+reads[1]+".fastq.gz")
    output:
        bam = os.path.join(outdir, "STARsolo/{sample}/{sample}.Aligned.sortedByCoord.out.bam"),
        raw_counts = os.path.join(outdir, "STARsolo/{sample}/{sample}.Solo.out/Gene/raw/matrix.mtx"),
        filtered_counts = os.path.join(outdir, "STARsolo/{sample}/{sample}.Solo.out/Gene/filtered/matrix.mtx"),
        filtered_bc = os.path.join(outdir, "STARsolo/{sample}/{sample}.Solo.out/Gene/filtered/barcodes.tsv"),
        tmpdir=temp(directory(os.path.join(outdir, "STARsolo/{sample}/{sample}._STARtmp"))),
        tmpgenome=temp(directory(os.path.join(outdir, "STARsolo/{sample}/{sample}._STARgenome")))
    params:
        star_bugfix='/hpc/hub_oudenaarden/vbhardwaj/programs/STAR-2.7.10a_alpha_220506/source',# use this unil version 2.7.11 is released, to avoid segfault
        gtf_file = gtf_file,
        index = star_index,
        prefix = "STARsolo/{sample}/{sample}.",
        samsort_memory = '2G',
        sample_dir = "STARsolo/{sample}",
        bclist = barcodes_cs2,
        UMIstart = STARsoloCoords[0],
        UMIlen = STARsoloCoords[1],
        CBstart = STARsoloCoords[2],
        CBlen = STARsoloCoords[3],
        outdir = "STARsolo/",
        tempDir = tempDir,
        sample = "{sample}"
    log: os.path.join(outdir, "logs/mapReads_{sample}.err")
    threads: 20
    resources:
        mem_mb=100000
    shell:
        """
    TMPDIR={params.tempDir}
    MYTEMP=$(mktemp -d ${{TMPDIR:-/tmp}}/snakepipes.XXXXXXXXXX);
    ( [ -d {params.sample_dir} ] || mkdir -p {params.sample_dir} )
    maxIntronLen=`expr $(awk '{{ print $3-$2 }}' {params.index}/sjdbList.fromgtf_file.out.tab | sort -n -r | head -1) + 1`
    {params.star_bugfix}/STAR --runThreadN {threads} \
    --sjdbOverhang 100 \
    --outSAMunmapped Within \
    --outSAMtype BAM SortedByCoordinate \
    --outBAMsortingBinsN 20 \
    --genomeDir {params.index} \
    --readFilesIn {input.read2} {input.read1}\
    --readFilesCommand zcat \
    --outFileNamePrefix {params.prefix} \
    --outSAMattributes NH HI AS nM MD jM jI MC ch CB UB GX GN \
    --outSAMstrandField intronMotif \
    --sjdbgtf_filefile {params.gtf_file} \
    --outFilterType BySJout \
    --outFilterIntronMotifs RemoveNoncanonical \
    --alignIntronMax ${{maxIntronLen}} \
    --alignSJoverhangMin 8 \
    --soloFeatures Gene Velocyto \
    --soloCBstart {params.CBstart} \
    --soloCBlen {params.CBlen} \
    --soloUMIstart {params.UMIstart} \
    --soloUMIlen {params.UMIlen} \
    --soloCBwhitelist {params.bclist} \
    --soloBarcodeReadLength 0 \
    --soloStrand Forward \
    --soloCBmatchWLtype Exact \
    --soloType CB_UMI_Simple \
    --soloUMIdedup NoDedup 2> {log}
    rm -rf $MYTEMP
        """

#    --soloUMIdedup Exact
rule filterBAMunique:
    input: os.path.join(outdir, "STARsolo/{sample}/{sample}.Aligned.sortedByCoord.out.bam"),
    output: temp(os.path.join(outdir, "STARsolo/{sample}.uniqueReads.sam"))
    threads: 10
    resources:
        mem_mb=50000
    shell:
        "samtools view -@ {threads} -h -d CB -F 4 -F 256 -F 2048 -q 255 {input} | grep -w -v 'CB:Z:-' > {output}"

rule dedupBAMunique:
    input: os.path.join(outdir, "STARsolo/{sample}.uniqueReads.sam")
    output:
        bam=os.path.join(outdir, "STARsolo/{sample}.uniqueReads.bam"),
        qc=os.path.join(outdir, "QC/umi_dedup/{sample}_per_umi_per_position.tsv")
    params:
        sample = "{sample}",
        tempdir= tempDir
    log: os.path.join(outdir, "logs/filterBAM.{sample}.log")
    threads: 1
    resources:
        mem_mb=80000
    shell:
        """
        umi_tools dedup --per-cell --cell-tag CB --umi-tag UB --extract-umi-method tag \
        --method unique --spliced-is-unique \
        --output-stats=QC/umi_dedup/{params.sample} \
        --temp-dir {params.tempdir} -L {log} -i -I {input} > {output.bam}
        """

rule indexBAM:
    input: "STARsolo/{sample}.uniqueReads.bam"
    output: "STARsolo/{sample}.uniqueReads.bam.bai"
    log: os.path.join(outdir, "logs/indexBAM.{sample}.log")
    threads: 5
    shell: 'samtools index -@ {threads} {input}'

# Get bigwig files of the trimmed fastq files check for A/T stretches
rule getBW:
    input:
        bam = os.path.join(outdir, "STARsolo/{sample}.uniqueReads.bam"),
        idx = os.path.join(outdir, "STARsolo/{sample}.uniqueReads.bam.bai")
    output: os.path.join(outdir, "bigwigs/{sample}.bw")
    log: os.path.join(outdir, "logs/getBW_{sample}.err")
    threads: 20
    shell:
        "bamCoverage -p {threads} --normalizeUsing CPM -b {input.bam} -o {output} > {log} 2>&1"

rule multiQC:
    input: expand(os.path.join(outdir, "STARsolo/{sample}.uniqueReads.bam"), sample = samples)
    output: os.path.join(outdir, "QC/multiqc_report.html")
    params:
        outdir = os.path.join(outdir, "QC")
    log: os.path.join(outdir, "logs/multiqc.log")
    threads: 1
    shell:
        "multiqc -f -o {params.outdir} . > {log} 2>&1"
