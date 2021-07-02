if protocol == "chic":
    rule FASTQ:
          input:
              r1=indir+"/{sample}"+reads[0]+ext,
              r2=indir+"/{sample}"+reads[1]+ext
          output:
              r1="FASTQ/{sample}"+reads[0]+".fastq.gz",
              r2="FASTQ/{sample}"+reads[1]+".fastq.gz"
          shell:
              "( [ -f {output.r1} ] || ln -s -r {input.r1} {output.r1} ) && \
               ( [ -f {output.r2} ] || ln -s -r {input.r2} {output.r2} )"
else:
    rule split_tChIC:
        input:
            r1=indir+"/{sample}"+reads[0]+ext,
            nla=barcode_list,
            cs2=barcodes_cs2
        output:
            cs2=temp("{sample}.CS2.txt"),
            nla=temp("{sample}.NLA.txt")
        params:
            script = os.path.join(workflow.basedir, "tools", "split_fastq.py"),
            prefix="{sample}"
        threads: 10
        shell:
            "{params.script} --ncpus={threads} --infile={input.r1} --prefix={params.prefix} \
            --nla_bc={input.nla} --celseq_bc={input.cs2}"

    rule split_fastq_dna:
        input:
            r1=indir+"/{sample}"+reads[0]+ext,
            r2=indir+"/{sample}"+reads[1]+ext,
            names="{sample}.NLA.txt"
        output:
            r1="FASTQ/{sample}"+reads[0]+".fastq.gz",
            r2="FASTQ/{sample}"+reads[1]+".fastq.gz"
        log: "logs/{sample}_dnaSplit.log"
        shell:
            "filterbyname.sh -Xmx50g include=t in={input.r1} in2={input.r2} \
             out={output.r1} out2={output.r2} names={input.names} > {log} 2>&1"

    rule split_fastq_rna:
        input:
            r1=indir+"/{sample}"+reads[0]+ext,
            r2=indir+"/{sample}"+reads[1]+ext,
            names="{sample}.CS2.txt"
        output:
            r1="FASTQ_RNA/{sample}"+reads[0]+".fastq.gz",
            r2="FASTQ_RNA/{sample}"+reads[1]+".fastq.gz"
        log: "logs/{sample}_rnaSplit.log"
        shell:
            "filterbyname.sh -Xmx50g include=t in={input.r1} in2={input.r2} \
             out={output.r1} out2={output.r2} names={input.names} > {log} 2>&1"
