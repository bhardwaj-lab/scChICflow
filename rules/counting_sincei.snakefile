if countRegions == "windows":
    rule count_regions:
        input:
            bam = expand("dedup_bam/{sample}.bam", sample = samples),
            bai = expand("dedup_bam/{sample}.bam.bai", sample = samples),
            barcodes = barcode_list,
            blk = blacklist_bed
        output: "counts/scCounts_"+binSize+"bp_bins.loom"
        params:
            bin = binSize,
            filters= "--minAlignedFraction 0.6 --GCcontentFilter '0.2,0.8'",
            prefix = "counts/scCounts_"+binSize+"bp_bins"
            #path = sincei_path if sincei_path else ""
        log: "logs/sincei_count_windows.err"
        threads: 10
        conda: CONDA_SHARED_ENV
        shell:
            "scCountReads bins \
            {params.filters} \
            --barcodes {input.barcodes} \
            -bl {input.blk} \
            -b {input.bam} \
            --smartLabels -p {threads} -bs {params.bin} \
            -o {params.prefix} > {log} 2>&1"

elif countRegions == "bed" or countRegions == "peaks":
    rule count_regions:
        input:
            bam = expand("dedup_bam/{sample}.bam", sample = samples),
            bai = expand("dedup_bam/{sample}.bam.bai", sample = samples),
            barcodes = barcode_list,
            blk = blacklist_bed,
            bed = lambda wildcards: bedFile if countRegions == "bed" else "macs2_peaks/peaks_union.bed"
        output: "counts/scCounts_"+countRegions+".loom"
        params:
            filters= "--minAlignedFraction 0.6 --GCcontentFilter '0.2,0.8'",
            prefix = "counts/scCounts_"+countRegions,
            #path = sincei_path if sincei_path else ""
        log: "logs/sincei_count_bed.err"
        threads: 10
        conda: CONDA_SHARED_ENV
        shell:
            "scCountReads features \
            --BED {input.bed} \
            {params.filters} \
            --barcodes {input.barcodes} \
            -bl {input.blk} \
            -b {input.bam} \
            --smartLabels -p {threads} \
            -o {params.prefix} > {log} 2>&1"

elif countRegions == "genes":
    print("Counting reads in genes per cell")
