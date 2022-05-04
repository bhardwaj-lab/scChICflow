if countRegions == "windows":
    rule count_regions:
        input:
            bam = expand("dedup_bam/{sample}.bam", sample = samples),
            bai = expand("dedup_bam/{sample}.bam.bai", sample = samples),
            barcodes = barcode_list,
            blk = blacklist_bed
        output:
            counts = "counts/scCounts_"+binSize+"bp_bins.counts.mtx",
            colnames = "counts/scCounts_"+binSize+"bp_bins.colnames.txt",
            rownames = "counts/scCounts_"+binSize+"bp_bins.rownames.txt",
        params:
            bin = binSize,
            prefix = "counts/scCounts_"+binSize+"bp_bins",
            path='~/programs/sincei/bin'
        log: "logs/sincei_count_windows.err"
        threads: 10
        conda: CONDA_SHARED_ENV
        shell:
            "{params.path}/scCountReads.py bins \
            --minAlignedFraction 0.6 --GCcontentFilter '0.2,0.8' \
            --barcodes {input.barcodes} \
            -bl {input.blk} \
            -b {input.bam} \
            --smartLabels -p {threads} -bs {params.bin} --outFileFormat mtx \
            -o {params.prefix} > {log} 2>&1"

elif countRegions == "bed" or countRegions == "peaks":
    rule count_regions:
        input:
            bam = expand("dedup_bam/{sample}.bam", sample = samples),
            bai = expand("dedup_bam/{sample}.bam.bai", sample = samples),
            barcodes = barcode_list,
            blk = blacklist_bed,
            bed = lambda wildcards: bedFile if countRegions == "bed" else "macs2_peaks/peaks_union.bed"
        output:
            counts = "counts/scCounts_"+countRegions+".counts.mtx",
            colnames = "counts/scCounts_"+countRegions+".colnames.txt",
            rownames = "counts/scCounts_"+countRegions+".rownames.txt",
        params:
            bin = binSize,
            prefix = "counts/scCounts_"+countRegions,
            path='~/programs/sincei/bin'
        log: "logs/sincei_count_bed.err"
        threads: 10
        conda: CONDA_SHARED_ENV
        shell:
            "{params.path}/scCountReads.py BED-file \
            --BED {input.bed} \
            --minAlignedFraction 0.6 --GCcontentFilter '0.2,0.8' \
            --barcodes {input.barcodes} \
            -bl {input.blk} \
            -b {input.bam} \
            --smartLabels -p {threads} -bs {params.bin} --outFileFormat mtx \
            -o {params.prefix} > {log} 2>&1"

elif countRegions == "genes":
    print("Counting reads in genes per cell")
