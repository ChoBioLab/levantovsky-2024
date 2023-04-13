#!/bin/bash
input="/mnt/cho_lab/disk1/jiayuzh/tmp/B_Cells_bam_loc.txt"
while IFS= read -r line
do
    barcodes=$(echo $line | cut -d " " -f 1)
    bam=$(echo $line | cut -d " " -f 2)
    sam_body=$(echo $line | cut -d " " -f 3)
    sam_header=$(echo $line | cut -d " " -f 4)
    filtered_sam=$(echo $line | cut -d " " -f 5)
    filtered_bam=$(echo $line | cut -d " " -f 6)
    samtools view -H $bam > $sam_header # Save the header lines
    samtools view $bam | LC_ALL=C grep -F -f $barcodes > $sam_body # Filter alignments using barcode_1. Use LC_ALL=C to set C locale instead of UTF-8
    cat $sam_header $sam_body > $filtered_sam # Combine header and body
    samtools view -b $filtered_sam > $filtered_bam # Convert filtered.sam to BAM format
    rm $sam_body $sam_header $filtered_sam
done < "$input"
samtools merge /mnt/cho_lab/disk1/jiayuzh/tmp/RL_Hu_B_Cells.bam /mnt/cho_lab/disk1/jiayuzh/tmp/*B_Cells_filtered.bam # merge subsetted bams from all samples
rm /mnt/cho_lab/disk1/jiayuzh/tmp/*_filtered.bam