#!/bin/bash

# ============================================================
# MFCGR - BATCH READ ALIGNMENT
# Align 120 trimmed FASTQ pairs to H37Rv
# ============================================================

PROJECT="/mnt/d/MFCGR_PROJECT"

REFERENCE="$PROJECT/01_H37Rv/ncbi_dataset/data/GCF_000195955.2/GCF_000195955.2_ASM19595v2_genomic.fna"

FASTQ_DIR="$PROJECT/08_RESULTS/trimmed_fastq"
ALIGN_DIR="$PROJECT/08_RESULTS/alignment"

THREADS=4

mkdir -p "$ALIGN_DIR"

TOTAL=0
SUCCESS=0
SKIPPED=0
FAILED=0

echo "============================================================"
echo "MFCGR BATCH ALIGNMENT"
echo "============================================================"
echo "Reference : $REFERENCE"
echo "FASTQ dir : $FASTQ_DIR"
echo "Output    : $ALIGN_DIR"
echo "Threads   : $THREADS"
echo

for R1 in "$FASTQ_DIR"/*_1.trimmed.fastq.gz
do

    R1_NAME=$(basename "$R1")
    ISOLATE="${R1_NAME%_1.trimmed.fastq.gz}"

    R2="$FASTQ_DIR/${ISOLATE}_2.trimmed.fastq.gz"

    TOTAL=$((TOTAL + 1))

    echo
    echo "------------------------------------------------------------"
    echo "Processing: $ISOLATE"
    echo "------------------------------------------------------------"

    if [ ! -f "$R2" ]; then
        echo "ERROR: R2 file missing"
        FAILED=$((FAILED + 1))
        continue
    fi

    SORTED_BAM="$ALIGN_DIR/${ISOLATE}.sorted.bam"
    BAM_INDEX="$ALIGN_DIR/${ISOLATE}.sorted.bam.bai"

    # Skip if BAM and index already exist
    if [ -f "$SORTED_BAM" ] && [ -f "$BAM_INDEX" ]; then
        echo "Already completed - skipping"
        SKIPPED=$((SKIPPED + 1))
        continue
    fi

    SAM="$ALIGN_DIR/${ISOLATE}.sam"

    echo "Step 1: BWA alignment..."

    bwa mem -t "$THREADS" \
        "$REFERENCE" \
        "$R1" \
        "$R2" \
        > "$SAM"

    if [ $? -ne 0 ]; then
        echo "ERROR: BWA alignment failed"
        FAILED=$((FAILED + 1))
        rm -f "$SAM"
        continue
    fi

    echo "Step 2: Sorting BAM..."

    samtools sort -@ "$THREADS" \
        -o "$SORTED_BAM" \
        "$SAM"

    if [ $? -ne 0 ]; then
        echo "ERROR: BAM sorting failed"
        FAILED=$((FAILED + 1))
        rm -f "$SAM"
        continue
    fi

    echo "Step 3: Indexing BAM..."

    samtools index "$SORTED_BAM"

    if [ $? -ne 0 ]; then
        echo "ERROR: BAM indexing failed"
        FAILED=$((FAILED + 1))
        rm -f "$SAM"
        continue
    fi

    rm -f "$SAM"

    echo "Alignment completed: $ISOLATE"

    SUCCESS=$((SUCCESS + 1))

done

echo
echo "============================================================"
echo "BATCH ALIGNMENT COMPLETE"
echo "============================================================"
echo "Total   : $TOTAL"
echo "Success : $SUCCESS"
echo "Skipped : $SKIPPED"
echo "Failed  : $FAILED"
echo
echo "Alignment directory:"
echo "$ALIGN_DIR"
echo "============================================================"