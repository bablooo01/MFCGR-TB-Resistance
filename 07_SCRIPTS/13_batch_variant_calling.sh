#!/bin/bash

# ============================================================
# MFCGR - BATCH VARIANT CALLING
# bcftools haploid variant calling for 120 TB isolates
# ============================================================

PROJECT="/mnt/d/MFCGR_PROJECT"

REFERENCE="$PROJECT/01_H37Rv/ncbi_dataset/data/GCF_000195955.2/GCF_000195955.2_ASM19595v2_genomic.fna"

ALIGN_DIR="$PROJECT/08_RESULTS/alignment"
VCF_DIR="$PROJECT/08_RESULTS/variants"

THREADS=4

mkdir -p "$VCF_DIR"

TOTAL=0
SUCCESS=0
SKIPPED=0
FAILED=0

echo "============================================================"
echo "MFCGR BATCH VARIANT CALLING"
echo "============================================================"
echo "Reference : $REFERENCE"
echo "BAM dir   : $ALIGN_DIR"
echo "VCF dir   : $VCF_DIR"
echo

for BAM in "$ALIGN_DIR"/*.sorted.bam
do

    BAM_NAME=$(basename "$BAM")
    ISOLATE="${BAM_NAME%.sorted.bam}"

    TOTAL=$((TOTAL + 1))

    echo
    echo "------------------------------------------------------------"
    echo "Processing: $ISOLATE"
    echo "------------------------------------------------------------"

    RAW_VCF="$VCF_DIR/${ISOLATE}.haploid.raw.vcf.gz"
    FILTERED_VCF="$VCF_DIR/${ISOLATE}.filtered.vcf.gz"

    # --------------------------------------------------------
    # Skip if filtered VCF already exists
    # --------------------------------------------------------

    if [ -f "$FILTERED_VCF" ] && [ -f "${FILTERED_VCF}.tbi" ]; then
        echo "Already completed - skipping"
        SKIPPED=$((SKIPPED + 1))
        continue
    fi

    # --------------------------------------------------------
    # Step 1: Haploid variant calling
    # --------------------------------------------------------

    if [ ! -f "$RAW_VCF" ]; then

        echo "Step 1: bcftools mpileup + haploid calling..."

        bcftools mpileup \
            -Ou \
            -f "$REFERENCE" \
            "$BAM" \
        | bcftools call \
            -mv \
            --ploidy 1 \
            -Oz \
            -o "$RAW_VCF"

        if [ $? -ne 0 ]; then
            echo "ERROR: Variant calling failed"
            FAILED=$((FAILED + 1))
            rm -f "$RAW_VCF"
            continue
        fi

    else

        echo "Raw VCF already exists - skipping calling"

    fi

    # --------------------------------------------------------
    # Step 2: Index raw VCF
    # --------------------------------------------------------

    if [ ! -f "${RAW_VCF}.tbi" ]; then

        echo "Indexing raw VCF..."

        bcftools index \
            -t \
            "$RAW_VCF"

        if [ $? -ne 0 ]; then
            echo "ERROR: Raw VCF indexing failed"
            FAILED=$((FAILED + 1))
            continue
        fi

    fi

    # --------------------------------------------------------
    # Step 3: Quality filtering
    #
    # Same criteria validated during pilot:
    # QUAL >= 30
    # INFO/DP >= 10
    # INFO/MQ >= 40
    # --------------------------------------------------------

    echo "Step 2: Quality filtering..."

    bcftools view \
        -i 'QUAL>=30 && INFO/DP>=10 && INFO/MQ>=40' \
        "$RAW_VCF" \
        -Oz \
        -o "$FILTERED_VCF"

    if [ $? -ne 0 ]; then
        echo "ERROR: Filtering failed"
        FAILED=$((FAILED + 1))
        rm -f "$FILTERED_VCF"
        continue
    fi

    # --------------------------------------------------------
    # Step 4: Index filtered VCF
    # --------------------------------------------------------

    echo "Indexing filtered VCF..."

    bcftools index \
        -t \
        "$FILTERED_VCF"

    if [ $? -ne 0 ]; then
        echo "ERROR: Filtered VCF indexing failed"
        FAILED=$((FAILED + 1))
        continue
    fi

    echo "Variant calling completed: $ISOLATE"

    SUCCESS=$((SUCCESS + 1))

done

echo
echo "============================================================"
echo "BATCH VARIANT CALLING COMPLETE"
echo "============================================================"
echo "Total   : $TOTAL"
echo "Success : $SUCCESS"
echo "Skipped : $SKIPPED"
echo "Failed  : $FAILED"
echo
echo "VCF directory:"
echo "$VCF_DIR"
echo "============================================================"