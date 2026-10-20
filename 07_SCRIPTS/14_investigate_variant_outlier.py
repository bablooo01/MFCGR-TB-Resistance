import subprocess
from pathlib import Path

PROJECT = Path("/mnt/d/MFCGR_PROJECT")

VCF_DIR = PROJECT / "08_RESULTS" / "variants"
ALIGN_DIR = PROJECT / "08_RESULTS" / "alignment"

ISOLATES = [
    "ERR4810741",
    "ERR038300",
]

for isolate in ISOLATES:

    print()
    print("=" * 70)
    print(f"VARIANT QC: {isolate}")
    print("=" * 70)

    vcf = VCF_DIR / f"{isolate}.filtered.vcf.gz"
    bam = ALIGN_DIR / f"{isolate}.sorted.bam"

    print("\nFiltered VCF:")
    result = subprocess.run(
        ["bcftools", "stats", str(vcf)],
        capture_output=True,
        text=True,
    )

    for line in result.stdout.splitlines():
        if line.startswith("SN") and (
            "number of records" in line
            or "number of SNPs" in line
            or "number of indels" in line
        ):
            print(line)

    print("\nBAM alignment:")
    result = subprocess.run(
        ["samtools", "flagstat", str(bam)],
        capture_output=True,
        text=True,
    )

    print(result.stdout)

print()
print("=" * 70)
print("OUTLIER INVESTIGATION COMPLETE")
print("=" * 70)