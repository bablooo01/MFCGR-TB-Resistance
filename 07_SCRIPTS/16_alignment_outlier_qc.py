import subprocess
from pathlib import Path

PROJECT = Path("/mnt/d/MFCGR_PROJECT")
ALIGN_DIR = PROJECT / "08_RESULTS" / "alignment"

ISOLATES = [
    "ERR4810741",
    "ERR038300",
]

for isolate in ISOLATES:

    bam = ALIGN_DIR / f"{isolate}.sorted.bam"

    print()
    print("=" * 70)
    print(f"ALIGNMENT QC: {isolate}")
    print("=" * 70)

    print("\n--- SAMTOOLS FLAGSTAT ---")
    result = subprocess.run(
        ["samtools", "flagstat", str(bam)],
        capture_output=True,
        text=True,
    )
    print(result.stdout)

    print("--- SAMTOOLS IDXSTATS ---")
    result = subprocess.run(
        ["samtools", "idxstats", str(bam)],
        capture_output=True,
        text=True,
    )
    print(result.stdout)

print()
print("=" * 70)
print("ALIGNMENT OUTLIER QC COMPLETE")
print("=" * 70)