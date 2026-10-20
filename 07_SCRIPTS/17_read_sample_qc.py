import gzip
from pathlib import Path
from collections import Counter

PROJECT = Path("/mnt/d/MFCGR_PROJECT")

FASTQ = (
    PROJECT
    / "08_RESULTS"
    / "trimmed_fastq"
    / "ERR4810741_1.trimmed.fastq.gz"
)

SAMPLE_READS = 10000

lengths = []
gc_values = []
base_counts = Counter()

print("=" * 70)
print("READ-LEVEL DIAGNOSTIC: ERR4810741")
print("=" * 70)

with gzip.open(FASTQ, "rt") as f:

    for i in range(SAMPLE_READS):

        header = f.readline()
        sequence = f.readline().strip()
        plus = f.readline()
        quality = f.readline()

        if not sequence:
            break

        lengths.append(len(sequence))

        counts = Counter(sequence.upper())
        base_counts.update(counts)

        gc = (
            (counts.get("G", 0) + counts.get("C", 0))
            / len(sequence)
            * 100
        )

        gc_values.append(gc)


print(f"\nReads sampled       : {len(lengths):,}")

print("\nREAD LENGTH")
print(f"Minimum             : {min(lengths)} bp")
print(f"Maximum             : {max(lengths)} bp")
print(f"Mean                : {sum(lengths) / len(lengths):.2f} bp")

print("\nGC CONTENT")
print(f"Minimum             : {min(gc_values):.2f}%")
print(f"Maximum             : {max(gc_values):.2f}%")
print(f"Mean                : {sum(gc_values) / len(gc_values):.2f}%")

print("\nBASE COMPOSITION")
total_bases = sum(base_counts.values())

for base in ["A", "C", "G", "T", "N"]:
    count = base_counts.get(base, 0)
    percentage = count / total_bases * 100
    print(f"{base}: {count:,} ({percentage:.2f}%)")

print()
print("=" * 70)
print("READ-LEVEL DIAGNOSTIC COMPLETE")
print("=" * 70)