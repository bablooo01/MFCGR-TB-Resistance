import csv
import subprocess
from pathlib import Path
from collections import Counter

ROOT = Path("/mnt/d/MFCGR_PROJECT")
ALIGN = ROOT / "08_RESULTS/alignment"
COHORT = ROOT / "08_RESULTS/cohort_accessions_for_methylation.csv"
OUT = ROOT / "08_RESULTS/cohort_top_variant_validation.csv"
REF = ROOT / "01_H37RV/ncbi_dataset/data/GCF_000195955.2/GCF_000195955.2_ASM19595v2_genomic.fna"

TARGETS = [
    ("Rv3263", 3643985, "A", "C"),
    ("Rv2024c", 2269545, "G", "A"),
]

def count_bases(bases, ref_base):
    counts = Counter()
    i = 0
    while i < len(bases):
        c = bases[i]
        if c == "^":
            i += 2  # skip read-start marker and mapping-quality character
            continue
        if c == "$":
            i += 1
            continue
        if c in "+-":
            i += 1
            while i < len(bases) and bases[i].isdigit():
                i += 1
            if i < len(bases):
                n = 0
                while i < len(bases) and n < 0:
                    i += 1
            continue
        if c in ".,":
            counts[ref_base] += 1
        elif c.upper() in "ACGT":
            counts[c.upper()] += 1
        i += 1
    return counts

with COHORT.open(newline="") as f:
    reader = csv.DictReader(f)
    rows = list(reader)

# Normalize header names and identify the accession/category columns.
def find_col(row, candidates):
    for key in row:
        if key and key.strip().lower() in candidates:
            return key
    return None

if not rows:
    raise SystemExit("ERROR: cohort CSV is empty.")

run_col = find_col(rows[0], {"run_accession", "run", "run accession"})
cat_col = find_col(rows[0], {"category", "resistance_category", "resistance category"})

if not run_col or not cat_col:
    raise SystemExit(f"ERROR: Could not identify accession/category columns. Headers: {list(rows[0])}")

results = []
missing = []

for idx, row in enumerate(rows, 1):
    accession = row[run_col].strip()
    category = row[cat_col].strip()
    bam = ALIGN / f"{accession}.sorted.bam"

    if not bam.exists() or not Path(str(bam) + ".bai").exists():
        missing.append(accession)
        continue

    for gene, pos, ref_base, alt_base in TARGETS:
        region = f"NC_000962.3:{pos}-{pos}"
        cmd = [
            "samtools", "mpileup", "-q", "20", "-Q", "20",
            "-f", str(REF), "-r", region, str(bam)
        ]
        proc = subprocess.run(cmd, capture_output=True, text=True)
        if proc.returncode != 0:
            raise RuntimeError(f"samtools failed for {accession} at {pos}: {proc.stderr}")

        lines = [x for x in proc.stdout.splitlines() if x.strip()]
        if not lines:
            depth, counts = 0, Counter()
        else:
            fields = lines[-1].split("\t")
            depth = int(fields[3])
            counts = count_bases(fields[4], fields[2])

        ref_count = counts[ref_base]
        alt_count = counts[alt_base]
        usable = sum(counts[b] for b in "ACGT")
        af = alt_count / usable if usable else None

        if usable < 10:
            call = "INSUFFICIENT_COVERAGE"
        elif af >= 0.8:
            call = "ALT_SUPPORTED"
        elif af <= 0.05:
            call = "REF_SUPPORTED"
        else:
            call = "MIXED_OR_LOW_FREQUENCY_ALT"

        results.append({
            "run_accession": accession,
            "CATEGORY": category,
            "gene": gene,
            "position": pos,
            "REF": ref_base,
            "ALT": alt_base,
            "mpileup_depth": depth,
            "usable_ACGT_depth": usable,
            "REF_reads": ref_count,
            "ALT_reads": alt_count,
            "ALT_frequency": "" if af is None else round(af, 5),
            "call": call,
        })

    if idx % 10 == 0:
        print(f"Processed {idx}/{len(rows)} isolates")

fields = [
    "run_accession", "CATEGORY", "gene", "position", "REF", "ALT",
    "mpileup_depth", "usable_ACGT_depth", "REF_reads", "ALT_reads",
    "ALT_frequency", "call"
]
with OUT.open("w", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=fields)
    writer.writeheader()
    writer.writerows(results)

print(f"\nCompleted: {len(results)} position checks")
print(f"Missing BAM/index pairs: {len(missing)}")
if missing:
    print("Missing accessions:", ", ".join(missing))
print(f"Output: {OUT}")
