import csv
import gzip
import re
from pathlib import Path
from collections import defaultdict

ROOT = Path(r"D:\MFCGR_PROJECT")
RESULTS = ROOT / "08_RESULTS"
VCF_DIR = RESULTS / "variants"
COHORT_FILE = RESULTS / "cohort_accessions_for_methylation.csv"
GFF_FILE = ROOT / "01_H37RV" / "ncbi_dataset" / "data" / "GCF_000195955.2" / "genomic.gff"

TARGETS = {
    "Rv3263": "mamA_candidate",
    "Rv2756c": "hsdM",
    "Rv2024c": "mamB_candidate",
    "Rv2761c": "hsdS",

}

def parse_attributes(text):
    result = {}
    for item in text.strip().split(";"):
        if "=" in item:
            k, v = item.split("=", 1)
            result[k] = v
    return result

# Load cohort and category information
with open(COHORT_FILE, newline="", encoding="utf-8-sig") as f:
    cohort = list(csv.DictReader(f))

category = {r["run_accession"]: r["CATEGORY"] for r in cohort}
expected = set(category)

# Get target gene coordinates from GFF
genes = {}
with open(GFF_FILE, encoding="utf-8") as f:
    for line in f:
        if line.startswith("#"):
            continue
        fields = line.rstrip("\n").split("\t")
        if len(fields) != 9 or fields[2] != "gene":
            continue
        attrs = parse_attributes(fields[8])
        locus = attrs.get("locus_tag", "")
        if locus in TARGETS:
            genes[locus] = {
                "name": TARGETS[locus],
                "contig": fields[0],
                "start": int(fields[3]),
                "end": int(fields[4]),
                "strand": fields[6],
            }

missing_genes = set(TARGETS) - set(genes)
if missing_genes:
    raise SystemExit(f"ERROR: target genes missing from GFF: {sorted(missing_genes)}")

print("Target gene coordinates:")
for locus, g in genes.items():
    print(f'  {locus} ({g["name"]}): {g["contig"]}:{g["start"]}-{g["end"]} ({g["strand"]})')

# Select one VCF per cohort accession from the variants folder
vcf_by_accession = {}
for path in VCF_DIR.glob("*.filtered.vcf.gz"):
    match = re.match(r"(.+)\.filtered\.vcf\.gz$", path.name)
    if match and match.group(1) in expected:
        vcf_by_accession[match.group(1)] = path

missing_vcfs = sorted(expected - set(vcf_by_accession))
print(f"\nCohort isolates: {len(expected)}")
print(f"Matching filtered VCFs in variants folder: {len(vcf_by_accession)}")
print(f"Missing VCFs: {len(missing_vcfs)}")

if missing_vcfs:
    print("Missing accessions:", ", ".join(missing_vcfs))

variant_rows = []
counts = defaultdict(int)
isolate_gene = defaultdict(int)
reference_contigs = set()

for accession, path in sorted(vcf_by_accession.items()):
    with gzip.open(path, "rt", encoding="utf-8") as f:
        for line in f:
            if line.startswith("##contig="):
                m = re.search(r"ID=([^,>]+)", line)
                if m:
                    reference_contigs.add(m.group(1))
                continue
            if line.startswith("#"):
                continue

            fields = line.rstrip("\n").split("\t")
            if len(fields) < 8:
                continue

            chrom = fields[0]
            pos = int(fields[1])
            ref = fields[3]
            alt = fields[4]
            qual = fields[5]
            filt = fields[6]
            info = fields[7]

            for locus, g in genes.items():
                if chrom != g["contig"] or not (g["start"] <= pos <= g["end"]):
                    continue

                dp = re.search(r"(?:^|;)DP=([^;]+)", info)
                mq = re.search(r"(?:^|;)MQ=([^;]+)", info)

                variant_rows.append({
                    "run_accession": accession,
                    "CATEGORY": category.get(accession, ""),
                    "gene": g["name"],
                    "locus_tag": locus,
                    "contig": chrom,
                    "position": pos,
                    "REF": ref,
                    "ALT": alt,
                    "QUAL": qual,
                    "FILTER": filt,
                    "DP": dp.group(1) if dp else "",
                    "MQ": mq.group(1) if mq else "",
                    "gene_start": g["start"],
                    "gene_end": g["end"],
                    "strand": g["strand"],
                    "interpretation": "Gene-overlapping variant; protein effect not yet annotated",
                })
                counts[locus] += 1
                isolate_gene[(accession, locus)] += 1

if reference_contigs:
    print("VCF contigs:", ", ".join(sorted(reference_contigs)))
    print("GFF contigs:", ", ".join(sorted({g["contig"] for g in genes.values()})))

variant_file = RESULTS / "methyltransferase_gene_variants_120.csv"
with open(variant_file, "w", newline="", encoding="utf-8-sig") as f:
    cols = [
        "run_accession", "CATEGORY", "gene", "locus_tag", "contig",
        "position", "REF", "ALT", "QUAL", "FILTER", "DP", "MQ",
        "gene_start", "gene_end", "strand", "interpretation"
    ]
    w = csv.DictWriter(f, fieldnames=cols)
    w.writeheader()
    w.writerows(variant_rows)

summary_file = RESULTS / "methyltransferase_gene_variant_summary.csv"
with open(summary_file, "w", newline="", encoding="utf-8-sig") as f:
    cols = ["run_accession", "CATEGORY", "gene", "locus_tag", "gene_variant_count"]
    w = csv.DictWriter(f, fieldnames=cols)
    w.writeheader()
    for accession in sorted(expected):
        for locus, g in genes.items():
            w.writerow({
                "run_accession": accession,
                "CATEGORY": category.get(accession, ""),
                "gene": g["name"],
                "locus_tag": locus,
                "gene_variant_count": isolate_gene.get((accession, locus), 0),
            })

print(f"\nVariants overlapping target genes: {len(variant_rows)}")
print("\nVariant counts by gene:")
for locus, g in genes.items():
    print(f'  {locus} ({g["name"]}): {counts[locus]} variant records')

print("\nSaved outputs:")
print(variant_file)
print(summary_file)
print("\nNOTE: zero means no variant record overlaps the gene in this VCF; it does not prove normal methylation function.")


