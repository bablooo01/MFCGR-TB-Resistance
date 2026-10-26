
import gzip
import pandas as pd
from pathlib import Path

ROOT = Path(r"D:\MFCGR_PROJECT")
VCF_DIR = ROOT / "08_RESULTS" / "variants"
COHORT_FILE = ROOT / "08_RESULTS" / "cohort_accessions_for_methylation.csv"
OUT = ROOT / "08_RESULTS" / "top_variant_vcf_validation.csv"

TARGETS = {
    ("NC_000962.3", 3643985): {
        "gene": "Rv3263",
        "ref": "A",
        "alt": "C",
    },
    ("NC_000962.3", 2269545): {
        "gene": "Rv2024c",
        "ref": "G",
        "alt": "A",
    },
}

def open_vcf(path):
    if path.suffix == ".gz":
        return gzip.open(path, "rt")
    return open(path, "r", encoding="utf-8")

def inspect_vcf(path):
    found = {}

    with open_vcf(path) as handle:
        for line in handle:
            if line.startswith("#"):
                continue

            fields = line.rstrip("\n").split("\t")
            if len(fields) < 8:
                continue

            chrom = fields[0]
            pos = int(fields[1])
            key = (chrom, pos)

            if key not in TARGETS:
                continue

            fmt = fields[8].split(":") if len(fields) > 8 else []
            sample = fields[9].split(":") if len(fields) > 9 else []
            values = dict(zip(fmt, sample))

            found[key] = {
                "VCF_REF": fields[3],
                "VCF_ALT": fields[4],
                "QUAL": fields[5],
                "FILTER": fields[6],
                "INFO": fields[7],
                "GT": values.get("GT", ""),
                "AD": values.get("AD", ""),
                "DP_sample": values.get("DP", ""),
                "GQ": values.get("GQ", ""),
            }

    return found

cohort = pd.read_csv(COHORT_FILE)
print("Cohort columns:", cohort.columns.tolist())

# Identify the run-accession column in the cohort file.
run_col = next(
    (c for c in cohort.columns if c.lower() == "run_accession"),
    None
)
if run_col is None:
    raise ValueError("Could not find run_accession column in cohort CSV")

accessions = cohort[run_col].dropna().astype(str).unique()
rows = []

for accession in accessions:
    path = VCF_DIR / f"{accession}.filtered.vcf.gz"

    if not path.exists():
        rows.append({
            "run_accession": accession,
            "gene": "",
            "position": "",
            "status": "VCF_MISSING",
        })
        continue

    calls = inspect_vcf(path)

    for key, target in TARGETS.items():
        record = calls.get(key)

        row = {
            "run_accession": accession,
            "gene": target["gene"],
            "position": key[1],
            "expected_REF": target["ref"],
            "expected_ALT": target["alt"],
        }

        if record is None:
            row.update({
                "status": "NO_VARIANT_RECORD",
                "GT": "",
                "AD": "",
                "DP_sample": "",
                "GQ": "",
                "VCF_REF": "",
                "VCF_ALT": "",
                "QUAL": "",
                "FILTER": "",
                "INFO": "",
            })
        else:
            row.update(record)
            if record["VCF_REF"] != target["ref"]:
                row["status"] = "REF_MISMATCH"
            elif target["alt"] not in record["VCF_ALT"].split(","):
                row["status"] = "ALT_MISMATCH"
            else:
                row["status"] = "VARIANT_RECORD_FOUND"

        rows.append(row)

result = pd.DataFrame(rows)
result.to_csv(OUT, index=False)

print("\nTotal cohort accessions:", len(accessions))
print("Validation rows:", len(result))
print("\nStatus counts:")
print(result["status"].value_counts(dropna=False).to_string())

print("\nResults by gene and status:")
print(
    result.groupby(["gene", "status"], dropna=False)
    .size().rename("records").reset_index().to_string(index=False)
)

print("\nSaved:", OUT)
print(
    "\nIMPORTANT: NO_VARIANT_RECORD does not automatically mean reference "
    "genotype. We must inspect callable coverage and the VCF FORMAT fields."
)
