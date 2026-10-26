import csv
from pathlib import Path
from collections import Counter

ROOT = Path("/mnt/d/MFCGR_PROJECT")
RES = ROOT / "08_RESULTS"
FASTA = ROOT / "01_H37RV/ncbi_dataset/data/GCF_000195955.2/GCF_000195955.2_ASM19595v2_genomic.fna"
GFF = ROOT / "01_H37RV/ncbi_dataset/data/GCF_000195955.2/genomic.gff"
INPUT = RES / "methyltransferase_variant_consequences.csv"
OUTPUT = RES / "consequence_independent_validation.csv"
SUMMARY = RES / "consequence_validation_summary.csv"

COMP = str.maketrans("ACGTN", "TGCAN")
CODON_TABLE = {
    a+b+c: aa
    for a, row in zip("TCAG", [
        ["FFLL", "SSSS", "YY**", "CC*W"],
        ["LLLL", "PPPP", "HHQQ", "RRRR"],
        ["IIIM", "TTTT", "NNKK", "SSRR"],
        ["VVVV", "AAAA", "DDEE", "GGGG"]
    ])
    for b, col in zip("TCAG", row)
    for c, aa in zip("TCAG", col)
}

# Use an explicit standard genetic code to avoid external packages.
bases = "TCAG"
table_rows = [
    "FFLLSSSSYY**CC*W", "LLLLPPPPHHQQRRRR",
    "IIIMTTTTNNKKSSRR", "VVVVAAAADDEEGGGG"
]
CODON = {}
for i, a in enumerate(bases):
    for j, b in enumerate(bases):
        for k, c in enumerate(bases):
            CODON[a+b+c] = table_rows[i][j*4+k]

def translate(codon):
    return CODON.get(codon, "X")

def attrs(s):
    d = {}
    for item in s.split(";"):
        if "=" in item:
            k, v = item.split("=", 1)
            d[k] = v
    return d

# Load reference FASTA.
parts = []
with FASTA.open() as f:
    for line in f:
        if not line.startswith(">"):
            parts.append(line.strip().upper())
genome = "".join(parts)

# Load candidate CDS features.
cds = {}
with GFF.open() as f:
    for line in f:
        if line.startswith("#"):
            continue
        p = line.rstrip().split("\t")
        if len(p) != 9 or p[0] != "NC_000962.3" or p[2] != "CDS":
            continue
        a = attrs(p[8])
        locus = a.get("locus_tag", "")
        if locus not in {"Rv2024c", "Rv2756c", "Rv2761c", "Rv3263"}:
            continue
        start, end, strand = int(p[3]), int(p[4]), p[6]
        seq = genome[start-1:end]
        if strand == "-":
            seq = seq.translate(COMP)[::-1]
        cds[locus] = {
            "start": start, "end": end, "strand": strand,
            "sequence": seq, "phase": p[7]
        }

with INPUT.open(newline="", encoding="utf-8-sig") as f:
    rows = list(csv.DictReader(f))

results = []
for r in rows:
    locus = r.get("locus_tag", "")
    pos = int(r["position"])
    ref, alt = r["REF"].upper(), r["ALT"].upper()
    status = ""
    consequence = protein_pos = protein_change = ""

    if pos < 1 or pos > len(genome):
        status = "POSITION_OUTSIDE_REFERENCE"
    elif genome[pos-1] != ref:
        status = "REFERENCE_MISMATCH"
    elif locus not in cds:
        status = "NO_CDS_MATCH"
    elif len(ref) != 1 or len(alt) != 1 or alt not in "ACGT":
        status = "NON_SNV_REQUIRES_SEPARATE_REVIEW"
    else:
        feat = cds[locus]
        idx = pos - feat["start"] if feat["strand"] == "+" else feat["end"] - pos
        tx = feat["sequence"]
        expected_ref = ref if feat["strand"] == "+" else ref.translate(COMP)
        alt_tx = alt if feat["strand"] == "+" else alt.translate(COMP)

        if idx < 0 or idx >= len(tx):
            status = "POSITION_OUTSIDE_CDS"
        elif tx[idx] != expected_ref:
            status = "CDS_REFERENCE_MISMATCH"
        else:
            cs = (idx // 3) * 3
            codon = tx[cs:cs+3]
            if len(codon) != 3:
                status = "INCOMPLETE_CODON"
            else:
                altered = list(codon)
                altered[idx % 3] = alt_tx
                altered = "".join(altered)
                aa1, aa2 = translate(codon), translate(altered)
                protein_pos = idx // 3 + 1
                if aa1 == aa2:
                    consequence = "synonymous"
                    protein_change = f"{aa1}{protein_pos}="
                elif aa2 == "*":
                    consequence = "stop_gained"
                    protein_change = f"{aa1}{protein_pos}*"
                elif aa1 == "*":
                    consequence = "stop_lost"
                    protein_change = f"*{protein_pos}{aa2}"
                else:
                    consequence = "missense"
                    protein_change = f"{aa1}{protein_pos}{aa2}"
                status = "OK"

    results.append({
        **r,
        "reference_base_in_FASTA": genome[pos-1] if 1 <= pos <= len(genome) else "",
        "calculated_consequence": consequence,
        "calculated_protein_position": protein_pos,
        "calculated_protein_change": protein_change,
        "existing_consequence_matches": (
            str(r.get("Consequence", "")).lower() == consequence.lower()
            if status == "OK" else ""
        ),
        "existing_position_matches": (
            str(r.get("ProteinPosition", "")) == str(protein_pos)
            if status == "OK" else ""
        ),
        "existing_change_matches": (
            str(r.get("ProteinChange", "")) == protein_change
            if status == "OK" else ""
        ),
        "validation_status": status
    })

fields = list(dict.fromkeys(k for r in results for k in r))
with OUTPUT.open("w", newline="", encoding="utf-8") as f:
    w = csv.DictWriter(f, fieldnames=fields)
    w.writeheader()
    w.writerows(results)

counts = Counter(r["validation_status"] for r in results)
summary = [{"validation_status": k, "records": v} for k, v in counts.items()]
for field in ["existing_consequence_matches", "existing_position_matches", "existing_change_matches"]:
    vals = [r[field] for r in results if r["validation_status"] == "OK"]
    summary.append({"validation_status": field + "_true", "records": sum(v is True for v in vals)})
    summary.append({"validation_status": field + "_false", "records": sum(v is False for v in vals)})

with SUMMARY.open("w", newline="", encoding="utf-8") as f:
    w = csv.DictWriter(f, fieldnames=["validation_status", "records"])
    w.writeheader()
    w.writerows(summary)

print("Validated records:", len(results))
print("Validation statuses:", dict(counts))
print("Consequence summary:", dict(Counter(r["calculated_consequence"] for r in results if r["validation_status"] == "OK")))
print("Created:", OUTPUT)
print("Created:", SUMMARY)
print("NOTE: This checks SNVs in single CDS features; split CDS/phase cases need additional review.")
