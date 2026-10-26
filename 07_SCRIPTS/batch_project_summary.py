import csv
from pathlib import Path
from collections import Counter, defaultdict

ROOT = Path("/mnt/d/MFCGR_PROJECT")
RES = ROOT / "08_RESULTS"
GFF = ROOT / "01_H37RV/ncbi_dataset/data/GCF_000195955.2/genomic.gff"

def read_csv(path):
    with path.open(newline="", encoding="utf-8-sig") as f:
        return list(csv.DictReader(f))

def write_csv(path, rows):
    if not rows:
        print("No rows to write:", path.name)
        return
    with path.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    print("Created:", path)

def attrs(text):
    out = {}
    for item in text.split(";"):
        if "=" in item:
            k, v = item.split("=", 1)
            out[k] = v
    return out

# 1. Read reference annotations.
features = []
if not GFF.exists():
    raise SystemExit(f"NCBI GFF not found: {GFF}")

with GFF.open(encoding="utf-8") as f:
    for line in f:
        if line.startswith("#"):
            continue
        p = line.rstrip().split("\t")
        if len(p) != 9 or p[0] != "NC_000962.3":
            continue
        a = attrs(p[8])
        if p[2] in {"gene", "CDS"}:
            features.append({
                "feature": p[2], "start": int(p[3]), "end": int(p[4]),
                "strand": p[6], "gene": a.get("gene", ""),
                "locus_tag": a.get("locus_tag", ""),
                "product": a.get("product", ""),
                "protein_id": a.get("protein_id", "")
            })

targets = [
    ("Rv2024c", 2269545, "G", "A"),
    ("Rv2756c", 3069167, "A", "G"),
    ("Rv2761c", 3072285, "A", "C"),
    ("Rv3263", 3643985, "A", "C"),
]

annotation_rows = []
for locus, pos, ref, alt in targets:
    matches = [
        x for x in features
        if x["start"] <= pos <= x["end"]
        and (x["locus_tag"] == locus or x["gene"] == locus)
    ]
    for x in matches:
        annotation_rows.append({
            "candidate_locus": locus, "tested_position": pos,
            "tested_REF": ref, "tested_ALT": alt,
            **x
        })
    if not matches:
        annotation_rows.append({
            "candidate_locus": locus, "tested_position": pos,
            "tested_REF": ref, "tested_ALT": alt,
            "feature": "NO_EXACT_MATCH", "start": "", "end": "",
            "strand": "", "gene": "", "locus_tag": "",
            "product": "", "protein_id": ""
        })

write_csv(RES / "candidate_annotation_check.csv", annotation_rows)

# 2. Summarize read-validated candidate alleles by resistance group.
calls_path = RES / "cohort_top_variant_validation.csv"
calls = read_csv(calls_path)
call_summary = []
for gene in sorted(set(r["gene"] for r in calls)):
    for group in ["Sensitive", "MDR", "Pre-XDR", "XDR"]:
        subset = [r for r in calls if r["gene"] == gene and r["CATEGORY"] == group]
        usable = [r for r in subset if r["call"] in {"ALT_SUPPORTED", "REF_SUPPORTED"}]
        alt = sum(r["call"] == "ALT_SUPPORTED" for r in usable)
        ref = sum(r["call"] == "REF_SUPPORTED" for r in usable)
        call_summary.append({
            "gene": gene, "category": group, "total_isolates": len(subset),
            "usable_calls": len(usable), "ALT_supported": alt,
            "REF_supported": ref,
            "ALT_fraction_among_usable": round(alt / len(usable), 4) if usable else ""
        })
write_csv(RES / "candidate_allele_group_summary.csv", call_summary)

# 3. Combine exploratory association statistics and consequence summaries.
tests = read_csv(RES / "top_variant_fisher_tests.csv")
consequence_counts = Counter()
consequence_path = RES / "methyltransferase_variant_consequences.csv"
if consequence_path.exists():
    for r in read_csv(consequence_path):
        gene = r.get("gene", "")
        consequence = (
            r.get("consequence") or r.get("predicted_consequence")
            or r.get("effect") or "UNSPECIFIED"
        )
        consequence_counts[(gene, consequence)] += 1

rank_rows = []
for gene in sorted(set(t["gene"] for t in tests)):
    gt = [t for t in tests if t["gene"] == gene]
    best = min(gt, key=lambda t: float(t["FDR_q_value"]))
    best_q = float(best["FDR_q_value"])
    best_or = float(best["odds_ratio"])
    significant = sum(float(t["FDR_q_value"]) < 0.05 for t in gt)
    allele_rows = [r for r in calls if r["gene"] == gene]
    sensitive = [r for r in allele_rows if r["CATEGORY"] == "Sensitive"]
    other = [r for r in allele_rows if r["CATEGORY"] != "Sensitive"]
    sens_alt = sum(r["call"] == "ALT_SUPPORTED" for r in sensitive)
    other_alt = sum(r["call"] == "ALT_SUPPORTED" for r in other)
    rank_rows.append({
        "gene": gene,
        "best_comparison": best["comparison"],
        "best_odds_ratio": best_or,
        "best_FDR_q_value": best_q,
        "comparisons_q_below_0.05": significant,
        "Sensitive_ALT_supported": sens_alt,
        "NonSensitive_ALT_supported": other_alt,
        "ranking_basis": "Exploratory association only; annotation and lineage validation pending"
    })

rank_rows.sort(key=lambda r: (r["comparisons_q_below_0.05"], -r["best_odds_ratio"], r["best_FDR_q_value"]), reverse=True)
for i, r in enumerate(rank_rows, 1):
    r["preliminary_rank"] = i
write_csv(RES / "preliminary_candidate_gene_ranking.csv", rank_rows)

# 4. Circos-ready tables: gene coordinates and allele group frequencies.
circos_rows = []
for locus, pos, ref, alt in targets:
    matches = [
        x for x in features
        if x["start"] <= pos <= x["end"]
        and (x["locus_tag"] == locus or x["gene"] == locus)
    ]
    f = matches[0] if matches else {}
    for group in ["Sensitive", "MDR", "Pre-XDR", "XDR"]:
        subset = [
            r for r in calls
            if r["gene"] == locus and r["CATEGORY"] == group
            and r["call"] in {"ALT_SUPPORTED", "REF_SUPPORTED"}
        ]
        alt_count = sum(r["call"] == "ALT_SUPPORTED" for r in subset)
        circos_rows.append({
            "contig": "NC_000962.3", "candidate_locus": locus,
            "gene_start": f.get("start", ""),
            "gene_end": f.get("end", ""),
            "strand": f.get("strand", ""),
            "tested_position": pos, "REF": ref, "ALT": alt,
            "category": group, "usable_isolates": len(subset),
            "ALT_supported": alt_count,
            "ALT_fraction": round(alt_count / len(subset), 4) if subset else ""
        })
write_csv(RES / "circos_candidate_data.csv", circos_rows)

print("\nBATCH SUMMARY")
print("Reference annotation check:", RES / "candidate_annotation_check.csv")
print("Allele group summary:", RES / "candidate_allele_group_summary.csv")
print("Preliminary ranking:", RES / "preliminary_candidate_gene_ranking.csv")
print("Circos input:", RES / "circos_candidate_data.csv")
print("\nIMPORTANT: The ranking is preliminary, not a final biological ranking.")
print("CGR/Circos figures should be generated after annotation and sequence validation.")
