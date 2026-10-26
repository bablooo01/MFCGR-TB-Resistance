from pathlib import Path
import shutil
import zipfile
import csv

ROOT = Path("/mnt/d/MFCGR_PROJECT")
RES = ROOT / "08_RESULTS"
FIG = RES / "figures"
OUT = ROOT / "09_FINAL_PACKAGE"
OUT.mkdir(parents=True, exist_ok=True)

tables = [
    "selected_120_isolates.csv",
    "methyltransferase_gene_variants_120.csv",
    "methyltransferase_variant_consequences_standardized.csv",
    "consequence_independent_validation.csv",
    "consequence_validation_summary.csv",
    "cohort_top_variant_validation.csv",
    "top_variant_fisher_tests.csv",
    "candidate_annotation_check.csv",
    "candidate_allele_group_summary.csv",
    "preliminary_candidate_gene_ranking.csv",
    "circos_candidate_data.csv",
    "candidate_gene_sequences.fasta",
    "candidate_sequence_summary.csv",
    "MFCGR_results_inventory.md",
]

figures = [
    "CGR_Rv2024c.png",
    "CGR_Rv2756c.png",
    "CGR_Rv2761c.png",
    "CGR_Rv3263.png",
    "candidate_loci_circular_overview.png",
]

for sub in ["tables", "figures"]:
    (OUT / sub).mkdir(exist_ok=True)

copied = []
missing = []

for name in tables:
    src = RES / name
    if src.exists():
        shutil.copy2(src, OUT / "tables" / name)
        copied.append(f"tables/{name}")
    else:
        missing.append(name)

for name in figures:
    src = FIG / name
    if src.exists():
        shutil.copy2(src, OUT / "figures" / name)
        copied.append(f"figures/{name}")
    else:
        missing.append(f"figures/{name}")

# Collect key counts for the summary.
cohort_path = RES / "selected_120_isolates.csv"
cohort_n = 120
if cohort_path.exists():
    with cohort_path.open(newline="", encoding="utf-8-sig") as f:
        cohort_n = sum(1 for _ in csv.DictReader(f))

report = """# MFCGR Project — Consolidated Preliminary Results

## Project title
MFCGR: Methylation-aware Fractal Chaos Game Representation for
Characterizing Drug-Resistant *Mycobacterium tuberculosis* Genomes

## Current computational progress
- A selected cohort of 120 isolates was prepared, with 30 isolates in each
  resistance category: Sensitive, MDR, Pre-XDR and XDR.
- Candidate-locus variants were extracted from the filtered VCF files.
- First-pass protein consequence predictions were produced and compared with
  a separate codon-level calculation.
- Two candidate SNV positions were assessed using read-level evidence across
  the cohort.
- Exploratory Fisher exact tests and FDR-adjusted q-values were calculated.
- Four reference-sequence CGR figures and one circular candidate-locus overview
  were generated.

## Candidate findings
The current exploratory analysis found associations between the tested alleles
at Rv2024c and Rv3263 and resistance categories in this cohort. These results
are hypothesis-generating and are not evidence of causation.

Rv2756c is annotated as a type I restriction/modification system DNA methylase
(HsdM). Rv2761c is annotated as the HsdS specificity determinant. Rv2024c is
annotated as a hypothetical protein in the reference annotation examined;
its function should not be assumed to be that of a methyltransferase.

## Important limitations
1. The candidate-gene ranking is preliminary and based on exploratory
   association statistics.
2. Lineage/population structure has not yet been controlled for in the
   association analysis.
3. The independent consequence check reproduced the consequence labels and
   protein positions, but full CDS phase, split-feature and annotation validation
   remains necessary. Matching summary counts alone do not establish correctness.
4. The read-level validation covered two candidate positions, not every variant
   in the cohort. The allele-calling parser and thresholds should be reviewed.
5. The CGR figures use reference CDS sequences, not isolate-specific sequences.
6. The circular overview is a preliminary Circos-style visualization, not a
   complete multi-track Circos analysis.
7. The current Illumina sequencing data do not directly measure DNA methylation.
   A true methylation-aware analysis requires compatible methylation measurements
   or a clearly defined and validated alternative.
8. The reported associations do not demonstrate that the variants cause drug
   resistance or change methylation.

## Recommended remaining work
- Validate gene annotations and variant consequences with a trusted annotation
  workflow that handles strand, CDS phase and split CDS features.
- Verify allele calls and coverage across the cohort.
- Reassess associations with lineage/population structure where possible.
- Verify candidate gene functions using authoritative sources.
- Obtain compatible methylation data or clearly narrow the methylation claims.
- Extend CGR to isolate-specific sequences if appropriate.
- Finalize the gene ranking only after these checks.
- Complete the methods, results, discussion and limitations sections of the report.

## Included files
See the `tables/` and `figures/` folders for the copied results.
"""

(OUT / "MFCGR_consolidated_summary.md").write_text(report, encoding="utf-8")

inventory = [
    "# Final package file inventory",
    "",
    f"Selected cohort rows reported by source file: {cohort_n}",
    f"Files copied: {len(copied)}",
    "",
    "## Included",
    *[f"- {x}" for x in copied],
    "",
    "## Missing source files",
    *([f"- {x}" for x in missing] if missing else ["- None"]),
    "",
    "The package contains preliminary results and explicitly documents remaining validation."
]
(OUT / "PACKAGE_INVENTORY.md").write_text("\n".join(inventory), encoding="utf-8")

zip_path = ROOT / "MFCGR_final_results_package.zip"
with zipfile.ZipFile(zip_path, "w", zipfile.ZIP_DEFLATED) as z:
    for p in OUT.rglob("*"):
        if p.is_file():
            z.write(p, p.relative_to(ROOT))

print("PACKAGE FOLDER:", OUT)
print("ZIP FILE:", zip_path)
print("Files copied:", len(copied))
print("Missing source files:", len(missing))
for x in missing:
    print("MISSING:", x)
print("Summary:", OUT / "MFCGR_consolidated_summary.md")
