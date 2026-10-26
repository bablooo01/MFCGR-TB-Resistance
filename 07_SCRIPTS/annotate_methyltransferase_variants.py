
import gzip
import pandas as pd
from pathlib import Path
from Bio import SeqIO
from Bio.Seq import Seq

# ============================================================
# PATHS
# ============================================================
ROOT = Path(r"D:\MFCGR_PROJECT")

FASTA = ROOT / "01_H37RV" / "ncbi_dataset" / "data" / "GCF_000195955.2" / "GCF_000195955.2_ASM19595v2_genomic.fna"
GFF = ROOT / "01_H37RV" / "ncbi_dataset" / "data" / "GCF_000195955.2" / "genomic.gff"
INPUT = ROOT / "08_RESULTS" / "methyltransferase_gene_variants_120.csv"
OUT = ROOT / "08_RESULTS" / "methyltransferase_variant_consequences.csv"
SUMMARY = ROOT / "08_RESULTS" / "methyltransferase_consequence_summary.csv"

# ============================================================
# LOAD REFERENCE
# ============================================================
print("Loading reference genome...")
records = list(SeqIO.parse(str(FASTA), "fasta"))

if len(records) != 1:
    raise ValueError(f"Expected one reference sequence; found {len(records)}")

reference_record = records[0]
reference = str(reference_record.seq).upper()
print("Reference:", reference_record.id)
print("Length:", len(reference))

# ============================================================
# PARSE GFF CDS FEATURES
# ============================================================
print("Reading CDS annotations...")
cds_features = []

with open(GFF, "r", encoding="utf-8") as handle:
    for line in handle:
        if line.startswith("#"):
            continue

        fields = line.rstrip("\n").split("\t")
        if len(fields) != 9 or fields[2] != "CDS":
            continue

        contig, _, _, start, end, _, strand, phase, attrs = fields

        attributes = {}
        for item in attrs.split(";"):
            if "=" in item:
                key, value = item.split("=", 1)
                attributes[key] = value

        gene_name = attributes.get("gene", "")
        locus_tag = attributes.get("locus_tag", "")
        product = attributes.get("product", "")

        if not locus_tag:
            continue

        cds_features.append({
            "contig": contig,
            "start": int(start),
            "end": int(end),
            "strand": strand,
            "phase": phase,
            "locus_tag": locus_tag,
            "gene": gene_name,
            "product": product
        })

print("CDS features loaded:", len(cds_features))

# ============================================================
# LOAD VARIANT CSV
# ============================================================
print("Loading extracted variants...")
df = pd.read_csv(INPUT)

# Standardize column names from the actual input CSV
df.columns = df.columns.str.strip()

df["gene"] = df["gene"].astype(str).str.strip()
df["locus_tag"] = df["locus_tag"].astype(str).str.strip()
df["position"] = pd.to_numeric(
    df["position"], errors="raise"
).astype(int)

required = {"gene", "locus_tag", "position", "REF", "ALT"}
missing = required - set(df.columns)

if missing:
    raise ValueError(
        f"Missing columns: {sorted(missing)}. "
        f"Available columns: {list(df.columns)}"
    )

# Candidate locus tags and their current working labels.
# Labels are hypotheses, not proof of methyltransferase function.
targets = {
    "Rv2024c": "mamB_candidate",
    "Rv2756c": "hsdM",
    "Rv2761c": "hsdS",
    "Rv3263": "mamA_candidate"
}

# Match GFF features using locus tag, gene name, or the input label.
def feature_matches(row, feature):
    input_gene = str(row["locus_tag"]).strip()
    return (
        input_gene == feature["locus_tag"]
        or input_gene == feature["gene"]
        or targets.get(feature["locus_tag"], "") == input_gene
    )

# ============================================================
# ANNOTATE SINGLE-NUCLEOTIDE VARIANTS
# ============================================================
def consequence_for_variant(row):
    pos = int(row["position"])
    ref = str(row["REF"]).upper()
    alt = str(row["ALT"]).upper()

    result = {
        "LocusTag": "",
        "GFF_gene": "",
        "Product": "",
        "Strand": "",
        "Consequence": "unannotated",
        "ProteinPosition": "",
        "ProteinChange": "",
        "AnnotationNote": ""
    }

    if pos < 1 or pos > len(reference):
        result["AnnotationNote"] = "Position outside reference"
        return result

    # This first-pass script intentionally handles only single-base
    # substitutions. Indels and complex alleles need separate handling.
    if len(ref) != 1 or len(alt) != 1 or alt not in "ACGT" or ref not in "ACGT":
        result["Consequence"] = "complex_or_indel"
        result["AnnotationNote"] = "Requires separate indel-aware annotation"
        return result

    actual_ref = reference[pos - 1]
    if actual_ref != ref:
        result["Consequence"] = "REF_mismatch"
        result["AnnotationNote"] = (
            f"VCF REF={ref}, reference base={actual_ref}"
        )
        return result

    matching = [
        feature for feature in cds_features
        if feature["contig"] == reference_record.id
        and feature["start"] <= pos <= feature["end"]
        and feature_matches(row, feature)
    ]

    if not matching:
        result["AnnotationNote"] = "No matching CDS feature in GFF"
        return result

    feature = matching[0]
    result["LocusTag"] = feature["locus_tag"]
    result["GFF_gene"] = feature["gene"]
    result["Product"] = feature["product"]
    result["Strand"] = feature["strand"]

    # Build the spliced CDS in coding orientation, including joined CDS
    # segments where the GFF provides multiple CDS records for a locus.
    segments = [
        item for item in cds_features
        if item["contig"] == feature["contig"]
        and item["locus_tag"] == feature["locus_tag"]
    ]
    segments.sort(key=lambda item: item["start"], reverse=(feature["strand"] == "-"))

    coding_seq = ""
    genomic_positions = []

    for segment in segments:
        positions = list(range(segment["start"], segment["end"] + 1))
        segment_seq = reference[segment["start"] - 1:segment["end"]]

        if feature["strand"] == "-":
            positions.reverse()
            segment_seq = str(Seq(segment_seq).reverse_complement())

        coding_seq += segment_seq
        genomic_positions.extend(positions)

    if pos not in genomic_positions:
        result["AnnotationNote"] = "Variant not found in constructed CDS"
        return result

    offset = genomic_positions.index(pos)
    coding_ref = coding_seq[offset]
    coding_alt = str(Seq(alt).complement()) if feature["strand"] == "-" else alt

    if coding_ref != coding_seq[offset] or coding_ref not in "ACGT":
        result["Consequence"] = "coding_base_check_failed"
        return result

    # GFF phase is needed to interpret partial CDS segments correctly.
    # Avoid asserting a protein consequence when phase is nonzero.
    if any(segment["phase"] not in ("0", ".") for segment in segments):
        result["Consequence"] = "CDS_phase_review_needed"
        result["AnnotationNote"] = "Nonzero CDS phase; verify transcript/CDS construction"
        return result

    codon_start = (offset // 3) * 3
    codon = coding_seq[codon_start:codon_start + 3]

    if len(codon) != 3 or any(base not in "ACGT" for base in codon):
        result["Consequence"] = "incomplete_codon"
        return result

    changed_codon_list = list(codon)
    changed_codon_list[offset % 3] = coding_alt
    changed_codon = "".join(changed_codon_list)

    aa_ref = str(Seq(codon).translate())
    aa_alt = str(Seq(changed_codon).translate())

    result["ProteinPosition"] = (offset // 3) + 1
    result["ProteinChange"] = f"{aa_ref}{(offset // 3) + 1}{aa_alt}"

    if aa_ref == aa_alt:
        result["Consequence"] = "synonymous"
    elif aa_alt == "*":
        result["Consequence"] = "stop_gained"
    elif aa_ref == "*":
        result["Consequence"] = "stop_lost"
    else:
        result["Consequence"] = "missense"

    return result

print("Annotating variants...")
annotations = df.apply(consequence_for_variant, axis=1, result_type="expand")
result = pd.concat([df.reset_index(drop=True), annotations], axis=1)

result.to_csv(OUT, index=False)

summary = (
    result.groupby(["gene", "Consequence"], dropna=False)
    .size()
    .reset_index(name="VariantRecords")
)
summary.to_csv(SUMMARY, index=False)

print("\nConsequence counts:")
print(result["Consequence"].value_counts(dropna=False).to_string())

print("\nGene-by-consequence summary:")
print(summary.to_string(index=False))

print("\nSaved:")
print(OUT)
print(SUMMARY)
print("\nIMPORTANT: Review REF mismatches, indels, CDS phase, and gene identity before ranking.")


