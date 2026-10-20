import pandas as pd
from pathlib import Path

# ============================================================
# MFCGR - WHO AMR ANNOTATION
# Pilot isolate: ERR038300
# ============================================================

PROJECT = Path(r"D:\MFCGR_PROJECT")

WHO_FILE = PROJECT / "06_WHO_CATALOGUE" / "WHO-UCN-TB-2023.7-eng.xlsx"

MATCH_FILE = (
    PROJECT / "08_RESULTS" / "alignment" / "ERR038300.WHO_matches.txt"
)

OUTPUT_FILE = (
    PROJECT / "08_RESULTS" / "alignment" / "ERR038300_WHO_AMR_annotation.csv"
)


# ------------------------------------------------------------
# 1. Load WHO catalogue
# ------------------------------------------------------------

print("Loading WHO catalogue...")

who = pd.read_excel(
    WHO_FILE,
    sheet_name="Catalogue_master_file",
    header=2,
    usecols=[
        "drug",
        "gene",
        "mutation",
        "variant",
        "tier",
        "effect",
        "Initial confidence grading WHO dataset",
        "FINAL CONFIDENCE GRADING",
        "Comment",
    ],
)

print(f"WHO catalogue rows loaded: {len(who):,}")


# ------------------------------------------------------------
# 2. Load WHO coordinate matches
# ------------------------------------------------------------

print("Loading WHO coordinate matches...")

matches = pd.read_csv(
    MATCH_FILE,
    sep=r"\s+",
    header=None,
    names=[
        "chromosome",
        "position",
        "REF",
        "ALT",
        "WHO_raw",
    ],
)

matches["WHO_VARIANT"] = (
    matches["WHO_raw"]
    .str.replace("WHO_VARIANT=", "", regex=False)
)

print(f"WHO coordinate matches: {len(matches)}")


# ------------------------------------------------------------
# 3. Join with WHO catalogue
# ------------------------------------------------------------

print("Joining variants with WHO catalogue...")

out = matches.merge(
    who,
    left_on="WHO_VARIANT",
    right_on="variant",
    how="left",
)


# ------------------------------------------------------------
# 4. Conservative AMR classification
# ------------------------------------------------------------

def classify_amr(confidence):

    if pd.isna(confidence):
        return "No WHO classification"

    confidence = str(confidence)

    if "Not assoc" in confidence:
        return "Not associated with resistance"

    if "Uncertain" in confidence:
        return "Uncertain significance"

    return "WHO resistance-associated / requires review"


out["AMR_CLASSIFICATION"] = (
    out["FINAL CONFIDENCE GRADING"]
    .apply(classify_amr)
)


# ------------------------------------------------------------
# 5. Select final columns
# ------------------------------------------------------------

final_columns = [
    "chromosome",
    "position",
    "REF",
    "ALT",
    "WHO_VARIANT",
    "drug",
    "gene",
    "mutation",
    "tier",
    "effect",
    "Initial confidence grading WHO dataset",
    "FINAL CONFIDENCE GRADING",
    "AMR_CLASSIFICATION",
    "Comment",
]

result = out[final_columns]


# ------------------------------------------------------------
# 6. Save result
# ------------------------------------------------------------

result.to_csv(
    OUTPUT_FILE,
    index=False,
)

print()
print("=" * 60)
print("WHO AMR ANNOTATION COMPLETE")
print("=" * 60)

print(f"Rows in final annotation : {len(result)}")
print()
print("AMR classification:")
print(
    result["AMR_CLASSIFICATION"]
    .value_counts(dropna=False)
    .to_string()
)

print()
print(f"Output file:")
print(OUTPUT_FILE)