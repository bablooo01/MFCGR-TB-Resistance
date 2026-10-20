import pandas as pd
from pathlib import Path

# ============================================================
# MFCGR PROJECT
# FINAL SELECTION OF 120 CRyPTIC ISOLATES
# 30 Sensitive + 30 MDR + 30 Pre-XDR + 30 XDR
# ============================================================

PROJECT_DIR = Path(r"D:\MFCGR_PROJECT")

INPUT_FILE = (
    PROJECT_DIR
    / "08_RESULTS"
    / "cryptic_classification_candidates.csv"
)

OUTPUT_FILE = (
    PROJECT_DIR
    / "08_RESULTS"
    / "selected_120_isolates.csv"
)

RANDOM_SEED = 42
N_PER_CATEGORY = 30

CATEGORIES = [
    "Sensitive",
    "MDR",
    "Pre-XDR",
    "XDR"
]

print("=" * 80)
print("MFCGR - FINAL 120 ISOLATE SELECTION")
print("=" * 80)

# ============================================================
# 1. LOAD CANDIDATE FILE
# ============================================================

print("\nLoading candidate file...")

df = pd.read_csv(INPUT_FILE)

print("Total candidate isolates:", len(df))

# ============================================================
# 2. BASIC CHECKS
# ============================================================

print("\nChecking UNIQUEID values...")

duplicate_ids = df["UNIQUEID"].duplicated().sum()

print("Duplicate UNIQUEIDs:", duplicate_ids)

if duplicate_ids > 0:
    print("\nWARNING: Duplicate UNIQUEIDs found.")
    df = df.drop_duplicates(subset=["UNIQUEID"])

print("Unique isolates:", df["UNIQUEID"].nunique())

# ============================================================
# 3. CHECK FASTQ INFORMATION
# ============================================================

print("\nChecking FASTQ availability...")

df["HAS_FASTQ"] = (
    df["fastq_ftp"].notna()
    &
    (df["fastq_ftp"].astype(str).str.strip() != "")
)

print(
    "Isolates with FASTQ:",
    df["HAS_FASTQ"].sum()
)

# ============================================================
# 4. CHECK PAIRED-END FASTQ
# ============================================================

def is_paired_fastq(value):

    if pd.isna(value):
        return False

    value = str(value).strip()

    if not value:
        return False

    # CRyPTIC/ENA may provide multiple FASTQ files
    # separated by semicolons.
    files = [
        x.strip()
        for x in value.split(";")
        if x.strip()
    ]

    return len(files) >= 2


df["PAIRED_FASTQ"] = df["fastq_ftp"].apply(
    is_paired_fastq
)

print(
    "Isolates with paired FASTQ:",
    df["PAIRED_FASTQ"].sum()
)

# ============================================================
# 5. KEEP ONLY SUITABLE ISOLATES
# ============================================================

eligible = df[
    df["HAS_FASTQ"]
    &
    df["PAIRED_FASTQ"]
    &
    df["CATEGORY"].isin(CATEGORIES)
].copy()

print("\nEligible isolates for selection:", len(eligible))

# ============================================================
# 6. CHECK CATEGORY COUNTS
# ============================================================

print("\n" + "=" * 80)
print("ELIGIBLE CATEGORY COUNTS")
print("=" * 80)

category_counts = (
    eligible["CATEGORY"]
    .value_counts()
    .reindex(CATEGORIES, fill_value=0)
)

print(category_counts)

# ============================================================
# 7. VERIFY THAT EACH CATEGORY HAS 30+
# ============================================================

for category in CATEGORIES:

    count = category_counts[category]

    if count < N_PER_CATEGORY:

        raise ValueError(
            f"\nERROR: Only {count} eligible isolates "
            f"available for {category}. "
            f"Need {N_PER_CATEGORY}."
        )

print("\nAll four categories have enough isolates.")

# ============================================================
# 8. RANDOMLY SELECT 30 FROM EACH CATEGORY
# ============================================================

print("\nSelecting isolates...")

selected_parts = []

for category in CATEGORIES:

    category_df = eligible[
        eligible["CATEGORY"] == category
    ].copy()

    selected = category_df.sample(
        n=N_PER_CATEGORY,
        random_state=RANDOM_SEED
    )

    selected_parts.append(selected)

# Combine all four groups
selected_120 = pd.concat(
    selected_parts,
    ignore_index=True
)

# ============================================================
# 9. SORT FOR EASY INSPECTION
# ============================================================

category_order = {
    "Sensitive": 1,
    "MDR": 2,
    "Pre-XDR": 3,
    "XDR": 4
}

selected_120["CATEGORY_ORDER"] = (
    selected_120["CATEGORY"]
    .map(category_order)
)

selected_120 = selected_120.sort_values(
    ["CATEGORY_ORDER", "UNIQUEID"]
).reset_index(drop=True)

selected_120 = selected_120.drop(
    columns=["CATEGORY_ORDER"]
)

# ============================================================
# 10. ADD SAMPLE NUMBER
# ============================================================

selected_120.insert(
    0,
    "SAMPLE_NUMBER",
    range(1, len(selected_120) + 1)
)

# ============================================================
# 11. FINAL VALIDATION
# ============================================================

print("\n" + "=" * 80)
print("FINAL SELECTION")
print("=" * 80)

print(
    selected_120["CATEGORY"].value_counts()
)

print("\nTotal selected:", len(selected_120))

# Check duplicate IDs
print(
    "Duplicate UNIQUEIDs:",
    selected_120["UNIQUEID"].duplicated().sum()
)

# Check FASTQ
print(
    "Missing FASTQ:",
    selected_120["fastq_ftp"].isna().sum()
)

# Check paired FASTQ
print(
    "Missing paired FASTQ:",
    (~selected_120["PAIRED_FASTQ"]).sum()
)

# ============================================================
# 12. FINAL SAFETY CHECK
# ============================================================

assert len(selected_120) == 120
assert selected_120["UNIQUEID"].nunique() == 120
assert selected_120["fastq_ftp"].notna().all()
assert selected_120["PAIRED_FASTQ"].all()

for category in CATEGORIES:

    assert (
        selected_120["CATEGORY"]
        .eq(category)
        .sum()
        == N_PER_CATEGORY
    )

# ============================================================
# 13. SAVE FINAL FILE
# ============================================================

OUTPUT_FILE.parent.mkdir(
    parents=True,
    exist_ok=True
)

selected_120.to_csv(
    OUTPUT_FILE,
    index=False
)

print("\n" + "=" * 80)
print("SUCCESS!")
print("=" * 80)

print("\nFinal dataset:")
print("Sensitive : 30")
print("MDR       : 30")
print("Pre-XDR   : 30")
print("XDR       : 30")
print("----------------")
print("TOTAL     : 120")

print("\nSaved to:")
print(OUTPUT_FILE)

print("\n" + "=" * 80)
print("DONE")
print("=" * 80)