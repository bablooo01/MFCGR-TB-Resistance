import pandas as pd
from pathlib import Path

# ============================================================
# MFCGR - VALIDATE FINAL 120 ISOLATES
# ============================================================

PROJECT_DIR = Path(r"D:\MFCGR_PROJECT")

INPUT_FILE = (
    PROJECT_DIR
    / "08_RESULTS"
    / "selected_120_isolates.csv"
)

print("=" * 80)
print("MFCGR - FINAL 120 ISOLATE VALIDATION")
print("=" * 80)

# ============================================================
# 1. LOAD
# ============================================================

df = pd.read_csv(
    INPUT_FILE,
    low_memory=False
)

print("\nTotal isolates:", len(df))

# ============================================================
# 2. BASIC VALIDATION
# ============================================================

print("\n" + "=" * 80)
print("BASIC CHECKS")
print("=" * 80)

print("Unique UNIQUEIDs:", df["UNIQUEID"].nunique())
print("Duplicate UNIQUEIDs:", df["UNIQUEID"].duplicated().sum())

print("\nCategory counts:")
print(df["CATEGORY"].value_counts())

# ============================================================
# 3. FASTQ CHECK
# ============================================================

print("\n" + "=" * 80)
print("FASTQ CHECK")
print("=" * 80)

def check_paired_fastq(value):

    if pd.isna(value):
        return False

    files = [
        x.strip()
        for x in str(value).split(";")
        if x.strip()
    ]

    return (
        len(files) >= 2
        and any("_1.fastq" in x for x in files)
        and any("_2.fastq" in x for x in files)
    )


df["PAIRED_FASTQ_CHECK"] = df["fastq_ftp"].apply(
    check_paired_fastq
)

print(
    "Paired FASTQ:",
    df["PAIRED_FASTQ_CHECK"].sum(),
    "/",
    len(df)
)

# ============================================================
# 4. WGS / SPECIES CHECK
# ============================================================

print("\n" + "=" * 80)
print("WGS / SPECIES CHECK")
print("=" * 80)

print(
    "M. tuberculosis:",
    (
        df["SPECIES"]
        .astype(str)
        .str.lower()
        .eq("m. tuberculosis")
    ).sum()
)

print(
    "has_genome = True:",
    (df["has_genome"] == True).sum()
)

print(
    "status values:"
)

print(
    df["status"].value_counts(dropna=False)
)

# ============================================================
# 5. COVERAGE / DEPTH
# ============================================================

print("\n" + "=" * 80)
print("SEQUENCING QUALITY INFORMATION")
print("=" * 80)

print("\nTB coverage:")
print(df["TB_COVERAGE"].describe())

print("\nTB depth:")
print(df["TB_DEPTH"].describe())

# ============================================================
# 6. LINEAGE
# ============================================================

print("\n" + "=" * 80)
print("LINEAGE DISTRIBUTION")
print("=" * 80)

print(
    df["LINEAGE"].value_counts(dropna=False)
)

print("\nSublineage distribution:")
print(
    df["SUBLINEAGE"].value_counts(dropna=False).head(20)
)

# ============================================================
# 7. RESISTANCE-RELATED INFORMATION
# ============================================================

print("\n" + "=" * 80)
print("RESISTANCE PROFILE")
print("=" * 80)

drugs = [
    "RIF",
    "INH",
    "LEV",
    "MXF",
    "BDQ",
    "LZD"
]

available_drugs = [
    drug
    for drug in drugs
    if drug in df.columns
]

print("\nAvailable resistance columns:")
print(available_drugs)

for drug in available_drugs:

    print(f"\n{drug}:")
    print(
        df[drug]
        .value_counts(dropna=False)
    )

# ============================================================
# 8. SHOW SELECTED ISOLATES
# ============================================================

print("\n" + "=" * 80)
print("SELECTED 120 ISOLATES")
print("=" * 80)

columns_to_show = [
    "SAMPLE_NUMBER",
    "UNIQUEID",
    "CATEGORY",
    "run_accession",
    "LINEAGE",
    "SUBLINEAGE",
    "TB_COVERAGE",
    "TB_DEPTH"
]

columns_to_show = [
    col for col in columns_to_show
    if col in df.columns
]

print(
    df[columns_to_show].to_string(index=False)
)

# ============================================================
# 9. FINAL VALIDATION
# ============================================================

print("\n" + "=" * 80)
print("FINAL VALIDATION")
print("=" * 80)

checks = {
    "Exactly 120 isolates":
        len(df) == 120,

    "120 unique UNIQUEIDs":
        df["UNIQUEID"].nunique() == 120,

    "30 Sensitive":
        (df["CATEGORY"] == "Sensitive").sum() == 30,

    "30 MDR":
        (df["CATEGORY"] == "MDR").sum() == 30,

    "30 Pre-XDR":
        (df["CATEGORY"] == "Pre-XDR").sum() == 30,

    "30 XDR":
        (df["CATEGORY"] == "XDR").sum() == 30,

    "All paired FASTQ":
        df["PAIRED_FASTQ_CHECK"].all(),

    "All M. tuberculosis":
        df["SPECIES"]
        .astype(str)
        .str.lower()
        .eq("m. tuberculosis")
        .all(),

    "All has_genome=True":
        (df["has_genome"] == True).all()
}

all_passed = True

for check, result in checks.items():

    status = "PASS" if result else "FAIL"

    print(
        f"{status:6} : {check}"
    )

    if not result:
        all_passed = False

# ============================================================
# 10. RESULT
# ============================================================

print("\n" + "=" * 80)

if all_passed:
    print("ALL VALIDATION CHECKS PASSED")
else:
    print("SOME VALIDATION CHECKS FAILED")

print("=" * 80)