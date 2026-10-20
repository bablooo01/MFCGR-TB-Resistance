import pandas as pd
from pathlib import Path

PROJECT_DIR = Path(r"D:\MFCGR_PROJECT")
CRYPTIC_DIR = PROJECT_DIR / "02_CRYPTIC"

dst_path = CRYPTIC_DIR / "DST_MEASUREMENTS.parquet"

print("=" * 70)
print("CRyPTIC DRUG / PHENOTYPE INSPECTION")
print("=" * 70)

# Load data
dst = pd.read_parquet(dst_path)

print("\nDATA SHAPE:")
print(dst.shape)

# ------------------------------------------------------------
# 1. DRUGS
# ------------------------------------------------------------

print("\n" + "=" * 70)
print("DRUGS PRESENT IN DATA")
print("=" * 70)

drugs = dst.index.get_level_values("DRUG")

print("\nNumber of unique drugs:", drugs.nunique())

print("\nDrug names:")
print(sorted(drugs.dropna().unique()))

# ------------------------------------------------------------
# 2. PHENOTYPE VALUES
# ------------------------------------------------------------

print("\n" + "=" * 70)
print("PHENOTYPE VALUES")
print("=" * 70)

print(dst["PHENOTYPE"].value_counts(dropna=False))

# ------------------------------------------------------------
# 3. DRUG × PHENOTYPE
# ------------------------------------------------------------

print("\n" + "=" * 70)
print("DRUG × PHENOTYPE COUNTS")
print("=" * 70)

drug_pheno = pd.crosstab(
    drugs,
    dst["PHENOTYPE"],
    dropna=False
)

print(drug_pheno)

# ------------------------------------------------------------
# 4. NUMBER OF DST RESULTS PER ISOLATE
# ------------------------------------------------------------

print("\n" + "=" * 70)
print("NUMBER OF DRUG RESULTS PER ISOLATE")
print("=" * 70)

uniqueids = dst.index.get_level_values("UNIQUEID")

drug_counts = uniqueids.value_counts()

print(drug_counts.describe())

print("\nDistribution of number of drug results:")
print(drug_counts.value_counts().sort_index())

# ------------------------------------------------------------
# 5. EXAMPLE ISOLATES
# ------------------------------------------------------------

print("\n" + "=" * 70)
print("EXAMPLE ISOLATES WITH DRUG RESULTS")
print("=" * 70)

for uid in uniqueids.unique()[:5]:

    print("\nIsolate:", uid)

    temp = dst.loc[uid]

    print(
        temp[["PHENOTYPE", "METHOD_MIC", "QUALITY"]]
        .to_string()
    )

print("\n" + "=" * 70)
print("INSPECTION COMPLETE")
print("=" * 70)