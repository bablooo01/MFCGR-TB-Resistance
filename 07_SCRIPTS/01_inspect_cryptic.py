import pandas as pd
from pathlib import Path

# ============================================================
# MFCGR PROJECT - CRyPTIC DATA INSPECTION
# ============================================================

# Project folders
PROJECT_DIR = Path(r"D:\MFCGR_PROJECT")
CRYPTIC_DIR = PROJECT_DIR / "02_CRYPTIC"

# Files
files = {
    "DST_MEASUREMENTS": CRYPTIC_DIR / "DST_MEASUREMENTS.parquet",
    "DST_SAMPLES": CRYPTIC_DIR / "DST_SAMPLES.parquet",
    "GENOMES": CRYPTIC_DIR / "GENOMES.parquet",
    "WGS_SAMPLES": CRYPTIC_DIR / "WGS_SAMPLES.parquet",
}

print("=" * 70)
print("MFCGR PROJECT - CRyPTIC DATA INSPECTION")
print("=" * 70)

for name, path in files.items():

    print("\n" + "=" * 70)
    print(f"{name}")
    print("=" * 70)

    if not path.exists():
        print(f"ERROR: File not found:")
        print(path)
        continue

    try:
        df = pd.read_parquet(path)

        print(f"File: {path.name}")
        print(f"Rows: {df.shape[0]:,}")
        print(f"Columns: {df.shape[1]}")

        print("\nCOLUMN NAMES:")
        for i, col in enumerate(df.columns, start=1):
            print(f"{i:3}. {col}")

        print("\nDATA TYPES:")
        print(df.dtypes)

        print("\nFIRST 5 ROWS:")
        print(df.head().to_string())

        print("\nMISSING VALUES:")
        print(df.isna().sum())

    except Exception as e:
        print(f"ERROR while reading {name}:")
        print(e)

print("\n" + "=" * 70)
print("INSPECTION COMPLETE")
print("=" * 70)