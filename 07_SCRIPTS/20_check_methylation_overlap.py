import os
import re
import glob
import pandas as pd

PROJECT = "/mnt/d/MFCGR_PROJECT"

SELECTED_FILE = "/mnt/d/MFCGR_PROJECT/08_RESULTS/selected_120_isolates.csv"

MODLIN_DIR = os.path.join(PROJECT, "04_MODLIN_METHYLATION")
PHELAN_DIR = os.path.join(PROJECT, "05_PHELAN_METHYLATION")
RESULTS_DIR = os.path.join(PROJECT, "08_RESULTS")

os.makedirs(RESULTS_DIR, exist_ok=True)


print("=" * 75)
print("MFCGR - CRYPTIC / METHYLOME METADATA MATCHING")
print("=" * 75)


# ============================================================
# NORMALIZATION
# ============================================================
def norm(x):

    if pd.isna(x):
        return ""

    x = str(x).strip().upper()

    x = re.sub(r"\s+", "", x)

    return x


# ============================================================
# LOAD CRYPTIC
# ============================================================
print("\n[1/4] Loading CRyPTIC selected isolates...")

selected = pd.read_csv(
    SELECTED_FILE,
    low_memory=False
)

print("CRyPTIC isolates:", len(selected))


important_cryptic_columns = [
    "SAMPLE_NUMBER",
    "UNIQUEID",
    "LINEAGE",
    "SUBLINEAGE",
    "run_accession",
    "sample_accession",
    "study_accession"
]

available_cryptic_columns = [
    c for c in important_cryptic_columns
    if c in selected.columns
]

print("\nCRyPTIC metadata columns used:")

for c in available_cryptic_columns:
    print("  -", c)


# ============================================================
# CREATE NORMALIZED CRYPTIC TABLE
# ============================================================
cryptic = selected[
    available_cryptic_columns
].copy()

for c in available_cryptic_columns:

    cryptic[c + "_NORM"] = cryptic[c].map(norm)


# ============================================================
# FIND METHYLATION FILES
# ============================================================
print("\n[2/4] Finding methylation files...")

methylation_files = []

for directory, source in [
    (MODLIN_DIR, "Modlin"),
    (PHELAN_DIR, "Phelan")
]:

    if not os.path.exists(directory):
        continue

    for ext in ["*.csv", "*.tsv", "*.txt"]:

        for f in glob.glob(
            os.path.join(directory, "**", ext),
            recursive=True
        ):

            methylation_files.append(
                (f, source)
            )


print(
    "Methylation files found:",
    len(methylation_files)
)


# ============================================================
# MATCHING
# ============================================================
print("\n[3/4] Comparing metadata...")

all_matches = []
lineage_matches = []


for filepath, source in methylation_files:

    filename = os.path.basename(filepath)

    print("\n" + "-" * 75)
    print(source, "|", filename)
    print("-" * 75)

    try:

        if filepath.lower().endswith(".tsv"):

            df = pd.read_csv(
                filepath,
                sep="\t",
                low_memory=False
            )

        else:

            df = pd.read_csv(
                filepath,
                low_memory=False
            )

    except Exception as e:

        print("Could not read:", e)
        continue


    print("Rows:", len(df))
    print("Columns:", list(df.columns))


    # --------------------------------------------------------
    # FIND POSSIBLE ISOLATE / STRAIN COLUMN
    # --------------------------------------------------------
    id_columns = []

    for col in df.columns:

        name = str(col).lower()

        if any(
            key in name
            for key in [
                "isolate",
                "strain",
                "sample",
                "accession",
                "id",
                "name"
            ]
        ):

            id_columns.append(col)


    # --------------------------------------------------------
    # EXACT VALUE MATCHES
    # --------------------------------------------------------
    for mcol in id_columns:

        methyl_values = (
            df[mcol]
            .dropna()
            .map(norm)
        )

        methyl_values = set(
            x for x in methyl_values
            if x
        )

        for ccol in available_cryptic_columns:

            cryptic_values = set(
                cryptic[ccol + "_NORM"]
                .dropna()
            )

            cryptic_values.discard("")

            overlap = (
                methyl_values &
                cryptic_values
            )

            if overlap:

                print(
                    "\nEXACT VALUE MATCH",
                    source,
                    filename,
                    mcol,
                    "<->",
                    ccol,
                    ":",
                    len(overlap)
                )

                for value in sorted(overlap):

                    cryptic_rows = cryptic[
                        cryptic[ccol + "_NORM"] == value
                    ]

                    for _, row in cryptic_rows.iterrows():

                        all_matches.append({
                            "source": source,
                            "file": filename,
                            "methylation_column": mcol,
                            "matched_value": value,
                            "cryptic_column": ccol,
                            "cryptic_run": row.get(
                                "run_accession",
                                ""
                            ),
                            "cryptic_uniqueid": row.get(
                                "UNIQUEID",
                                ""
                            ),
                            "cryptic_sample": row.get(
                                "SAMPLE_NUMBER",
                                ""
                            )
                        })


    # --------------------------------------------------------
    # LINEAGE COMPARISON
    # --------------------------------------------------------
    if "LINEAGE" in df.columns:

        methyl_lineages = set(
            df["LINEAGE"]
            .dropna()
            .map(norm)
        )

        cryptic_lineages = set(
            cryptic["LINEAGE_NORM"]
            .dropna()
        )

        common_lineages = (
            methyl_lineages &
            cryptic_lineages
        )

        if common_lineages:

            print(
                "\nCOMMON LINEAGES:",
                ", ".join(sorted(common_lineages))
            )

            for lineage in sorted(common_lineages):

                cryptic_subset = cryptic[
                    cryptic["LINEAGE_NORM"] == lineage
                ]

                methyl_count = (
                    df["LINEAGE"]
                    .map(norm)
                    .eq(lineage)
                    .sum()
                )

                lineage_matches.append({
                    "source": source,
                    "file": filename,
                    "lineage": lineage,
                    "cryptic_isolates": len(
                        cryptic_subset
                    ),
                    "methylation_rows": methyl_count
                })


# ============================================================
# SAVE RESULTS
# ============================================================
print("\n[4/4] Saving results...")

exact_file = os.path.join(
    RESULTS_DIR,
    "cryptic_methylation_exact_metadata_matches.csv"
)

lineage_file = os.path.join(
    RESULTS_DIR,
    "cryptic_methylation_lineage_matches.csv"
)


pd.DataFrame(all_matches).to_csv(
    exact_file,
    index=False
)

pd.DataFrame(lineage_matches).to_csv(
    lineage_file,
    index=False
)


# ============================================================
# FINAL SUMMARY
# ============================================================
print("\n" + "=" * 75)
print("FINAL MATCHING SUMMARY")
print("=" * 75)

print(
    "CRyPTIC isolates:",
    len(selected)
)

print(
    "Exact metadata matches:",
    len(all_matches)
)

print(
    "Lineage-level overlaps:",
    len(lineage_matches)
)


if all_matches:

    print("\n*** EXACT MATCHES FOUND ***")

    for m in all_matches:

        print(
            m["source"],
            "|",
            m["file"],
            "|",
            m["methylation_column"],
            "=",
            m["matched_value"],
            "| CRyPTIC",
            m["cryptic_run"]
        )

else:

    print(
        "\nNo exact isolate/sample/accession matches found."
    )


if lineage_matches:

    print("\n*** COMMON LINEAGES ***")

    for m in lineage_matches:

        print(
            m["source"],
            "|",
            m["file"],
            "| lineage =",
            m["lineage"],
            "| CRyPTIC isolates =",
            m["cryptic_isolates"],
            "| methylation rows =",
            m["methylation_rows"]
        )

else:

    print(
        "\nNo lineage overlap found."
    )


print("\nReports saved:")

print(exact_file)
print(lineage_file)

print("\n" + "=" * 75)
print("COMPLETE")
print("=" * 75)