from pathlib import Path
import pandas as pd

# ============================================================
# MFCGR PROJECT
# Batch WHO AMR Annotation for 120 Isolates
# ============================================================

BASE = Path("/mnt/d/MFCGR_PROJECT")

MATCH_DIR = BASE / "08_RESULTS" / "alignment"
WHO_XLSX = BASE / "06_WHO_CATALOGUE" / "WHO-UCN-TB-2023.7-eng.xlsx"
OUTPUT = BASE / "08_RESULTS" / "WHO_AMR_all_120.csv"


print("=" * 70)
print("MFCGR - BATCH WHO AMR ANNOTATION")
print("=" * 70)


# ============================================================
# STEP 1: LOAD WHO CATALOGUE
# ============================================================

print("\n[1/4] Loading WHO catalogue...")

who = pd.read_excel(
    WHO_XLSX,
    sheet_name="Catalogue_master_file",
    header=2
)

# Clean column names
who.columns = [str(col).strip() for col in who.columns]

print(f"WHO catalogue rows loaded: {len(who)}")


# Check required column
if "variant" not in who.columns:
    raise ValueError(
        "ERROR: 'variant' column not found in WHO catalogue."
    )

# Clean variant column
who["variant"] = (
    who["variant"]
    .astype(str)
    .str.strip()
)


# ============================================================
# STEP 2: SELECT WHO INFORMATION WE NEED
# ============================================================

wanted_columns = [
    "drug",
    "gene",
    "mutation",
    "variant",
    "tier",
    "effect",
    "genomic position",
    "Initial confidence grading WHO dataset",
    "FINAL CONFIDENCE GRADING",
    "Comment"
]

available_columns = [
    col for col in wanted_columns
    if col in who.columns
]

who_lookup = who[available_columns].copy()

# Remove empty variant entries
who_lookup = who_lookup[
    (who_lookup["variant"] != "") &
    (who_lookup["variant"] != "nan")
]

print(
    f"WHO catalogue entries available for annotation: "
    f"{len(who_lookup)}"
)


# ============================================================
# STEP 3: FIND ALL WHO MATCH FILES
# ============================================================

print("\n[2/4] Finding WHO coordinate-match files...")

match_files = sorted(
    MATCH_DIR.glob("*.WHO_matches.txt")
)

print(
    f"WHO match files found: {len(match_files)}"
)

if len(match_files) == 0:
    raise FileNotFoundError(
        "No *.WHO_matches.txt files found."
    )


# ============================================================
# STEP 4: ANNOTATE ALL ISOLATES
# ============================================================

print("\n[3/4] Annotating isolates...\n")

all_results = []
failed = []

for index, match_file in enumerate(match_files, start=1):

    isolate = match_file.name.replace(
        ".WHO_matches.txt",
        ""
    )

    try:

        rows = []

        with open(
            match_file,
            "r",
            encoding="utf-8"
        ) as file:

            for line in file:

                line = line.strip()

                if not line:
                    continue

                parts = line.split("\t")

                if len(parts) < 5:
                    continue

                chromosome = parts[0]
                position = parts[1]
                ref = parts[2]
                alt = parts[3]
                who_field = parts[4]

                if "WHO_VARIANT=" not in who_field:
                    continue

                who_variant = (
                    who_field
                    .split("WHO_VARIANT=", 1)[1]
                    .strip()
                )

                rows.append({
                    "isolate": isolate,
                    "chromosome": chromosome,
                    "position": position,
                    "vcf_ref": ref,
                    "vcf_alt": alt,
                    "WHO_VARIANT": who_variant
                })

        if len(rows) == 0:

            print(
                f"[{index}/{len(match_files)}] "
                f"{isolate}: 0 matches"
            )

            continue


        # ----------------------------------------------------
        # Join our variants with WHO catalogue
        # ----------------------------------------------------

        matches_df = pd.DataFrame(rows)

        annotated = matches_df.merge(
            who_lookup,
            left_on="WHO_VARIANT",
            right_on="variant",
            how="left"
        )


        # ----------------------------------------------------
        # Conservative WHO interpretation
        # ----------------------------------------------------

        def classify(row):

            text = " ".join(
                str(row.get(column, ""))
                for column in [
                    "effect",
                    "Comment",
                    "FINAL CONFIDENCE GRADING"
                ]
            ).lower()

            if "not assoc" in text:
                return "Not associated with resistance"

            elif "uncertain" in text:
                return "Uncertain significance"

            else:
                return (
                    "WHO resistance-associated / "
                    "requires review"
                )


        annotated["AMR_interpretation"] = (
            annotated.apply(
                classify,
                axis=1
            )
        )


        all_results.append(annotated)

        print(
            f"[{index}/{len(match_files)}] "
            f"{isolate}: "
            f"{len(annotated)} annotation rows"
        )


    except Exception as error:

        failed.append(
            (isolate, str(error))
        )

        print(
            f"[{index}/{len(match_files)}] "
            f"{isolate}: FAILED"
        )

        print(f"    Reason: {error}")


# ============================================================
# STEP 5: COMBINE ALL RESULTS
# ============================================================

print("\n[4/4] Combining annotation results...")

if len(all_results) == 0:

    raise RuntimeError(
        "No WHO annotation results were generated."
    )


final_df = pd.concat(
    all_results,
    ignore_index=True
)


# Remove exact duplicate rows
final_df = final_df.drop_duplicates()


# Save master CSV
final_df.to_csv(
    OUTPUT,
    index=False
)


# ============================================================
# FINAL SUMMARY
# ============================================================

print("\n")
print("=" * 70)
print("WHO AMR ANNOTATION COMPLETE")
print("=" * 70)

print(
    f"WHO match files found : {len(match_files)}"
)

print(
    f"Isolates annotated    : "
    f"{final_df['isolate'].nunique()}"
)

print(
    f"Total annotation rows : "
    f"{len(final_df)}"
)

print(
    f"Failed isolates       : "
    f"{len(failed)}"
)


# ------------------------------------------------------------
# AMR interpretation summary
# ------------------------------------------------------------

print("\nAMR interpretation counts:")

print(
    final_df[
        "AMR_interpretation"
    ].value_counts()
)


# ------------------------------------------------------------
# Drug summary
# ------------------------------------------------------------

if "drug" in final_df.columns:

    print("\nTop drugs represented:")

    print(
        final_df["drug"]
        .value_counts()
        .head(15)
    )


# ------------------------------------------------------------
# Gene summary
# ------------------------------------------------------------

if "gene" in final_df.columns:

    print("\nTop genes represented:")

    print(
        final_df["gene"]
        .value_counts()
        .head(15)
    )


# ------------------------------------------------------------
# Failed isolates
# ------------------------------------------------------------

if failed:

    print("\nFailed isolates:")

    for isolate, error in failed:
        print(
            f"  {isolate}: {error}"
        )


print("\nOutput file:")
print(OUTPUT)

print("=" * 70)