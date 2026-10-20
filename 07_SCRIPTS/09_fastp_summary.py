# ============================================================
# 
# MFCGR PROJECT
# fastp Summary + Statistics + Plots
# ============================================================

import json
import re
from pathlib import Path

import pandas as pd

import matplotlib
matplotlib.use("Agg")

import matplotlib.pyplot as plt


# ============================================================
# PATHS
# ============================================================

PROJECT_DIR = Path(r"D:\MFCGR_PROJECT")

REPORT_DIR = PROJECT_DIR / "08_RESULTS" / "fastp_reports"
RESULT_DIR = PROJECT_DIR / "08_RESULTS"

CSV_OUTPUT = RESULT_DIR / "fastp_summary.csv"
ERROR_OUTPUT = RESULT_DIR / "fastp_summary_errors.txt"


# ============================================================
# SETTINGS
# ============================================================

EXPECTED_ISOLATES = 120


# ============================================================
# HEADER
# ============================================================

print("=" * 60)
print("FASTP SUMMARY")
print("=" * 60)


# ============================================================
# CHECK REPORT DIRECTORY
# ============================================================

if not REPORT_DIR.exists():
    raise FileNotFoundError(
        f"fastp report directory not found:\n{REPORT_DIR}"
    )


json_files = sorted(REPORT_DIR.glob("*.json"))

print(f"JSON reports found: {len(json_files)}")

if len(json_files) != EXPECTED_ISOLATES:
    print(
        f"WARNING: Expected {EXPECTED_ISOLATES} JSON reports "
        f"but found {len(json_files)}."
    )


# ============================================================
# JSON REPAIR FUNCTION
# ============================================================

def repair_invalid_backslashes(text):
    """
    fastp JSON files generated on Windows can contain Windows-style
    backslashes such as D:\MFCGR_PROJECT\...

    Python's json parser expects valid JSON escape sequences.

    This function escapes only backslashes that are NOT already
    valid JSON escapes.
    """

    pattern = r'\\(?!["\\/bfnrtu])'

    return re.sub(pattern, r'\\\\', text)


# ============================================================
# SAFE JSON LOADER
# ============================================================

def load_fastp_json(json_file):

    with open(
        json_file,
        "r",
        encoding="utf-8",
        errors="replace"
    ) as f:

        text = f.read()

    try:
        return json.loads(text)

    except json.JSONDecodeError:

        repaired_text = repair_invalid_backslashes(text)

        return json.loads(repaired_text)


# ============================================================
# HELPER FUNCTION
# ============================================================

def get_value(dictionary, key, default=0):

    value = dictionary.get(key, default)

    if value is None:
        return default

    return value


# ============================================================
# PARSE REPORTS
# ============================================================

records = []
errors = []

for json_file in json_files:

    isolate = json_file.stem

    try:

        data = load_fastp_json(json_file)

        # ----------------------------------------------------
        # SUMMARY
        # ----------------------------------------------------

        summary = data.get("summary", {})

        before = summary.get("before_filtering", {})
        after = summary.get("after_filtering", {})

        # ----------------------------------------------------
        # FILTERING RESULT
        # ----------------------------------------------------

        filtering = data.get("filtering_result", {})

        # ----------------------------------------------------
        # ADAPTER CUTTING
        # ----------------------------------------------------
        # IMPORTANT:
        # fastp stores adapter statistics under:
        #
        # adapter_cutting
        #
        # Example:
        # "adapter_trimmed_reads": 5592
        # "adapter_trimmed_bases": 109766
        # ----------------------------------------------------

        adapter = data.get("adapter_cutting", {})

        adapter_trimmed_reads = get_value(
            adapter,
            "adapter_trimmed_reads",
            0
        )

        adapter_trimmed_bases = get_value(
            adapter,
            "adapter_trimmed_bases",
            0
        )

        # ----------------------------------------------------
        # DUPLICATION
        # ----------------------------------------------------

        duplication = data.get("duplication", {})

        duplication_rate = get_value(
            duplication,
            "rate",
            0
        )

        # fastp JSON stores duplication rate as a fraction.
        # Convert to percentage.
        duplication_percent = float(duplication_rate) * 100

        # ----------------------------------------------------
        # INSERT SIZE
        # ----------------------------------------------------

        insert_size = data.get("insert_size", {})

        insert_peak = get_value(
            insert_size,
            "peak",
            0
        )

        # ----------------------------------------------------
        # READ COUNTS
        # ----------------------------------------------------

        before_reads = get_value(
            before,
            "total_reads",
            0
        )

        after_reads = get_value(
            after,
            "total_reads",
            0
        )

        # ----------------------------------------------------
        # READ RETENTION
        # ----------------------------------------------------

        if float(before_reads) > 0:

            read_retention = (
                float(after_reads)
                / float(before_reads)
            ) * 100

        else:

            read_retention = 0

        # ----------------------------------------------------
        # QUALITY
        # ----------------------------------------------------

        q20_before = (
            float(get_value(before, "q20_rate", 0)) * 100
        )

        q20_after = (
            float(get_value(after, "q20_rate", 0)) * 100
        )

        q30_before = (
            float(get_value(before, "q30_rate", 0)) * 100
        )

        q30_after = (
            float(get_value(after, "q30_rate", 0)) * 100
        )

        # ----------------------------------------------------
        # GC CONTENT
        # ----------------------------------------------------

        gc_before = (
            float(get_value(before, "gc_content", 0)) * 100
        )

        gc_after = (
            float(get_value(after, "gc_content", 0)) * 100
        )

        # ----------------------------------------------------
        # FAILED READS
        # ----------------------------------------------------

        low_quality_reads = get_value(
            filtering,
            "low_quality_reads",
            0
        )

        too_many_n_reads = get_value(
            filtering,
            "too_many_N_reads",
            0
        )

        too_short_reads = get_value(
            filtering,
            "too_short_reads",
            0
        )

        adapter_dimer_reads = get_value(
            filtering,
            "adapter_dimer_reads",
            0
        )

        # ----------------------------------------------------
        # RECORD
        # ----------------------------------------------------

        record = {

            "Isolate": isolate,

            "Reads_Before": before_reads,
            "Reads_After": after_reads,

            "Read_Retention_Percent": read_retention,

            "Q20_Before_Percent": q20_before,
            "Q20_After_Percent": q20_after,

            "Q30_Before_Percent": q30_before,
            "Q30_After_Percent": q30_after,

            "GC_Before_Percent": gc_before,
            "GC_After_Percent": gc_after,

            "Low_Quality_Reads": low_quality_reads,
            "Too_Many_N_Reads": too_many_n_reads,
            "Too_Short_Reads": too_short_reads,
            "Adapter_Dimer_Reads": adapter_dimer_reads,

            "Adapter_Trimmed_Reads": adapter_trimmed_reads,
            "Adapter_Trimmed_Bases": adapter_trimmed_bases,

            "Duplication_Percent": duplication_percent,

            "Insert_Size_Peak": insert_peak
        }

        records.append(record)

    except Exception as e:

        errors.append(
            f"{isolate}: {type(e).__name__}: {e}"
        )


# ============================================================
# REPORT PARSING STATUS
# ============================================================

if errors:

    print()
    print("WARNING: Some reports could not be parsed.")
    print(f"Successful reports : {len(records)}")
    print(f"Failed reports     : {len(errors)}")

    with open(
        ERROR_OUTPUT,
        "w",
        encoding="utf-8"
    ) as f:

        for error in errors:
            f.write(error + "\n")

    print()
    print(f"Error log: {ERROR_OUTPUT}")

else:

    print("All JSON reports parsed successfully.")


# ============================================================
# CREATE DATAFRAME
# ============================================================

if not records:

    raise RuntimeError(
        "No fastp reports were successfully parsed."
    )


df = pd.DataFrame(records)


# ============================================================
# SORT BY ISOLATE
# ============================================================

df = df.sort_values(
    by="Isolate"
).reset_index(drop=True)


# ============================================================
# ROUND VALUES
# ============================================================

numeric_columns = [

    "Read_Retention_Percent",

    "Q20_Before_Percent",
    "Q20_After_Percent",

    "Q30_Before_Percent",
    "Q30_After_Percent",

    "GC_Before_Percent",
    "GC_After_Percent",

    "Duplication_Percent",

    "Insert_Size_Peak"
]


for column in numeric_columns:

    if column in df.columns:

        df[column] = pd.to_numeric(
            df[column],
            errors="coerce"
        )


# ============================================================
# SAVE CSV
# ============================================================

df.to_csv(
    CSV_OUTPUT,
    index=False
)


print()
print("=" * 60)
print("SUMMARY CREATED")
print("=" * 60)

print(
    f"Isolates summarized : {len(df)}"
)

print(
    f"CSV file            : {CSV_OUTPUT}"
)


# ============================================================
# BASIC STATISTICS
# ============================================================

print()
print("BASIC STATISTICS")
print("-" * 60)


print(
    f"Mean read retention : "
    f"{df['Read_Retention_Percent'].mean():.2f}%"
)

print(
    f"Mean Q20 before     : "
    f"{df['Q20_Before_Percent'].mean():.2f}%"
)

print(
    f"Mean Q20 after      : "
    f"{df['Q20_After_Percent'].mean():.2f}%"
)

print(
    f"Mean Q30 before     : "
    f"{df['Q30_Before_Percent'].mean():.2f}%"
)

print(
    f"Mean Q30 after      : "
    f"{df['Q30_After_Percent'].mean():.2f}%"
)

print(
    f"Mean GC before      : "
    f"{df['GC_Before_Percent'].mean():.2f}%"
)

print(
    f"Mean GC after       : "
    f"{df['GC_After_Percent'].mean():.2f}%"
)

print(
    f"Mean duplication    : "
    f"{df['Duplication_Percent'].mean():.2f}%"
)

print(
    f"Total adapter-trimmed reads : "
    f"{df['Adapter_Trimmed_Reads'].sum():,}"
)

print(
    f"Total adapter-trimmed bases : "
    f"{df['Adapter_Trimmed_Bases'].sum():,}"
)

print(
    f"Mean insert-size peak       : "
    f"{df['Insert_Size_Peak'].mean():.2f} bp"
)


# ============================================================
# PLOT DIRECTORY
# ============================================================

PLOT_DIR = RESULT_DIR / "fastp_plots"

PLOT_DIR.mkdir(
    parents=True,
    exist_ok=True
)


# ============================================================
# PLOT 1
# READ RETENTION DISTRIBUTION
# ============================================================

plt.figure(figsize=(8, 5))

plt.hist(
    df["Read_Retention_Percent"],
    bins=15
)

plt.xlabel("Read Retention (%)")
plt.ylabel("Number of Isolates")
plt.title("FASTP Read Retention Distribution")

plt.tight_layout()

plt.savefig(
    PLOT_DIR / "01_read_retention_distribution.png",
    dpi=300
)

plt.close()


# ============================================================
# PLOT 2
# Q20 BEFORE VS AFTER
# ============================================================

plt.figure(figsize=(7, 5))

plt.boxplot(
    [
        df["Q20_Before_Percent"],
        df["Q20_After_Percent"]
    ],
    tick_labels=["Before", "After"]
)

plt.ylabel("Q20 (%)")
plt.title("Q20 Quality Before and After FASTP")

plt.tight_layout()

plt.savefig(
    PLOT_DIR / "02_q20_before_after.png",
    dpi=300
)

plt.close()


# ============================================================
# PLOT 3
# Q30 BEFORE VS AFTER
# ============================================================

plt.figure(figsize=(7, 5))

plt.boxplot(
    [
        df["Q30_Before_Percent"],
        df["Q30_After_Percent"]
    ],
    tick_labels=["Before", "After"]
)

plt.ylabel("Q30 (%)")
plt.title("Q30 Quality Before and After FASTP")

plt.tight_layout()

plt.savefig(
    PLOT_DIR / "03_q30_before_after.png",
    dpi=300
)

plt.close()


# ============================================================
# PLOT 4
# ADAPTER TRIMMED READS
# ============================================================

plt.figure(figsize=(10, 5))

plt.bar(
    range(len(df)),
    df["Adapter_Trimmed_Reads"]
)

plt.xlabel("Isolate Index")
plt.ylabel("Adapter-Trimmed Reads")
plt.title("Adapter-Trimmed Reads Across Isolates")

plt.tight_layout()

plt.savefig(
    PLOT_DIR / "04_adapter_trimmed_reads.png",
    dpi=300
)

plt.close()


# ============================================================
# PLOT 5
# GC CONTENT BEFORE VS AFTER
# ============================================================

plt.figure(figsize=(7, 5))

plt.boxplot(
    [
        df["GC_Before_Percent"],
        df["GC_After_Percent"]
    ],
    tick_labels=["Before", "After"]
)

plt.ylabel("GC Content (%)")
plt.title("GC Content Before and After FASTP")

plt.tight_layout()

plt.savefig(
    PLOT_DIR / "05_gc_before_after.png",
    dpi=300
)

plt.close()


# ============================================================
# ADDITIONAL PLOT
# DUPLICATION DISTRIBUTION
# ============================================================

plt.figure(figsize=(8, 5))

plt.hist(
    df["Duplication_Percent"],
    bins=15
)

plt.xlabel("Duplication (%)")
plt.ylabel("Number of Isolates")
plt.title("FASTP Duplication Rate Distribution")

plt.tight_layout()

plt.savefig(
    PLOT_DIR / "06_duplication_distribution.png",
    dpi=300
)

plt.close()


# ============================================================
# ADDITIONAL PLOT
# INSERT SIZE DISTRIBUTION
# ============================================================

plt.figure(figsize=(8, 5))

plt.hist(
    df["Insert_Size_Peak"],
    bins=15
)

plt.xlabel("Insert Size Peak (bp)")
plt.ylabel("Number of Isolates")
plt.title("Insert Size Peak Distribution")

plt.tight_layout()

plt.savefig(
    PLOT_DIR / "07_insert_size_distribution.png",
    dpi=300
)

plt.close()


# ============================================================
# FINAL CHECK
# ============================================================

print()
print("=" * 60)
print("PLOTS CREATED")
print("=" * 60)

print(
    f"Plot directory: {PLOT_DIR}"
)

for plot_file in sorted(PLOT_DIR.glob("*.png")):

    print(
        f"  {plot_file.name}"
    )


# ============================================================
# FINAL SUMMARY
# ============================================================

print()
print("=" * 60)
print("FASTP SUMMARY COMPLETE")
print("=" * 60)

print(
    f"Reports processed : {len(df)} / {len(json_files)}"
)

print(
    f"CSV               : {CSV_OUTPUT}"
)

print(
    f"Plots             : {PLOT_DIR}"
)

if errors:

    print(
        f"Warnings          : {len(errors)} report(s) failed"
    )

    print(
        f"Error log         : {ERROR_OUTPUT}"
    )

else:

    print(
        "Status            : SUCCESS"
    )

print("=" * 60)