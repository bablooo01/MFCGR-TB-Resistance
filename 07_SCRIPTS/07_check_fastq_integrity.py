from pathlib import Path
from urllib.parse import urlparse
import hashlib
import pandas as pd

ROOT = Path(r"D:\MFCGR_PROJECT")
CSV_PATH = ROOT / "08_RESULTS" / "selected_120_isolates.csv"
FASTQ_DIR = ROOT / "03_GENOMES"
OUT_DIR = ROOT / "08_RESULTS"

def split_values(value):
    """Split semicolon-separated metadata; return an empty list if missing."""
    if pd.isna(value) or str(value).strip() == "":
        return []
    return [x.strip() for x in str(value).split(";") if x.strip()]

def get_column(df, name):
    """Find a column regardless of capitalization or surrounding spaces."""
    for col in df.columns:
        if str(col).strip().lower() == name.lower():
            return col
    return None

def md5_file(path):
    digest = hashlib.md5()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(8 * 1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()

def filename_from_url(url):
    # The source may omit the ftp:// prefix.
    if url.startswith("ftp.sra.ebi.ac.uk/"):
        url = "https://" + url
    elif url.startswith("ftp://"):
        url = "https://" + url[len("ftp://"):]

    return Path(urlparse(url).path).name

def find_file(filename):
    matches = list(FASTQ_DIR.rglob(filename))
    return matches[0] if matches else None

df = pd.read_csv(CSV_PATH, low_memory=False)

url_col = get_column(df, "fastq_ftp")
bytes_col = get_column(df, "fastq_bytes")
md5_col = get_column(df, "fastq_md5")
id_col = get_column(df, "UNIQUEID")

if url_col is None:
    raise ValueError("Could not find the fastq_ftp column in the CSV.")

results = []

for _, row in df.iterrows():
    isolate_id = str(row[id_col]) if id_col else "Unknown isolate"

    urls = split_values(row[url_col])
    expected_sizes = split_values(row[bytes_col]) if bytes_col else []
    expected_md5s = split_values(row[md5_col]) if md5_col else []

    for i, url in enumerate(urls):
        filename = filename_from_url(url)
        path = find_file(filename)

        expected_size = None
        if i < len(expected_sizes):
            try:
                expected_size = int(float(expected_sizes[i]))
            except ValueError:
                pass

        expected_md5 = (
            expected_md5s[i].lower()
            if i < len(expected_md5s)
            else None
        )

        status = ""
        actual_size = None
        actual_md5 = None

        if path is None:
            status = "MISSING"
        else:
            actual_size = path.stat().st_size

            if expected_size is not None and actual_size != expected_size:
                status = "SIZE_MISMATCH"
            elif expected_md5:
                actual_md5 = md5_file(path).lower()
                status = "OK" if actual_md5 == expected_md5 else "MD5_MISMATCH"
            elif expected_size is not None:
                status = "OK_SIZE_ONLY"
            else:
                status = "PRESENT_NO_REFERENCE"

        results.append({
            "UNIQUEID": isolate_id,
            "filename": filename,
            "status": status,
            "expected_bytes": expected_size,
            "actual_bytes": actual_size,
            "expected_md5": expected_md5,
            "actual_md5": actual_md5,
            "local_path": str(path) if path else ""
        })

result_df = pd.DataFrame(results)
OUT_DIR.mkdir(parents=True, exist_ok=True)

file_report = OUT_DIR / "fastq_integrity_report.csv"
result_df.to_csv(file_report, index=False)

print("\nFILE CHECK SUMMARY")
print(result_df["status"].value_counts(dropna=False).to_string())
print(f"\nFile report saved to:\n{file_report}")

# Summarize whether every expected file for each isolate passed.
good_statuses = {"OK", "OK_SIZE_ONLY"}

isolate_report = (
    result_df.assign(passed=result_df["status"].isin(good_statuses))
    .groupby("UNIQUEID", as_index=False)
    .agg(expected_files=("filename", "count"),
         passed_files=("passed", "sum"))
)

isolate_report["all_files_passed"] = (
    isolate_report["expected_files"] == isolate_report["passed_files"]
)

isolate_report_path = OUT_DIR / "fastq_isolate_completeness.csv"
isolate_report.to_csv(isolate_report_path, index=False)

print("\nISOLATE SUMMARY")
print(f"Isolates checked: {len(isolate_report)}")
print(f"Isolates with all files passing: "
      f"{isolate_report['all_files_passed'].sum()}")
print(f"Isolate report saved to:\n{isolate_report_path}")