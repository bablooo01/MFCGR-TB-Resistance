
import csv
import json
import time
import urllib.parse
import urllib.request
from pathlib import Path

ROOT = Path(r"D:\MFCGR_PROJECT")
INPUT = ROOT / "08_RESULTS" / "selected_120_isolates.csv"
OUTPUT = ROOT / "08_RESULTS" / "all_runs_for_120.csv"

def fetch(url):
    req = urllib.request.Request(
        url, headers={"User-Agent": "MFCGR-research/1.0"}
    )
    with urllib.request.urlopen(req, timeout=30) as response:
        return response.read().decode("utf-8")

with INPUT.open(newline="", encoding="utf-8-sig") as f:
    isolates = list(csv.DictReader(f))

results = []
seen = set()

for i, isolate in enumerate(isolates, 1):
    sample = (isolate.get("sample_accession") or "").strip()
    if not sample:
        continue

    query = f'sample_accession="{sample}"'
    url = (
        "https://www.ebi.ac.uk/ena/portal/api/search?"
        + urllib.parse.urlencode({
            "result": "read_run",
            "query": query,
            "fields": (
                "run_accession,sample_accession,secondary_sample_accession,"
                "experiment_accession,study_accession,instrument_platform,"
                "instrument_model,library_strategy,library_source,"
                "fastq_ftp,submitted_ftp"
            ),
            "format": "json",
            "limit": "0",
        })
    )

    try:
        records = json.loads(fetch(url))
        for record in records:
            run = record.get("run_accession", "")
            key = (sample, run)
            if key in seen:
                continue
            seen.add(key)
            results.append({
                "UNIQUEID": isolate.get("UNIQUEID", ""),
                "CATEGORY": isolate.get("CATEGORY", ""),
                "original_sample_accession": sample,
                **record,
            })
        print(f"{i}/{len(isolates)} checked: {sample}")
    except Exception as e:
        print(f"Could not query {sample}: {e}")

    time.sleep(0.15)

fields = [
    "UNIQUEID", "CATEGORY", "original_sample_accession",
    "run_accession", "sample_accession", "secondary_sample_accession",
    "experiment_accession", "study_accession", "instrument_platform",
    "instrument_model", "library_strategy", "library_source",
    "fastq_ftp", "submitted_ftp",
]

with OUTPUT.open("w", newline="", encoding="utf-8-sig") as f:
    writer = csv.DictWriter(f, fieldnames=fields, extrasaction="ignore")
    writer.writeheader()
    writer.writerows(results)

print(f"\nFinished. Linked runs found: {len(results)}")
print(f"Results saved to: {OUTPUT}")
print("Inspect platforms for PacBio/ONT and files for methylation calls.")
