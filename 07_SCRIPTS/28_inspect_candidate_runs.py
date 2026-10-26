
import csv
import io
import time
import urllib.parse
import urllib.request
from pathlib import Path

BASE = Path(r"D:\MFCGR_PROJECT")
SOURCE = BASE / "08_RESULTS" / "candidate_run_details.csv"
OUTPUT = BASE / "08_RESULTS" / "candidate_run_metadata.csv"

def fetch(url):
    request = urllib.request.Request(
        url, headers={"User-Agent": "MFCGR-research/1.0"}
    )
    with urllib.request.urlopen(request, timeout=45) as response:
        return response.read().decode("utf-8", errors="replace")

with SOURCE.open(newline="", encoding="utf-8-sig") as f:
    candidates = list(csv.DictReader(f))

results = []

for i, row in enumerate(candidates, 1):
    run = (row.get("run_accession") or "").strip()
    print(f"{i}/{len(candidates)} Checking {run or 'missing run'}")

    result = {
        "sample_accession": row.get("sample_accession", ""),
        "experiment": row.get("experiment", ""),
        "run_accession": run,
        "organism": "",
        "study_accession": "",
        "sample_name": "",
        "platform": "",
        "instrument": "",
        "library_strategy": "",
        "spots": "",
        "bases": "",
        "download_path": "",
        "status": ""
    }

    if not run:
        result["status"] = "Missing run accession"
        results.append(result)
        continue

    try:
        url = (
            "https://trace.ncbi.nlm.nih.gov/Traces/sra-db-be/runinfo?"
            + urllib.parse.urlencode({"acc": run})
        )
        content = fetch(url)
        records = list(csv.DictReader(io.StringIO(content)))

        if not records:
            result["status"] = (
                "No RunInfo row; inspect public archive record manually"
            )
        else:
            r = records[0]
            for output_key, possible_keys in {
                "organism": ["ScientificName", "scientific_name"],
                "study_accession": ["BioProject", "Study", "study_accession"],
                "sample_name": ["SampleName", "sample_name"],
                "platform": ["Platform", "platform"],
                "instrument": ["Model", "instrument"],
                "library_strategy": ["LibraryStrategy", "library_strategy"],
                "spots": ["spots"],
                "bases": ["bases"],
                "download_path": ["download_path"]
            }.items():
                for key in possible_keys:
                    if r.get(key):
                        result[output_key] = r[key]
                        break

            result["status"] = "RunInfo retrieved"

    except Exception as e:
        result["status"] = "ERROR: " + str(e)

    results.append(result)
    time.sleep(0.5)

with OUTPUT.open("w", newline="", encoding="utf-8-sig") as f:
    writer = csv.DictWriter(f, fieldnames=list(results[0].keys()))
    writer.writeheader()
    writer.writerows(results)

print("\nSaved:", OUTPUT)
print("Records processed:", len(results))
