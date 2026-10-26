
import csv
import time
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET
from pathlib import Path

BASE = Path(r"D:\MFCGR_PROJECT")
SOURCE = BASE / "08_RESULTS" / "linked_sra_experiments.csv"
OUTPUT = BASE / "08_RESULTS" / "candidate_run_details.csv"

def fetch(url):
    req = urllib.request.Request(
        url, headers={"User-Agent": "MFCGR-research/1.0"}
    )
    with urllib.request.urlopen(req, timeout=40) as response:
        return response.read()

with SOURCE.open(newline="", encoding="utf-8-sig") as f:
    rows = list(csv.DictReader(f))

candidates = [
    r for r in rows
    if (r.get("instrument") or "").lower()
       in {"sequel", "sequel ii", "minion", "gridion"}
]

results = []

for i, row in enumerate(candidates, 1):
    exp = row["experiment"].strip()
    print(f"{i}/{len(candidates)} Checking {exp}")

    try:
        url = (
            "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
            "efetch.fcgi?"
            + urllib.parse.urlencode({
                "db": "sra",
                "id": exp,
                "retmode": "xml"
            })
        )
        root = ET.fromstring(fetch(url))
        runs = root.findall(".//RUN")

        if not runs:
            results.append({
                "sample_accession": row["sample_accession"],
                "experiment": exp,
                "run_accession": "",
                "instrument": row["instrument"],
                "library_strategy": row["library_strategy"],
                "notes": "No RUN element returned; needs manual record review"
            })

        for run in runs:
            results.append({
                "sample_accession": row["sample_accession"],
                "experiment": exp,
                "run_accession": run.get("accession", ""),
                "instrument": row["instrument"],
                "library_strategy": row["library_strategy"],
                "notes": ""
            })

    except Exception as e:
        results.append({
            "sample_accession": row["sample_accession"],
            "experiment": exp,
            "run_accession": "",
            "instrument": row["instrument"],
            "library_strategy": row["library_strategy"],
            "notes": "ERROR: " + str(e)
        })

    time.sleep(0.4)

with OUTPUT.open("w", newline="", encoding="utf-8-sig") as f:
    fields = [
        "sample_accession", "experiment", "run_accession",
        "instrument", "library_strategy", "notes"
    ]
    writer = csv.DictWriter(f, fieldnames=fields)
    writer.writeheader()
    writer.writerows(results)

print("\nSaved:", OUTPUT)
print("Records:", len(results))
