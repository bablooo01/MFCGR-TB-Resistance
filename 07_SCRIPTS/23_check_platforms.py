import csv
import json
import time
from pathlib import Path
from urllib.request import urlopen, Request
from urllib.parse import urlencode

root = Path("/mnt/d/MFCGR_PROJECT")
source = root / "08_RESULTS/selected_120_isolates.csv"
output = root / "08_RESULTS/sequencing_platforms_120.csv"

with source.open(newline="", encoding="utf-8-sig") as f:
    samples = list(csv.DictReader(f))

fields = [
    "run_accession", "sample_accession",
    "instrument_platform", "instrument_model",
    "library_strategy", "status"
]

with output.open("w", newline="", encoding="utf-8") as f:
    writer = csv.DictWriter(f, fieldnames=fields)
    writer.writeheader()

    for i, sample in enumerate(samples, 1):
        run = sample["run_accession"].strip()
        params = {
            "accession": run,
            "result": "read_run",
            "fields": ",".join(fields[:-1]),
            "format": "json"
        }
        url = "https://www.ebi.ac.uk/ena/portal/api/filereport?" + urlencode(params)
        row = {key: "" for key in fields}
        row["run_accession"] = run
        row["sample_accession"] = sample.get("sample_accession", "")
        row["status"] = "NOT_FOUND"

        try:
            req = Request(url, headers={"User-Agent": "MFCGR-research/1.0"})
            with urlopen(req, timeout=20) as response:
                data = json.load(response)
            if data:
                for key in fields[:-1]:
                    if data[0].get(key):
                        row[key] = data[0][key]
                row["status"] = "OK"
        except Exception as e:
            row["status"] = "ERROR: " + str(e)[:80]

        writer.writerow(row)
        f.flush()
        print(f"[{i}/120] {run}: {row['status']}")
        time.sleep(0.1)

print("\nReport saved to:", output)
