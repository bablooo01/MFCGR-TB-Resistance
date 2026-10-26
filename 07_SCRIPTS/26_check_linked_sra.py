
import csv
import json
import time
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET
from pathlib import Path

BASE = Path(r"D:\MFCGR_PROJECT")
INPUT = BASE / "08_RESULTS" / "ncbi_biosample_check.csv"
OUTPUT = BASE / "08_RESULTS" / "linked_sra_experiments.csv"

def get(url):
    req = urllib.request.Request(
        url, headers={"User-Agent": "MFCGR-research/1.0"}
    )
    with urllib.request.urlopen(req, timeout=40) as r:
        return r.read()

def api(endpoint, params):
    return get(
        "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
        + endpoint + "?" + urllib.parse.urlencode(params)
    )

with INPUT.open(newline="", encoding="utf-8-sig") as f:
    samples = list(csv.DictReader(f))

results = []

for i, sample in enumerate(samples, 1):
    accession = sample["accession"]
    print(f"{i}/{len(samples)} {accession}")

    try:
        search = json.loads(api("esearch.fcgi", {
            "db": "biosample",
            "term": f'"{accession}"[Accession]',
            "retmode": "json"
        }).decode())

        ids = search["esearchresult"]["idlist"]
        if not ids:
            continue

        linked = ET.fromstring(api("elink.fcgi", {
            "dbfrom": "biosample",
            "db": "sra",
            "id": ",".join(ids),
            "retmode": "xml"
        }))

        sra_ids = list(dict.fromkeys(
            e.text for e in linked.iter("Id") if e.text
        ))

        if not sra_ids:
            results.append({
                "sample_accession": accession,
                "sra_record_id": "",
                "experiment": "",
                "run": "",
                "title": "",
                "platform": "",
                "instrument": "",
                "library_strategy": "",
                "library_source": "",
                "notes": "No linked SRA record returned"
            })
            continue

        xml = ET.fromstring(api("efetch.fcgi", {
            "db": "sra",
            "id": ",".join(sra_ids),
            "retmode": "xml"
        }))

        for exp in xml.findall(".//EXPERIMENT"):
            parent = next(
                (p for p in xml.iter()
                 if exp in list(p)), None
            )
            run = exp.find(".//RUN_REF")
            platform = exp.findtext(
                ".//PLATFORM/*"
            ) or ""
            instrument = exp.findtext(
                ".//INSTRUMENT_MODEL"
            ) or ""
            strategy = exp.findtext(
                ".//LIBRARY_STRATEGY"
            ) or ""
            source = exp.findtext(
                ".//LIBRARY_SOURCE"
            ) or ""
            title = exp.findtext(".//TITLE") or ""

            run_accession = (
                run.get("accession", "") if run is not None else ""
            )

            results.append({
                "sample_accession": accession,
                "sra_record_id": "",
                "experiment": exp.get("accession", ""),
                "run": run_accession,
                "title": title,
                "platform": platform,
                "instrument": instrument,
                "library_strategy": strategy,
                "library_source": source,
                "notes": ""
            })

        time.sleep(0.4)

    except Exception as e:
        results.append({
            "sample_accession": accession,
            "sra_record_id": "",
            "experiment": "",
            "run": "",
            "title": "",
            "platform": "",
            "instrument": "",
            "library_strategy": "",
            "library_source": "",
            "notes": "ERROR: " + str(e)
        })

with OUTPUT.open("w", newline="", encoding="utf-8-sig") as f:
    fields = [
        "sample_accession", "sra_record_id", "experiment", "run",
        "title", "platform", "instrument", "library_strategy",
        "library_source", "notes"
    ]
    writer = csv.DictWriter(f, fieldnames=fields)
    writer.writeheader()
    writer.writerows(results)

print("\nFinished. Records:", len(results))
print("Saved to:", OUTPUT)
