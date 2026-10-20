import csv
import time
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET
from pathlib import Path

PROJECT = Path("/mnt/d/MFCGR_PROJECT")
OUT = PROJECT / "08_RESULTS"
OUT.mkdir(parents=True, exist_ok=True)

SAMPLES = [
    ("UTB-20", "SAMN50582507"),
    ("UTB-21", "SAMN50582508"),
    ("UTB-28", "SAMN50582512"),
    ("UTB-29", "SAMN50582513"),
    ("UTB-42", "SAMN50582516"),
    ("UTB-65", "SAMN50582519"),
    ("UTB-82", "SAMN50582522"),
    ("UTB-94", "SAMN50582524"),
    ("UTB-100", "SAMN50582528"),
    ("UTB-127", "SAMN50582530"),
]

BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"


def fetch(url):
    request = urllib.request.Request(
        url,
        headers={"User-Agent": "MFCGR-research/1.0"}
    )
    with urllib.request.urlopen(request, timeout=45) as response:
        return response.read()


def get_run_records(biosample):
    query = urllib.parse.urlencode({
        "db": "sra",
        "term": biosample,
        "retmax": 100,
        "retmode": "json",
    })

    data = fetch(BASE + "esearch.fcgi?" + query)
    import json
    ids = json.loads(data)["esearchresult"]["idlist"]

    if not ids:
        return []

    query = urllib.parse.urlencode({
        "db": "sra",
        "id": ",".join(ids),
        "rettype": "full",
        "retmode": "xml",
    })

    xml_data = fetch(BASE + "efetch.fcgi?" + query)
    root = ET.fromstring(xml_data)

    records = []

    for package in root.iter():
        if package.tag.split("}")[-1] != "EXPERIMENT_PACKAGE":
            continue

        sample_text = ET.tostring(
            package, encoding="unicode"
        )

        if biosample not in sample_text:
            continue

        sample_name = ""
        for node in package.iter():
            tag = node.tag.split("}")[-1]
            if tag in ("SAMPLE_NAME", "TITLE") and node.text:
                if not sample_name:
                    sample_name = node.text.strip()

        platform = ""
        instrument = ""
        for node in package.iter():
            tag = node.tag.split("}")[-1]
            if tag == "INSTRUMENT_MODEL" and node.text:
                instrument = node.text.strip()
            if tag in ("PACBIO_SMRT", "ILLUMINA"):
                platform = tag

        runs = []
        for node in package.iter():
            if node.tag.split("}")[-1] == "RUN":
                runs.append({
                    "run_accession": node.attrib.get("accession", ""),
                    "total_bases": node.attrib.get("total_bases", ""),
                    "size": node.attrib.get("size", ""),
                })

        if not runs:
            runs = [{}]

        for run in runs:
            records.append({
                "biosample": biosample,
                "sample_name": sample_name,
                "platform": platform,
                "instrument": instrument,
                **run,
                "record_found": "YES",
            })

    return records


print("=" * 72)
print("MFCGR - PACBIO REVIO PUBLIC DATA AVAILABILITY CHECK")
print("=" * 72)

results = []

for i, (isolate, biosample) in enumerate(SAMPLES, start=1):
    print(f"\n[{i}/{len(SAMPLES)}] {isolate} | {biosample}")

    try:
        records = get_run_records(biosample)

        if not records:
            print("No SRA run record returned for this BioSample.")
            results.append({
                "isolate": isolate,
                "biosample": biosample,
                "record_found": "NO SRA RUN RETURNED",
            })
        else:
            for record in records:
                record["isolate"] = isolate
                results.append(record)
                print(
                    "Run:", record.get("run_accession", ""),
                    "| Platform:", record.get("platform", ""),
                    "| Instrument:", record.get("instrument", ""),
                    "| Bases:", record.get("total_bases", ""),
                    "| Size:", record.get("size", "")
                )

    except Exception as exc:
        print("NCBI query error:", str(exc))
        results.append({
            "isolate": isolate,
            "biosample": biosample,
            "record_found": "QUERY ERROR",
            "error": str(exc),
        })

    time.sleep(0.4)

report = OUT / "south_africa_pacbio_sra_inventory.csv"

columns = [
    "isolate", "biosample", "record_found", "sample_name",
    "platform", "instrument", "run_accession",
    "total_bases", "size", "error"
]

with report.open("w", newline="", encoding="utf-8") as handle:
    writer = csv.DictWriter(handle, fieldnames=columns)
    writer.writeheader()
    writer.writerows(results)

print("\n" + "=" * 72)
print("CHECK COMPLETE")
print("BioSamples checked:", len(SAMPLES))
print("Run records returned:", sum(
    1 for r in results if r.get("run_accession")
))
print("Report:", report)
print("=" * 72)
print(
    "\nNote: an SRA run record confirms sequencing data metadata, "
    "not necessarily that methylation-call files are deposited."
)
