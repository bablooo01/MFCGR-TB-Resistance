
import csv
import time
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET
from pathlib import Path

BASE = Path(r"D:\MFCGR_PROJECT")
INPUT = BASE / "08_RESULTS" / "all_runs_for_120.csv"
OUTPUT = BASE / "08_RESULTS" / "ncbi_biosample_check.csv"

KEYWORDS = (
    "methylation", "modified base", "modification",
    "methylome", "kinetics", "supplementary file",
    "modbam", "mm ml tag", "pacbio", "nanopore",
    "oxford nanopore"
)

def fetch(url):
    request = urllib.request.Request(
        url, headers={"User-Agent": "MFCGR-research/1.0"}
    )
    with urllib.request.urlopen(request, timeout=30) as response:
        return response.read()

with INPUT.open(newline="", encoding="utf-8-sig") as f:
    rows = list(csv.DictReader(f))

# Collect sample identifiers already present in the existing results.
accessions = set()
for row in rows:
    for field in ("sample_accession", "secondary_sample_accession",
                  "original_sample_accession"):
        value = (row.get(field) or "").strip()
        if value and value.upper().startswith(("SAMN", "SAMEA")):
            accessions.add(value)

results = []

for i, accession in enumerate(sorted(accessions), 1):
    print(f"{i}/{len(accessions)} Checking {accession}")

    result = {
        "accession": accession,
        "ncbi_record_found": "No",
        "biosample_title": "",
        "matching_details": "",
        "record_url": ""
    }

    try:
        # Search NCBI BioSample using the exact accession.
        params = urllib.parse.urlencode({
            "db": "biosample",
            "term": f'"{accession}"[Accession]',
            "retmode": "json"
        })
        search_url = (
            "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?"
            + params
        )
        data = fetch(search_url)
        import json
        search = json.loads(data.decode("utf-8"))
        ids = search.get("esearchresult", {}).get("idlist", [])

        if ids:
            result["ncbi_record_found"] = "Yes"
            result["record_url"] = (
                "https://www.ncbi.nlm.nih.gov/biosample/?term=" + accession
            )

            efetch_url = (
                "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
                "efetch.fcgi?"
                + urllib.parse.urlencode({
                    "db": "biosample",
                    "id": ",".join(ids),
                    "retmode": "xml"
                })
            )
            root = ET.fromstring(fetch(efetch_url))
            result["biosample_title"] = " ".join(
                x.text or "" for x in root.iter("Title")
            )

            found = []
            for element in root.iter():
                text = " ".join(element.itertext()).strip()
                if text and any(k in text.lower() for k in KEYWORDS):
                    found.append(" ".join(text.split())[:500])

            result["matching_details"] = " | ".join(dict.fromkeys(found))

    except Exception as e:
        result["matching_details"] = "CHECK_ERROR: " + str(e)

    results.append(result)
    time.sleep(0.4)

with OUTPUT.open("w", newline="", encoding="utf-8-sig") as f:
    writer = csv.DictWriter(f, fieldnames=[
        "accession", "ncbi_record_found", "biosample_title",
        "matching_details", "record_url"
    ])
    writer.writeheader()
    writer.writerows(results)

print("\nFinished.")
print("Accessions checked:", len(accessions))
print("Results saved to:", OUTPUT)
