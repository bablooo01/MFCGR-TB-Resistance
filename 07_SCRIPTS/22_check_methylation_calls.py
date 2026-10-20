import csv
import urllib.request
from pathlib import Path

BASE = Path("/mnt/d/MFCGR_PROJECT")
OUT = BASE / "08_RESULTS"
OUT.mkdir(parents=True, exist_ok=True)

biosamples = {
    "SAMN50582507": "UTB-20",
    "SAMN50582508": "UTB-21",
    "SAMN50582512": "UTB-28",
    "SAMN50582513": "UTB-29",
    "SAMN50582516": "UTB-42",
    "SAMN50582519": "UTB-65",
    "SAMN50582522": "UTB-82",
    "SAMN50582524": "UTB-94",
    "SAMN50582528": "UTB-100",
    "SAMN50582530": "UTB-127",
}

url = (
    "https://ftp.ncbi.nlm.nih.gov/"
    "pub/supplementary_data/basemodification.csv"
)

print("=" * 72)
print("MFCGR - CHECKING PUBLIC PACBIO METHYLATION-CALL FILES")
print("=" * 72)
print("Fetching NCBI's public base-modification catalogue...")

try:
    request = urllib.request.Request(
        url, headers={"User-Agent": "MFCGR-research/1.0"}
    )
    with urllib.request.urlopen(request, timeout=60) as response:
        content = response.read().decode("utf-8-sig")

except Exception as error:
    print("\nCould not retrieve the NCBI catalogue.")
    print("Error:", error)
    print("No conclusion can be drawn about file availability.")
    raise SystemExit(1)

reader = csv.DictReader(content.splitlines())
catalogue = list(reader)

print("Catalogue rows:", len(catalogue))
if catalogue:
    print("Catalogue columns:", list(catalogue[0].keys()))

matches = []

for row in catalogue:
    values = {str(v).strip() for v in row.values() if v is not None}

    for accession, isolate in biosamples.items():
        if accession in values:
            record = dict(row)
            record["BioSample_checked"] = accession
            record["Isolate_label"] = isolate
            matches.append(record)

output = OUT / "south_africa_methylation_call_matches.csv"

if matches:
    fieldnames = list(dict.fromkeys(
        key for row in matches for key in row.keys()
    ))
    with output.open("w", newline="", encoding="utf-8") as file:
        writer = csv.DictWriter(file, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(matches)

    print("\nEXACT MATCHES FOUND:", len(matches))
    for row in matches:
        print(
            row["Isolate_label"],
            row["BioSample_checked"],
            "| File type:", row.get("File Type", ""),
            "| URI:", row.get("URI", "")
        )
    print("\nMatches saved to:", output)

else:
    print("\nNo exact BioSample matches found in this catalogue.")
    print("This does not prove that no methylation files exist elsewhere.")
    print("Check the study's supplementary files and BioSample records too.")

print("\nCHECK COMPLETE")
