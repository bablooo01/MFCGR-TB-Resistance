import json
import re
from pathlib import Path

PROJECT = Path("/mnt/d/MFCGR_PROJECT")
REPORT_DIR = PROJECT / "08_RESULTS" / "fastp_reports"

ISOLATES = [
    "ERR4810741",
    "ERR038300",
]


def repair_json(text):
    pattern = r'\\(?!["\\/bfnrtu])'
    return re.sub(pattern, r"\\\\", text)


def load_report(isolate):
    path = REPORT_DIR / f"{isolate}.json"

    with open(path, "r", encoding="utf-8") as f:
        text = f.read()

    try:
        return json.loads(text)
    except json.JSONDecodeError:
        print(f"JSON repair required: {isolate}")
        repaired = repair_json(text)
        return json.loads(repaired)


for isolate in ISOLATES:

    print()
    print("=" * 70)
    print(f"FASTP QC: {isolate}")
    print("=" * 70)

    data = load_report(isolate)

    before = data.get("summary", {}).get("before_filtering", {})
    after = data.get("summary", {}).get("after_filtering", {})
    filtering = data.get("filtering_result", {})
    adapter = data.get("adapter_cutting", {})

    print("\nBEFORE FILTERING")
    print(f"Reads:       {before.get('total_reads')}")
    print(f"Bases:       {before.get('total_bases')}")
    print(f"Q20:         {before.get('q20_rate')}")
    print(f"Q30:         {before.get('q30_rate')}")
    print(f"GC:          {before.get('gc_content')}")

    print("\nAFTER FILTERING")
    print(f"Reads:       {after.get('total_reads')}")
    print(f"Bases:       {after.get('total_bases')}")
    print(f"Q20:         {after.get('q20_rate')}")
    print(f"Q30:         {after.get('q30_rate')}")
    print(f"GC:          {after.get('gc_content')}")

    print("\nFILTERING")
    print(f"Passed reads:       {filtering.get('passed_filter_reads')}")
    print(f"Low-quality reads:  {filtering.get('low_quality_reads')}")
    print(f"Too many N reads:   {filtering.get('too_many_N_reads')}")
    print(f"Too short reads:    {filtering.get('too_short_reads')}")

    print("\nADAPTER")
    print(f"Adapter-trimmed reads: {adapter.get('adapter_trimmed_reads')}")
    print(f"Adapter-trimmed bases: {adapter.get('adapter_trimmed_bases')}")


print()
print("=" * 70)
print("FASTP OUTLIER INVESTIGATION COMPLETE")
print("=" * 70)