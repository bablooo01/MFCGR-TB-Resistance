from pathlib import Path
import subprocess
import csv

BASE = Path("/mnt/d/MFCGR_PROJECT")
VCF_DIR = BASE / "08_RESULTS/variants"
OUT_DIR = BASE / "08_RESULTS/alignment"
WHO_FILE = BASE / "06_WHO_CATALOGUE/WHO-UCN-TB-2023.7-eng_genomic_coordinates.txt"

# Load WHO variants using exact chromosome, position, REF and ALT matching.
who_index = {}
with WHO_FILE.open(encoding="utf-8-sig") as f:
    reader = csv.DictReader(f, delimiter="\t")
    for row in reader:
        chrom = row["chromosome"].strip()
        pos = row["position"].strip()
        ref = row["reference_nucleotide"].strip()
        alt = row["alternative_nucleotide"].strip()
        variant = row["variant"].strip()
        who_index.setdefault((chrom, pos, ref, alt), set()).add(variant)

vcfs = sorted(VCF_DIR.glob("*.filtered.vcf.gz"))
print(f"WHO catalogue entries loaded: {sum(map(len, who_index.values()))}")
print(f"Filtered VCF files found: {len(vcfs)}")

success = 0
failed = []

for i, vcf in enumerate(vcfs, 1):
    sample = vcf.name.replace(".filtered.vcf.gz", "")
    output = OUT_DIR / f"{sample}.WHO_matches.txt"

    try:
        result = subprocess.run(
            ["bcftools", "query", "-f",
             "%CHROM\\t%POS\\t%REF\\t%ALT\\n", str(vcf)],
            capture_output=True, text=True, check=True
        )

        matches = set()
        for line in result.stdout.splitlines():
            fields = line.split("\t")
            if len(fields) != 4:
                continue
            chrom, pos, ref, alts = fields
            for alt in alts.split(","):
                key = (chrom, pos, ref, alt)
                for variant in who_index.get(key, ()):
                    matches.add(
                        f"{chrom}\t{pos}\t{ref}\t{alt}\tWHO_VARIANT={variant}"
                    )

        with output.open("w", encoding="utf-8") as f:
            for line in sorted(matches):
                f.write(line + "\n")

        success += 1
        print(f"[{i}/{len(vcfs)}] {sample}: {len(matches)} WHO matches")

    except Exception as exc:
        failed.append((sample, str(exc)))
        print(f"[{i}/{len(vcfs)}] FAILED: {sample}: {exc}")

print("\n===== BATCH WHO MATCHING SUMMARY =====")
print(f"VCFs found : {len(vcfs)}")
print(f"Successful : {success}")
print(f"Failed     : {len(failed)}")
if failed:
    for sample, error in failed:
        print(f"FAILED {sample}: {error}")
print(f"Output directory: {OUT_DIR}")
