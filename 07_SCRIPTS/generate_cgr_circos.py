import csv
import math
from pathlib import Path
from collections import Counter

ROOT = Path("/mnt/d/MFCGR_PROJECT")
RES = ROOT / "08_RESULTS"
FIG = RES / "figures"
FIG.mkdir(parents=True, exist_ok=True)

REF = ROOT / "01_H37RV/ncbi_dataset/data/GCF_000195955.2/GCF_000195955.2_ASM19595v2_genomic.fna"
GFF = ROOT / "01_H37RV/ncbi_dataset/data/GCF_000195955.2/genomic.gff"

try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np
except ImportError as e:
    raise SystemExit(f"Missing plotting dependency: {e}. Tell me the exact error.")

TARGETS = {"Rv2024c", "Rv2756c", "Rv2761c", "Rv3263"}

def read_fasta(path):
    seq = []
    with path.open() as f:
        for line in f:
            if line.startswith(">"):
                continue
            seq.append(line.strip().upper())
    return "".join(seq)

def parse_attrs(s):
    d = {}
    for item in s.split(";"):
        if "=" in item:
            k, v = item.split("=", 1)
            d[k] = v
    return d

def revcomp(seq):
    return seq.translate(str.maketrans("ACGTN", "TGCAN"))[::-1]

def cgr_matrix(seq, k=6):
    size = 2**k
    mat = np.zeros((size, size), dtype=np.float64)
    x = y = 0.5
    corners = {
        "A": (0.0, 0.0), "C": (0.0, 1.0),
        "G": (1.0, 1.0), "T": (1.0, 0.0)
    }
    points = 0
    for base in seq:
        if base not in corners:
            continue
        cx, cy = corners[base]
        x = (x + cx) / 2
        y = (y + cy) / 2
        ix = min(size - 1, max(0, int(x * size)))
        iy = min(size - 1, max(0, int(y * size)))
        mat[iy, ix] += 1
        points += 1
    if points:
        mat /= points
    return mat, points

if not REF.exists() or not GFF.exists():
    raise SystemExit("Reference FASTA or GFF is missing. Check the NCBI file paths.")

genome = read_fasta(REF)
genes = []
with GFF.open() as f:
    for line in f:
        if line.startswith("#"):
            continue
        p = line.rstrip().split("\t")
        if len(p) != 9 or p[0] != "NC_000962.3" or p[2] != "CDS":
            continue
        a = parse_attrs(p[8])
        locus = a.get("locus_tag", "")
        if locus in TARGETS:
            start, end = int(p[3]), int(p[4])
            genes.append({
                "locus": locus, "start": start, "end": end,
                "strand": p[6], "product": a.get("product", ""),
                "sequence": genome[start-1:end]
            })

found = {g["locus"] for g in genes}
missing = TARGETS - found
if missing:
    print("WARNING: no CDS annotation found for:", sorted(missing))

if not genes:
    raise SystemExit("No candidate CDS sequences found.")

# Save extracted sequences.
fasta_out = RES / "candidate_gene_sequences.fasta"
with fasta_out.open("w") as f:
    for g in genes:
        seq = g["sequence"]
        if g["strand"] == "-":
            seq = revcomp(seq)
        f.write(f'>{g["locus"]} {g["start"]}-{g["end"]} strand={g["strand"]}\n')
        for i in range(0, len(seq), 70):
            f.write(seq[i:i+70] + "\n")
print("Created:", fasta_out)

# CGR plots: one per candidate CDS.
summary = []
for g in genes:
    seq = g["sequence"]
    if g["strand"] == "-":
        seq = revcomp(seq)
    matrix, n = cgr_matrix(seq, k=6)
    fig, ax = plt.subplots(figsize=(6, 5))
    im = ax.imshow(matrix, origin="lower", interpolation="nearest", aspect="equal")
    ax.set_title(f'CGR: {g["locus"]} ({n:,} valid bases)')
    ax.set_xlabel("CGR x-coordinate")
    ax.set_ylabel("CGR y-coordinate")
    fig.colorbar(im, ax=ax, label="Fraction of sequence points")
    fig.tight_layout()
    out = FIG / f'CGR_{g["locus"]}.png'
    fig.savefig(out, dpi=220)
    plt.close(fig)
    summary.append({
        "locus": g["locus"], "start": g["start"], "end": g["end"],
        "strand": g["strand"], "length_bp": len(seq),
        "valid_ACGT_bases": n, "product": g["product"],
        "CGR_file": str(out)
    })
    print("Created:", out)

with (RES / "candidate_sequence_summary.csv").open("w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=summary[0].keys())
    w.writeheader()
    w.writerows(summary)

# Circular genome overview: candidate gene arcs and tested positions.
genome_len = len(genome)
fig, ax = plt.subplots(figsize=(9, 9), subplot_kw={"projection": "polar"})
ax.set_theta_direction(-1)
ax.set_theta_offset(math.pi / 2)
ax.set_ylim(0, 1.2)
ax.set_yticks([])
ax.set_xticks([])
ax.spines["polar"].set_visible(False)

theta = np.linspace(0, 2*math.pi, 600)
ax.plot(theta, np.full_like(theta, 0.88), linewidth=2)

for g in genes:
    mid = (g["start"] + g["end"]) / 2
    angle = 2 * math.pi * mid / genome_len
    width = 2 * math.pi * (g["end"] - g["start"] + 1) / genome_len
    ax.bar(angle, 0.10, width=max(width, 0.008), bottom=0.88,
           align="center", alpha=0.85)
    ax.text(angle, 1.04, g["locus"], ha="center", va="center", fontsize=9)

# Add tested positions as inner markers.
targets = [
    ("Rv2024c", 2269545), ("Rv2756c", 3069167),
    ("Rv2761c", 3072285), ("Rv3263", 3643985)
]
for locus, pos in targets:
    angle = 2 * math.pi * pos / genome_len
    ax.scatter([angle], [0.72], s=28)
    ax.text(angle, 0.63, f"{pos:,}", ha="center", va="center", fontsize=7)

ax.set_title("M. tuberculosis H37Rv: candidate loci and tested positions",
             pad=24)
out = FIG / "candidate_loci_circular_overview.png"
fig.savefig(out, dpi=240, bbox_inches="tight")
plt.close(fig)
print("Created:", out)

print("\nDONE")
print("Sequence summary:", RES / "candidate_sequence_summary.csv")
print("Figures folder:", FIG)
print("Note: these are sequence-based CGR and a Circos-style overview.")
print("They do not represent measured DNA methylation.")
