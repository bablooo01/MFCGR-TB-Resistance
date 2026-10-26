import gzip
import sys

bam = sys.argv[1]

with gzip.open(bam, "rb") as f:
    magic = f.read(4)
    if magic != b"BAM\x01":
        print("Unexpected BAM signature")
        sys.exit(1)

    header_length = int.from_bytes(f.read(4), "little", signed=True)
    header = f.read(header_length).decode("utf-8", errors="replace")

print("BAM header for:", bam)
print(header[:6000])
