import pandas as pd
import requests
from pathlib import Path
import time


# ============================================================
# MFCGR PROJECT
# Download FASTQ files for selected 120 CRyPTIC isolates
# ============================================================

PROJECT_DIR = Path(r"D:\MFCGR_PROJECT")

# Input file containing the selected 120 isolates
INPUT_FILE = PROJECT_DIR / "08_RESULTS" / "selected_120_isolates.csv"

# Folder where FASTQ files will be stored
OUTPUT_DIR = PROJECT_DIR / "03_GENOMES"

# Download log
LOG_FILE = PROJECT_DIR / "08_RESULTS" / "fastq_download_log.csv"


# ============================================================
# CREATE REQUIRED DIRECTORIES
# ============================================================

OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
LOG_FILE.parent.mkdir(parents=True, exist_ok=True)


print("=" * 80)
print("MFCGR - FASTQ DOWNLOAD")
print("=" * 80)

print("\nProject directory:")
print(PROJECT_DIR)

print("\nInput file:")
print(INPUT_FILE)

print("\nFASTQ output directory:")
print(OUTPUT_DIR)


# ============================================================
# CHECK INPUT FILE
# ============================================================

if not INPUT_FILE.exists():
    print("\nERROR: Input file not found!")
    print(INPUT_FILE)
    print("\nPlease make sure selected_120_isolates.csv exists in:")
    print(PROJECT_DIR / "08_RESULTS")
    raise SystemExit(1)


# ============================================================
# LOAD SELECTED ISOLATES
# ============================================================

print("\nLoading selected isolates...")

df = pd.read_csv(INPUT_FILE)

print("Number of selected isolates:", len(df))

print("\nColumns available:")
print(df.columns.tolist())


# ============================================================
# CHECK REQUIRED COLUMNS
# ============================================================

if "UNIQUEID" not in df.columns:
    print("\nERROR: UNIQUEID column not found!")
    raise SystemExit(1)

if "fastq_ftp" not in df.columns:
    print("\nERROR: fastq_ftp column not found!")
    print("\nAvailable columns:")
    print(df.columns.tolist())
    raise SystemExit(1)


# ============================================================
# DOWNLOAD FUNCTION
# ============================================================

def convert_to_https(url):
    """
    Convert CRyPTIC/ENA FASTQ path into a downloadable HTTPS URL.
    """

    if pd.isna(url):
        return None

    url = str(url).strip()

    if url == "":
        return None

    # Already HTTPS
    if url.startswith("https://"):
        return url

    # Already HTTP
    if url.startswith("http://"):
        return url

    # FTP URL with scheme
    if url.startswith("ftp://"):
        return url.replace("ftp://", "https://", 1)

    # CRyPTIC format:
    # ftp.sra.ebi.ac.uk/vol1/fastq/...
    if url.startswith("ftp.sra.ebi.ac.uk/"):
        return "https://" + url

    # If the URL somehow starts with www
    if url.startswith("www."):
        return "https://" + url

    # Unknown format
    print("\nWARNING: Unknown URL format:")
    print(url)

    return None


# ============================================================
# DOWNLOAD ONE FILE
# ============================================================

def download_file(url, output_file, chunk_size=1024 * 1024):

    try:

        print("\nDownloading:")
        print(url)

        print("Saving to:")
        print(output_file)

        # ----------------------------------------------------
        # Request file
        # ----------------------------------------------------

        response = requests.get(
            url,
            stream=True,
            timeout=(30, 300)
        )

        # Raise error for HTTP problems
        response.raise_for_status()

        # ----------------------------------------------------
        # Get expected file size
        # ----------------------------------------------------

        total_size = response.headers.get("content-length")

        if total_size is not None:
            total_size = int(total_size)

        downloaded = 0

        # ----------------------------------------------------
        # Write file
        # ----------------------------------------------------

        with open(output_file, "wb") as f:

            for chunk in response.iter_content(
                chunk_size=chunk_size
            ):

                if chunk:

                    f.write(chunk)

                    downloaded += len(chunk)

                    # Print progress
                    if total_size:

                        percentage = (
                            downloaded / total_size
                        ) * 100

                        print(
                            f"\rProgress: "
                            f"{percentage:6.2f}% "
                            f"({downloaded / (1024**2):.1f} MB)",
                            end=""
                        )

                    else:

                        print(
                            f"\rDownloaded: "
                            f"{downloaded / (1024**2):.1f} MB",
                            end=""
                        )

        print()

        # ----------------------------------------------------
        # Verify file exists
        # ----------------------------------------------------

        if output_file.exists():

            file_size = output_file.stat().st_size

            if file_size > 0:

                print("DONE")

                print(
                    f"File size: "
                    f"{file_size / (1024**2):.2f} MB"
                )

                return True

        print("FAILED: Empty file")

        return False

    except Exception as e:

        print("\nFAILED")

        print("Error:", e)

        # Delete incomplete file
        if output_file.exists():

            try:
                output_file.unlink()
            except:
                pass

        return False


# ============================================================
# PREPARE DOWNLOAD LOG
# ============================================================

download_log = []


# ============================================================
# PROCESS ALL SELECTED ISOLATES
# ============================================================

total_isolates = len(df)

successful_files = 0
failed_files = 0
skipped_files = 0


print("\n")
print("=" * 80)
print("STARTING DOWNLOAD")
print("=" * 80)


for index, row in df.iterrows():

    uniqueid = str(row["UNIQUEID"])

    fastq_value = row["fastq_ftp"]

    print("\n")
    print("-" * 80)

    print(
        f"Isolate {index + 1} / {total_isolates}"
    )

    print("UNIQUEID:")
    print(uniqueid)


    # ========================================================
    # CHECK FASTQ URL
    # ========================================================

    if pd.isna(fastq_value):

        print("No FASTQ URL available.")

        download_log.append({
            "UNIQUEID": uniqueid,
            "URL": "",
            "LOCAL_FILE": "",
            "STATUS": "NO_URL"
        })

        failed_files += 1

        continue


    # ========================================================
    # CRyPTIC may contain multiple FASTQ URLs separated by ;
    # ========================================================

    urls = str(fastq_value).split(";")

    print(
        f"FASTQ files for this isolate: {len(urls)}"
    )


    # ========================================================
    # DOWNLOAD EACH FASTQ FILE
    # ========================================================

    for raw_url in urls:

        raw_url = raw_url.strip()

        if raw_url == "":
            continue


        # ----------------------------------------------------
        # Convert URL to HTTPS
        # ----------------------------------------------------

        url = convert_to_https(raw_url)

        if url is None:

            print("\nCould not convert URL:")
            print(raw_url)

            download_log.append({
                "UNIQUEID": uniqueid,
                "URL": raw_url,
                "LOCAL_FILE": "",
                "STATUS": "INVALID_URL"
            })

            failed_files += 1

            continue


        # ----------------------------------------------------
        # Extract filename
        # ----------------------------------------------------

        filename = url.split("/")[-1]

        if filename == "":

            print("ERROR: Could not determine filename.")

            download_log.append({
                "UNIQUEID": uniqueid,
                "URL": url,
                "LOCAL_FILE": "",
                "STATUS": "INVALID_FILENAME"
            })

            failed_files += 1

            continue


        # ----------------------------------------------------
        # Output path
        # ----------------------------------------------------

        output_file = OUTPUT_DIR / filename


        # ----------------------------------------------------
        # If file already exists
        # ----------------------------------------------------

        if output_file.exists():

            existing_size = output_file.stat().st_size

            if existing_size > 0:

                print("\nAlready exists:")
                print(output_file)

                print(
                    f"Existing size: "
                    f"{existing_size / (1024**2):.2f} MB"
                )

                download_log.append({
                    "UNIQUEID": uniqueid,
                    "URL": url,
                    "LOCAL_FILE": str(output_file),
                    "STATUS": "ALREADY_EXISTS"
                })

                skipped_files += 1

                continue

            else:

                # Delete zero-byte file
                output_file.unlink()


        # ----------------------------------------------------
        # Download
        # ----------------------------------------------------

        success = download_file(
            url,
            output_file
        )


        # ----------------------------------------------------
        # Save status
        # ----------------------------------------------------

        if success:

            successful_files += 1

            status = "SUCCESS"

        else:

            failed_files += 1

            status = "FAILED"


        download_log.append({
            "UNIQUEID": uniqueid,
            "URL": url,
            "LOCAL_FILE": str(output_file),
            "STATUS": status
        })


        # ----------------------------------------------------
        # Small delay
        # ----------------------------------------------------

        time.sleep(0.5)


# ============================================================
# SAVE DOWNLOAD LOG
# ============================================================

log_df = pd.DataFrame(download_log)

log_df.to_csv(
    LOG_FILE,
    index=False
)


# ============================================================
# FINAL SUMMARY
# ============================================================

print("\n")
print("=" * 80)
print("DOWNLOAD SUMMARY")
print("=" * 80)

print(
    f"Total selected isolates : {total_isolates}"
)

print(
    f"Successful files        : {successful_files}"
)

print(
    f"Skipped existing files  : {skipped_files}"
)

print(
    f"Failed files            : {failed_files}"
)

print("\nFASTQ directory:")
print(OUTPUT_DIR)

print("\nDownload log:")
print(LOG_FILE)

print("\n")


# ============================================================
# LIST DOWNLOADED FILES
# ============================================================

fastq_files = list(
    OUTPUT_DIR.glob("*.fastq.gz")
)

print("=" * 80)
print("FASTQ FILES PRESENT")
print("=" * 80)

print(
    "Number of FASTQ files:",
    len(fastq_files)
)


if len(fastq_files) > 0:

    print("\nFirst 20 files:")

    for file in sorted(fastq_files)[:20]:

        size_mb = file.stat().st_size / (1024 ** 2)

        print(
            f"{file.name:45} "
            f"{size_mb:10.2f} MB"
        )

else:

    print("\nWARNING:")
    print("No FASTQ files were downloaded.")


# ============================================================
# DONE
# ============================================================

print("\n")
print("=" * 80)
print("FASTQ DOWNLOAD PROCESS COMPLETED")
print("=" * 80)