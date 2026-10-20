import os
import time
import hashlib
import pandas as pd
import requests

# ============================================================
# MFCGR PROJECT - FASTQ VERIFICATION + RETRY DOWNLOADER
# ============================================================

BASE = r"D:\MFCGR_PROJECT"

INPUT_CSV = os.path.join(
    BASE, "08_RESULTS", "selected_120_isolates.csv"
)

FASTQ_DIR = os.path.join(
    BASE, "03_GENOMES"
)

RESULTS_DIR = os.path.join(
    BASE, "08_RESULTS"
)

FILE_STATUS_CSV = os.path.join(
    RESULTS_DIR, "fastq_retry_status.csv"
)

ISOLATE_STATUS_CSV = os.path.join(
    RESULTS_DIR, "fastq_isolate_completeness.csv"
)

os.makedirs(FASTQ_DIR, exist_ok=True)
os.makedirs(RESULTS_DIR, exist_ok=True)


# ------------------------------------------------------------
# SETTINGS
# ------------------------------------------------------------

MAX_RETRIES = 5
CHUNK_SIZE = 1024 * 1024       # 1 MB
TIMEOUT = (30, 300)            # connect, read
RETRY_WAIT = 10


# ------------------------------------------------------------
# URL CONVERSION
# ------------------------------------------------------------

def make_https(url):

    url = str(url).strip()

    if not url:
        return None

    if url.startswith("ftp://"):
        return url.replace("ftp://", "https://", 1)

    if url.startswith("ftp.sra.ebi.ac.uk/"):
        return "https://" + url

    if url.startswith("https://"):
        return url

    if url.startswith("http://"):
        return url.replace("http://", "https://", 1)

    return url


# ------------------------------------------------------------
# MD5
# ------------------------------------------------------------

def calculate_md5(filepath):

    md5 = hashlib.md5()

    with open(filepath, "rb") as f:

        while True:

            chunk = f.read(8 * 1024 * 1024)

            if not chunk:
                break

            md5.update(chunk)

    return md5.hexdigest()


# ------------------------------------------------------------
# PARSE SEMICOLON-SEPARATED VALUES
# ------------------------------------------------------------

def split_values(value):

    if pd.isna(value):
        return []

    value = str(value).strip()

    if not value:
        return []

    return [x.strip() for x in value.split(";") if x.strip()]


# ------------------------------------------------------------
# DOWNLOAD ONE FILE
# ------------------------------------------------------------

def download_file(url, output_path, expected_size=None, expected_md5=None):

    url = make_https(url)

    filename = os.path.basename(output_path)

    # --------------------------------------------------------
    # CHECK EXISTING FILE
    # --------------------------------------------------------

    if os.path.exists(output_path):

        current_size = os.path.getsize(output_path)

        # Exact expected size
        if expected_size is not None:

            try:
                expected_size_int = int(float(expected_size))

                if current_size == expected_size_int:

                    if expected_md5:

                        print(f"Checking MD5: {filename}")

                        actual_md5 = calculate_md5(output_path)

                        if actual_md5.lower() == expected_md5.lower():

                            print(f"VALID: {filename}")
                            return True

                        else:

                            print(
                                f"MD5 mismatch: {filename}"
                            )

                            os.remove(output_path)

                    else:

                        print(f"VALID BY SIZE: {filename}")
                        return True

            except Exception:
                pass

        # If no metadata, at least don't redownload
        if expected_size is None and current_size > 0:

            print(f"EXISTING FILE: {filename}")
            return True

    # --------------------------------------------------------
    # RETRIES
    # --------------------------------------------------------

    for attempt in range(1, MAX_RETRIES + 1):

        print()
        print(
            f"Attempt {attempt}/{MAX_RETRIES}: {filename}"
        )

        try:

            # ------------------------------------------------
            # RESUME FROM EXISTING PARTIAL FILE
            # ------------------------------------------------

            start_byte = 0

            if os.path.exists(output_path):

                start_byte = os.path.getsize(output_path)

            headers = {}

            if start_byte > 0:

                headers["Range"] = f"bytes={start_byte}-"

                print(
                    f"Resuming from {start_byte / (1024**2):.2f} MB"
                )

            # ------------------------------------------------
            # REQUEST
            # ------------------------------------------------

            response = requests.get(
                url,
                headers=headers,
                stream=True,
                timeout=TIMEOUT
            )

            # ------------------------------------------------
            # SERVER RESPONSE
            # ------------------------------------------------

            if response.status_code == 416:

                print("Range not accepted; restarting file.")

                if os.path.exists(output_path):
                    os.remove(output_path)

                continue

            response.raise_for_status()

            # ------------------------------------------------
            # DETERMINE APPEND / WRITE MODE
            # ------------------------------------------------

            if start_byte > 0 and response.status_code == 206:

                mode = "ab"

            else:

                # Server ignored Range request.
                # Start again from zero.
                if start_byte > 0:

                    print(
                        "Server did not provide partial content."
                    )

                    print("Restarting from beginning.")

                mode = "wb"
                start_byte = 0

            # ------------------------------------------------
            # DOWNLOAD
            # ------------------------------------------------

            downloaded = start_byte

            with open(output_path, mode) as f:

                for chunk in response.iter_content(
                    chunk_size=CHUNK_SIZE
                ):

                    if chunk:

                        f.write(chunk)

                        downloaded += len(chunk)

                        if expected_size:

                            try:

                                percent = (
                                    downloaded /
                                    int(float(expected_size))
                                ) * 100

                                print(
                                    f"\rProgress: {percent:.2f}% "
                                    f"({downloaded/(1024**2):.1f} MB)",
                                    end=""
                                )

                            except Exception:
                                pass

            print()

            # ------------------------------------------------
            # VERIFY SIZE
            # ------------------------------------------------

            final_size = os.path.getsize(output_path)

            if expected_size is not None:

                try:

                    expected_size_int = int(float(expected_size))

                    if final_size != expected_size_int:

                        print(
                            f"SIZE MISMATCH: "
                            f"{final_size} vs "
                            f"{expected_size_int}"
                        )

                        print("Retrying...")

                        time.sleep(RETRY_WAIT)
                        continue

                except Exception:
                    pass

            # ------------------------------------------------
            # VERIFY MD5
            # ------------------------------------------------

            if expected_md5:

                print("Calculating MD5...")

                actual_md5 = calculate_md5(output_path)

                if actual_md5.lower() != expected_md5.lower():

                    print("MD5 MISMATCH!")
                    print("Retrying...")

                    try:
                        os.remove(output_path)
                    except:
                        pass

                    time.sleep(RETRY_WAIT)
                    continue

                print("MD5 VERIFIED.")

            print(
                f"SUCCESS: {filename} "
                f"({final_size/(1024**2):.2f} MB)"
            )

            return True

        except Exception as e:

            print()
            print(f"FAILED: {filename}")
            print(f"Error: {e}")

            if attempt < MAX_RETRIES:

                wait = RETRY_WAIT * attempt

                print(
                    f"Waiting {wait} seconds before retry..."
                )

                time.sleep(wait)

    return False


# ------------------------------------------------------------
# LOAD SELECTED 120
# ------------------------------------------------------------

print("=" * 70)
print("MFCGR FASTQ VERIFICATION + RETRY")
print("=" * 70)

print()
print("Reading:")
print(INPUT_CSV)

df = pd.read_csv(INPUT_CSV)

print()
print(f"Selected isolates: {len(df)}")


# ------------------------------------------------------------
# CHECK REQUIRED COLUMN
# ------------------------------------------------------------

if "fastq_ftp" not in df.columns:

    raise ValueError(
        "Column 'fastq_ftp' was not found "
        "in selected_120_isolates.csv"
    )


# ------------------------------------------------------------
# FILE-LEVEL STATUS
# ------------------------------------------------------------

file_results = []

total_files = 0
valid_before = 0
downloaded_now = 0
failed = 0


# ------------------------------------------------------------
# PROCESS EACH ISOLATE
# ------------------------------------------------------------

for index, row in df.iterrows():

    isolate_number = index + 1

    uniqueid = row.get(
        "UNIQUEID",
        f"ISOLATE_{isolate_number}"
    )

    category = row.get(
        "CATEGORY",
        ""
    )

    urls = split_values(
        row["fastq_ftp"]
    )

    sizes = split_values(
        row.get("fastq_bytes", "")
    )

    md5s = split_values(
        row.get("fastq_md5", "")
    )

    print()
    print("=" * 70)
    print(
        f"Isolate {isolate_number} / {len(df)}"
    )
    print(f"UNIQUEID: {uniqueid}")
    print(f"CATEGORY: {category}")
    print(f"FASTQ files: {len(urls)}")
    print("=" * 70)

    for file_index, original_url in enumerate(urls):

        total_files += 1

        url = make_https(original_url)

        filename = os.path.basename(
            url.split("?")[0]
        )

        output_path = os.path.join(
            FASTQ_DIR,
            filename
        )

        expected_size = None

        if file_index < len(sizes):

            try:
                expected_size = int(
                    float(sizes[file_index])
                )
            except:
                expected_size = None

        expected_md5 = None

        if file_index < len(md5s):

            expected_md5 = md5s[file_index]

        print()
        print("-" * 70)
        print(f"File: {filename}")
        print(f"URL: {url}")

        # ----------------------------------------------------
        # VERIFY EXISTING FILE
        # ----------------------------------------------------

        existing_valid = False

        if os.path.exists(output_path):

            current_size = os.path.getsize(
                output_path
            )

            if expected_size is not None:

                if current_size == expected_size:

                    if expected_md5:

                        print(
                            "Existing file has expected size."
                        )

                        print("Checking MD5...")

                        actual_md5 = calculate_md5(
                            output_path
                        )

                        if (
                            actual_md5.lower()
                            ==
                            expected_md5.lower()
                        ):

                            existing_valid = True

                    else:

                        existing_valid = True

            else:

                if current_size > 0:

                    existing_valid = True

        if existing_valid:

            print(
                f"ALREADY VALID: {filename}"
            )

            valid_before += 1

            file_results.append({
                "UNIQUEID": uniqueid,
                "CATEGORY": category,
                "FILENAME": filename,
                "URL": url,
                "EXPECTED_SIZE": expected_size,
                "ACTUAL_SIZE": os.path.getsize(
                    output_path
                ),
                "EXPECTED_MD5": expected_md5,
                "STATUS": "VALID_EXISTING"
            })

            continue

        # ----------------------------------------------------
        # DOWNLOAD / RETRY
        # ----------------------------------------------------

        success = download_file(
            url,
            output_path,
            expected_size,
            expected_md5
        )

        if success:

            downloaded_now += 1
            status = "DOWNLOADED"

        else:

            failed += 1
            status = "FAILED"

        actual_size = None

        if os.path.exists(output_path):

            actual_size = os.path.getsize(
                output_path
            )

        file_results.append({
            "UNIQUEID": uniqueid,
            "CATEGORY": category,
            "FILENAME": filename,
            "URL": url,
            "EXPECTED_SIZE": expected_size,
            "ACTUAL_SIZE": actual_size,
            "EXPECTED_MD5": expected_md5,
            "STATUS": status
        })


# ------------------------------------------------------------
# SAVE FILE STATUS
# ------------------------------------------------------------

status_df = pd.DataFrame(file_results)

status_df.to_csv(
    FILE_STATUS_CSV,
    index=False
)


# ------------------------------------------------------------
# BUILD ISOLATE COMPLETENESS REPORT
# ------------------------------------------------------------

isolate_results = []

for uniqueid, group in status_df.groupby(
    "UNIQUEID",
    sort=False
):

    total = len(group)

    valid = group[
        group["STATUS"].isin(
            ["VALID_EXISTING", "DOWNLOADED"]
        )
    ].shape[0]

    failed_count = total - valid

    isolate_row = group.iloc[0]

    if failed_count == 0:

        status = "COMPLETE"

    else:

        status = "INCOMPLETE"

    isolate_results.append({
        "UNIQUEID": uniqueid,
        "CATEGORY": isolate_row["CATEGORY"],
        "EXPECTED_FASTQ_FILES": total,
        "VALID_FASTQ_FILES": valid,
        "FAILED_FASTQ_FILES": failed_count,
        "STATUS": status
    })


isolate_df = pd.DataFrame(
    isolate_results
)

isolate_df.to_csv(
    ISOLATE_STATUS_CSV,
    index=False
)


# ------------------------------------------------------------
# FINAL SUMMARY
# ------------------------------------------------------------

complete_isolates = (
    isolate_df["STATUS"] == "COMPLETE"
).sum()

incomplete_isolates = (
    isolate_df["STATUS"] == "INCOMPLETE"
).sum()


print()
print()
print("=" * 70)
print("FINAL FASTQ RETRY SUMMARY")
print("=" * 70)

print(
    f"Selected isolates       : {len(df)}"
)

print(
    f"Expected FASTQ files    : {total_files}"
)

print(
    f"Already valid           : {valid_before}"
)

print(
    f"Downloaded in retry     : {downloaded_now}"
)

print(
    f"Failed after retries    : {failed}"
)

print()
print(
    f"COMPLETE isolates       : {complete_isolates}"
)

print(
    f"INCOMPLETE isolates     : {incomplete_isolates}"
)

print()
print("FASTQ directory:")
print(FASTQ_DIR)

print()
print("File status:")
print(FILE_STATUS_CSV)

print()
print("Isolate completeness:")
print(ISOLATE_STATUS_CSV)

print()
print("=" * 70)
print("RETRY PROCESS COMPLETED")
print("=" * 70)