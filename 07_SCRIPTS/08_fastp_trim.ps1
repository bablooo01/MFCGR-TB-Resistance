# ============================================================
# MFCGR PROJECT
# FASTQ Adapter Trimming and Quality Filtering
# Tool: fastp v1.3.3
# Input: 120 paired-end isolates / 240 FASTQ files
# ============================================================

$ErrorActionPreference = "Stop"

# -----------------------------
# PATHS
# -----------------------------

$fastp = "D:\fastp-1.3.3-windows-ucrt64\fastp-1.3.3-windows-ucrt64\fastp.exe"

$rawDir = "D:\MFCGR_PROJECT\03_GENOMES"

$trimmedDir = "D:\MFCGR_PROJECT\08_RESULTS\trimmed_fastq"

$reportDir = "D:\MFCGR_PROJECT\08_RESULTS\fastp_reports"

$logFile = "D:\MFCGR_PROJECT\08_RESULTS\fastp_batch_log.txt"


# -----------------------------
# CHECK FASTP
# -----------------------------

if (-not (Test-Path $fastp)) {
    Write-Host "ERROR: fastp.exe not found:" -ForegroundColor Red
    Write-Host $fastp
    exit 1
}

Write-Host ""
Write-Host "MFCGR FASTP TRIMMING PIPELINE" -ForegroundColor Cyan
Write-Host "========================================" -ForegroundColor Cyan

$version = & $fastp --version
Write-Host "Using: $version" -ForegroundColor Green


# -----------------------------
# CREATE OUTPUT DIRECTORIES
# -----------------------------

New-Item -ItemType Directory -Force $trimmedDir | Out-Null
New-Item -ItemType Directory -Force $reportDir | Out-Null


# -----------------------------
# INITIALIZE LOG
# -----------------------------

"============================================================" | Out-File $logFile
"MFCGR FASTP BATCH TRIMMING" | Out-File $logFile -Append
"Started: $(Get-Date)" | Out-File $logFile -Append
"============================================================" | Out-File $logFile -Append


# -----------------------------
# FIND R1 FILES
# -----------------------------

$r1Files = Get-ChildItem $rawDir -Filter "*_1.fastq.gz" | Sort-Object Name

$total = $r1Files.Count

Write-Host ""
Write-Host "R1 files found: $total" -ForegroundColor Yellow

if ($total -ne 120) {
    Write-Host "WARNING: Expected 120 R1 files, but found $total." -ForegroundColor Yellow
}


# -----------------------------
# COUNTERS
# -----------------------------

$processed = 0
$success = 0
$skipped = 0
$failed = 0


# -----------------------------
# PROCESS EACH ISOLATE
# -----------------------------

foreach ($r1 in $r1Files) {

    $processed++

    # Isolate ID
    $id = $r1.Name -replace "_1\.fastq\.gz$", ""

    # Input R2
    $r2 = Join-Path $rawDir ($id + "_2.fastq.gz")

    # Output files
    $out1 = Join-Path $trimmedDir ($id + "_1.trimmed.fastq.gz")
    $out2 = Join-Path $trimmedDir ($id + "_2.trimmed.fastq.gz")

    # Reports
    $html = Join-Path $reportDir ($id + ".html")
    $json = Join-Path $reportDir ($id + ".json")


    Write-Host ""
    Write-Host "------------------------------------------------------------" -ForegroundColor DarkGray
    Write-Host "[$processed / $total] $id" -ForegroundColor Cyan
    Write-Host "------------------------------------------------------------"


    # -------------------------
    # CHECK R2
    # -------------------------

    if (-not (Test-Path $r2)) {

        Write-Host "FAILED: R2 file missing" -ForegroundColor Red

        "FAILED`t$id`tR2 file missing" | Out-File $logFile -Append

        $failed++
        continue
    }


    # -------------------------
    # SKIP IF COMPLETE OUTPUT EXISTS
    # -------------------------

    if (
        (Test-Path $out1) -and
        (Test-Path $out2) -and
        (Test-Path $html) -and
        (Test-Path $json)
    ) {

        Write-Host "SKIPPED: Output already exists" -ForegroundColor Yellow

        "SKIPPED`t$id`tOutput already exists" | Out-File $logFile -Append

        $skipped++
        continue
    }


    # -------------------------
    # RUN FASTP
    # -------------------------

    Write-Host "Running fastp..." -ForegroundColor White

    try {

        & $fastp `
            -i $r1.FullName `
            -I $r2 `
            -o $out1 `
            -O $out2 `
            -h $html `
            -j $json `
            -w 4

        if (
            ($LASTEXITCODE -eq 0) -and
            (Test-Path $out1) -and
            (Test-Path $out2) -and
            (Test-Path $html) -and
            (Test-Path $json)
        ) {

            Write-Host "SUCCESS: $id" -ForegroundColor Green

            "SUCCESS`t$id" | Out-File $logFile -Append

            $success++
        }
        else {

            Write-Host "FAILED: fastp did not produce all expected outputs" -ForegroundColor Red

            "FAILED`t$id`tMissing output or non-zero exit code" | Out-File $logFile -Append

            $failed++
        }

    }
    catch {

        Write-Host "FAILED: $id" -ForegroundColor Red
        Write-Host $_.Exception.Message -ForegroundColor Red

        "FAILED`t$id`t$($_.Exception.Message)" | Out-File $logFile -Append

        $failed++
    }
}


# -----------------------------
# FINAL SUMMARY
# -----------------------------

Write-Host ""
Write-Host "============================================================" -ForegroundColor Cyan
Write-Host "FASTP BATCH TRIMMING COMPLETE" -ForegroundColor Cyan
Write-Host "============================================================" -ForegroundColor Cyan

Write-Host "Total isolates : $total"
Write-Host "Successful     : $success" -ForegroundColor Green
Write-Host "Skipped        : $skipped" -ForegroundColor Yellow
Write-Host "Failed         : $failed" -ForegroundColor Red

Write-Host ""
Write-Host "Trimmed FASTQ directory:"
Write-Host $trimmedDir

Write-Host ""
Write-Host "fastp reports directory:"
Write-Host $reportDir

Write-Host ""
Write-Host "Log file:"
Write-Host $logFile

Write-Host ""
Write-Host "Finished: $(Get-Date)"

"============================================================" | Out-File $logFile -Append
"Finished: $(Get-Date)" | Out-File $logFile -Append
"Total: $total | Successful: $success | Skipped: $skipped | Failed: $failed" | Out-File $logFile -Append