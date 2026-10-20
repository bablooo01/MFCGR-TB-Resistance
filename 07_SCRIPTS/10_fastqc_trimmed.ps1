# ============================================================
# FastQC on trimmed FASTQ files
# ============================================================

$ErrorActionPreference = "Stop"

# ============================================================
# PATHS
# ============================================================

$ProjectDir = "D:\MFCGR_PROJECT"

$FastQCDir = "D:\fastqc_v0.12.1\FastQC"

$InputDir = Join-Path $ProjectDir "08_RESULTS\trimmed_fastq"

$OutputDir = Join-Path $ProjectDir "08_RESULTS\fastqc_trimmed"


# ============================================================
# HEADER
# ============================================================

Write-Host ""
Write-Host "============================================================"
Write-Host "FASTQC - TRIMMED FASTQ"
Write-Host "============================================================"
Write-Host ""


# ============================================================
# CHECK DIRECTORIES
# ============================================================

if (!(Test-Path $FastQCDir)) {
    throw "FastQC directory not found: $FastQCDir"
}

if (!(Test-Path $InputDir)) {
    throw "Trimmed FASTQ directory not found: $InputDir"
}

New-Item -ItemType Directory -Force -Path $OutputDir | Out-Null


# ============================================================
# FIND TRIMMED FASTQ FILES
# ============================================================

$FastqFiles = Get-ChildItem `
    -Path $InputDir `
    -Filter "*.trimmed.fastq.gz" `
    -File |
    Sort-Object Name


Write-Host "Trimmed FASTQ files found: $($FastqFiles.Count)"
Write-Host ""


# ============================================================
# EXPECTED FILE COUNT
# ============================================================

$ExpectedFiles = 240

if ($FastqFiles.Count -ne $ExpectedFiles) {

    Write-Host "WARNING:"
    Write-Host "Expected files : $ExpectedFiles"
    Write-Host "Found files    : $($FastqFiles.Count)"
    Write-Host ""

}


# ============================================================
# MOVE TO FASTQC DIRECTORY
# ============================================================

Set-Location $FastQCDir


# ============================================================
# PROCESS FILES
# ============================================================

$Success = 0
$Failed = 0
$Skipped = 0

$Total = $FastqFiles.Count
$Current = 0


foreach ($FastqFile in $FastqFiles) {

    $Current++

    Write-Host ""
    Write-Host "------------------------------------------------------------"
    Write-Host "Processing $Current / $Total"
    Write-Host "File: $($FastqFile.Name)"
    Write-Host "------------------------------------------------------------"


    # Expected FastQC HTML filename
    $BaseName = $FastqFile.Name -replace "\.fastq\.gz$", ""

    $ExpectedHTML = Join-Path `
        $OutputDir `
        "${BaseName}_fastqc.html"


    # --------------------------------------------------------
    # SKIP IF ALREADY COMPLETED
    # --------------------------------------------------------

    if (Test-Path $ExpectedHTML) {

        Write-Host "Already exists - skipping."

        $Skipped++

        continue
    }


    # --------------------------------------------------------
    # RUN FASTQC
    # --------------------------------------------------------

    try {

        & java `
            "-Dfastqc.output_dir=$OutputDir" `
            "-Xmx2g" `
            "-cp" `
            ".;sam-1.103.jar;jbzip2-0.9.jar" `
            "uk.ac.babraham.FastQC.FastQCApplication" `
            "$($FastqFile.FullName)"


        if ($LASTEXITCODE -eq 0) {

            Write-Host "SUCCESS"

            $Success++

        }
        else {

            Write-Host "FAILED - exit code $LASTEXITCODE"

            $Failed++
        }

    }
    catch {

        Write-Host "FAILED"
        Write-Host $_

        $Failed++
    }
}


# ============================================================
# FINAL SUMMARY
# ============================================================

Write-Host ""
Write-Host "============================================================"
Write-Host "FASTQC TRIMMED COMPLETE"
Write-Host "============================================================"

Write-Host "Total FASTQ files : $Total"
Write-Host "Successful        : $Success"
Write-Host "Skipped           : $Skipped"
Write-Host "Failed            : $Failed"

Write-Host ""
Write-Host "Output directory:"
Write-Host $OutputDir

Write-Host ""
Write-Host "============================================================"