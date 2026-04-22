#!/bin/bash

# =============================================================================
# download_data.sh
# Downloads processed RNA-seq data for the CARA experiment (OSD-120 / GLDS-120)
# from the NASA Open Science Data Repository (OSDR).
#
# Research Question:
#   Do different Arabidopsis thaliana genotypes (Col-0, WS, phyD) show
#   different transcriptional responses to spaceflight microgravity, and
#   does light environment (light vs. dark) modulate these differences?
#
# Data Source:
#   NASA OSDR Study OSD-120 (CARA - Characterizing Arabidopsis Root Attraction)
#   https://osdr.nasa.gov/bio/repo/data/studies/OSD-120
#
# How it works:
#   Rather than hardcoding file names (which break when NASA updates files),
#   this script queries the OSDR API to discover available files dynamically,
#   then downloads only the files relevant to this analysis.
#
# Usage:
#   bash data/download_data.sh
#   Run from the ROOT of the project directory (not from inside data/)
#
# Requirements:
#   Python 3 (checks automatically)
#
# Output:
#   CSV files saved to the data/ directory.
#   Large files (>50MB) are listed in .gitignore and not tracked by git.
# =============================================================================

echo "=============================================="
echo " NASA OSDR Data Downloader - OSD-120 (CARA)"
echo "=============================================="
echo ""

# --- Check Python is available -----------------------------------------------
if ! command -v python3 &> /dev/null && ! command -v python &> /dev/null; then
    echo "[ERROR] Python is required but was not found."
    echo "  Please install Python from https://www.python.org/downloads/"
    exit 1
fi

# Use whichever Python command is available
PYTHON=$(command -v python3 || command -v python)
echo "[INFO] Using Python: $PYTHON"
echo ""

# --- Run the download logic in Python ----------------------------------------
# Python handles JSON parsing and file downloads cleanly and cross-platform.
# The script queries the OSDR API to get real file names and URLs,
# then filters for files relevant to differential expression analysis.

# Resolve data/ dir from bash and pass to Python via environment variable
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
export OSDR_DATA_DIR="$SCRIPT_DIR"
echo "[INFO] Data will be saved to: $SCRIPT_DIR"
echo ""

$PYTHON - <<'PYEOF'
import urllib.request
import urllib.error
import json
import os
import sys

# ── Configuration ─────────────────────────────────────────────────────────────

STUDY_ID    = "120"
GLDS_KEY    = f"OSD-{STUDY_ID}"   # API returns keys as OSD-XXX, not GLDS-XXX
API_URL     = f"https://osdr.nasa.gov/osdr/data/osd/files/{STUDY_ID}"
OSDR_BASE   = "https://osdr.nasa.gov"

# DATA_DIR is passed in from bash as an environment variable.
# This avoids the __file__ issue that occurs when Python runs inline from bash.
DATA_DIR = os.environ.get("OSDR_DATA_DIR")
if not DATA_DIR or not os.path.isdir(DATA_DIR):
    print(f"[ERROR] Could not locate data/ directory (got: {DATA_DIR})")
    print("  Make sure you run this script from the project root:")
    print("  bash data/download_data.sh")
    sys.exit(1)
print(f"[INFO] Saving files to: {DATA_DIR}")

# We want files whose names contain ANY of these keywords.
# This avoids hardcoding exact file names while still filtering out
# irrelevant files (raw FASTQ, alignment files, etc.)
KEYWORDS = [
    "differential_expression",   # main DE results
    "contrasts",                 # table of comparisons made
    "SampleTable",               # sample metadata
]

# ── Step 1: Query the OSDR API ────────────────────────────────────────────────

print(f"Querying NASA OSDR API for study OSD-{STUDY_ID}...")
print(f"  API endpoint: {API_URL}")
print("")

try:
    req = urllib.request.Request(
        API_URL,
        headers={"User-Agent": "Mozilla/5.0 (research download script)"}
    )
    with urllib.request.urlopen(req, timeout=30) as response:
        raw = response.read()
        data = json.loads(raw)
except urllib.error.URLError as e:
    print(f"[ERROR] Could not reach NASA OSDR API: {e}")
    print("  Check your internet connection and try again.")
    sys.exit(1)
except json.JSONDecodeError:
    print("[ERROR] API returned unexpected data (not valid JSON).")
    sys.exit(1)

# ── Step 2: Parse the file listing ───────────────────────────────────────────

studies = data.get("studies", {})
if GLDS_KEY not in studies:
    print(f"[ERROR] Study key {GLDS_KEY} not found in API response.")
    print(f"  Keys found: {list(studies.keys())}")
    sys.exit(1)

all_files = studies[GLDS_KEY].get("study_files", [])
print(f"[INFO] Found {len(all_files)} total files in OSD-{STUDY_ID}.")

# ── Step 3: Filter for relevant files ────────────────────────────────────────

target_files = []
for f in all_files:
    name = f.get("file_name", "")
    url  = f.get("remote_url", "")
    size = f.get("file_size", 0)

    # Only keep CSV files that match our keywords
    if name.endswith(".csv") and any(kw in name for kw in KEYWORDS):
        size_mb = size / (1024 * 1024)
        target_files.append({
            "name":    name,
            "url":     OSDR_BASE + url,
            "size_mb": size_mb,
        })

if not target_files:
    print("[ERROR] No matching files found. The API structure may have changed.")
    print("  Visit https://osdr.nasa.gov/bio/repo/data/studies/OSD-120 to download manually.")
    sys.exit(1)

print(f"[INFO] Identified {len(target_files)} files to download:\n")
for f in target_files:
    flag = " <- preferred DE file" if "rRNArm" in f["name"] else ""
    print(f"  * {f['name']} ({f['size_mb']:.1f} MB){flag}")
print("")

# ── Step 4: Download each file ────────────────────────────────────────────────

def download_file(name, url, size_mb):
    out_path = os.path.join(DATA_DIR, name)

    if os.path.isfile(out_path):
        print(f"[SKIP] Already exists: {name}")
        return True

    print(f"[DOWNLOAD] {name} ({size_mb:.1f} MB)...")

    def progress(block_count, block_size, total_size):
        if total_size > 0:
            pct = min(block_count * block_size / total_size * 100, 100)
            print(f"\r  Progress: {pct:.1f}%", end="", flush=True)

    try:
        tmp_path = out_path + ".tmp"
        urllib.request.urlretrieve(url, tmp_path, reporthook=progress)
        os.rename(tmp_path, out_path)
        print(f"\r[OK] Saved: {name}        ")
        return True
    except Exception as e:
        print(f"\r[ERROR] Failed: {name}: {e}")
        if os.path.exists(tmp_path):
            os.remove(tmp_path)
        return False

errors = []
for f in target_files:
    success = download_file(f["name"], f["url"], f["size_mb"])
    if not success:
        errors.append(f["name"])

# ── Step 5: Summary ───────────────────────────────────────────────────────────

print("")
print("==============================================")
if not errors:
    print(" All files downloaded successfully!")
    print(" Check your data/ folder.")
    print("")
    print(" NOTE: Large files (>50 MB) are listed in")
    print(" .gitignore and will not be pushed to GitHub.")
    print(" The download script makes your project")
    print(" reproducible without committing large files.")
else:
    print(" Download complete with some errors.")
    print(" The following files failed:")
    for e in errors:
        print(f"   - {e}")
    print("")
    print(" Try running again, or download manually from:")
    print(" https://osdr.nasa.gov/bio/repo/data/studies/OSD-120")
print("==============================================")

PYEOF
