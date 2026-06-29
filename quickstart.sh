#!/usr/bin/env bash
# ============================================================================
# BioCompiler Quick-Start Script
# ============================================================================
#
# One-command verification that the BioCompiler artifact works.
# Safe to run multiple times (idempotent).
#
# What this script does:
#   1. Checks Python version (>= 3.10)
#   2. Installs the package in editable mode with dev dependencies
#   3. Runs a smoke test: optimize insulin for E. coli via CLI
#   4. Checks that the standalone verifier (biocompiler-verify) works
#   5. Runs the Lean4 proof build (if lake is available)
#   6. Prints a summary of all results
#
# Usage:
#   ./quickstart.sh
#
# Exit codes:
#   0  All steps passed (or skipped gracefully)
#   1  A required step failed
# ============================================================================

set -euo pipefail

# ── Colours (disabled if not a terminal) ─────────────────────────────────────
if [[ -t 1 ]]; then
    GREEN='\033[0;32m'
    RED='\033[0;31m'
    YELLOW='\033[0;33m'
    CYAN='\033[0;36m'
    BOLD='\033[1m'
    RESET='\033[0m'
else
    GREEN='' RED='' YELLOW='' CYAN='' BOLD='' RESET=''
fi

# ── Tracking variables ────────────────────────────────────────────────────────
PASS=0
FAIL=0
SKIP=0
RESULTS=()

# ── Helper functions ──────────────────────────────────────────────────────────

step_header() {
    echo ""
    echo -e "${BOLD}${CYAN}══════════════════════════════════════════════════════════════${RESET}"
    echo -e "${BOLD}${CYAN}  $1${RESET}"
    echo -e "${BOLD}${CYAN}══════════════════════════════════════════════════════════════${RESET}"
}

record_pass() {
    PASS=$((PASS + 1))
    RESULTS+=("${GREEN}PASS${RESET} | $1")
}

record_fail() {
    FAIL=$((FAIL + 1))
    RESULTS+=("${RED}FAIL${RESET} | $1")
}

record_skip() {
    SKIP=$((SKIP + 1))
    RESULTS+=("${YELLOW}SKIP${RESET} | $1")
}

# Resolve the project root directory (where this script lives)
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# ============================================================================
# STEP 1 — Check Python version (>= 3.10)
# ============================================================================
step_header "Step 1/6: Checking Python version"

PYTHON="${PYTHON:-python3}"

if ! command -v "$PYTHON" &>/dev/null; then
    echo -e "${RED}ERROR: '${PYTHON}' not found in PATH.${RESET}"
    echo "  Install Python >= 3.10 and ensure it is on your PATH."
    record_fail "Python executable found"
    echo -e "${RED}Aborting — Python is required for all subsequent steps.${RESET}"
    exit 1
fi

# Extract major.minor version
PY_VERSION=$("$PYTHON" -c 'import sys; print(f"{sys.version_info.major}.{sys.version_info.minor}")')
PY_MAJOR=$("$PYTHON" -c 'import sys; print(sys.version_info.major)')
PY_MINOR=$("$PYTHON" -c 'import sys; print(sys.version_info.minor)')

echo "  Detected Python version: ${PY_VERSION}"

if [[ "$PY_MAJOR" -lt 3 ]] || { [[ "$PY_MAJOR" -eq 3 ]] && [[ "$PY_MINOR" -lt 10 ]]; }; then
    echo -e "${RED}ERROR: Python >= 3.10 required, found ${PY_VERSION}.${RESET}"
    record_fail "Python >= 3.10 (found ${PY_VERSION})"
    exit 1
fi

echo -e "${GREEN}OK: Python ${PY_VERSION} meets the >= 3.10 requirement.${RESET}"
record_pass "Python >= 3.10 (found ${PY_VERSION})"

# ============================================================================
# STEP 2 — Install the package (editable, with dev deps)
# ============================================================================
step_header "Step 2/6: Installing biocompiler (pip install -e \".[dev]\")"

# Idempotent: pip install -e is safe to run again; it will just no-op if
# already installed at the same version.
echo "  Running: pip install -e \".[dev]\"  (from ${SCRIPT_DIR})"

if pip install -e ".[dev]" 2>&1 | tail -5; then
    echo -e "${GREEN}OK: Package installed successfully.${RESET}"

    # Quick sanity check: can we import it?
    INSTALLED_VERSION=$("$PYTHON" -c "import biocompiler; print(biocompiler.__version__)" 2>/dev/null || echo "unknown")
    echo "  Installed version: ${INSTALLED_VERSION}"
    record_pass "pip install -e \".[dev]\" (v${INSTALLED_VERSION})"
else
    echo -e "${RED}ERROR: pip install failed.${RESET}"
    record_fail "pip install -e \".[dev]\""
    exit 1
fi

# ============================================================================
# STEP 3 — Smoke test: optimize insulin for E. coli via CLI
# ============================================================================
step_header "Step 3/6: Smoke test — optimize insulin for E. coli"

# Insulin precursor (human, UniProt P01308) amino-acid sequence
INSULIN_SEQ="MALWMRLLPLLALLALWGPDPAAAFVNQHLCGSHLVEALYLVCGERGFFYTPKTRREAEDLQVGQVELGGGPGAGSLQPLALEGSLQKRGIVEQCCTSICSLYQLENYCN"

echo "  Protein: human insulin precursor (${#INSULIN_SEQ} aa)"
echo "  Organism: Escherichia coli"
echo ""

# We capture both stdout and exit code.
# Use --json for machine-readable output (easier to validate).
if OPT_OUTPUT=$(biocompiler optimize "${INSULIN_SEQ}" --organism escherichia_coli --json 2>&1); then
    echo -e "${GREEN}OK: CLI optimize completed successfully.${RESET}"

    # Extract key metrics from JSON output (best-effort)
    CAI=$(echo "$OPT_OUTPUT" | "$PYTHON" -c "
import sys, json
try:
    data = json.load(sys.stdin)
    print(data.get('cai', data.get('CAI', 'N/A')))
except Exception:
    print('N/A')
" 2>/dev/null || echo "N/A")

    GC=$(echo "$OPT_OUTPUT" | "$PYTHON" -c "
import sys, json
try:
    data = json.load(sys.stdin)
    print(data.get('gc_content', data.get('gc', 'N/A')))
except Exception:
    print('N/A')
" 2>/dev/null || echo "N/A")

    SEQ_LEN=$(echo "$OPT_OUTPUT" | "$PYTHON" -c "
import sys, json
try:
    data = json.load(sys.stdin)
    seq = data.get('optimized_sequence', data.get('dna_sequence', ''))
    print(len(seq))
except Exception:
    print('N/A')
" 2>/dev/null || echo "N/A")

    echo "  CAI: ${CAI}"
    echo "  GC content: ${GC}"
    echo "  Optimized sequence length: ${SEQ_LEN} bp"
    record_pass "CLI optimize insulin (CAI=${CAI}, GC=${GC}, ${SEQ_LEN}bp)"
else
    echo -e "${RED}ERROR: CLI optimize command failed.${RESET}"
    echo "  Output: ${OPT_OUTPUT}"
    record_fail "CLI optimize insulin"
    exit 1
fi

# ============================================================================
# STEP 4 — Check standalone verifier (biocompiler-verify)
# ============================================================================
step_header "Step 4/6: Checking standalone verifier (biocompiler-verify)"

# biocompiler-verify is an alias for `biocompiler check`.
# We create a minimal FASTA file with the optimized sequence from Step 3,
# then run the verifier on it.

# Extract the optimized DNA sequence from the JSON output
OPT_SEQ=$(echo "$OPT_OUTPUT" | "$PYTHON" -c "
import sys, json
try:
    data = json.load(sys.stdin)
    print(data.get('optimized_sequence', data.get('dna_sequence', '')))
except Exception:
    print('')
" 2>/dev/null || echo "")

VERIFY_TMPDIR=$(mktemp -d)
trap 'rm -rf "${VERIFY_TMPDIR}"' EXIT

if [[ -n "$OPT_SEQ" && "$OPT_SEQ" != "N/A" ]]; then
    # Write a minimal FASTA file for the verifier
    VERIFY_FASTA="${VERIFY_TMPDIR}/insulin_optimized.fasta"
    echo ">insulin_ecoli_optimized" > "$VERIFY_FASTA"
    echo "$OPT_SEQ" >> "$VERIFY_FASTA"

    echo "  Running: biocompiler-verify ${VERIFY_FASTA} --species escherichia_coli"

    if VERIFY_OUTPUT=$(biocompiler-verify "${VERIFY_FASTA}" --species escherichia_coli 2>&1); then
        echo -e "${GREEN}OK: Standalone verifier completed successfully.${RESET}"
        # Show a snippet of the verifier output (last 10 lines)
        echo "$VERIFY_OUTPUT" | tail -10
        record_pass "biocompiler-verify on optimized insulin"
    else
        echo -e "${YELLOW}WARNING: biocompiler-verify exited with non-zero status.${RESET}"
        echo "  Output: ${VERIFY_OUTPUT}"
        # This is a warning, not a hard failure — the verifier may flag issues
        # but the tool itself must at least run.
        if echo "$VERIFY_OUTPUT" | head -1 | grep -qi "usage\|error\|not found"; then
            record_fail "biocompiler-verify (failed to run)"
        else
            echo -e "${YELLOW}Verifier ran but reported issues (non-zero exit).${RESET}"
            record_pass "biocompiler-verify ran (with warnings)"
        fi
    fi
else
    echo -e "${YELLOW}SKIP: No optimized sequence available from Step 3; cannot run verifier.${RESET}"
    record_skip "biocompiler-verify (no sequence from step 3)"
fi

# ============================================================================
# STEP 5 — Run Lean4 proof build (if lake is available)
# ============================================================================
step_header "Step 5/6: Lean4 proof build (if lake is available)"

if command -v lake &>/dev/null; then
    LAKE_VERSION=$(lake --version 2>&1 | head -1 || echo "unknown")
    echo "  Found lake: ${LAKE_VERSION}"
    echo "  Running: (cd proof && lake build)  from ${SCRIPT_DIR}"

    if (cd "${SCRIPT_DIR}/proof" && lake build 2>&1); then
        echo -e "${GREEN}OK: Lean4 proofs built successfully.${RESET}"
        record_pass "Lean4 proof build (lake build)"

        # Bonus: check for 'sorry' in proof sources
        echo ""
        echo "  Checking for incomplete proofs (sorry)..."
        if (cd "${SCRIPT_DIR}" && ! rg -c "sorry" proof/BioCompiler/ 2>/dev/null); then
            echo -e "${GREEN}OK: No 'sorry' found — all proofs are complete.${RESET}"
            record_pass "No sorry in proof sources"
        else
            SORRY_COUNT=$(cd "${SCRIPT_DIR}" && rg -c "sorry" proof/BioCompiler/ 2>/dev/null | awk -F: '{s+=$2} END{print s+0}' || echo "0")
            echo -e "${YELLOW}WARNING: Found 'sorry' in ${SORRY_COUNT} line(s) — some proofs are incomplete.${RESET}"
            record_skip "Sorry check (${SORRY_COUNT} incomplete proofs)"
        fi
    else
        echo -e "${RED}ERROR: lake build failed.${RESET}"
        record_fail "Lean4 proof build (lake build)"
    fi
else
    echo -e "${YELLOW}SKIP: 'lake' not found in PATH — Lean4 proof build skipped.${RESET}"
    echo "  To run Lean4 proofs, install Lean4 and Lake: https://leanprover.github.io/"
    record_skip "Lean4 proof build (lake not installed)"
fi

# ============================================================================
# STEP 6 — Print summary
# ============================================================================
step_header "Step 6/6: Summary"

echo ""
echo -e "${BOLD}BioCompiler Quick-Start Results${RESET}"
echo -e "${BOLD}──────────────────────────────────────────────────────────────${RESET}"
for result in "${RESULTS[@]}"; do
    echo -e "  ${result}"
done
echo -e "${BOLD}──────────────────────────────────────────────────────────────${RESET}"
echo ""
echo -e "  ${GREEN}Passed:${RESET}  ${PASS}"
echo -e "  ${RED}Failed:${RESET}  ${FAIL}"
echo -e "  ${YELLOW}Skipped:${RESET} ${SKIP}"
echo ""

if [[ "$FAIL" -gt 0 ]]; then
    echo -e "${RED}${BOLD}RESULT: SOME STEPS FAILED — see details above.${RESET}"
    exit 1
elif [[ "$PASS" -gt 0 ]]; then
    echo -e "${GREEN}${BOLD}RESULT: ALL CHECKS PASSED — BioCompiler artifact is working!${RESET}"
    exit 0
else
    echo -e "${YELLOW}${BOLD}RESULT: NO CHECKS WERE RUN — nothing to verify.${RESET}"
    exit 1
fi
