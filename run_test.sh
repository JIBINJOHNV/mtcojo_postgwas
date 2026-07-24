#!/usr/bin/env bash
# =============================================================================
#  run_test.sh — End-to-end test using bundled test data
#  Runs mtCOJO on chr1:1Mb-2Mb subset of SCZ and BIP GWAS VCF files.
#  All input data is self-contained inside the package.
# =============================================================================

set -e
PACKAGE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
TEST_DATA="${PACKAGE_DIR}/mtcojo_postgwas/test_data"
OUT_DIR="${PACKAGE_DIR}/test_output"

echo ""
echo "============================================================"
echo "  mtcojo_postgwas — Running Bundled Test"
echo "  Package dir : ${PACKAGE_DIR}"
echo "  Test data   : ${TEST_DATA}"
echo "  Output dir  : ${OUT_DIR}"
echo "============================================================"
echo ""

mkdir -p "${OUT_DIR}"

# ── Resolve manifest (fill in absolute PACKAGE_DIR path) ───────────────────
RESOLVED_MANIFEST="${OUT_DIR}/resolved_manifest.csv"
sed "s|PACKAGE_DIR|${PACKAGE_DIR}|g" \
    "${TEST_DATA}/test_manifest.csv" > "${RESOLVED_MANIFEST}"

echo "[→] Resolved manifest:"
cat "${RESOLVED_MANIFEST}"
echo ""

# ── Detect bcftools ─────────────────────────────────────────────────────────
if command -v bcftools &>/dev/null; then
    BCFTOOLS_BIN="$(command -v bcftools)"
elif [[ -f "${HOME}/miniconda3/bin/bcftools" ]]; then
    BCFTOOLS_BIN="${HOME}/miniconda3/bin/bcftools"
elif [[ -f "${HOME}/miniforge3/envs/mtcojo_postgwas/bin/bcftools" ]]; then
    BCFTOOLS_BIN="${HOME}/miniforge3/envs/mtcojo_postgwas/bin/bcftools"
else
    echo "ERROR: bcftools not found. Activate the conda environment or pass --bcftools."
    exit 1
fi

# ── Detect gcta64 ───────────────────────────────────────────────────────────
if command -v gcta64 &>/dev/null; then
    GCTA64_BIN="$(command -v gcta64)"
elif [[ -f "${HOME}/miniforge3/envs/mtcojo_postgwas/bin/gcta64" ]]; then
    GCTA64_BIN="${HOME}/miniforge3/envs/mtcojo_postgwas/bin/gcta64"
else
    # Try to find gcta64 in common conda environment locations.
    GCTA64_BIN=$(find "${HOME}/miniconda3/envs" "${HOME}/miniforge3/envs" -name "gcta64" -type f 2>/dev/null | head -1)
    if [[ -z "${GCTA64_BIN}" ]]; then
        echo "ERROR: gcta64 not found. Activate the conda environment first:"
        echo "  conda activate mtcojo_postgwas"
        exit 1
    fi
fi

echo "[✓] bcftools : ${BCFTOOLS_BIN}"
echo "[✓] gcta64   : ${GCTA64_BIN}"
echo ""

# ── Run mtcojo-postgwas ─────────────────────────────────────────────────────
mtcojo-postgwas \
    -m "${RESOLVED_MANIFEST}" \
    -d "${OUT_DIR}" \
    -o scz_bip_test \
    -b "${TEST_DATA}/ref_panel_chr1_mini" \
    -l "${TEST_DATA}/ld_scores" \
    --bcftools "${BCFTOOLS_BIN}" \
    --gcta64   "${GCTA64_BIN}" \
    --gwas-thresh 1e-5

echo ""
echo "============================================================"
echo "  Test complete! Output files in: ${OUT_DIR}"
echo "  Key output files:"
echo "    SCZ.ma              — Target summary statistics (GCTA .ma format)"
echo "    BIP.ma              — Covariate summary statistics (GCTA .ma format)"
echo "    scz_bip_test.mtcojo.list — Trait manifest list"
echo "    scz_bip_test.mtcojo.cma  — Final conditional summary statistics"
echo "============================================================"
echo ""
