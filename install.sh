#!/usr/bin/env bash
# =============================================================================
#  install.sh — Auto-installer for mtcojo_postgwas
#
#  Steps performed:
#    1. Detect OS (macOS or Linux/Ubuntu)
#    2. Install Miniforge3 (mamba) if neither mamba nor conda is found
#    3. Create/update conda environment (includes bcftools, gcta, Python deps)
#    4. Install mtcojo_postgwas Python package in editable mode
#    5. Install CBIIT/ldsc (Python 3 fork) via pip + clone scripts
#    6. Verify all tools
# =============================================================================

set -euo pipefail

PACKAGE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ENV_NAME="mtcojo_postgwas"
ENV_FILE="${PACKAGE_DIR}/environment.yml"
LDSC_DIR="${PACKAGE_DIR}/tools/ldsc"          # where ldsc scripts are cloned

# Colour helpers
RED='\033[0;31m'; GRN='\033[0;32m'; YLW='\033[1;33m'
CYN='\033[0;36m'; BLD='\033[1m'; RST='\033[0m'

info()  { echo -e "${CYN}[→]${RST} $*"; }
pass()  { echo -e "${GRN}[✓]${RST} $*"; }
warn()  { echo -e "${YLW}[!]${RST} $*"; }
error() { echo -e "${RED}[✗]${RST} $*" >&2; exit 1; }

echo ""
echo -e "${BLD}============================================================${RST}"
echo -e "${BLD}  mtcojo_postgwas Installer${RST}"
echo -e "${BLD}============================================================${RST}"
echo ""

# ─────────────────────────────────────────────────────────────────────────────
# 1. Detect operating system
# ─────────────────────────────────────────────────────────────────────────────
OS="$(uname -s)"
ARCH="$(uname -m)"
info "Detected OS : ${OS} / ${ARCH}"

if [[ "$OS" == "Darwin" ]]; then
    MINIFORGE_URL="https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-MacOSX-${ARCH}.sh"
    MINIFORGE_SCRIPT="Miniforge3-MacOSX-${ARCH}.sh"
elif [[ "$OS" == "Linux" ]]; then
    MINIFORGE_URL="https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-${ARCH}.sh"
    MINIFORGE_SCRIPT="Miniforge3-Linux-${ARCH}.sh"
else
    error "Unsupported operating system: $OS  (supported: Darwin, Linux)"
fi

# ─────────────────────────────────────────────────────────────────────────────
# 2. Locate mamba/conda; install Miniforge3 if neither found
# ─────────────────────────────────────────────────────────────────────────────
if command -v mamba &>/dev/null; then
    CONDA_BIN="mamba"
    pass "Found mamba: $(command -v mamba)"
elif command -v conda &>/dev/null; then
    CONDA_BIN="conda"
    pass "Found conda: $(command -v conda)"
else
    warn "Neither mamba nor conda found — installing Miniforge3..."
    curl -fsSL "${MINIFORGE_URL}" -o "/tmp/${MINIFORGE_SCRIPT}"
    bash "/tmp/${MINIFORGE_SCRIPT}" -b -p "${HOME}/miniforge3"
    export PATH="${HOME}/miniforge3/bin:${PATH}"
    CONDA_BIN="mamba"
    pass "Miniforge3 installed at ${HOME}/miniforge3"
fi

# Upgrade to mamba for faster solves if only conda was found
if [[ "$CONDA_BIN" == "conda" ]]; then
    info "Installing mamba into base environment for faster dependency resolution..."
    conda install -n base -c conda-forge mamba -y
    CONDA_BIN="mamba"
fi

# ─────────────────────────────────────────────────────────────────────────────
# 3. Create / update conda environment
# ─────────────────────────────────────────────────────────────────────────────
echo ""
echo -e "${BLD}[Step 3/6] Conda environment: ${ENV_NAME}${RST}"
if conda env list | grep -q "^${ENV_NAME}"; then
    info "Environment '${ENV_NAME}' exists — updating..."
    ${CONDA_BIN} env update -n "${ENV_NAME}" -f "${ENV_FILE}" --prune
else
    info "Creating environment '${ENV_NAME}' (bcftools, gcta, Python, ldsc deps)..."
    ${CONDA_BIN} env create -f "${ENV_FILE}"
fi
pass "Conda environment ready: ${ENV_NAME}"

# Resolve the environment's Python/pip paths
CONDA_PREFIX=$(conda run -n "${ENV_NAME}" python -c "import sys; print(sys.prefix)")
ENV_PYTHON="${CONDA_PREFIX}/bin/python"
ENV_PIP="${CONDA_PREFIX}/bin/pip"

# ─────────────────────────────────────────────────────────────────────────────
# 4. Install mtcojo_postgwas Python package in editable mode
# ─────────────────────────────────────────────────────────────────────────────
echo ""
echo -e "${BLD}[Step 4/6] Install mtcojo_postgwas package${RST}"
info "Installing mtcojo_postgwas into '${ENV_NAME}'..."
"${ENV_PIP}" install -e "${PACKAGE_DIR}"
pass "mtcojo_postgwas installed"

# ─────────────────────────────────────────────────────────────────────────────
# 5. Install CBIIT/ldsc (Python 3 compatible fork)
#
#    Strategy:
#      a) Install 'ldsc' PyPI package (provides Python library + pip entrypoints)
#      b) Clone https://github.com/CBIIT/ldsc into tools/ldsc for the
#         standalone scripts (munge_sumstats.py, ldsc.py) that the pipeline
#         calls directly via subprocess.
# ─────────────────────────────────────────────────────────────────────────────
echo ""
echo -e "${BLD}[Step 5/6] Install CBIIT/ldsc (Python 3 fork)${RST}"

# 5a. PyPI package (ldsc Python library + deps already in environment.yml,
#     but install explicitly to be sure)
info "Installing ldsc via pip..."
"${ENV_PIP}" install ldsc
pass "ldsc Python package installed"

info "Ensuring Python report plotting packages are installed..."
conda run -n "${ENV_NAME}" python -c "import matplotlib, matplotlib_venn" >/dev/null 2>&1 || \
    "${ENV_PIP}" install matplotlib matplotlib-venn
pass "Python report plotting packages available"

# 5b. Clone the CBIIT/ldsc repository for standalone scripts
mkdir -p "${PACKAGE_DIR}/tools"
if [[ -d "${LDSC_DIR}/.git" ]]; then
    info "CBIIT/ldsc repo already cloned — pulling latest..."
    git -C "${LDSC_DIR}" pull --ff-only
else
    info "Cloning https://github.com/CBIIT/ldsc → ${LDSC_DIR}"
    git clone https://github.com/CBIIT/ldsc.git "${LDSC_DIR}"
fi

# Make scripts executable
chmod +x "${LDSC_DIR}/ldsc.py" "${LDSC_DIR}/munge_sumstats.py" 2>/dev/null || true
pass "CBIIT/ldsc scripts ready: ${LDSC_DIR}"

# 5c. Install R plotting packages (CMplot, rMVP) if Rscript is available
if command -v Rscript &>/dev/null; then
    info "Rscript found — installing R packages (data.table, CMplot, rMVP)..."
    conda run -n "${ENV_NAME}" Rscript -e "if (!requireNamespace('data.table', quietly=TRUE)) install.packages('data.table', repos='https://cloud.r-project.org')" || true
    conda run -n "${ENV_NAME}" Rscript -e "if (!requireNamespace('CMplot', quietly=TRUE)) install.packages('CMplot', repos='https://cloud.r-project.org')" || true
    conda run -n "${ENV_NAME}" Rscript -e "if (!requireNamespace('rMVP', quietly=TRUE)) install.packages('rMVP', repos='https://cloud.r-project.org')" || true
fi

# Write the ldsc script paths to a config file that the pipeline reads
LDSC_CFG="${PACKAGE_DIR}/tools/ldsc_paths.txt"
cat > "${LDSC_CFG}" << EOF2
LDSC_DIR=${LDSC_DIR}
LDSC_SCRIPT=${LDSC_DIR}/ldsc.py
MUNGE_SCRIPT=${LDSC_DIR}/munge_sumstats.py
LDSC_PYTHON=${ENV_PYTHON}
EOF2
pass "ldsc paths written to: ${LDSC_CFG}"

# ─────────────────────────────────────────────────────────────────────────────
# 6. Verify installation
# ─────────────────────────────────────────────────────────────────────────────
echo ""
echo -e "${BLD}[Step 6/6] Verifying tools${RST}"

# mtcojo-postgwas CLI
conda run -n "${ENV_NAME}" mtcojo-postgwas --help > /dev/null 2>&1 && \
    pass "mtcojo-postgwas CLI available" || \
    warn "mtcojo-postgwas CLI check failed — check pip install"

# bcftools
BCFTOOLS_VER=$(conda run -n "${ENV_NAME}" bcftools --version | head -1)
pass "bcftools : ${BCFTOOLS_VER}"

# gcta64
GCTA_VER=$(conda run -n "${ENV_NAME}" gcta64 --version 2>&1 | head -1 || echo "gcta64 OK")
pass "gcta64   : ${GCTA_VER}"

# ldsc scripts
if "${ENV_PYTHON}" "${LDSC_DIR}/munge_sumstats.py" --help > /dev/null 2>&1; then
    pass "munge_sumstats.py : OK"
else
    warn "munge_sumstats.py check failed — check ${LDSC_DIR}"
fi
if "${ENV_PYTHON}" "${LDSC_DIR}/ldsc.py" --help > /dev/null 2>&1; then
    pass "ldsc.py           : OK"
else
    warn "ldsc.py check failed — check ${LDSC_DIR}"
fi

if conda run -n "${ENV_NAME}" python -c "import matplotlib, matplotlib_venn" > /dev/null 2>&1; then
    pass "Python plotting  : matplotlib + matplotlib-venn OK"
else
    warn "Python plotting packages missing — run: conda install -n ${ENV_NAME} matplotlib matplotlib-venn"
fi

echo ""
echo -e "${BLD}============================================================${RST}"
echo -e "${GRN}${BLD}  Installation complete!${RST}"
echo ""
echo -e "  Activate environment:"
echo -e "    ${CYN}conda activate ${ENV_NAME}${RST}"
echo ""
echo -e "  Run bundled test:"
echo -e "    ${CYN}bash ${PACKAGE_DIR}/run_test.sh${RST}"
echo ""
echo -e "  ldsc scripts location:"
echo -e "    ${CYN}${LDSC_DIR}/${RST}"
echo -e "${BLD}============================================================${RST}"
echo ""
