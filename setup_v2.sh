#!/usr/bin/env bash
set -euo pipefail

echo "=========================================="
echo "  PacBio HDV Pipeline Installation v2.0"
echo "=========================================="

# ---------- helpers ----------
command_exists() { command -v "$1" >/dev/null 2>&1; }

add_path_once() {
  local line='export PATH="$HOME/bin:$PATH"'
  grep -qxF "$line" "$HOME/.bashrc" 2>/dev/null || echo "$line" >> "$HOME/.bashrc"
  export PATH="$HOME/bin:$PATH"
}

# ---------- 1) conda/mamba detection ----------
if command_exists micromamba; then
  CONDA_CMD="micromamba"
elif command_exists mamba; then
  CONDA_CMD="mamba"
elif command_exists conda; then
  CONDA_CMD="conda"
else
  echo "[INFO] No conda/mamba found. Installing Miniconda to \$HOME/miniconda3 ..."
  wget -qO miniconda.sh https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
  bash miniconda.sh -b -p "$HOME/miniconda3"
  rm -f miniconda.sh
  export PATH="$HOME/miniconda3/bin:$PATH"
  CONDA_CMD="conda"
  echo "[INFO] Miniconda installed. Please open a new shell or 'source ~/.bashrc' if needed."
fi

# ---------- 2) Nextflow ----------
if ! command_exists nextflow; then
  echo "[INFO] Installing Nextflow ..."
  curl -s https://get.nextflow.io | bash
  mkdir -p "$HOME/bin"
  mv nextflow "$HOME/bin/"
  add_path_once
fi

# ---------- 3) Main environment (single env for PacBio/ONT) ----------
ENV_FILE="environment_v2.yml"
ENV_NAME="pacbio-hdv-v2"

if [[ ! -f "$ENV_FILE" ]]; then
  echo "[ERROR] $ENV_FILE not found in $(pwd)"
  exit 1
fi

if $CONDA_CMD info --envs | awk '{print $1}' | grep -qx "$ENV_NAME"; then
  echo "[INFO] Updating environment '$ENV_NAME' from $ENV_FILE ..."
  if [[ "$CONDA_CMD" == "micromamba" ]]; then
    $CONDA_CMD env update -n "$ENV_NAME" -f "$ENV_FILE" --prune
  else
    $CONDA_CMD env update -n "$ENV_NAME" -f "$ENV_FILE" --prune
  fi
else
  echo "[INFO] Creating environment '$ENV_NAME' from $ENV_FILE ..."
  $CONDA_CMD env create -f "$ENV_FILE"
fi

# ---------- 4) RVHaplo code ----------
mkdir -p tools
if [[ -d tools/RVHaplo ]]; then
  echo "[INFO] Found local tools/RVHaplo — using that."
else
  echo "[INFO] Cloning RVHaplo into tools/RVHaplo ..."
  (
    cd tools
    git clone https://github.com/dhcai21/RVHaplo.git
  )
fi
chmod +x tools/RVHaplo/rvhaplo.sh || true

echo "------------------------------------------"
echo " Setup complete!"
echo " Activate with:  conda activate $ENV_NAME"
echo " Run QC:         nextflow run qc_only.nf --reads \"test_data/*.fastq.gz\""
echo " Run HDV v2:     nextflow run main_v2.nf --input_dir test_data --outdir results --polisher racon"
echo "------------------------------------------"
