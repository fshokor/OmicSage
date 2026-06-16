#!/bin/bash
# ============================================================
# OmicSage — Docker Run Wrapper
# Convenience script: handles volume mounts, env vars, and GPU.
#
# Usage:
#   ./scripts/run_docker.sh scrna --config config/runs/gse166635_hcc.yaml
#   ./scripts/run_docker.sh spatial --config config/runs/kuppe_heart.yaml --step ingest
#   ./scripts/run_docker.sh cite --config config/runs/gse194122_cite.yaml --from-step normalize
#   ./scripts/run_docker.sh multiome --config config/runs/gse194122_multiome.yaml
#   ./scripts/run_docker.sh dev jupyter
#   ./scripts/run_docker.sh dev pytest tests/ -q
#   ./scripts/run_docker.sh dev bash
#
# GPU mode (requires nvidia-docker):
#   OMICSAGE_GPU=1 ./scripts/run_docker.sh spatial --config ...
# ============================================================

set -e

# ------------------------------------------------------------
# Resolve repo root (script can be called from anywhere)
# ------------------------------------------------------------
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"

# ------------------------------------------------------------
# Validate arguments
# ------------------------------------------------------------
MODALITY="${1:-}"
if [ -z "${MODALITY}" ]; then
    echo ""
    echo "OmicSage Docker Runner"
    echo "======================"
    echo ""
    echo "Usage: $(basename $0) <modality> [pipeline args...]"
    echo ""
    echo "Modalities:"
    echo "  scrna      scRNA-seq pipeline"
    echo "  cite       CITE-seq pipeline"
    echo "  multiome   Multiome (RNA + ATAC) pipeline"
    echo "  spatial    Spatial transcriptomics pipeline"
    echo "  dev        Development environment (Jupyter / pytest / bash)"
    echo ""
    echo "Examples:"
    echo "  $(basename $0) scrna --config config/runs/gse166635_hcc.yaml"
    echo "  $(basename $0) spatial --config config/runs/kuppe_heart.yaml --step ingest"
    echo "  $(basename $0) cite --config config/runs/gse194122_cite.yaml --from-step normalize"
    echo "  $(basename $0) dev jupyter"
    echo "  $(basename $0) dev pytest tests/ -q"
    echo "  $(basename $0) dev bash"
    echo ""
    echo "Environment variables:"
    echo "  OMICSAGE_GPU=1          Enable GPU (requires nvidia-docker)"
    echo "  ANTHROPIC_API_KEY=...   Enable Anthropic Claude AI layer"
    echo "  OPENAI_API_KEY=...      Enable OpenAI AI layer"
    echo "  OLLAMA_HOST=...         Ollama endpoint (default: host.docker.internal:11434)"
    echo ""
    exit 0
fi

shift  # remaining args passed through to pipeline

# ------------------------------------------------------------
# Image selection
# ------------------------------------------------------------
if [ "${MODALITY}" = "dev" ]; then
    IMAGE="omicsage:dev"
    EXTRA_FLAGS="-it -p 8888:8888"
else
    IMAGE="omicsage:latest"
    EXTRA_FLAGS=""
fi

# ------------------------------------------------------------
# GPU support
# ------------------------------------------------------------
GPU_FLAGS=""
if [ "${OMICSAGE_GPU:-0}" = "1" ]; then
    echo "[OmicSage] GPU mode enabled"
    GPU_FLAGS="--gpus all"
fi

# ------------------------------------------------------------
# Volume mounts
# ------------------------------------------------------------
VOLUMES=(
    "-v ${REPO_ROOT}/data:/app/data"
    "-v ${REPO_ROOT}/reports/output:/app/reports/output"
    "-v ${REPO_ROOT}/config:/app/config:ro"
    "-v ${REPO_ROOT}/pipeline:/app/pipeline:ro"
    "-v ${REPO_ROOT}/ai:/app/ai:ro"
)

# Dev image gets full read-write source mount for hot-reload
if [ "${MODALITY}" = "dev" ]; then
    VOLUMES=(
        "-v ${REPO_ROOT}/data:/app/data"
        "-v ${REPO_ROOT}/reports:/app/reports"
        "-v ${REPO_ROOT}/config:/app/config"
        "-v ${REPO_ROOT}/pipeline:/app/pipeline"
        "-v ${REPO_ROOT}/ai:/app/ai"
        "-v ${REPO_ROOT}/tests:/app/tests"
        "-v ${REPO_ROOT}/run_scrna_pipeline.py:/app/run_scrna_pipeline.py"
        "-v ${REPO_ROOT}/run_cite_pipeline.py:/app/run_cite_pipeline.py"
        "-v ${REPO_ROOT}/run_multiome_pipeline.py:/app/run_multiome_pipeline.py"
        "-v ${REPO_ROOT}/run_spatial_pipeline.py:/app/run_spatial_pipeline.py"
    )
fi

# ------------------------------------------------------------
# Environment variables
# ------------------------------------------------------------
ENV_FLAGS=(
    "-e OMICSAGE_GPU=${OMICSAGE_GPU:-0}"
    "-e PYTHONUNBUFFERED=1"
    "-e HDF5_USE_FILE_LOCKING=FALSE"
    "-e OLLAMA_HOST=${OLLAMA_HOST:-http://host.docker.internal:11434}"
)

# Pass API keys only if set in host shell (never hardcode)
if [ -n "${ANTHROPIC_API_KEY:-}" ]; then
    ENV_FLAGS+=("-e ANTHROPIC_API_KEY=${ANTHROPIC_API_KEY}")
fi
if [ -n "${OPENAI_API_KEY:-}" ]; then
    ENV_FLAGS+=("-e OPENAI_API_KEY=${OPENAI_API_KEY}")
fi

# ------------------------------------------------------------
# Build the docker run command
# ------------------------------------------------------------
CMD="docker run --rm \
    ${GPU_FLAGS} \
    ${EXTRA_FLAGS} \
    ${VOLUMES[*]} \
    ${ENV_FLAGS[*]} \
    ${IMAGE} \
    ${MODALITY} ${@}"

echo "[OmicSage] Running: ${MODALITY}"
echo "[OmicSage] Image:   ${IMAGE}"
echo "[OmicSage] Args:    ${@}"
echo ""

# Execute
eval ${CMD}
