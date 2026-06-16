#!/bin/bash
# ============================================================
# OmicSage — Container Entrypoint
# Routes to the correct pipeline runner based on:
#   1. OMICSAGE_MODALITY env var (set by docker-compose service)
#   2. First argument if it's a known modality name
#   3. Special commands: jupyter, pytest, bash, shell
# ============================================================

set -e

# Map modality names to runner scripts
declare -A RUNNERS=(
    ["scrna"]="run_scrna_pipeline.py"
    ["cite"]="run_cite_pipeline.py"
    ["multiome"]="run_multiome_pipeline.py"
    ["spatial"]="run_spatial_pipeline.py"
)

PYTHON="/opt/conda/envs/omicsage/bin/python"

# ------------------------------------------------------------
# GPU awareness
# ------------------------------------------------------------
if [ "${OMICSAGE_GPU}" = "1" ]; then
    echo "[OmicSage] GPU mode enabled"
else
    # Prevent scVI/PyTorch from consuming GPU if not wanted
    export CUDA_VISIBLE_DEVICES=""
fi

# ------------------------------------------------------------
# Command routing
# ------------------------------------------------------------
CMD_ARG="${1:-}"

case "${CMD_ARG}" in

    # Launch Jupyter Lab (dev image only)
    jupyter)
        echo "[OmicSage] Starting Jupyter Lab on port 8888..."
        exec /opt/conda/envs/omicsage/bin/jupyter lab \
            --ip=0.0.0.0 \
            --port=8888 \
            --no-browser \
            --allow-root \
            --NotebookApp.token='' \
            --NotebookApp.password=''
        ;;

    # Run pytest (dev image only)
    pytest)
        shift
        echo "[OmicSage] Running tests..."
        exec ${PYTHON} -m pytest "${@}"
        ;;

    # Interactive shell
    bash|shell)
        exec /bin/bash
        ;;

    # Explicit modality override: docker run omicsage:latest scrna --config ...
    scrna|cite|multiome|spatial)
        MODALITY="${CMD_ARG}"
        RUNNER="${RUNNERS[$MODALITY]}"
        shift
        echo "[OmicSage] Running ${MODALITY} pipeline → ${RUNNER}"
        exec ${PYTHON} "${RUNNER}" "${@}"
        ;;

    # Default: use OMICSAGE_MODALITY env var
    *)
        MODALITY="${OMICSAGE_MODALITY:-scrna}"
        RUNNER="${RUNNERS[$MODALITY]}"

        if [ -z "${RUNNER}" ]; then
            echo "[OmicSage] ERROR: Unknown modality '${MODALITY}'"
            echo "  Valid options: scrna | cite | multiome | spatial"
            echo "  Set via: -e OMICSAGE_MODALITY=<modality>"
            exit 1
        fi

        echo "[OmicSage] Modality: ${MODALITY} → ${RUNNER}"
        exec ${PYTHON} "${RUNNER}" "${@}"
        ;;
esac
