#!/bin/bash
# Run all 5 stability simulations sequentially: 4 R4 peptides + native NRF2 control
#
# GTX 1080 CUDA estimate: ~30-50 ns/day → ~1-1.5 days per run → ~5-7.5 days total
# M2 OpenCL estimate: ~7.6 ns/day → ~6.5 days per run → ~33 days total
#
# Usage:
#   bash run_all_5.sh           # default 50 ns each
#   bash run_all_5.sh 20        # 20 ns each (faster, still informative)
#   bash run_all_5.sh resume    # resume all from checkpoints

set -e
cd "$(dirname "$0")"

# Allow custom duration or resume mode
NS=${1:-50}
RESUME=""
if [ "$1" = "resume" ]; then
    NS=50
    RESUME="--resume"
fi

INPUT_DIR="input_structures"

echo "============================================================"
echo "SEQUENTIAL STABILITY SIMULATIONS — 5 x ${NS} ns"
echo "  4 R4 peptide candidates + 1 native NRF2 control"
echo "Started: $(date)"
echo "============================================================"
echo ""
echo "  To PAUSE: press Ctrl+C (saves checkpoint, can resume later)"
echo "  To STOP between runs: create a file called STOP"
echo "      touch STOP"
echo "  To RESUME all: bash run_all_5.sh resume"
echo ""

check_stop() {
    if [ -f "STOP" ]; then
        echo ""
        echo ">>> STOP file detected — halting between runs."
        echo ">>> Completed runs are saved. To continue:"
        echo ">>>   rm STOP"
        echo ">>>   bash run_all_5.sh resume"
        exit 0
    fi
}

# Run 1: Native NRF2-KEAP1 (2FLU crystal structure — CONTROL)
if [ -f "results/NRF2_native/NRF2_native_final.pdb" ]; then
    echo "[1/5] NRF2_native (CONTROL) — already complete, skipping"
else
    echo "[1/5] NRF2_native (CONTROL) — starting $(date)"
    python run_stability.py \
        --pdb "$INPUT_DIR/2FLU_KEAP1_NRF2_native.pdb" \
        --name NRF2_native \
        --ns $NS $RESUME
    echo "[1/5] NRF2_native (CONTROL) — finished $(date)"
fi
echo ""
check_stop

# Run 2: R4_DPETGEDF (HADDOCK -106.1, Z-score -2.0, largest cluster)
if [ -f "results/R4_DPETGEDF/R4_DPETGEDF_final.pdb" ]; then
    echo "[2/5] R4_DPETGEDF — already complete, skipping"
else
    echo "[2/5] R4_DPETGEDF — starting $(date)"
    python run_stability.py \
        --pdb "$INPUT_DIR/R4_DPETGEDF_KEAP1_best.pdb" \
        --name R4_DPETGEDF \
        --ns $NS $RESUME
    echo "[2/5] R4_DPETGEDF — finished $(date)"
fi
echo ""
check_stop

# Run 3: R4_LDEETGEWY (HADDOCK -106.3, highest raw score)
if [ -f "results/R4_LDEETGEWY/R4_LDEETGEWY_final.pdb" ]; then
    echo "[3/5] R4_LDEETGEWY — already complete, skipping"
else
    echo "[3/5] R4_LDEETGEWY — starting $(date)"
    python run_stability.py \
        --pdb "$INPUT_DIR/R4_LDEETGEWY_KEAP1_best.pdb" \
        --name R4_LDEETGEWY \
        --ns $NS $RESUME
    echo "[3/5] R4_LDEETGEWY — finished $(date)"
fi
echo ""
check_stop

# Run 4: R4_DPETGEWY (HADDOCK -103.4)
if [ -f "results/R4_DPETGEWY/R4_DPETGEWY_final.pdb" ]; then
    echo "[4/5] R4_DPETGEWY — already complete, skipping"
else
    echo "[4/5] R4_DPETGEWY — starting $(date)"
    python run_stability.py \
        --pdb "$INPUT_DIR/R4_DPETGEWY_KEAP1_best.pdb" \
        --name R4_DPETGEWY \
        --ns $NS $RESUME
    echo "[4/5] R4_DPETGEWY — finished $(date)"
fi
echo ""
check_stop

# Run 5: R4_DPETGEFF (HADDOCK -102.7)
if [ -f "results/R4_DPETGEFF/R4_DPETGEFF_final.pdb" ]; then
    echo "[5/5] R4_DPETGEFF — already complete, skipping"
else
    echo "[5/5] R4_DPETGEFF — starting $(date)"
    python run_stability.py \
        --pdb "$INPUT_DIR/R4_DPETGEFF_KEAP1_best.pdb" \
        --name R4_DPETGEFF \
        --ns $NS $RESUME
    echo "[5/5] R4_DPETGEFF — finished $(date)"
fi
echo ""

echo "============================================================"
echo "ALL 5 SIMULATIONS COMPLETE"
echo "Finished: $(date)"
echo "============================================================"
echo ""
echo "Results in:"
echo "  results/NRF2_native/    (CONTROL — native NRF2 peptide)"
echo "  results/R4_DPETGEDF/"
echo "  results/R4_LDEETGEWY/"
echo "  results/R4_DPETGEWY/"
echo "  results/R4_DPETGEFF/"
echo ""
echo "Compare RMSD of each R4 peptide vs native to assess binding stability."
