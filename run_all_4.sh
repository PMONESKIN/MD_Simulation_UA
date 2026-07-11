#!/bin/bash
# Run all 4 R4 peptide stability simulations sequentially (50 ns each)
#
# GTX 1080 CUDA estimate: ~30-50 ns/day → ~1-1.5 days per run → ~4-6 days total
# M2 OpenCL estimate: ~7.6 ns/day → ~6.5 days per run → ~26 days total
#
# Usage:
#   bash run_all_4.sh           # default 50 ns each
#   bash run_all_4.sh 20        # 20 ns each (faster, still informative)
#   bash run_all_4.sh resume    # resume all from checkpoints

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
echo "SEQUENTIAL STABILITY SIMULATIONS — 4 x ${NS} ns"
echo "Started: $(date)"
echo "============================================================"
echo ""

# Run 1: R4_DPETGEDF (HADDOCK -106.1, Z-score -2.0, largest cluster)
if [ -f "results/R4_DPETGEDF/R4_DPETGEDF_final.pdb" ]; then
    echo "[1/4] R4_DPETGEDF — already complete, skipping"
else
    echo "[1/4] R4_DPETGEDF — starting $(date)"
    python run_stability.py \
        --pdb "$INPUT_DIR/R4_DPETGEDF_KEAP1_best.pdb" \
        --name R4_DPETGEDF \
        --ns $NS $RESUME
    echo "[1/4] R4_DPETGEDF — finished $(date)"
fi
echo ""

# Run 2: R4_LDEETGEWY (HADDOCK -106.3, highest raw score)
if [ -f "results/R4_LDEETGEWY/R4_LDEETGEWY_final.pdb" ]; then
    echo "[2/4] R4_LDEETGEWY — already complete, skipping"
else
    echo "[2/4] R4_LDEETGEWY — starting $(date)"
    python run_stability.py \
        --pdb "$INPUT_DIR/R4_LDEETGEWY_KEAP1_best.pdb" \
        --name R4_LDEETGEWY \
        --ns $NS $RESUME
    echo "[2/4] R4_LDEETGEWY — finished $(date)"
fi
echo ""

# Run 3: R4_DPETGEWY (HADDOCK -103.4)
if [ -f "results/R4_DPETGEWY/R4_DPETGEWY_final.pdb" ]; then
    echo "[3/4] R4_DPETGEWY — already complete, skipping"
else
    echo "[3/4] R4_DPETGEWY — starting $(date)"
    python run_stability.py \
        --pdb "$INPUT_DIR/R4_DPETGEWY_KEAP1_best.pdb" \
        --name R4_DPETGEWY \
        --ns $NS $RESUME
    echo "[3/4] R4_DPETGEWY — finished $(date)"
fi
echo ""

# Run 4: R4_DPETGEFF (HADDOCK -102.7)
if [ -f "results/R4_DPETGEFF/R4_DPETGEFF_final.pdb" ]; then
    echo "[4/4] R4_DPETGEFF — already complete, skipping"
else
    echo "[4/4] R4_DPETGEFF — starting $(date)"
    python run_stability.py \
        --pdb "$INPUT_DIR/R4_DPETGEFF_KEAP1_best.pdb" \
        --name R4_DPETGEFF \
        --ns $NS $RESUME
    echo "[4/4] R4_DPETGEFF — finished $(date)"
fi
echo ""

echo "============================================================"
echo "ALL 4 SIMULATIONS COMPLETE"
echo "Finished: $(date)"
echo "============================================================"
echo ""
echo "Results in:"
echo "  results/R4_DPETGEDF/"
echo "  results/R4_LDEETGEWY/"
echo "  results/R4_DPETGEWY/"
echo "  results/R4_DPETGEFF/"
echo ""
echo "View in VMD:"
echo "  vmd results/R4_DPETGEDF/R4_DPETGEDF_solvated.pdb results/R4_DPETGEDF/R4_DPETGEDF_trajectory.dcd"
