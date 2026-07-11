# MD_Simulation_UA

OpenMM molecular dynamics stability simulations for R4-tagged KEAP1 peptide candidates + native NRF2 control.

Tests whether HADDOCK-docked peptides remain stably bound to KEAP1 under explicit-solvent MD at 310 K (37°C). The native NRF2-KEAP1 complex (PDB 2FLU) runs as a control for direct comparison.

## Peptides (5 simulations)

| Name | Sequence | HADDOCK Score | Z-Score | Notes |
|---|---|---|---|---|
| NRF2_native | AFFAQLQLDEETGEFL | — | — | CONTROL (2FLU crystal) |
| R4_DPETGEDF | RRRR-DPETGEDF | -106.1 ± 1.0 | -2.0 | Best Z-score, largest cluster |
| R4_LDEETGEWY | RRRR-LDEETGEWY | -106.3 ± 4.0 | -1.5 | Highest raw HADDOCK score |
| R4_DPETGEWY | RRRR-DPETGEWY | -103.4 ± 1.4 | -1.8 | |
| R4_DPETGEFF | RRRR-DPETGEFF | -102.7 ± 3.8 | -1.7 | |

## Setup

```bash
conda env create -f environment.yml
conda activate openmm-md
```

## Run

```bash
# All 5 sequentially (50 ns each — native control first)
bash run_all_5.sh

# Custom duration (20 ns each)
bash run_all_5.sh 20

# Single run
python run_stability.py --pdb input_structures/R4_DPETGEDF_KEAP1_best.pdb --name R4_DPETGEDF --ns 50

# Resume from checkpoint after interruption
python run_stability.py --name R4_DPETGEDF --resume
bash run_all_5.sh resume
```

## Performance

| GPU | Speed | Time per 50 ns | All 5 |
|---|---|---|---|
| GTX 1080 (CUDA) | ~30-50 ns/day | ~1-1.5 days | ~5-7.5 days |
| Apple M2 (OpenCL) | ~7.6 ns/day | ~6.5 days | ~33 days |

## Output

Each run produces in `results/<name>/`:
- `*_solvated.pdb` — solvated system topology
- `*_trajectory.dcd` — MD trajectory (50 ps frames)
- `*_energy.csv` — energy, temperature, volume over time
- `*_checkpoint.chk` — restart checkpoint (every 500 ps)
- `*_final.pdb` — final structure after simulation

View in VMD:
```
vmd results/R4_DPETGEDF/R4_DPETGEDF_solvated.pdb results/R4_DPETGEDF/R4_DPETGEDF_trajectory.dcd
```

## Analysis

Compare peptide RMSD, H-bonds, and contact distance vs the NRF2_native control to determine which R4 candidates maintain stable KEAP1 binding.
