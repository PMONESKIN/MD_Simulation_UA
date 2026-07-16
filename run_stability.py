#!/usr/bin/env python3
"""
OpenMM Binding Stability Simulation
====================================
Tests whether a HADDOCK-docked peptide stays bound to KEAP1 over time.

Usage:
    python run_stability.py --pdb input_structures/R4_DPETGEDF_KEAP1_best.pdb --name R4_DPETGEDF --ns 50
    python run_stability.py --pdb input_structures/R4_DPETGEFF_KEAP1_best.pdb --name R4_DPETGEFF --ns 50
    python run_stability.py --resume --name R4_DPETGEDF  # resume from checkpoint

OneSkin Technologies — Summer 2026
"""

import argparse
import os
import signal
import sys
import time
from pathlib import Path

# Global flag for graceful pause
_pause_requested = False

def _handle_pause(signum, frame):
    """Handle Ctrl+C: set pause flag instead of crashing."""
    global _pause_requested
    _pause_requested = True
    print("\n\n>>> PAUSE REQUESTED — will save checkpoint after current batch and stop.")
    print(">>> To resume later: python run_stability.py --name <name> --resume\n")

def run_simulation(pdb_path, name, sim_ns, output_dir, resume=False):
    """Set up and run an OpenMM stability simulation."""

    from openmm.app import (
        PDBFile, ForceField, Modeller, PME, HBonds,
        Simulation, DCDReporter, StateDataReporter,
        CheckpointReporter
    )
    from openmm import (
        LangevinMiddleIntegrator, MonteCarloBarostat, Platform
    )
    from openmm.unit import (
        nanometers, kelvin, atmospheres, picoseconds,
        picosecond, femtoseconds, kilojoules_per_mole, molar
    )
    from pdbfixer import PDBFixer

    print(f"{'='*60}")
    print(f"OPENMM BINDING STABILITY SIMULATION")
    print(f"Complex: {name}")
    print(f"PDB: {pdb_path}")
    print(f"Duration: {sim_ns} ns")
    print(f"Output: {output_dir}/")
    print(f"Resume: {resume}")
    print(f"{'='*60}\n")

    os.makedirs(output_dir, exist_ok=True)

    # Output file paths
    solvated_pdb = os.path.join(output_dir, f"{name}_solvated.pdb")
    min_pdb = os.path.join(output_dir, f"{name}_minimized.pdb")
    traj_file = os.path.join(output_dir, f"{name}_trajectory.dcd")
    log_file = os.path.join(output_dir, f"{name}_energy.csv")
    checkpoint_file = os.path.join(output_dir, f"{name}_checkpoint.chk")
    final_pdb = os.path.join(output_dir, f"{name}_final.pdb")

    # CHARMM36 for protein + peptide, TIP3P for water
    forcefield = ForceField('charmm36.xml', 'charmm36/water.xml')

    # ============================================================
    # Step 1-2: Load or build solvated system
    # ============================================================
    # When resuming, load the saved solvated PDB to guarantee same atom count
    if resume and os.path.exists(solvated_pdb):
        print("[1-2/6] Loading saved solvated system for resume...")
        pdb_loaded = PDBFile(solvated_pdb)
        topology = pdb_loaded.topology
        positions = pdb_loaded.positions
        print(f"  Loaded: {solvated_pdb}")
        print(f"  Total atoms: {topology.getNumAtoms()}")
    else:
        print("[1/6] Fixing PDB structure...")

        fixer = PDBFixer(filename=str(pdb_path))
        fixer.findMissingResidues()
        fixer.findMissingAtoms()
        fixer.addMissingAtoms()
        fixer.addMissingHydrogens(7.4)  # pH 7.4

        # Remove water from HADDOCK (we'll add our own)
        fixer.removeHeterogens(keepWater=False)

        print(f"  Missing residues found: {len(fixer.missingResidues)}")
        print(f"  Structure fixed")

        print("[2/6] Setting up force field and solvation...")

        modeller = Modeller(fixer.topology, fixer.positions)

        # Add water box with 1.0 nm padding and neutralizing ions
        modeller.addSolvent(
            forcefield,
            model='tip3p',
            padding=1.0 * nanometers,
            ionicStrength=0.15 * molar,  # 150 mM NaCl (physiological)
        )

        topology = modeller.topology
        positions = modeller.positions

        n_atoms = topology.getNumAtoms()
        print(f"  Force field: CHARMM36")
        print(f"  Water model: TIP3P")
        print(f"  Ionic strength: 150 mM NaCl")
        print(f"  Total atoms: {n_atoms}")

        # Save solvated system
        with open(solvated_pdb, 'w') as f:
            PDBFile.writeFile(topology, positions, f)
        print(f"  Saved: {solvated_pdb}")

    # ============================================================
    # Step 3: Create simulation system
    # ============================================================
    print("[3/6] Creating simulation system...")

    system = forcefield.createSystem(
        topology,
        nonbondedMethod=PME,
        nonbondedCutoff=1.0 * nanometers,
        constraints=HBonds,
    )

    # Add pressure coupling (NPT ensemble)
    system.addForce(MonteCarloBarostat(1 * atmospheres, 310.15 * kelvin))

    # Langevin integrator: 310.15 K (37°C body temp), 2 fs timestep
    integrator = LangevinMiddleIntegrator(
        310.15 * kelvin,    # temperature
        1.0 / picosecond,   # friction coefficient
        2.0 * femtoseconds  # timestep
    )

    # Select best available platform (CUDA > OpenCL > Metal > CPU)
    platform = None
    for plat_name in ['CUDA', 'OpenCL', 'Metal', 'CPU']:
        try:
            platform = Platform.getPlatformByName(plat_name)
            print(f"  Platform: {plat_name}")
            break
        except Exception:
            continue

    simulation = Simulation(topology, system, integrator, platform)
    simulation.context.setPositions(positions)

    print(f"  Temperature: 310.15 K (37°C)")
    print(f"  Pressure: 1 atm (NPT)")
    print(f"  Timestep: 2 fs")

    # ============================================================
    # Resume from checkpoint if requested
    # ============================================================
    total_steps = int(sim_ns * 1000 / 0.002)  # ns → ps → steps (2 fs/step)
    report_interval = 25000     # Report energy every 50 ps
    traj_interval = 25000       # Save frame every 50 ps
    checkpoint_interval = 250000 # Checkpoint every 500 ps

    if resume and os.path.exists(checkpoint_file):
        print(f"\n[RESUME] Loading checkpoint: {checkpoint_file}")
        simulation.loadCheckpoint(checkpoint_file)
        current_step = simulation.context.getState().getStepCount()
        remaining_steps = total_steps - current_step
        current_ns = current_step * 0.002 / 1000
        print(f"  Resumed at step {current_step} ({current_ns:.2f} ns)")
        print(f"  Remaining: {remaining_steps} steps ({remaining_steps * 0.002 / 1e6:.1f} ns)")

        if remaining_steps <= 0:
            print("  Simulation already complete!")
            return

        # Set up reporters — new files for resumed segment to avoid corruption
        resume_traj = os.path.join(output_dir, f"{name}_trajectory_resumed_{current_step}.dcd")
        resume_log = os.path.join(output_dir, f"{name}_energy_resumed_{current_step}.csv")
        simulation.reporters.append(DCDReporter(resume_traj, traj_interval))
        simulation.reporters.append(StateDataReporter(
            resume_log, report_interval,
            step=True, time=True, potentialEnergy=True, kineticEnergy=True,
            totalEnergy=True, temperature=True, volume=True, speed=True,
            separator=','
        ))
        simulation.reporters.append(StateDataReporter(
            sys.stdout, report_interval * 4,
            step=True, time=True, temperature=True, speed=True,
            remainingTime=True, totalSteps=total_steps
        ))
        simulation.reporters.append(CheckpointReporter(checkpoint_file, checkpoint_interval))

        # Run remaining steps in chunks (allows graceful pause with Ctrl+C)
        signal.signal(signal.SIGINT, _handle_pause)
        chunk_size = checkpoint_interval
        steps_done = 0
        start_time = time.time()

        while steps_done < remaining_steps:
            batch = min(chunk_size, remaining_steps - steps_done)
            simulation.step(batch)
            steps_done += batch

            # Check for pause request (Ctrl+C) or STOP file
            stop_file = os.path.join(Path(__file__).parent, "STOP")
            if _pause_requested or os.path.exists(stop_file):
                elapsed = time.time() - start_time
                total_done = current_step + steps_done
                ns_done = total_done * 0.002 / 1000
                print(f"\n>>> PAUSED at step {total_done} ({ns_done:.2f} ns)")
                print(f">>> Checkpoint saved: {checkpoint_file}")
                print(f">>> To resume: python run_stability.py --name {name} --resume")
                if os.path.exists(stop_file):
                    print(f">>> (Remove STOP file before resuming)")
                paused_pdb = os.path.join(output_dir, f"{name}_paused.pdb")
                positions = simulation.context.getState(getPositions=True).getPositions()
                with open(paused_pdb, 'w') as f:
                    PDBFile.writeFile(simulation.topology, positions, f)
                sys.exit(0)

        elapsed = time.time() - start_time
        signal.signal(signal.SIGINT, signal.SIG_DFL)

    else:
        # ============================================================
        # Step 4: Energy minimization
        # ============================================================
        print("[4/6] Energy minimization...")

        e_before = simulation.context.getState(getEnergy=True).getPotentialEnergy()
        print(f"  Before: {e_before}")

        simulation.minimizeEnergy(maxIterations=1000)

        e_after = simulation.context.getState(getEnergy=True).getPotentialEnergy()
        print(f"  After:  {e_after}")

        # Save minimized structure
        positions = simulation.context.getState(getPositions=True).getPositions()
        with open(min_pdb, 'w') as f:
            PDBFile.writeFile(simulation.topology, positions, f)

        # ============================================================
        # Step 5: Equilibration (100 ps NVT + 100 ps NPT)
        # ============================================================
        print("[5/6] Equilibration...")

        # Set velocities
        simulation.context.setVelocitiesToTemperature(310.15 * kelvin)

        # NVT equilibration (50,000 steps = 100 ps)
        print(f"  NVT: 100 ps...")
        simulation.step(50000)

        # NPT equilibration (50,000 steps = 100 ps)
        print(f"  NPT: 100 ps...")
        simulation.step(50000)

        print(f"  Equilibration complete")

        # ============================================================
        # Step 6: Production run
        # ============================================================
        print(f"[6/6] Production run: {sim_ns} ns ({total_steps} steps)...")
        print(f"  Trajectory frame every 50 ps")
        print(f"  Energy report every 50 ps")
        print(f"  Checkpoint every 500 ps")
        print()

        # Set up reporters
        simulation.reporters.append(DCDReporter(traj_file, traj_interval))
        simulation.reporters.append(StateDataReporter(
            log_file, report_interval,
            step=True, time=True, potentialEnergy=True, kineticEnergy=True,
            totalEnergy=True, temperature=True, volume=True, speed=True,
            separator=','
        ))
        simulation.reporters.append(StateDataReporter(
            sys.stdout, report_interval * 4,
            step=True, time=True, temperature=True, speed=True,
            remainingTime=True, totalSteps=total_steps
        ))
        simulation.reporters.append(CheckpointReporter(checkpoint_file, checkpoint_interval))

        # Run production in chunks (allows graceful pause with Ctrl+C)
        signal.signal(signal.SIGINT, _handle_pause)
        chunk_size = checkpoint_interval  # run in 500 ps chunks
        steps_done = 0
        start_time = time.time()

        while steps_done < total_steps:
            batch = min(chunk_size, total_steps - steps_done)
            simulation.step(batch)
            steps_done += batch

            # Check for pause request (Ctrl+C) or STOP file
            stop_file = os.path.join(Path(__file__).parent, "STOP")
            if _pause_requested or os.path.exists(stop_file):
                elapsed = time.time() - start_time
                ns_done = steps_done * 0.002 / 1000
                print(f"\n>>> PAUSED at step {steps_done} ({ns_done:.2f} ns)")
                print(f">>> Checkpoint saved: {checkpoint_file}")
                print(f">>> To resume: python run_stability.py --name {name} --resume")
                if os.path.exists(stop_file):
                    print(f">>> (Remove STOP file before resuming)")
                # Save current structure
                paused_pdb = os.path.join(output_dir, f"{name}_paused.pdb")
                positions = simulation.context.getState(getPositions=True).getPositions()
                with open(paused_pdb, 'w') as f:
                    PDBFile.writeFile(simulation.topology, positions, f)
                sys.exit(0)

        elapsed = time.time() - start_time
        signal.signal(signal.SIGINT, signal.SIG_DFL)

    # Save final structure
    positions = simulation.context.getState(getPositions=True).getPositions()
    with open(final_pdb, 'w') as f:
        PDBFile.writeFile(simulation.topology, positions, f)

    print(f"\n{'='*60}")
    print(f"SIMULATION COMPLETE: {name}")
    print(f"{'='*60}")
    print(f"  Duration: {sim_ns} ns")
    print(f"  Wall time: {elapsed/3600:.1f} hours")
    print(f"  Speed: {sim_ns / (elapsed/86400):.1f} ns/day")
    print(f"  Trajectory: {traj_file}")
    print(f"  Energy log: {log_file}")
    print(f"  Final structure: {final_pdb}")
    print(f"  Solvated PDB: {solvated_pdb}")
    print(f"\n  View in VMD:")
    print(f"    vmd {solvated_pdb} {traj_file}")


def main():
    parser = argparse.ArgumentParser(description="OpenMM binding stability simulation")
    parser.add_argument("--pdb", help="HADDOCK docked complex PDB")
    parser.add_argument("--name", required=True, help="Run name")
    parser.add_argument("--ns", type=float, default=50, help="Simulation length in ns (default: 50)")
    parser.add_argument("--output", default=None, help="Output directory (default: results/<name>)")
    parser.add_argument("--resume", action="store_true", help="Resume from checkpoint")

    args = parser.parse_args()

    if not args.resume and not args.pdb:
        parser.error("--pdb is required unless --resume is set")

    if args.resume and not args.pdb:
        # When resuming, we still need the PDB for topology
        default_pdb = Path(__file__).parent / "input_structures" / f"{args.name}_KEAP1_best.pdb"
        if default_pdb.exists():
            args.pdb = str(default_pdb)
        else:
            parser.error(f"Cannot find PDB for resume: {default_pdb}")

    pdb_path = Path(args.pdb)
    if not pdb_path.exists():
        # Check in input_structures folder
        alt = Path(__file__).parent / "input_structures" / pdb_path.name
        if alt.exists():
            pdb_path = alt
        else:
            print(f"ERROR: PDB not found: {args.pdb}")
            sys.exit(1)

    output_dir = args.output or str(Path(__file__).parent / "results" / args.name)

    run_simulation(pdb_path, args.name, args.ns, output_dir, resume=args.resume)


if __name__ == "__main__":
    main()
