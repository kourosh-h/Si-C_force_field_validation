Force-field energies of the scan structures

The force-field energy of octylsilanetriol, Si(OH)3-C8H17, at every structure of the ORCA bond, angle and
dihedral scans in ../ORCA_results/scans, computed with GROMACS as single points (gmx mdrun -rerun), and the minimum
of the force-field energy of the molecule, from the relaxed ORCA structure. No molecular dynamics is run: every
energy is the energy of one structure.
The topology has the atom types, charges and parameters of the force field of the alkylsilane monolayers: CHARMM36
(GROMACS port of July 2022) for the alkyl chain and the OH groups on Si, and the terms at the Si-C bond. The model
molecule has OH groups on all three oxygens, so with these charges it carries a net charge of +0.19 e.

Files
octylsilanetriol.top     topology of the molecule with the parameters it needs (input of every run)
rerun.mdp                GROMACS settings of the single points: all atom pairs counted, no constraints
minimize.mdp             GROMACS settings of the minimization
make_frames.py           writes the structures of each scan as <scan>/frames.g96 (positions in nm), and the relaxed
                         ORCA structure as minimized/start.g96
run_energies.sh          runs GROMACS for the eight scans and the minimization
g96_to_xyz.py            writes the steps of the minimization as an .xyz file (used by run_energies.sh)

In each scan folder (Si-C, C4-C5, Si-C-C, C-C-C, O-Si-C, Si-C-C-C, C-C-C-C, O-Si-C-C):
frames.g96               the structures of the scan, in the order of the ORCA scan (input)
rerun.tpr                the run input of GROMACS, made from frames.g96, ../octylsilanetriol.top and ../rerun.mdp
mdout.mdp                all GROMACS settings of the run, as written by gmx grompp
rerun.log                the log of gmx mdrun
rerun.edr                all energy terms of every structure
potential.xvg            structure number (from 0) and force-field energy (kJ/mol) of every structure

In minimized:
start.g96                the relaxed ORCA structure (input)
minimized.tpr            the run input of GROMACS, made from start.g96, ../octylsilanetriol.top and ../minimize.mdp
mdout.mdp                all GROMACS settings of the run, as written by gmx grompp
minimized.log            the log of gmx mdrun
minimized.edr            the energy terms of the minimization steps
minimized.trr            every step of the minimization
minimized.g96            the minimum of the force-field energy (positions in nm)
minimized_trj.xyz        every step of the minimization as an .xyz file (Angstrom), shown in the notebook

The lines with the name of the computer and of the working folder were removed from the .log and mdout.mdp files.

Scans: Si-C and C4-C5 bonds, Si-C-C, C-C-C and O-Si-C angles (50 structures each), Si-C-C-C, C-C-C-C and
O-Si-C-C dihedrals (37 structures each).

Run (GROMACS 2025, gmx in the PATH)
python3 make_frames.py
bash run_energies.sh
