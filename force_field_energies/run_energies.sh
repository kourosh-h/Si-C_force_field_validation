#!/bin/bash
# Force-field energy of every structure of the eight scans (GROMACS, gmx in the PATH)
for scan in Si-C C4-C5 Si-C-C O-Si-C C-C-C O-Si-C-C Si-C-C-C C-C-C-C
do
    cd $scan
    gmx grompp -f ../rerun.mdp -c frames.g96 -p ../octylsilanetriol.top -o rerun.tpr
    gmx mdrun -s rerun.tpr -rerun frames.g96 -deffnm rerun -nt 1
    echo Potential | gmx energy -f rerun.edr -o potential.xvg -xvg none
    cd ..
done

# the minimum of the force-field energy of the molecule, from the relaxed ORCA structure
cd minimized
gmx grompp -f ../minimize.mdp -c start.g96 -p ../octylsilanetriol.top -o minimized.tpr
gmx mdrun -s minimized.tpr -deffnm minimized -c minimized.g96 -nt 1
echo 0 | gmx trjconv -f minimized.trr -s minimized.tpr -o minimized_trj.g96
python3 ../g96_to_xyz.py minimized_trj.g96 ../../ORCA_results/relaxation/octylsilanetriol.xyz minimized_trj.xyz
rm minimized_trj.g96
cd ..
