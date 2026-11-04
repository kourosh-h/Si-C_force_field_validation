GROMACS force field of alkylsilane monolayers on silica in water

alkylsilane_silica.ff is the force field of alkylsilane monolayers on silica in water in one GROMACS force-field
folder. It has the 14 atom types of these systems and every parameter between them, with the values of the files in
../ForceField:
CHARMM36 (GROMACS port of July 2022)    alkyl chain: CG321, CG331, HGA2, HGA3
                                        OH groups on the silicon of the head group: SI, OSIH, HSIO
                                        water: OT, HT
INTERFACE force field of silica         ISC4, IOC23, IOC24, IOC25 (the oxygen of the bridge to the head group), IHOY
the terms at the Si-C bond and at the bridge to the silica
CHARMM TIP3P water                      the water molecule

Files
alkylsilane_silica.ff/forcefield.itp     [ defaults ] (Lennard-Jones with the Lorentz-Berthelot rule, 1-4 pairs
                                         generated, fudgeLJ and fudgeQQ 1.0) and the two files below
alkylsilane_silica.ff/ffnonbonded.itp    atom types and 1-4 pairs
alkylsilane_silica.ff/ffbonded.itp       bonds, angles with their Urey-Bradley terms, and dihedrals
alkylsilane_silica.ff/tip3p.itp          the TIP3P water molecule (TIP3), rigid (SETTLE)
make_force_field.py                      writes alkylsilane_silica.ff from the files in ../ForceField

Use
Copy alkylsilane_silica.ff into the folder of the topology and include it in the .top file:

#include "alkylsilane_silica.ff/forcefield.itp"
#include "monolayer.itp"
#include "alkylsilane_silica.ff/tip3p.itp"

monolayer.itp is the molecule topology of the silica slab with its chains: the atoms with their types and charges,
the bonds, 1-4 pairs, angles and dihedrals. These terms take their parameters from alkylsilane_silica.ff by atom
type. The charge column of [ atomtypes ] has the charge of each atom type in the monolayers and in TIP3P water;
GROMACS uses the charges in [ atoms ] of the molecule topology.

Check
For a monolayer in water (80,903 atoms), the files in ../ForceField and alkylsilane_silica.ff give the same energy
terms (GROMACS 2025.3, double precision, single point).

Run
python3 make_force_field.py
