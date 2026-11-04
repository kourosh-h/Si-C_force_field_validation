Relaxed scans of octylsilanetriol

Relaxed scans of octylsilanetriol, Si(OH)3-C8H17, at the HF/6-31G* level with ORCA. Every scan starts from the relaxed
molecule (octylsilanetriol.xyz, from ../relaxation). At each point the scanned coordinate is fixed and the rest of the
molecule is optimized. The atoms are counted from 0, as in the ORCA files: Si 0, O 5, C1 7, C2 10, C3 13, C4 16,
C5 19, C6 22.

Folder     Coordinate (atoms)                    Points
Si-C       Si-C1 bond (0 7)                      50, from 1.80 to 1.94 Angstrom
C4-C5      C4-C5 bond (16 19)                    50, from 1.45 to 1.59 Angstrom
Si-C-C     Si-C1-C2 angle (0 7 10)               50, from 100 to 130 degrees
C-C-C      C3-C4-C5 angle (13 16 19)             50, from 98 to 128 degrees
O-Si-C     O-Si-C1 angle (5 0 7)                 50, from 97 to 127 degrees
Si-C-C-C   Si-C1-C2-C3 dihedral (0 7 10 13)      37, from -180 to 180 degrees
C-C-C-C    C3-C4-C5-C6 dihedral (13 16 19 22)    37, from -180 to 180 degrees
O-Si-C-C   O-Si-C1-C2 dihedral (5 0 7 10)        37, from -185 to 175 degrees

Files of each folder (<name> is the name of the scan, for example si_c_scan)
<name>.inp                  ORCA input
octylsilanetriol.xyz        start structure, the relaxed molecule
<name>.out                  ORCA output
<name>.relaxscanact.dat     the coordinate and the energy (Hartree) of every point
<name>.relaxscanscf.dat     the coordinate and the SCF energy (Hartree) of every point
<name>.allxyz               the structures of all points
<name>.001.xyz ...          the structure of each point
<name>.xyz                  the structure of the last point
<name>_trj.xyz              the structures of all optimization steps
<name>_property.txt         the properties of every point, as written by ORCA
<name>.opt                  the data of the optimizer, as written by ORCA
<name>.densities            the electron densities, as written by ORCA
<name>.err                  the error messages of ORCA (empty)

The .gbw files, the binary wavefunctions of ORCA (about 1 MB for every point), are not included.

Run, in each folder
orca <name>.inp > <name>.out
