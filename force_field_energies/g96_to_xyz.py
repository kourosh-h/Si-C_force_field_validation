"""
Writes the structures of a GROMACS .g96 file (nm; POSITION or POSITIONRED blocks) as an .xyz file (Angstrom), with
the elements of an .xyz file of the same atoms in the same order.

Run:  python3 g96_to_xyz.py structures.g96 elements.xyz structures.xyz
"""
import sys
from pathlib import Path

g96, elements_xyz, output = sys.argv[1], sys.argv[2], sys.argv[3]
lines = Path(elements_xyz).read_text().split('\n')
elements = []
for line in lines[2:2 + int(lines[0])]:
    elements.append(line.split()[0])
out = []
step = 0
lines = Path(g96).read_text().split('\n')
for i in range(len(lines)):
    if lines[i] in ('POSITION', 'POSITIONRED'):
        out.append(str(len(elements)))
        out.append(f'step {step}')
        for k in range(len(elements)):
            line = lines[i + 1 + k]
            # POSITION lines start with residue and atom names, POSITIONRED lines with the coordinates
            if lines[i] == 'POSITION':
                line = line[24:]
            x, y, z = float(line[0:15]) * 10, float(line[15:30]) * 10, float(line[30:45]) * 10
            out.append(f'{elements[k]:<2s} {x:14.6f} {y:14.6f} {z:14.6f}')
        step += 1
Path(output).write_text('\n'.join(out) + '\n')
