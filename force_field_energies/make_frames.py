from pathlib import Path

import numpy as np


SCANS = Path('../ORCA_results/scans')
FILES = {'Si-C': 'si_c_scan', 'C4-C5': 'c4_c5_scan', 'Si-C-C': 'si_c_c_scan', 'O-Si-C': 'o_si_c_scan', 'C-C-C': 'c_c_c_scan',
         'O-Si-C-C': 'o_si_c_c_scan', 'Si-C-C-C': 'si_c_c_c_scan', 'C-C-C-C': 'c_c_c_c_scan'}
BOX = 10.0


def read_frames(path):
    lines = Path(path).read_text().split('\n')
    frames = []
    i = 0
    while i < len(lines):
        if lines[i].strip().isdigit():
            n = int(lines[i])
            positions = []
            for line in lines[i + 2:i + 2 + n]:
                parts = line.split()
                positions.append([float(parts[1]), float(parts[2]), float(parts[3])])
            frames.append(np.array(positions))
            i = i + n + 2
        else:
            i = i + 1
    return frames


def atom_names(path):
    names = []
    atoms = Path(path).read_text().split('[ atoms ]')[1].split('[ bonds ]')[0]
    for line in atoms.split('\n'):
        parts = line.split()
        if parts and parts[0].isdigit():
            names.append(parts[4])
    return names


def write_g96(frames, names, path):
    lines = []
    for step in range(len(frames)):
        xyz = frames[step] / 10
        xyz = xyz - xyz.mean(axis=0) + BOX / 2
        lines += ['TITLE', f'scan point {step + 1}', 'END', 'TIMESTEP', f'{step:15d}{float(step):15.6f}', 'END', 'POSITION']
        for i in range(len(xyz)):
            lines.append(f'{1:5d} OSL   {names[i]:<5s}{i + 1:7d}{xyz[i, 0]:15.9f}{xyz[i, 1]:15.9f}{xyz[i, 2]:15.9f}')
        lines += ['END', 'BOX', f'{BOX:15.9f}{BOX:15.9f}{BOX:15.9f}', 'END']
    Path(path).write_text('\n'.join(lines) + '\n')


names = atom_names('octylsilanetriol.top')
frames = read_frames('../ORCA_results/relaxation/octylsilanetriol.xyz')    # the relaxed molecule, start of the minimization
Path('minimized').mkdir(exist_ok=True)
write_g96(frames, names, Path('minimized') / 'start.g96')
for scan in FILES:
    frames = read_frames(SCANS / scan / f'{FILES[scan]}.allxyz')
    Path(scan).mkdir(exist_ok=True)
    write_g96(frames, names, Path(scan) / 'frames.g96')
    print(f'{scan}: {len(frames)} structures')
