"""
The force field of alkylsilane monolayers on silica in water as one GROMACS force-field folder, alkylsilane_silica.ff,
written from the force-field files in ../ForceField: CHARMM36 (GROMACS port of July 2022), the INTERFACE force field of
silica, the terms at the Si-C bond and at the bridge to the silica, and CHARMM TIP3P water. It keeps the 14 atom types
of these systems and every parameter between them. Where two files define the same term, the file read later in the
order of the simulations is kept (CHARMM36, INTERFACE, the terms at the Si-C bond), as GROMACS does.

The charge column of [ atomtypes ] holds the charge of each atom type in the monolayers and in TIP3P water; the
molecule topologies give the charges in [ atoms ].

Run:  python3 make_force_field.py
"""
from pathlib import Path

SOURCE = Path('../ForceField')
OUT = Path('alkylsilane_silica.ff')
# the atom types: alkyl chain and head group (CHARMM36), silica (INTERFACE), water (CHARMM TIP3P)
TYPES = ['CG321', 'CG331', 'HGA2', 'HGA3', 'SI', 'OSIH', 'HSIO', 'ISC4', 'IOC23', 'IOC24', 'IOC25', 'IHOY', 'OT', 'HT']
# their charges (e) in the monolayers and in TIP3P water
CHARGES = {'CG321': -0.18, 'CG331': -0.27, 'HGA2': 0.09, 'HGA3': 0.09, 'SI': 0.76, 'OSIH': -0.51, 'HSIO': 0.32,
           'ISC4': 1.1, 'IOC23': -0.55, 'IOC24': -0.675, 'IOC25': -0.655, 'IHOY': 0.4, 'OT': -0.834, 'HT': 0.417}
FILES = [('CHARMM36', [SOURCE / 'charmm36-jul2022.ff' / 'ffnonbonded.itp', SOURCE / 'charmm36-jul2022.ff' / 'ffbonded.itp']),
         ('INTERFACE silica', [SOURCE / 'interface_silica_emami.ff' / 'ff_interface_emami.itp']),
         ('Si-C bond and bridge', [SOURCE / 'charmm_sam_connection_point.ff' / 'ff_charmm_connection_point.itp'])]
# number of atom types in an entry of each section
SIZES = {'atomtypes': 1, 'pairtypes': 2, 'bondtypes': 2, 'angletypes': 3, 'dihedraltypes': 4}
HEADERS = {'atomtypes': '; name  at.num        mass    charge  ptype          sigma        epsilon',
           'pairtypes': '; i      j      func        sigma14      epsilon14',
           'bondtypes': '; i      j      func           b0             kb',
           'angletypes': '; i      j      k      func         theta0        ktheta             s0            kub',
           'dihedraltypes': '; i      j      k      l      func           phi0           kphi  mult'}
DESCRIPTIONS = {'CHARMM36': 'CHARMM36: the alkyl chain, the silanol groups of the head group and TIP3P water',
                'INTERFACE silica': 'INTERFACE force field of silica',
                'Si-C bond and bridge': 'the terms at the Si-C bond and at the bridge to the silica'}


def read_entries(path):
    """The entries of a force-field file: (section, words), with the #ifdef blocks of undefined names left out."""
    entries = []
    section = None
    active = []
    for line in path.read_text().split('\n'):
        line = line.split(';')[0].strip()
        if line.startswith('#ifdef'):
            active.append(False)
        elif line.startswith('#ifndef'):
            active.append(True)
        elif line.startswith('#else'):
            active[-1] = not active[-1]
        elif line.startswith('#endif'):
            active.pop()
        elif not line or line.startswith('#') or not all(active):
            continue
        elif line.startswith('['):
            section = line.strip('[] ')
        elif section in SIZES:
            entries.append((section, line.split()))
    return entries


def key_of(section, words):
    """The atom types and function of an entry, the same for both directions of the types."""
    size = SIZES[section]
    types = tuple(words[:size])
    if section == 'atomtypes':
        return types
    return min(types, types[::-1]) + (words[size],)


# every entry between the 14 atom types; a term that a later file defines again is replaced
terms = {}
replaced = []
for source, paths in FILES:
    found = {}
    for path in paths:
        for section, words in read_entries(path):
            if all(name in TYPES for name in words[:SIZES[section]]):
                found.setdefault((section, key_of(section, words)), []).append(words)
    for key, rows in found.items():
        if key in terms and terms[key][1] != rows:
            replaced.append((key, terms[key][0], source))
        terms[key] = (source, rows)


def row_text(section, words):
    """An entry as a line, with the values as written in the source file (the charge of an atom type from CHARGES)."""
    size = SIZES[section]
    names = ' '.join(f'{name:<6s}' for name in words[:size])
    if section == 'atomtypes':
        number, mass, charge, ptype, sigma, epsilon = words[1:7]
        return f'{words[0]:<6s} {number:>5s} {mass:>11s} {CHARGES[words[0]]:9.4f}  {ptype:>5s} {sigma:>14s} {epsilon:>14s}'
    values = ' '.join(f'{value:>14s}' for value in words[size + 1:])
    return f'{names} {words[size]:>5s} {values}'


def write_sections(path, sections, title):
    lines = [f'; {title}', '; written by make_force_field.py from ../ForceField', '']
    for section in sections:
        lines += [f'[ {section} ]', HEADERS[section]]
        for source, paths in FILES:
            rows = []
            for (kind, key), (origin, entries) in terms.items():
                if kind == section and origin == source:
                    rows += entries
            if section == 'atomtypes':
                rows.sort(key=lambda words: TYPES.index(words[0]))
            if rows:
                lines.append(f'; {DESCRIPTIONS[source]}')
                lines += [row_text(section, words) for words in rows]
        lines.append('')
    path.write_text('\n'.join(lines))


defaults = None
for line in (SOURCE / 'charmm36-jul2022.ff' / 'forcefield.itp').read_text().split('\n'):
    if defaults == 'next' and line.strip() and not line.startswith(';'):
        defaults = line.split()
    if line.strip() == '[ defaults ]':
        defaults = 'next'
OUT.mkdir(exist_ok=True)
(OUT / 'forcefield.itp').write_text('\n'.join([
    '; The force field of alkylsilane monolayers on silica in water: CHARMM36, the INTERFACE force field of silica,',
    '; the terms at the Si-C bond and at the bridge to the silica, and CHARMM TIP3P water (tip3p.itp).',
    '; written by make_force_field.py from ../ForceField',
    '', '[ defaults ]', '; nbfunc  comb-rule  gen-pairs  fudgeLJ  fudgeQQ',
    f'  {defaults[0]:>6s} {defaults[1]:>10s} {defaults[2]:>10s} {defaults[3]:>8s} {defaults[4]:>8s}', '',
    '#include "ffnonbonded.itp"', '#include "ffbonded.itp"', '']))
write_sections(OUT / 'ffnonbonded.itp', ['atomtypes', 'pairtypes'], 'Nonbonded parameters: atom types and 1-4 pairs')
write_sections(OUT / 'ffbonded.itp', ['bondtypes', 'angletypes', 'dihedraltypes'],
               'Bonded parameters: bonds, angles with their Urey-Bradley terms, and dihedrals')

# TIP3P water: the molecule of ../ForceField/charmm_water_TIP3.ff/TIP3.itp
water = {}
section = None
for line in (SOURCE / 'charmm_water_TIP3.ff' / 'TIP3.itp').read_text().split('\n'):
    line = line.split(';')[0].strip()
    if line.startswith('['):
        section = line.strip('[] ')
    elif line:
        water.setdefault(section, []).append(line.split())
atoms = water['atoms']
settles = water['settles'][0]
lines = ['; CHARMM TIP3P water, rigid (SETTLE); the atom types OT and HT are in ffnonbonded.itp',
         '; written by make_force_field.py from ../ForceField', '',
         '[ moleculetype ]', '; name  nrexcl', f'{water["moleculetype"][0][0]}  {water["moleculetype"][0][1]}', '',
         '[ atoms ]', ';  nr  type  resnr  residue  atom  cgnr     charge      mass']
for nr, kind, resnr, residue, atom, cgnr, charge, mass in atoms:
    lines.append(f'{int(nr):5d}  {kind:<4s} {int(resnr):5d}  {residue:<7s}  {atom:<4s} {int(cgnr):4d} {float(charge):10.6f} '
                 f'{float(mass):9.4f}')
lines += ['', '[ settles ]', '; OW  funct      doh      dhh',
          f'{int(settles[0]):4d} {int(settles[1]):6d} {float(settles[2]):9.6f} {float(settles[3]):9.6f}', '',
          '[ exclusions ]'] + ['  '.join(row) for row in water['exclusions']] + ['']
(OUT / 'tip3p.itp').write_text('\n'.join(lines))

counts = {}
for (section, key), (source, rows) in terms.items():
    counts[(source, section)] = counts.get((source, section), 0) + len(rows)
for (source, section), count in counts.items():
    print(f'{source:22s} {section:15s} {count:4d}')
for key, first, second in replaced:
    print(f'{key} of {first} replaced by {second}')
