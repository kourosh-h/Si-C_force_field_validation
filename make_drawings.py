"""
The schematic drawings of the notebook, written to doc/ as PDF files and as PNG files, which the notebook shows: the
model molecule, the scanned bonds, angles and dihedrals, and the scans with an oxygen. They use molecule_drawing.py
and the structures in ORCA_results.

Run: python3 make_drawings.py
"""
from pathlib import Path

import matplotlib.pyplot as plt

from molecule_drawing import MoleculeDrawing

RELAXED = Path('ORCA_results/relaxation')
DOC = Path('doc')
OXYGEN = 5                 # the oxygen of the scans with an oxygen: atom 5 of the ORCA files, counted from 0
# colors: the coordinate at the Si-C bond, the same coordinate in the chain, with an oxygen
SI_C_COLOR = '#B7371B'
CHAIN_COLOR = '#0D4F8B'
OXYGEN_COLOR = '#756BB1'
# the look of the notebook: Arial, and the text sizes and lines of the drawings
plt.rcParams.update({'font.family': 'sans-serif', 'font.sans-serif': ['Arial'], 'pdf.fonttype': 42,
                     'mathtext.default': 'it', 'axes.unicode_minus': False, 'font.size': 8.5,
                     'axes.labelsize': 8.5, 'xtick.labelsize': 8.0, 'ytick.labelsize': 8.0, 'legend.fontsize': 7.5,
                     'axes.linewidth': 1.2, 'axes.spines.top': False, 'axes.spines.right': False,
                     'axes.labelpad': 2.5, 'xtick.direction': 'in', 'ytick.direction': 'in',
                     'xtick.major.size': 4.2, 'ytick.major.size': 4.2, 'xtick.major.width': 1.0,
                     'ytick.major.width': 1.0, 'xtick.minor.size': 2.0, 'ytick.minor.size': 2.0,
                     'xtick.minor.width': 0.7, 'ytick.minor.width': 0.7, 'xtick.minor.visible': True,
                     'ytick.minor.visible': True, 'xtick.major.pad': 3.0, 'ytick.major.pad': 3.0,
                     'lines.linewidth': 1.45, 'lines.markersize': 4.25, 'lines.dashed_pattern': [4.0, 2.2],
                     'lines.scale_dashes': False, 'legend.frameon': False, 'legend.handlelength': 1.8,
                     'legend.handletextpad': 0.4, 'legend.columnspacing': 1.6, 'savefig.dpi': 300})
DOC.mkdir(exist_ok=True)


def save(fig, name):
    """The drawing as doc/<name>.pdf and doc/<name>.png."""
    fig.savefig(DOC / f'{name}.pdf')
    fig.savefig(DOC / f'{name}.png', dpi=200)
    plt.close(fig)


# Figure 1: the model molecule
molecule = MoleculeDrawing(RELAXED / 'octylsilanetriol.xyz', OXYGEN)
fig, ax = molecule.new_figure(2.3)
molecule.draw(ax)
molecule.mark_bonds(ax, ['Si', 'C1'], SI_C_COLOR)
molecule.label_carbons(ax)
molecule.label_atom(ax, 'Si', r'Si$_{\rm anc}$')
molecule.label_atom(ax, 'O', r'O$_{\rm silanol}$')
molecule.label_atom(ax, 6, r'H$_{\rm silanol}$')
molecule.label_atom(ax, 24, r'H$_{\rm CH_2}$')
molecule.label_atom(ax, 31, r'H$_{\rm CH_3}$')
middle = 0.5 * (molecule.position('Si') + molecule.position('C1'))
molecule.name_mark(ax, middle, r'Si$_{\rm anc}$–C$_1$ = 0.186 nm', SI_C_COLOR)
molecule.place_name(ax, '$L$ = 8', molecule.position('C8'), [1.2, 1.5, 1.8], first=300)
molecule.atom_key(ax)
save(fig, 'molecule')

# Figure 2: the scanned bonds
fig, ax = molecule.new_figure(2.0)
molecule.draw(ax)
molecule.mark_bonds(ax, ['Si', 'C1'], SI_C_COLOR)
molecule.mark_bonds(ax, ['C4', 'C5'], CHAIN_COLOR)
molecule.label_carbons(ax)
molecule.label_atom(ax, 'Si', r'Si$_{\rm anc}$')
molecule.name_mark(ax, 0.5 * (molecule.position('Si') + molecule.position('C1')), r'Si$_{\rm anc}$–C$_1$', SI_C_COLOR)
molecule.name_mark(ax, 0.5 * (molecule.position('C4') + molecule.position('C5')), r'C$_4$–C$_5$', CHAIN_COLOR)
save(fig, 'scanned_bonds')

# Figure 4: the scanned angles
fig, ax = molecule.new_figure(2.0)
molecule.draw(ax)
arc_1 = molecule.mark_angle(ax, ['Si', 'C1', 'C2'], SI_C_COLOR)
arc_2 = molecule.mark_angle(ax, ['C3', 'C4', 'C5'], CHAIN_COLOR)
molecule.label_carbons(ax)
molecule.label_atom(ax, 'Si', r'Si$_{\rm anc}$')
molecule.name_mark(ax, arc_1, r'Si$_{\rm anc}$–C$_1$–C$_2$', SI_C_COLOR)
molecule.name_mark(ax, arc_2, r'C$_3$–C$_4$–C$_5$', CHAIN_COLOR)
save(fig, 'scanned_angles')

# Figure 6: the scanned dihedrals
fig, ax = molecule.new_figure(2.0)
molecule.draw(ax)
head_1 = molecule.mark_dihedral(ax, ['Si', 'C1', 'C2', 'C3'], SI_C_COLOR)
head_2 = molecule.mark_dihedral(ax, ['C3', 'C4', 'C5', 'C6'], CHAIN_COLOR)
molecule.label_carbons(ax)
molecule.label_atom(ax, 'Si', r'Si$_{\rm anc}$')
molecule.name_mark(ax, head_1, r'Si$_{\rm anc}$–C$_1$–C$_2$–C$_3$', SI_C_COLOR)
molecule.name_mark(ax, head_2, r'C$_3$–C$_4$–C$_5$–C$_6$', CHAIN_COLOR)
save(fig, 'scanned_dihedrals')

# Figure 8: the scans with an oxygen
fig, ax = molecule.new_figure(2.0)
molecule.draw(ax)
arc = molecule.mark_angle(ax, ['O', 'Si', 'C1'], OXYGEN_COLOR)
head = molecule.mark_dihedral(ax, ['O', 'Si', 'C1', 'C2'], OXYGEN_COLOR)
molecule.label_carbons(ax)
molecule.label_atom(ax, 'Si', r'Si$_{\rm anc}$')
molecule.label_atom(ax, 'O', 'O')
molecule.name_mark(ax, arc, r'O–Si$_{\rm anc}$–C$_1$', OXYGEN_COLOR)
molecule.name_mark(ax, head, r'O–Si$_{\rm anc}$–C$_1$–C$_2$', OXYGEN_COLOR)
save(fig, 'scanned_oxygen')
