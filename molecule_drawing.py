"""
Drawings of the model molecule from its relaxed ORCA structure, used by the notebook: the atoms as shaded spheres,
the bonds as gray lines, the names of the atoms, and the bonds, angles and dihedrals of the scans in color.

One drawing unit is 0.1 nm, and every drawing is 7.0 inches wide, with the same scale. The molecule lies along x,
seen from the side of its carbon zigzag and turned by VIEW_TURN about the vertical axis, so that the atoms behind
each other separate. Every name is placed where it overlaps no atom, bond, mark or other name.
"""
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import to_rgb
from matplotlib.patches import Circle, Polygon

WIDTH = 7.0                 # inches
UNITS_PER_INCH = 3.7066     # drawing units (0.1 nm) per inch
FONT_SIZE = 8.5             # pt, all texts
ATOM_SCALE = 0.75           # size of the atoms
BOND_WIDTH = 1.3            # pt, the bonds
MARK_WIDTH = 2.4            # pt, the bonds of a scanned coordinate
VIEW_TURN = 30.0            # degrees
GAP = 0.18                  # drawing units kept free around every name
# color and radius (drawing units) of the atoms
ELEMENTS = {'Si': ('#E3B23C', 0.55), 'O': ('#D2473A', 0.45), 'C': ('#5F6B73', 0.45), 'H': ('#F4F4F4', 0.27)}


def read_xyz(path):
    """Elements and positions (0.1 nm) of an xyz file."""
    lines = Path(path).read_text().split('\n')
    elements = []
    positions = []
    for line in lines[2:2 + int(lines[0])]:
        parts = line.split()
        elements.append(parts[0])
        positions.append([float(parts[1]), float(parts[2]), float(parts[3])])
    return elements, np.array(positions)


def find_bonds(elements, positions):
    """Pairs of bonded atoms: closer than 0.12 nm for bonds to H, 0.20 nm for bonds to Si and 0.17 nm for C-C."""
    bonds = []
    for i in range(len(elements)):
        for j in range(i + 1, len(elements)):
            pair = [elements[i], elements[j]]
            if pair == ['H', 'H']:
                continue
            if 'H' in pair:
                limit = 1.2
            elif 'Si' in pair:
                limit = 2.0
            else:
                limit = 1.7
            if np.linalg.norm(positions[i] - positions[j]) < limit:
                bonds.append((i, j))
    return bonds


def neighbors_of(atom, bonds):
    """The atoms bonded to an atom."""
    neighbors = []
    for i, j in bonds:
        if i == atom:
            neighbors.append(j)
        if j == atom:
            neighbors.append(i)
    return neighbors


def chain_of(elements, bonds):
    """The silicon and the carbons of the chain in their order: Si_anc, C1, ..., C8."""
    chain = [elements.index('Si')]
    while True:
        following = []
        for atom in neighbors_of(chain[-1], bonds):
            if elements[atom] == 'C' and atom not in chain:
                following.append(atom)
        if len(following) == 0:
            return chain
        chain.append(following[0])


def atom_radius(element):
    """Radius of a drawn atom (drawing units)."""
    radius = ELEMENTS[element][1]
    if element == 'H':
        return radius * ATOM_SCALE * 0.7
    return radius * ATOM_SCALE * 0.62


def draw_atom(ax, position, element, zorder, radius=None):
    """Atom as a shaded sphere: a darker disk with a dark rim and nine lighter disks shifted to the upper left."""
    if radius is None:
        radius = atom_radius(element)
    base = np.array(to_rgb(ELEMENTS[element][0]))
    ax.add_patch(Circle(position, radius, facecolor=base * 0.78, edgecolor='0.15', lw=0.5, zorder=zorder))
    for k in range(1, 10):
        fraction = k / 9
        shade = base * 0.78 + (np.minimum(base + 0.35, 1) - base * 0.78) * fraction
        center = (position[0] - 0.28 * radius * fraction, position[1] + 0.28 * radius * fraction)
        ax.add_patch(Circle(center, radius * (1 - 0.62 * fraction), facecolor=shade, lw=0, zorder=zorder))


def unit(vector):
    """The vector divided by its length."""
    return vector / np.linalg.norm(vector)


def point_to_box(point, box):
    """Distance of a point from a box (x0, y0, x1, y1); 0 inside it."""
    dx = max(box[0] - point[0], 0, point[0] - box[2])
    dy = max(box[1] - point[1], 0, point[1] - box[3])
    return np.hypot(dx, dy)


def segment_to_box(a, b, box):
    """Smallest distance between the segment a-b and a box, from 25 points along the segment."""
    smallest = 100.0
    for t in np.linspace(0, 1, 25):
        smallest = min(smallest, point_to_box(a + t * (b - a), box))
    return smallest


def arc_points(tail, head, rad):
    """Nine points along the curve of a matplotlib arc3 connection from tail to head (a quadratic Bezier curve)."""
    control = 0.5 * (tail + head) + rad * np.array([head[1] - tail[1], tail[0] - head[0]])
    points = []
    for t in np.linspace(0, 1, 9):
        points.append((1 - t) ** 2 * tail + 2 * (1 - t) * t * control + t ** 2 * head)
    return points


def segments_cross(a, b, c, d):
    """True if the segments a-b and c-d cross."""
    def side(p, q, r):
        return (q[0] - p[0]) * (r[1] - p[1]) - (q[1] - p[1]) * (r[0] - p[0])
    return side(a, b, c) * side(a, b, d) < 0 and side(c, d, a) * side(c, d, b) < 0


def segment_to_point(a, b, point):
    """Distance of a point from the segment a-b."""
    t = np.clip((point - a) @ (b - a) / ((b - a) @ (b - a)), 0, 1)
    return np.linalg.norm(point - (a + t * (b - a)))


class MoleculeDrawing:
    """The relaxed molecule in the plane of the drawing, with the names of its atoms: Si, C1 ... C8 and O."""

    def __init__(self, path, oxygen, normal=None):
        """
        The molecule of an .xyz file; oxygen: the number of the oxygen named O. normal: for a cluster on the silica,
        the outward normal of the surface, drawn upwards (otherwise the molecule lies along x).
        """
        self.elements, positions = read_xyz(path)
        self.bonds = find_bonds(self.elements, positions)
        self.chain = chain_of(self.elements, self.bonds)
        if normal is None:
            self.xy, self.depth = self.drawing_positions(positions)
        else:
            self.xy, self.depth = self.side_view(positions, np.array(normal, dtype=float))
        # Si, C1 ... C8, and O, the oxygen given (counted from 0 in the xyz file)
        self.atoms = {'Si': self.chain[0], 'O': oxygen}
        for k in range(1, len(self.chain)):
            self.atoms[f'C{k}'] = self.chain[k]

    def drawing_positions(self, positions):
        """
        Positions in the plane of the drawing and their depth: x along the chain from C1 to C8 and y in the plane of
        the carbon zigzag, with the molecule turned by VIEW_TURN about the vertical axis y.
        """
        chain = self.chain
        x_axis = unit(positions[chain[-1]] - positions[chain[1]])
        zigzag = positions[chain[2]] - positions[chain[1]]
        y_axis = unit(zigzag - (zigzag @ x_axis) * x_axis)
        z_axis = np.cross(x_axis, y_axis)
        turn = np.radians(VIEW_TURN)
        x_turned = np.cos(turn) * x_axis + np.sin(turn) * z_axis
        z_turned = np.cross(x_turned, y_axis)
        relative = positions - positions[chain[0]]
        return relative @ np.array([x_turned, y_axis]).T, relative @ z_turned

    def side_view(self, positions, normal):
        """
        Positions in the plane of the drawing and their depth for a cluster on the silica: y along the surface normal,
        x along the part of the Si-C1 bond parallel to the surface, turned by VIEW_TURN about y.
        """
        y_axis = unit(normal)
        bond = positions[self.chain[1]] - positions[self.chain[0]]
        x_axis = unit(bond - (bond @ y_axis) * y_axis)
        z_axis = np.cross(x_axis, y_axis)
        turn = np.radians(VIEW_TURN)
        x_turned = np.cos(turn) * x_axis + np.sin(turn) * z_axis
        z_turned = np.cross(x_turned, y_axis)
        relative = positions - positions[self.chain[0]]
        return relative @ np.array([x_turned, y_axis]).T, relative @ z_turned

    def index(self, atom):
        """Number of an atom in the xyz file, from its name in the drawing (Si, C1 ... C8, O) or its number."""
        if atom in self.atoms:
            return self.atoms[atom]
        return atom

    def position(self, atom):
        """Position of an atom in the drawing."""
        return self.xy[self.index(atom)]

    def new_figure(self, height, shift=0.0):
        """
        Figure of the full width and the given height (inches), at the scale of all drawings, with the molecule in
        the middle (moved to the left by shift drawing units). Starts the lists of the drawn shapes that the names
        must keep clear of: the atoms (circles), the bonds and other lines (segments) and the names (boxes).
        """
        fig = plt.figure(figsize=(WIDTH, height))
        ax = fig.add_axes([0, 0, 1, 1])
        ax.set_aspect('equal')
        ax.axis('off')
        center_x = 0.5 * (self.xy[:, 0].min() + self.xy[:, 0].max()) + shift
        center_y = 0.5 * (self.xy[:, 1].min() + self.xy[:, 1].max())
        half_width = 0.5 * WIDTH * UNITS_PER_INCH
        half_height = 0.5 * height * UNITS_PER_INCH
        ax.set_xlim(center_x - half_width, center_x + half_width)
        ax.set_ylim(center_y - half_height, center_y + half_height)
        self.circles = []
        self.segments = []
        self.boxes = []
        for i in range(len(self.elements)):
            self.circles.append((self.xy[i], atom_radius(self.elements[i])))
        for i, j in self.bonds:
            self.segments.append((self.xy[i], self.xy[j]))
        return fig, ax

    def draw(self, ax):
        """All atoms, the farthest first, and the bonds as gray lines just below the lower of their two atoms."""
        order = np.argsort(self.depth)
        zorder = {}
        for rank in range(len(order)):
            zorder[order[rank]] = 3 + 0.02 * rank
        for i, j in self.bonds:
            ax.plot([self.xy[i, 0], self.xy[j, 0]], [self.xy[i, 1], self.xy[j, 1]], color='0.45', lw=BOND_WIDTH,
                    solid_capstyle='round', zorder=min(zorder[i], zorder[j]) - 0.01)
        for i in range(len(self.elements)):
            draw_atom(ax, self.xy[i], self.elements[i], zorder[i])

    def room(self, box, own_atoms=()):
        """Free room around a box: its distance to the nearest atom, line or name (the atoms own_atoms excepted)."""
        room = 100.0
        for k in range(len(self.circles)):
            if k in own_atoms:
                continue
            center, radius = self.circles[k]
            room = min(room, point_to_box(center, box) - radius)
        for a, b in self.segments:
            room = min(room, segment_to_box(a, b, box))
        for other in self.boxes:
            gap_x = max(other[0] - box[2], box[0] - other[2], 0)
            gap_y = max(other[1] - box[3], box[1] - other[3], 0)
            room = min(room, max(gap_x, gap_y))
        return room

    def text_box(self, ax, text, x, y):
        """Box (x0, y0, x1, y1) of a name centered at (x, y), in drawing units, from its printed size."""
        renderer = ax.figure.canvas.get_renderer()
        extent = text.get_window_extent(renderer)
        to_data = ax.transData.inverted()
        corner_low = to_data.transform((extent.x0, extent.y0))
        corner_high = to_data.transform((extent.x1, extent.y1))
        half_width = 0.5 * (corner_high[0] - corner_low[0])
        half_height = 0.5 * (corner_high[1] - corner_low[1])
        return (x - half_width, y - half_height, x + half_width, y + half_height)

    def inside(self, ax, box):
        """True if the box lies inside the drawing."""
        x0, x1 = ax.get_xlim()
        y0, y1 = ax.get_ylim()
        return box[0] > x0 + 0.2 and box[2] < x1 - 0.2 and box[1] > y0 + 0.1 and box[3] < y1 - 0.1

    def place_name(self, ax, text, anchor, distances, color='black', own_atoms=(), leader=False, first=None):
        """
        Write a name near the point anchor: at the first of the places tried (the given distances, every 30 degrees,
        starting at the direction first, in degrees, if given) where it keeps GAP from every atom, line and name;
        otherwise at the place with the most room. leader=True draws a thin line from the anchor to the name, which
        must not cross an atom or another line. Returns the place.
        """
        label = ax.text(0, 0, text, color=color, ha='center', va='center', fontsize=FONT_SIZE, zorder=12)
        if first is None:
            first = 0
        # the direction first, then the directions next to it on both sides
        directions = []
        for step in range(0, 181, 30):
            for degrees in [(first + step) % 360, (first - step) % 360]:
                if degrees not in directions:
                    directions.append(degrees)
        best = None
        for distance in distances:
            for degrees in directions:
                direction = np.array([np.cos(np.radians(degrees)), np.sin(np.radians(degrees))])
                place = anchor + distance * direction
                box = self.text_box(ax, label, place[0], place[1])
                if not self.inside(ax, box):
                    continue
                room = self.room(box, own_atoms)
                if leader:
                    end = self.leader_end(anchor, box)
                    room = min(room, self.leader_room(anchor, end, own_atoms))
                if best is None or room > best[0]:
                    best = (room, place, box)
                if room >= GAP:
                    break
            if best is not None and best[0] >= GAP:
                break
        room, place, box = best
        label.set_position(place)
        self.boxes.append(box)
        if leader:
            end = self.leader_end(anchor, box)
            ax.plot([anchor[0], end[0]], [anchor[1], end[1]], color=color, lw=0.7, zorder=8.4)
            self.segments.append((anchor, end))
        return place

    def leader_end(self, anchor, box):
        """The point of a box nearest to the anchor, a little outside the box."""
        x = np.clip(anchor[0], box[0] - 0.08, box[2] + 0.08)
        y = np.clip(anchor[1], box[1] - 0.08, box[3] + 0.08)
        return np.array([x, y])

    def leader_room(self, anchor, end, own_atoms):
        """Free room along a leader from the anchor to end, the first 0.4 drawing units excepted."""
        length = np.linalg.norm(end - anchor)
        if length < 0.5:
            return 100.0
        start = anchor + 0.4 * (end - anchor) / length
        room = 100.0
        for k in range(len(self.circles)):
            if k in own_atoms:
                continue
            center, radius = self.circles[k]
            room = min(room, segment_to_point(start, end, center) - radius)
        for other in self.boxes:
            room = min(room, segment_to_box(start, end, other))
        # a leader crosses no bond and no other line
        for a, b in self.segments:
            if segments_cross(start, end, a, b):
                room = -1.0
        return room

    def label_atom(self, ax, atom, text):
        """Name of an atom just outside it, first tried on the side away from the atoms bonded to it."""
        k = self.index(atom)
        radius = atom_radius(self.elements[k])
        neighbors = neighbors_of(k, self.bonds)
        outward = self.xy[k] - np.mean(self.xy[neighbors], axis=0)
        first = round(np.degrees(np.arctan2(outward[1], outward[0])) / 30) * 30
        self.place_name(ax, text, self.xy[k], [radius + 0.35, radius + 0.5, radius + 0.7, radius + 0.95],
                        own_atoms=(k,), first=first)

    def label_carbons(self, ax):
        """C1 ... C8 next to their atoms, first tried on the inner side of the zigzag."""
        for k in range(1, len(self.chain)):
            before = self.xy[self.chain[k - 1], 1]
            if k + 1 < len(self.chain):
                after = self.xy[self.chain[k + 1], 1]
            else:
                after = before
            if self.xy[self.chain[k], 1] > 0.5 * (before + after):
                first = 270
            else:
                first = 90
            radius = atom_radius('C')
            self.place_name(ax, f'C$_{k}$', self.xy[self.chain[k]],
                            [radius + 0.4, radius + 0.55, radius + 0.75, radius + 1.0, radius + 1.3],
                            own_atoms=(self.chain[k],), first=first)

    def mark_bonds(self, ax, atoms, color):
        """The bonds between the given atoms, in their order, as thick lines in the color of the coordinate."""
        for k in range(len(atoms) - 1):
            a, b = self.position(atoms[k]), self.position(atoms[k + 1])
            ax.plot([a[0], b[0]], [a[1], b[1]], color=color, lw=MARK_WIDTH, solid_capstyle='round', zorder=8.5)

    def mark_angle(self, ax, atoms, color, radius=0.75):
        """An angle a-b-c: its two bonds and a shaded wedge with an arc at b. Returns the middle of the arc."""
        self.mark_bonds(ax, atoms, color)
        a, b, c = self.position(atoms[0]), self.position(atoms[1]), self.position(atoms[2])
        first = np.arctan2(a[1] - b[1], a[0] - b[0])
        second = np.arctan2(c[1] - b[1], c[0] - b[0])
        # the wedge runs the short way from one bond to the other
        sweep = (second - first) % (2 * np.pi)
        if sweep > np.pi:
            first, sweep = second, 2 * np.pi - sweep
        angles = first + np.linspace(0, sweep, 40)
        rim = b + radius * np.column_stack([np.cos(angles), np.sin(angles)])
        ax.add_patch(Polygon(np.vstack([b, rim]), closed=True, facecolor=color, alpha=0.35, lw=0, zorder=8.3))
        ax.plot(rim[:, 0], rim[:, 1], color=color, lw=1.1, zorder=8.6)
        for k in range(len(rim) - 1):
            self.segments.append((rim[k], rim[k + 1]))
        return rim[20]

    def mark_dihedral(self, ax, atoms, color, half=0.62):
        """
        A dihedral a-b-c-d: its three bonds and a curved arrow across the middle bond b-c, the rotation about it.
        Returns the head of the arrow.
        """
        self.mark_bonds(ax, atoms, color)
        p, q = self.position(atoms[1]), self.position(atoms[2])
        middle = 0.5 * (p + q)
        along = unit(q - p)
        normal = np.array([-along[1], along[0]])
        tail = middle - half * normal + 0.22 * along
        head = middle + half * normal - 0.22 * along
        ax.annotate('', xy=head, xytext=tail, zorder=10,
                    arrowprops=dict(arrowstyle='->', lw=1.1, color=color, connectionstyle='arc3,rad=0.6',
                                    mutation_scale=8, shrinkA=0, shrinkB=0))
        # the arc of the arrow as short lines, for the names placed later
        points = arc_points(tail, head, 0.6)
        for k in range(len(points) - 1):
            self.segments.append((points[k], points[k + 1]))
        return head

    def name_mark(self, ax, anchor, text, color, distances=(1.2, 1.5, 1.8, 2.2, 2.6, 3.0, 3.5)):
        """Name of a scanned coordinate in its color, near the point anchor of its mark, with a thin leader."""
        self.place_name(ax, text, anchor, distances, color=color, leader=True)

    def atom_key(self, ax, step=1.0):
        """Si, O, C and H in one column right of the molecule, with the names of the elements."""
        x = self.xy[:, 0].max() + 1.6
        y = self.xy[:, 1].max()
        for k, element in enumerate(['Si', 'O', 'C', 'H']):
            radius = ELEMENTS[element][1] * ATOM_SCALE * 0.62
            if element == 'H':
                radius = radius * 1.13
            draw_atom(ax, (x, y - step * k), element, 12, radius)
            ax.text(x + 0.4, y - step * k, element, ha='left', va='center', color='0.25', fontsize=FONT_SIZE,
                    zorder=12)
            self.boxes.append((x - 0.3, y - step * k - 0.3, x + 1.2, y - step * k + 0.3))
