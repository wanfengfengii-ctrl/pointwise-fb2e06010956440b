from dataclasses import dataclass
from .periodic_graph import bonds_in_identity_order, resolve_bond_sites
from .rational_matrix import add3, subtract3, dot3


@dataclass(frozen=True)
class BondGeometry:
    id: str
    u: str
    v: str
    offset: tuple
    displacement: tuple
    length_squared: object


def bond_displacement(lattice, u, v, offset):
    fractional = subtract3(add3(v.fractional, offset), u.fractional)
    return lattice.cartesian(fractional)


def bond_length_squared(displacement):
    return dot3(displacement, displacement)


def compute_bond_geometry(lattice, graph):
    rows = []
    for bond in bonds_in_identity_order(graph):
        u, v = resolve_bond_sites(graph, bond)
        displacement = bond_displacement(lattice, u, v, bond.offset)
        rows.append(BondGeometry(
            bond.id, bond.u, bond.v, bond.offset, displacement,
            bond_length_squared(displacement)))
    return tuple(rows)
