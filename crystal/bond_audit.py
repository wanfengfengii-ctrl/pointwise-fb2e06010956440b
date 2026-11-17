from fractions import Fraction
from .periodic_graph import bonds_in_identity_order, resolve_bond_sites
from .bond_geometry import bond_displacement, bond_length_squared


def ledger_identity_map(rows):
    result = {}
    for row in rows:
        if row.id in result:
            raise ValueError("duplicate ledger identity")
        result[row.id] = row
    return result


def audit_site_coverage(lattice, sites, positions):
    indexed = {}
    for identity, position in positions:
        if identity in indexed or len(position) != 3:
            raise ValueError("duplicate or incomplete position")
        if any(not isinstance(x, Fraction) for x in position):
            raise ValueError("position must be exact")
        indexed[identity] = position
    if set(indexed) != {site.id for site in sites}:
        raise ValueError("site coverage differs from source")
    for site in sites:
        if indexed[site.id] != lattice.cartesian(site.fractional):
            raise ValueError("site position differs from source")


def audit_bond_ledger(lattice, graph, rows):
    indexed = ledger_identity_map(rows)
    if set(indexed) != set(graph.bonds):
        raise ValueError("bond coverage differs from source")
    for bond in bonds_in_identity_order(graph):
        row = indexed[bond.id]
        if (row.u, row.v, row.offset) != (bond.u, bond.v, bond.offset):
            raise ValueError("bond source differs")
        u, v = resolve_bond_sites(graph, bond)
        expected = bond_displacement(lattice, u, v, bond.offset)
        if len(row.displacement) != 3 or any(
                not isinstance(x, Fraction) for x in row.displacement):
            raise ValueError("displacement must be exact")
        if not isinstance(row.length_squared, Fraction):
            raise ValueError("length squared must be exact")
        if row.displacement != expected:
            raise ValueError("bond displacement differs from source")
        if row.length_squared != bond_length_squared(expected):
            raise ValueError("bond length squared differs from displacement")
