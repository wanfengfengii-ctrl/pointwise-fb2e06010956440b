"""Finite component structure of the periodic bond graph.

A bond (id,u,v,t) connects images (u,n) and (v,n+t) for every integer
translation n. Within one connected component of the finite undirected
site graph, a spanning tree assigns every site a base translation pi(w);
each non-tree bond then closes a fundamental cycle whose integer offset
is pi(u) + t - pi(v). Closed walks at a root realise exactly the integer
lattice spanned by these chord offsets, so query reachability reduces to
one exact integer system per component plus tree-path bookkeeping.
"""
from .rational_matrix import add3, subtract3


class ComponentBasis:
    """Spanning-tree and chord data for one connected component."""

    __slots__ = ("root", "parent", "step", "pi", "chords", "chord_offsets")

    def __init__(self, root, parent, step, pi, chords, chord_offsets):
        self.root = root                # site id where pi is (0,0,0)
        self.parent = parent            # site id -> parent site id (root: None)
        self.step = step                # site id -> (bond id, sign) to parent
        self.pi = pi                    # site id -> tree offset from root
        self.chords = chords            # non-tree bonds, identity order
        self.chord_offsets = chord_offsets  # offset of each chord cycle

    def path_steps(self, site):
        """(bond id, sign) pairs from site up to the root."""
        steps = []
        current = site
        while current != self.root:
            bond_id, sign = self.step[current]
            steps.append((bond_id, sign))
            current = self.parent[current]
        return steps


def _sorted_bonds(bonds):
    return tuple(sorted(bonds, key=lambda bond: bond.id))


def build_component_basis(start, bonds):
    """Build the basis of the component containing site start.

    bonds is the revision's bond tuple; endpoints are already validated.
    Deterministic: sites are reached by identity-ordered bonds, chords keep
    identity order.
    """
    ordered = _sorted_bonds(bonds)
    parent = {start: None}
    step = {}
    pi = {start: (0, 0, 0)}
    tree_bond_ids = set()
    queue = [start]
    while queue:
        site = queue.pop(0)
        for bond in ordered:
            if bond.u == bond.v:
                continue
            if bond.u == site and bond.v not in parent:
                child, sign = bond.v, 1
            elif bond.v == site and bond.u not in parent:
                child, sign = bond.u, -1
            else:
                continue
            parent[child] = site
            step[child] = (bond.id, sign)
            offset = bond.offset if sign > 0 else tuple(-t for t in bond.offset)
            pi[child] = add3(pi[site], offset)
            tree_bond_ids.add(bond.id)
            queue.append(child)
    chords = tuple(bond for bond in ordered
                   if bond.id not in tree_bond_ids and bond.u in parent)
    chord_offsets = tuple(
        subtract3(add3(pi[bond.u], bond.offset), pi[bond.v])
        for bond in chords)
    return ComponentBasis(start, parent, step, pi, chords, chord_offsets)
