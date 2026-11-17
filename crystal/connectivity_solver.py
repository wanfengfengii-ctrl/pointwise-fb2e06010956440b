"""Exact integer reachability between periodic images.

For a query (u, v, translation) the tree path from u to v contributes the
base offset sigma = pi(v) - pi(u); the remaining translation - sigma must
be an integer combination of the component's chord cycle offsets. That
condition is decided exactly by the Smith normal form, and the solution
coefficients are unfolded back onto original bond ids: tree path steps
plus, for each used chord, its fundamental closed walk. Walking is
bidirectional; a negative coefficient records net backward use of the
same original bond.
"""
from .connectivity_graph import build_component_basis
from .connectivity_result import ConnectivityAnswer
from .rational_matrix import subtract3
from .smith_normal_form import prepare_integer_system, solve_prepared

ZERO = (0, 0, 0)


def _chord_system(basis):
    """Chord cycle offsets as a prepared 3 x len(chords) integer system."""
    matrix = [[offset[axis] for offset in basis.chord_offsets]
              for axis in range(3)]
    return prepare_integer_system(matrix)


def _accumulate(coefficients, bond_id, delta):
    coefficients[bond_id] = coefficients.get(bond_id, 0) + delta


def _answer_coefficients(basis, u, v, chord_uses):
    """Net forward uses per original bond id for one reachable query."""
    coefficients = {}
    for bond_id, sign in basis.path_steps(v):
        _accumulate(coefficients, bond_id, sign)
    for bond_id, sign in basis.path_steps(u):
        _accumulate(coefficients, bond_id, -sign)
    for bond, count in zip(basis.chords, chord_uses):
        if not count:
            continue
        _accumulate(coefficients, bond.id, count)
        for bond_id, sign in basis.path_steps(bond.u):
            _accumulate(coefficients, bond_id, count * sign)
        for bond_id, sign in basis.path_steps(bond.v):
            _accumulate(coefficients, bond_id, -count * sign)
    return tuple(sorted((bond_id, coefficient)
                        for bond_id, coefficient in coefficients.items()
                        if coefficient))


def solve_connectivity(cell, queries):
    """Solve canonical (u, v, translation) queries against the fixed cell.

    Returns a tuple of ConnectivityAnswer in input order. Every answer is
    exact: reachability is decided by integer elimination, never by real
    spans, floating point tolerances or a fixed expansion window.
    """
    bases = {}
    systems = {}
    answers = []
    for u, v, translation in queries:
        basis = bases.get(u)
        if basis is None:
            basis = build_component_basis(u, cell.bonds)
            for member in basis.pi:
                bases[member] = basis
        if v not in basis.pi:
            answers.append(ConnectivityAnswer(u, v, translation, False, None))
            continue
        sigma = subtract3(basis.pi[v], basis.pi[u])
        remaining = subtract3(translation, sigma)
        if basis.chords:
            system = systems.get(basis.root)
            if system is None:
                system = _chord_system(basis)
                systems[basis.root] = system
            chord_uses = solve_prepared(system, list(remaining))
        else:
            chord_uses = () if remaining == ZERO else None
        if chord_uses is None:
            answers.append(ConnectivityAnswer(u, v, translation, False, None))
            continue
        answers.append(ConnectivityAnswer(
            u, v, translation, True,
            _answer_coefficients(basis, u, v, chord_uses)))
    return tuple(answers)
