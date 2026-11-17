"""Independent audit of persisted integer-chain evidence.

This checker never re-solves reachability. It re-derives the finite
component from the source cell and verifies every reachable answer's
certificate: original bond ids of that component only, nonzero strict
integer coefficients, per-atom inflow-minus-outflow equal to the
endpoint indicators, and coefficient-weighted offsets summing exactly to
the requested translation. Unreachable answers must carry no
coefficients.
"""


def _components(cell):
    adjacency = {site.id: set() for site in cell.sites}
    for bond in cell.bonds:
        adjacency[bond.u].add(bond.v)
        adjacency[bond.v].add(bond.u)
    return adjacency


def _component_of(adjacency, start):
    component = set()
    pending = [start]
    while pending:
        site = pending.pop()
        if site not in component:
            component.add(site)
            pending.extend(adjacency[site] - component)
    return component


def audit_connectivity_answers(cell, answers):
    """Raise ValueError unless every answer's evidence holds for cell."""
    adjacency = _components(cell)
    site_ids = set(adjacency)
    bonds = {bond.id: bond for bond in cell.bonds}
    components = {}
    for answer in answers:
        if answer.u not in site_ids or answer.v not in site_ids:
            raise ValueError("answer endpoint is absent from the source")
        if not answer.reachable:
            if answer.coefficients is not None:
                raise ValueError("unreachable answer carries coefficients")
            continue
        component = components.get(answer.u)
        if component is None:
            component = _component_of(adjacency, answer.u)
            components[answer.u] = component
        if answer.v not in component:
            raise ValueError("reachable answer crosses components")
        coefficients = dict(answer.coefficients)
        balance = dict.fromkeys(site_ids, 0)
        displacement = [0, 0, 0]
        for bond_id, coefficient in coefficients.items():
            bond = bonds.get(bond_id)
            if bond is None:
                raise ValueError("coefficient references an unknown bond")
            if type(coefficient) is not int or coefficient == 0:
                raise ValueError("coefficients are nonzero strict integers")
            if bond.u not in component or bond.v not in component:
                raise ValueError("coefficient borrows another component")
            balance[bond.u] -= coefficient
            balance[bond.v] += coefficient
            for axis in range(3):
                displacement[axis] += coefficient * bond.offset[axis]
        expected = {site: int(site == answer.v) - int(site == answer.u)
                    for site in site_ids}
        if balance != expected:
            raise ValueError("coefficient endpoint balance mismatch")
        if displacement != list(answer.translation):
            raise ValueError("coefficient offset sum mismatch")
    return answers
