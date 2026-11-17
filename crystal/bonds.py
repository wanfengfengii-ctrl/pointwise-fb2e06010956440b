from dataclasses import dataclass
from .exact_numbers import exact_keys, read_identity, read_integer


def read_integer_offset(value):
    if type(value) not in (list, tuple) or len(value) != 3:
        raise ValueError("offset must have three integers")
    return tuple(read_integer(x, "offset component") for x in value)


@dataclass(frozen=True)
class PeriodicBond:
    id: str
    u: str
    v: str
    offset: tuple

    def __post_init__(self):
        for value in (self.id, self.u, self.v):
            read_identity(value)
        object.__setattr__(self, "offset", read_integer_offset(self.offset))


def validate_bond_endpoints(bond, site_index):
    if bond.u not in site_index or bond.v not in site_index:
        raise ValueError("bond endpoint is absent from this revision")


def read_bonds(value, site_index):
    if type(value) not in (list, tuple):
        raise ValueError("bonds must be a sequence")
    result = []
    identities = set()
    for item in value:
        exact_keys(item, ("id", "u", "v", "offset"))
        bond = PeriodicBond(item["id"], item["u"], item["v"], item["offset"])
        validate_bond_endpoints(bond, site_index)
        if bond.id in identities:
            raise ValueError("duplicate bond identity")
        identities.add(bond.id)
        result.append(bond)
    return tuple(sorted(result, key=lambda bond: bond.id))
