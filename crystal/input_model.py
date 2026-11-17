from dataclasses import dataclass
from fractions import Fraction
from .exact_numbers import exact_keys, encode_rational, json_dump, json_load


def snapshot_raw_input(value):
    def copy(item):
        if isinstance(item, Fraction):
            return encode_rational(item)
        if type(item) in (int, str):
            return item
        if type(item) in (list, tuple):
            return [copy(x) for x in item]
        if type(item) is dict and all(type(k) is str for k in item):
            return {k: copy(v) for k, v in item.items()}
        raise ValueError("unsupported input value")
    snapshot = copy(value)
    read_cell_fields(snapshot)
    return snapshot


def read_cell_fields(value):
    exact_keys(value, ("unit", "lattice", "sites", "bonds"))
    return value["unit"], value["lattice"], value["sites"], value["bonds"]


@dataclass(frozen=True)
class CellInput:
    lattice: object
    sites: tuple
    bonds: tuple
    raw_json: str

    def raw_mapping(self):
        return json_load(self.raw_json)


def cell_to_mapping(cell):
    return {
        "unit": "nm",
        "lattice": [[encode_rational(x) for x in row]
                    for row in cell.lattice.rows],
        "sites": [{"id": site.id,
                   "fractional": [encode_rational(x) for x in site.fractional]}
                  for site in sorted(cell.sites, key=lambda s: s.id)],
        "bonds": [{"id": bond.id, "u": bond.u, "v": bond.v,
                   "offset": list(bond.offset)}
                  for bond in sorted(cell.bonds, key=lambda b: b.id)],
    }
