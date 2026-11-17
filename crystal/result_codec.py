from dataclasses import dataclass
from fractions import Fraction
from .exact_numbers import (
    exact_keys, read_identity, read_integer, encode_rational,
    decode_rational, json_dump, json_load)
from .bonds import read_integer_offset
from .bond_geometry import BondGeometry

RESULT_FORMAT_VERSION = 1
CALCULATION = "primitive-v1"


def source_to_mapping(source):
    product, revision, digest, calculation = source
    read_identity(product)
    read_integer(revision, "revision")
    if revision < 1 or type(digest) is not str or len(digest) != 64:
        raise ValueError("invalid source binding")
    if any(c not in "0123456789abcdef" for c in digest):
        raise ValueError("invalid source digest")
    if calculation != CALCULATION:
        raise ValueError("unsupported calculation")
    return {"product_id": product, "revision": revision,
            "digest": digest, "calculation": calculation}


def read_source(value):
    exact_keys(value, ("product_id", "revision", "digest", "calculation"))
    source = tuple(value[k] for k in
                   ("product_id", "revision", "digest", "calculation"))
    source_to_mapping(source)
    return source


@dataclass(frozen=True)
class GeometryResult:
    source: tuple
    volume: Fraction
    positions: tuple
    bonds: tuple


def result_to_mapping(result):
    return {
        "source": source_to_mapping(result.source),
        "units": {"length": "nm", "length_squared": "nm2", "volume": "nm3"},
        "volume": result.volume,
        "positions": [{"id": identity, "cartesian": list(position)}
                      for identity, position in result.positions],
        "bonds": [{"id": row.id, "u": row.u, "v": row.v,
                   "offset": list(row.offset),
                   "displacement": list(row.displacement),
                   "length_squared": row.length_squared}
                  for row in result.bonds],
    }


def encode_result(result):
    value = result_to_mapping(result)
    value["volume"] = encode_rational(value["volume"])
    for site in value["positions"]:
        site["cartesian"] = [encode_rational(x) for x in site["cartesian"]]
    for bond in value["bonds"]:
        bond["displacement"] = [encode_rational(x) for x in bond["displacement"]]
        bond["length_squared"] = encode_rational(bond["length_squared"])
    text = json_dump({"version": RESULT_FORMAT_VERSION, "result": value})
    decode_result(text)
    return text


def _vector(value):
    if type(value) is not list or len(value) != 3:
        raise ValueError("incomplete result vector")
    return tuple(decode_rational(x) for x in value)


def decode_result(text):
    payload = json_load(text)
    exact_keys(payload, ("version", "result"))
    if type(payload["version"]) is not int or payload["version"] != 1:
        raise ValueError("unsupported result format")
    value = exact_keys(payload["result"],
                       ("source", "units", "volume", "positions", "bonds"))
    source = read_source(value["source"])
    if value["units"] != {"length": "nm", "length_squared": "nm2", "volume": "nm3"}:
        raise ValueError("incorrect result units")
    volume = decode_rational(value["volume"])
    if volume <= 0:
        raise ValueError("nonpositive volume")
    if type(value["positions"]) is not list or type(value["bonds"]) is not list:
        raise ValueError("incomplete result ledger")
    positions, bonds = [], []
    site_ids, bond_ids = set(), set()
    for site in value["positions"]:
        exact_keys(site, ("id", "cartesian"))
        identity = read_identity(site["id"])
        if identity in site_ids:
            raise ValueError("duplicate result site")
        site_ids.add(identity)
        positions.append((identity, _vector(site["cartesian"])))
    for row in value["bonds"]:
        exact_keys(row, ("id", "u", "v", "offset",
                         "displacement", "length_squared"))
        identity = read_identity(row["id"])
        u, v = read_identity(row["u"]), read_identity(row["v"])
        if identity in bond_ids or u not in site_ids or v not in site_ids:
            raise ValueError("invalid result bond identity or endpoint")
        bond_ids.add(identity)
        displacement = _vector(row["displacement"])
        length = decode_rational(row["length_squared"])
        if length != sum((x*x for x in displacement), Fraction(0)):
            raise ValueError("inconsistent result length")
        bonds.append(BondGeometry(identity, u, v,
                                  read_integer_offset(row["offset"]),
                                  displacement, length))
    return GeometryResult(source, volume,
                          tuple(sorted(positions)),
                          tuple(sorted(bonds, key=lambda row: row.id)))
