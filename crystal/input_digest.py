from hashlib import sha256
from .input_model import cell_to_mapping
from .exact_numbers import json_dump


def canonical_cell_bytes(cell):
    return json_dump(cell_to_mapping(cell)).encode("ascii")


def cell_digest(cell):
    return sha256(canonical_cell_bytes(cell)).hexdigest()


def verify_cell_digest(cell, expected):
    actual = cell_digest(cell)
    if type(expected) is not str or actual != expected:
        raise ValueError("canonical source digest mismatch")
    return actual
