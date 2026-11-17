from fractions import Fraction
from types import MappingProxyType
from .rational_matrix import read_matrix3

CANONICAL_LENGTH_UNIT = "nm"
LENGTH_SCALES = MappingProxyType({
    "nm": Fraction(1), "angstrom": Fraction(1, 10),
    "pm": Fraction(1, 1000),
})


def length_scale(unit):
    if type(unit) is not str or unit not in LENGTH_SCALES:
        raise ValueError("unsupported length unit")
    return LENGTH_SCALES[unit]


def normalize_lattice_lengths(matrix, unit):
    scale = length_scale(unit)
    return tuple(tuple(x * scale for x in row)
                 for row in read_matrix3(matrix))
