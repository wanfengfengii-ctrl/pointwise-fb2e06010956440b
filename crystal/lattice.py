from dataclasses import dataclass
from .rational_matrix import read_matrix3, determinant3, matvec3


@dataclass(frozen=True)
class Lattice:
    rows: tuple

    def __post_init__(self):
        object.__setattr__(self, "rows", read_matrix3(self.rows))
        validate_lattice(self)

    @property
    def signed_volume(self):
        return determinant3(self.rows)

    def cartesian(self, fractional):
        return matvec3(self.rows, fractional)


def validate_lattice(lattice):
    if determinant3(lattice.rows) == 0:
        raise ValueError("singular lattice")
    return lattice


def read_lattice(value):
    return Lattice(read_matrix3(value))
