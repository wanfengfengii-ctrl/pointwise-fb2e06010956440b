from dataclasses import dataclass
from .exact_numbers import exact_keys, read_identity
from .rational_matrix import read_vector3


def validate_fractional_position(value):
    position = read_vector3(value)
    if any(x < 0 or x >= 1 for x in position):
        raise ValueError("fractional coordinates must lie in [0,1)")
    return position


@dataclass(frozen=True)
class Site:
    id: str
    fractional: tuple

    def __post_init__(self):
        read_identity(self.id)
        object.__setattr__(self, "fractional",
                           validate_fractional_position(self.fractional))


def index_sites(sites):
    result = {}
    for site in sites:
        if site.id in result:
            raise ValueError("duplicate site identity")
        result[site.id] = site
    return result


def read_sites(value):
    if type(value) not in (list, tuple):
        raise ValueError("sites must be a sequence")
    sites = []
    for item in value:
        exact_keys(item, ("id", "fractional"))
        sites.append(Site(item["id"], item["fractional"]))
    index_sites(sites)
    return tuple(sorted(sites, key=lambda site: site.id))
