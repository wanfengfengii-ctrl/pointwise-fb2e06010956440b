"""Query input boundaries for the connectivity calculation.

Only this calculation enforces the 16 site / 32 bond / 16 query limits;
intake and the V1 geometry chain stay unrestricted. Every rejection is a
ValueError raised before any database write, so an illegal batch leaves
no task and no result.
"""
from .bonds import read_integer_offset
from .exact_numbers import exact_keys, read_identity

MAX_SITES = 16
MAX_BONDS = 32
MAX_QUERIES = 16


def check_connectivity_limits(cell):
    if len(cell.sites) > MAX_SITES or len(cell.bonds) > MAX_BONDS:
        raise ValueError(
            "connectivity accepts at most 16 sites and 32 bonds")


def read_connectivity_queries(value, site_ids):
    """Validate the query batch against the fixed revision's site ids.

    Returns a tuple of canonical (u, v, translation) triples in input
    order; translation is a tuple of three strict integers. Duplicates are
    preserved. Any defect rejects the whole batch with ValueError.
    """
    if type(value) is not list or len(value) > MAX_QUERIES:
        raise ValueError("queries must be a list of at most 16 items")
    queries = []
    for item in value:
        exact_keys(item, ("u", "v", "translation"))
        u = read_identity(item["u"])
        v = read_identity(item["v"])
        if u not in site_ids or v not in site_ids:
            raise ValueError("query endpoint is absent from this revision")
        translation = read_integer_offset(item["translation"])
        queries.append((u, v, translation))
    return tuple(queries)
