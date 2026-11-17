from .periodic_graph import build_periodic_graph
from .cell_geometry import compute_cell_geometry
from .bond_geometry import compute_bond_geometry
from .bond_audit import audit_site_coverage, audit_bond_ledger
from .result_codec import GeometryResult, encode_result, decode_result
from .input_digest import verify_cell_digest


def validate_geometry_result(cell, result):
    verify_cell_digest(cell, result.source[2])
    volume, _ = compute_cell_geometry(cell)
    if result.volume != volume:
        raise ValueError("result volume differs from source")
    graph = build_periodic_graph(cell.sites, cell.bonds)
    audit_site_coverage(cell.lattice, cell.sites, result.positions)
    audit_bond_ledger(cell.lattice, graph, result.bonds)
    return result


def execute_geometry(cell, source):
    verify_cell_digest(cell, source[2])
    graph = build_periodic_graph(cell.sites, cell.bonds)
    volume, positions = compute_cell_geometry(cell)
    result = GeometryResult(source, volume, positions,
                            compute_bond_geometry(cell.lattice, graph))
    validate_geometry_result(cell, result)
    return decode_result(encode_result(result))
