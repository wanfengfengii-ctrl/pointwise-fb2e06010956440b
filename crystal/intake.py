from .input_model import CellInput, snapshot_raw_input, read_cell_fields
from .exact_numbers import read_identity, json_dump
from .units import normalize_lattice_lengths
from .lattice import read_lattice
from .sites import read_sites
from .periodic_graph import read_periodic_graph, bonds_in_identity_order
from .input_digest import cell_digest


def validate_product_id(value):
    return read_identity(value)


def accept_cell(value):
    raw = snapshot_raw_input(value)
    unit, matrix, site_records, bond_records = read_cell_fields(raw)
    lattice = read_lattice(normalize_lattice_lengths(matrix, unit))
    sites = read_sites(site_records)
    graph = read_periodic_graph(sites, bond_records)
    cell = CellInput(lattice, sites, bonds_in_identity_order(graph),
                     json_dump(raw))
    return cell, cell_digest(cell)
