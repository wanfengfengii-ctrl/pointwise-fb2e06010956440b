from dataclasses import dataclass
from types import MappingProxyType
from .sites import index_sites
from .bonds import read_bonds, validate_bond_endpoints


@dataclass(frozen=True)
class PeriodicGraph:
    sites: object
    bonds: object


def build_periodic_graph(sites, bonds):
    site_index = index_sites(sites)
    bond_index = {}
    for bond in bonds:
        validate_bond_endpoints(bond, site_index)
        if bond.id in bond_index:
            raise ValueError("duplicate bond identity")
        bond_index[bond.id] = bond
    return PeriodicGraph(MappingProxyType(site_index),
                         MappingProxyType(bond_index))


def read_periodic_graph(sites, records):
    return build_periodic_graph(sites, read_bonds(records, index_sites(sites)))


def bonds_in_identity_order(graph):
    return tuple(graph.bonds[key] for key in sorted(graph.bonds))


def resolve_bond_sites(graph, bond):
    if graph.bonds.get(bond.id) != bond:
        raise ValueError("bond is not owned by this graph")
    return graph.sites[bond.u], graph.sites[bond.v]
