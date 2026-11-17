def cell_volume(lattice):
    return abs(lattice.signed_volume)


def cartesian_positions(lattice, sites):
    return tuple((site.id, lattice.cartesian(site.fractional))
                 for site in sorted(sites, key=lambda site: site.id))


def compute_cell_geometry(cell):
    return cell_volume(cell.lattice), cartesian_positions(cell.lattice, cell.sites)
