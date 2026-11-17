import copy
from fractions import Fraction as F
import unittest
from crystal import Archive
from crystal.intake import accept_cell
from crystal.input_digest import cell_digest


def unit_cell():
    return {"unit": "nm", "lattice": [[1,0,0],[0,1,0],[0,0,1]],
            "sites": [{"id": "u", "fractional": [0,0,0]},
                      {"id": "v", "fractional": [[1,2],[1,2],0]}],
            "bonds": [{"id": "uv", "u": "u", "v": "v",
                       "offset": [0,1,0]}]}


def skew_cell():
    cell = unit_cell()
    cell["lattice"] = [[1,[1,2],0],[0,1,0],[0,0,2]]
    return cell


class TestIntake(unittest.TestCase):
    def test_unit_cell(self):
        cell = unit_cell()
        with Archive() as archive:
            revision = archive.put("sample", cell)
            cell["sites"][0]["fractional"][0] = 1
            cell["bonds"].clear()
            read = archive.read_revision("sample", revision)
            self.assertEqual(read["sites"][0]["fractional"], [[0,1]]*3)
            self.assertEqual(read["bonds"][0]["offset"], [0,1,0])
            read["lattice"][0][0][0] = 99
            self.assertEqual(archive.read_revision("sample", revision)
                             ["lattice"][0][0], [1,1])
            self.assertEqual(archive.current_revision("sample"), revision)
            second = archive.put("sample", unit_cell())
            self.assertEqual(second, revision+1)
            self.assertEqual(archive.revisions.read_revision("sample", revision)
                             .digest,
                             archive.revisions.read_revision("sample", second)
                             .digest)

    def test_unit_scaling(self):
        baseline, digest = accept_cell(skew_cell())
        for unit, scale in (("nm",1), ("angstrom",10), ("pm",1000)):
            cell = skew_cell()
            cell["unit"] = unit
            cell["lattice"] = [[F(x[0],x[1])*scale
                                 if isinstance(x,list) else x*scale
                                 for x in row] for row in cell["lattice"]]
            accepted, actual_digest = accept_cell(cell)
            self.assertEqual(actual_digest, digest)
            self.assertEqual(accepted.lattice, baseline.lattice)
            self.assertEqual(accepted.sites, baseline.sites)
            self.assertEqual(accepted.bonds, baseline.bonds)
            self.assertEqual(accepted.raw_mapping()["unit"], unit)

    def test_site_domain_and_types(self):
        bad = []
        cell = unit_cell(); cell["sites"].append(copy.deepcopy(cell["sites"][0]))
        bad.append(cell)
        for coordinate in (-1, 1, [3,2], True, 0.5, [1,0], None):
            cell = unit_cell(); cell["sites"][0]["fractional"][0] = coordinate
            bad.append(cell)
        cell = unit_cell(); del cell["sites"]; bad.append(cell)
        cell = unit_cell(); cell["sites"] = None; bad.append(cell)
        cell = unit_cell(); cell["sites"][0]["id"] = ""; bad.append(cell)
        cell = unit_cell(); cell["lattice"] = [[1,0],[0,1]]; bad.append(cell)
        with Archive() as archive:
            for cell in bad:
                with self.subTest(cell=cell), self.assertRaises(ValueError):
                    archive.put("p", cell)
            with self.assertRaises(ValueError):
                archive.put("", unit_cell())
            self.assertEqual(archive.connection.execute(
                "SELECT COUNT(*) FROM products").fetchone()[0], 0)
            empty = unit_cell(); empty["sites"] = []; empty["bonds"] = []
            self.assertEqual(archive.put("empty", empty), 1)

    def test_reject_invalid_lattice_and_unit(self):
        with Archive() as archive:
            old = archive.put("p", unit_cell())
            for matrix, unit in (([[1,0,0],[1,0,0],[0,0,1]], "nm"),
                                 (unit_cell()["lattice"], "meter")):
                cell = unit_cell(); cell["lattice"] = matrix; cell["unit"] = unit
                with self.assertRaises(ValueError):
                    archive.put("p", cell)
                self.assertEqual(archive.current_revision("p"), old)
            self.assertEqual(archive.connection.execute(
                "SELECT COUNT(*) FROM revisions").fetchone()[0], 1)

    def test_digest_identity_and_order(self):
        cell = unit_cell()
        accepted, digest = accept_cell(cell)
        cell["sites"].reverse()
        self.assertEqual(accept_cell(cell)[1], digest)
        cell["sites"][0]["id"] = "V"
        cell["bonds"][0]["v"] = "V"
        self.assertNotEqual(accept_cell(cell)[1], digest)
        self.assertEqual(cell_digest(accepted), digest)
