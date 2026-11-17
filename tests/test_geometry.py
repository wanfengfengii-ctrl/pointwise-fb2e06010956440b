from dataclasses import replace
from fractions import Fraction as F
import json
import random
import unittest
from crystal import Archive, Tasks, Results
from crystal.executor import execute_geometry, validate_geometry_result
from tests.test_intake import skew_cell, unit_cell


class TestGeometry(unittest.TestCase):
    def test_skew_cell(self):
        with Archive() as archive:
            revision = archive.put("skew", skew_cell())
            result_id = Tasks(archive).run_geometry("skew", revision)
            result = Results(archive).read(result_id)
            positions = {p["id"]: p["cartesian"] for p in result["positions"]}
            self.assertEqual(result["volume"], F(2))
            self.assertEqual(positions["u"], [F(0)]*3)
            self.assertEqual(positions["v"], [F(3,4), F(1,2), F(0)])
            bond = result["bonds"][0]
            self.assertEqual(bond["displacement"], [F(5,4),F(3,2),F(0)])
            self.assertEqual(bond["length_squared"], F(61,16))
            packet = json.loads(Results(archive).export(result_id))
            self.assertEqual(packet["result"]["bonds"][0]
                             ["length_squared"], [61,16])
            print("POINTWISE_SCENARIO " + json.dumps({
                "case_id": "N2",
                "input": skew_cell(),
                "expected": {"volume":[2,1], "v":[[3,4],[1,2],[0,1]],
                             "d":[[5,4],[3,2],[0,1]], "l2":[61,16]},
                "actual": {"volume":packet["result"]["volume"],
                           "v":packet["result"]["positions"][1]["cartesian"],
                           "d":packet["result"]["bonds"][0]["displacement"],
                           "l2":packet["result"]["bonds"][0]["length_squared"]},
                "downstream_expected": [61,16],
                "downstream_actual": packet["result"]["bonds"][0]["length_squared"],
                "code_paths": ["crystal/tasks.py","crystal/cell_geometry.py",
                               "crystal/bond_geometry.py","crystal/export.py"],
                "oracle": "independent hand arithmetic: 25/16+9/4=61/16"}))

    def test_directed_offset_and_left_handed(self):
        cell = skew_cell()
        cell["bonds"] += [
            {"id":"reverse","u":"v","v":"u","offset":[0,-1,0]},
            {"id":"negative","u":"u","v":"v","offset":[0,-1,0]}]
        with Archive() as archive:
            r = archive.put("p", cell)
            value = Results(archive).read(Tasks(archive).run_geometry("p",r))
            rows = {b["id"]:b for b in value["bonds"]}
            self.assertEqual(rows["reverse"]["displacement"],
                             [F(-5,4), F(-3,2), F(0)])
            self.assertEqual(rows["negative"]["displacement"],
                             [F(1,4), F(-1,2), F(0)])
            self.assertEqual(rows["negative"]["length_squared"], F(5,16))
            cell["lattice"][2][2] = -2
            r = archive.put("p",cell)
            value = Results(archive).read(Tasks(archive).run_geometry("p",r))
            self.assertEqual(value["volume"], F(2))

    def test_seeded_cells_independent_arithmetic(self):
        rng = random.Random(1729)
        with Archive() as archive:
            for number in range(24):
                a,b,c = [rng.randint(1,5) for _ in range(3)]
                shear = rng.randint(-3,3)
                cell = unit_cell()
                cell["lattice"] = [[a,shear,0],[0,b,0],[0,0,c]]
                u = [F(rng.randrange(6),7) for _ in range(3)]
                v = [F(rng.randrange(6),7) for _ in range(3)]
                t = [rng.randint(-2,2) for _ in range(3)]
                cell["sites"][0]["fractional"] = u
                cell["sites"][1]["fractional"] = v
                cell["bonds"][0]["offset"] = t
                r = archive.put(str(number),cell)
                value = Results(archive).read(
                    Tasks(archive).run_geometry(str(number),r))
                dx,dy,dz = [v[i]+t[i]-u[i] for i in range(3)]
                expected = [a*dx+shear*dy,b*dy,c*dz]
                self.assertEqual(value["volume"], F(a*b*c))
                self.assertEqual(value["bonds"][0]["displacement"], expected)
                self.assertEqual(value["bonds"][0]["length_squared"],
                                 sum(x*x for x in expected))

    def test_audit_rejects_wrong_coverage_and_values(self):
        with Archive() as archive:
            r = archive.put("p",skew_cell())
            fixed = archive.revisions.read_revision("p",r)
            value = execute_geometry(fixed.cell,("p",r,fixed.digest,"primitive-v1"))
            bad = [
                replace(value, positions=value.positions[:1]),
                replace(value, positions=value.positions+value.positions[:1]),
                replace(value, bonds=()),
                replace(value, bonds=value.bonds+value.bonds),
                replace(value, bonds=(replace(value.bonds[0], offset=(1,0,0)),)),
                replace(value, bonds=(replace(value.bonds[0],
                                              length_squared=F(1)),)),
                replace(value, volume=F(3))]
            for result in bad:
                with self.assertRaises(ValueError):
                    validate_geometry_result(fixed.cell,result)
