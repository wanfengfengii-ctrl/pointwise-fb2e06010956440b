import copy
from fractions import Fraction as F
import unittest
from crystal import Archive, Tasks, Results
from tests.test_intake import unit_cell


class TestBonds(unittest.TestCase):
    def test_reject_invalid_reference(self):
        with Archive() as archive:
            old = archive.put("p",unit_cell())
            variants = []
            for endpoint in ("u","v"):
                cell = unit_cell(); cell["bonds"][0][endpoint] = "absent"
                variants.append(cell)
            for offset in ([1,2], [1,2,3,4], [0,True,0],
                           [0,[1,2],0], [0,F(1),0], [0,0.0,0], None):
                cell = unit_cell(); cell["bonds"][0]["offset"] = offset
                variants.append(cell)
            cell = unit_cell(); cell["bonds"].append(copy.deepcopy(cell["bonds"][0]))
            variants.append(cell)
            for cell in variants:
                with self.subTest(cell=cell), self.assertRaises(ValueError):
                    archive.put("p",cell)
                self.assertEqual(archive.current_revision("p"), old)
                self.assertEqual(archive.read_revision("p",old)["bonds"][0]["id"],"uv")
            with self.assertRaises(KeyError):
                archive.read_revision("p",old+1)

    def test_identity_and_empty_graph(self):
        cell = unit_cell()
        cell["sites"][1]["fractional"] = [0,0,0]
        cell["bonds"] = [
            {"id":"a","u":"u","v":"v","offset":[0,0,0]},
            {"id":"b","u":"u","v":"v","offset":[0,0,0]},
            {"id":"self","u":"u","v":"u","offset":[0,1,0]},
            {"id":"back","u":"u","v":"u","offset":[0,-1,0]}]
        with Archive() as archive:
            r = archive.put("p",cell)
            value = Results(archive).read(Tasks(archive).run_geometry("p",r))
            self.assertEqual([p["id"] for p in value["positions"]], ["u","v"])
            rows = {b["id"]:b for b in value["bonds"]}
            self.assertEqual(set(rows), {"a","b","self","back"})
            self.assertEqual(rows["a"]["length_squared"], F(0))
            self.assertEqual(rows["b"]["length_squared"], F(0))
            self.assertEqual(rows["self"]["displacement"], [F(0),F(1),F(0)])
            self.assertEqual(rows["back"]["displacement"], [F(0),F(-1),F(0)])
            cell["bonds"] = []
            r = archive.put("p",cell)
            value = Results(archive).read(Tasks(archive).run_geometry("p",r))
            self.assertEqual(value["bonds"], [])
            self.assertEqual(len(value["positions"]), 2)
            cell["sites"] = []
            r = archive.put("empty",cell)
            value = Results(archive).read(Tasks(archive).run_geometry("empty",r))
            self.assertEqual(value["positions"], [])
            self.assertEqual(value["volume"], F(1))
