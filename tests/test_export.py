import copy
from dataclasses import replace
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path
import tempfile
import unittest
from crystal import Archive, Tasks, Results
from crystal.export import build_export_packet
from crystal.input_model import CellInput
from tests.test_intake import skew_cell


def independent_recalculate(packet):
    def q(pair):
        assert type(pair) is list and len(pair) == 2
        assert all(type(x) is int for x in pair) and pair[1] > 0
        value = F(*pair)
        assert [value.numerator,value.denominator] == pair
        return value

    assert packet["format"] == "crystal-primitive-evidence"
    assert packet["version"] == 1
    source, result = packet["source"], packet["result"]
    cell = source["cell"]
    assert cell["unit"] == "nm"
    canonical = json.dumps(cell,sort_keys=True,separators=(",",":"),
                           ensure_ascii=True).encode("ascii")
    assert hashlib.sha256(canonical).hexdigest() == source["digest"]
    for key in ("product_id","revision","digest"):
        assert result["source"][key] == source[key]
    assert result["source"]["calculation"] == "primitive-v1"
    assert result["units"] == {
        "length":"nm","length_squared":"nm2","volume":"nm3"}
    a = [[q(x) for x in row] for row in cell["lattice"]]
    determinant = (
        a[0][0]*a[1][1]*a[2][2]+a[0][1]*a[1][2]*a[2][0]
        +a[0][2]*a[1][0]*a[2][1]-a[0][2]*a[1][1]*a[2][0]
        -a[0][1]*a[1][0]*a[2][2]-a[0][0]*a[1][2]*a[2][1])
    assert q(result["volume"]) == abs(determinant)
    sites = {s["id"]:[q(x) for x in s["fractional"]] for s in cell["sites"]}
    positions = {s["id"]:[q(x) for x in s["cartesian"]]
                 for s in result["positions"]}
    assert len(sites) == len(cell["sites"]) == len(positions)
    assert len(positions) == len(result["positions"])
    assert set(sites) == set(positions)
    for identity, f in sites.items():
        assert positions[identity] == [
            sum(a[i][j]*f[j] for j in range(3)) for i in range(3)]
    rows = {b["id"]:b for b in result["bonds"]}
    assert len(rows) == len(result["bonds"]) == len(cell["bonds"])
    assert set(rows) == {b["id"] for b in cell["bonds"]}
    for bond in cell["bonds"]:
        row = rows[bond["id"]]
        assert (row["u"],row["v"],row["offset"]) == (
            bond["u"],bond["v"],bond["offset"])
        assert len(bond["offset"]) == 3
        assert all(type(t) is int for t in bond["offset"])
        delta = [sites[bond["v"]][j]+bond["offset"][j]-sites[bond["u"]][j]
                 for j in range(3)]
        displacement = [sum(a[i][j]*delta[j] for j in range(3))
                        for i in range(3)]
        assert [q(x) for x in row["displacement"]] == displacement
        assert q(row["length_squared"]) == sum(x*x for x in displacement)
    raw = source["raw_cell"]
    scale = {"nm":F(1),"angstrom":F(1,10),"pm":F(1,1000)}[raw["unit"]]
    def raw_q(value):
        return F(*value) if isinstance(value,list) else F(value)
    assert [[raw_q(x)*scale for x in row] for row in raw["lattice"]] == a
    assert {s["id"]:[raw_q(x) for x in s["fractional"]]
            for s in raw["sites"]} == sites
    assert sorted(raw["bonds"],key=lambda b:b["id"]) == cell["bonds"]
    return {identity:F(*row["length_squared"]) for identity,row in rows.items()}


class TestExport(unittest.TestCase):
    def test_independent_recalculation(self):
        with tempfile.TemporaryDirectory(prefix="crystal-test-") as directory:
            path = Path(directory)/"archive.sqlite"
            with Archive(path) as archive:
                cell = skew_cell()
                cell["bonds"] += [
                    {"id":"duplicate","u":"u","v":"v","offset":[0,1,0]},
                    {"id":"self","u":"u","v":"u","offset":[1,0,0]}]
                r1 = archive.put("p",cell)
                rid = Tasks(archive).run_geometry("p",r1)
                original = Results(archive).read(rid)
                cell["lattice"][0][0] = 2
                r2 = archive.put("p",cell)
                rid2 = Tasks(archive).run_geometry("p",r2)
                expected_new = Results(archive).read(rid2)
            for _ in range(2):
                with Archive(path) as archive:
                    self.assertEqual(archive.current_revision("p"),r2)
                    self.assertEqual(Results(archive).read(rid),original)
                    self.assertEqual(Results(archive).read(rid2),expected_new)
                    packet = Results(archive).export(rid)
            values = independent_recalculate(json.loads(packet))
            self.assertEqual(values,{"duplicate":F(61,16),
                                     "self":F(1),"uv":F(61,16)})

    def test_exact_json_large_numbers(self):
        with tempfile.TemporaryDirectory(prefix="crystal-test-") as directory:
            path = Path(directory)/"archive.sqlite"
            with Archive(path) as archive:
                r1 = archive.put("p",skew_cell())
                old = Tasks(archive).run_geometry("p",r1)
                cell = skew_cell()
                cell["unit"] = "angstrom"
                cell["lattice"] = [[F(9007199254740993,7),F(1,3),0],
                                   [0,F(10,7),0],[0,0,20]]
                cell["sites"][1]["fractional"] = [F(1,3),F(2,5),F(1,11)]
                r2 = archive.put("p",cell)
                new = Tasks(archive).run_geometry("p",r2)
                new_value = Results(archive).read(new)
            with Archive(path) as archive:
                self.assertEqual(Results(archive).read(new),new_value)
                packets = [json.loads(Results(archive).export(rid))
                           for rid in (old,new)]
            for packet in packets:
                independent_recalculate(packet)
            self.assertEqual(F(*packets[1]["source"]["cell"]["lattice"][0][0]),
                             F(9007199254740993,70))
            self.assertEqual(packets[0]["source"]["revision"],r1)
            self.assertEqual(packets[1]["source"]["revision"],r2)

    def test_evidence_rejects_incomplete_or_wrong_source(self):
        with Archive() as archive:
            r = archive.put("p",skew_cell())
            rid = Tasks(archive).run_geometry("p",r)
            fixed = archive.revisions.read_revision("p",r)
            result = Results(archive).store.read_bundle(rid)
            raw = fixed.cell.raw_mapping()
            del raw["bonds"][0]["offset"]
            bad_cell = replace(fixed.cell,raw_json=json.dumps(raw))
            with self.assertRaises(ValueError):
                build_export_packet(replace(fixed,cell=bad_cell),result)
            with self.assertRaises(ValueError):
                build_export_packet(replace(fixed,revision=r+1),result)
            raw = fixed.cell.raw_mapping()
            raw["lattice"][0][0] = 2
            with self.assertRaises(ValueError):
                build_export_packet(
                    replace(fixed,cell=replace(fixed.cell,raw_json=json.dumps(raw))),
                    result)
            packet = json.loads(Results(archive).export(rid))
            bad = copy.deepcopy(packet)
            del bad["source"]["cell"]["bonds"][0]["offset"]
            with self.assertRaises((KeyError,AssertionError)):
                independent_recalculate(bad)
            self.assertEqual(independent_recalculate(packet),{"uv":F(61,16)})
