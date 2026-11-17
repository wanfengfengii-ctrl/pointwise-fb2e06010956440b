import sys
sys.dont_write_bytecode = True
import argparse
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path
import tempfile
from crystal import Archive, Tasks, Results


def offline_values(packet):
    cell = packet["source"]["cell"]
    a = [[F(*x) for x in row] for row in cell["lattice"]]
    f = {s["id"]:[F(*x) for x in s["fractional"]] for s in cell["sites"]}
    determinant = (
        a[0][0]*(a[1][1]*a[2][2]-a[1][2]*a[2][1])
        -a[0][1]*(a[1][0]*a[2][2]-a[1][2]*a[2][0])
        +a[0][2]*(a[1][0]*a[2][1]-a[1][1]*a[2][0]))
    rows = {}
    for bond in cell["bonds"]:
        delta = [f[bond["v"]][j]+bond["offset"][j]-f[bond["u"]][j]
                 for j in range(3)]
        d = [sum(a[i][j]*delta[j] for j in range(3)) for i in range(3)]
        rows[bond["id"]] = {"displacement":d,
                           "length_squared":sum(x*x for x in d)}
    return abs(determinant),rows


def offline_connectivity_audit(packet):
    """Off-database check of a connectivity evidence packet.

    Uses only json and hashlib on the exported text: packet shape, source
    digest recomputation, and the integer-chain certificate (component,
    endpoint balance, offset sum) for every reachable answer.
    """
    assert set(packet) == {"format", "version", "source", "result"}
    assert packet["format"] == "crystal-connectivity-evidence"
    assert packet["version"] == 1
    source, result = packet["source"], packet["result"]
    assert set(source) == {"product_id", "revision", "digest",
                           "cell", "raw_cell"}
    cell = source["cell"]
    canonical = json.dumps(cell, ensure_ascii=True, sort_keys=True,
                           separators=(",", ":")).encode("ascii")
    assert hashlib.sha256(canonical).hexdigest() == source["digest"]
    for key in ("product_id", "revision", "digest"):
        assert result["source"][key] == source[key]
    assert result["source"]["calculation"] == "connectivity-v1"
    ids = {site["id"] for site in cell["sites"]}
    edges = {bond["id"]: bond for bond in cell["bonds"]}
    adjacency = {identity: set() for identity in ids}
    for bond in cell["bonds"]:
        adjacency[bond["u"]].add(bond["v"])
        adjacency[bond["v"]].add(bond["u"])
    for answer in result["answers"]:
        q = answer["query"]
        coefficients = answer["coefficients"]
        if not answer["reachable"]:
            assert coefficients is None
            continue
        component, pending = set(), [q["u"]]
        while pending:
            site = pending.pop()
            if site not in component:
                component.add(site)
                pending.extend(adjacency[site] - component)
        assert q["v"] in component
        balance = dict.fromkeys(ids, 0)
        displacement = [0, 0, 0]
        for key, coefficient in coefficients.items():
            edge = edges[key]
            assert type(coefficient) is int and coefficient != 0
            assert edge["u"] in component and edge["v"] in component
            balance[edge["u"]] -= coefficient
            balance[edge["v"]] += coefficient
            for axis in range(3):
                displacement[axis] += coefficient * edge["offset"][axis]
        assert balance == {
            site: int(site == q["v"]) - int(site == q["u"]) for site in ids}
        assert displacement == q["translation"]
    return packet


def demonstrate_connectivity(path):
    old_cell = {"unit": "nm", "lattice": [[1, 0, 0], [0, 1, 0], [0, 0, 1]],
                "sites": [{"id": "u", "fractional": [0, 0, 0]}],
                "bonds": [{"id": "a", "u": "u", "v": "u",
                           "offset": [2, 0, 0]}]}
    queries = [{"u": "u", "v": "u", "translation": [1, 0, 0]}]
    with Archive(path) as archive:
        old_revision = archive.put("specimen", old_cell)
        old_id = Tasks(archive).run_connectivity(
            "specimen", old_revision, queries)
        old = Results(archive).read(old_id)
        assert old["answers"][0]["reachable"] is False
        assert old["answers"][0]["coefficients"] is None
        new_cell = json.loads(json.dumps(old_cell))
        new_cell["bonds"].append(
            {"id": "b", "u": "u", "v": "u", "offset": [3, 0, 0]})
        new_revision = archive.put("specimen", new_cell)
        new_id = Tasks(archive).run_connectivity(
            "specimen", new_revision, queries)
        new = Results(archive).read(new_id)
        assert new["answers"][0]["reachable"] is True
    with Archive(path) as archive:
        # Reopening changes nothing: both revisions stay readable and fixed.
        assert Results(archive).read(old_id) == old
        assert Results(archive).read(new_id) == new
        old_packet = json.loads(Results(archive).export(old_id))
        new_packet = json.loads(Results(archive).export(new_id))
    # These audits run after the database is closed.
    offline_connectivity_audit(old_packet)
    offline_connectivity_audit(new_packet)
    print("POINTWISE_SCENARIO " + json.dumps({
        "case_id": "periodic-connectivity",
        "old_revision": old_revision, "new_revision": new_revision,
        "expected": [False, True],
        "actual": [old["answers"][0]["reachable"],
                   new["answers"][0]["reachable"]],
        "integer_chain": new["answers"][0]["coefficients"],
        "offline_source_digest": new_packet["source"]["digest"],
        "code_paths": ["crystal/tasks.py", "crystal/connectivity_solver.py",
                       "crystal/smith_normal_form.py",
                       "crystal/connectivity_graph.py",
                       "crystal/result_store.py", "crystal/results.py",
                       "crystal/export.py"],
        "oracle": "json+hashlib certificate audit after archive close"},
        sort_keys=True))
    return 0


def demonstrate(path, product, revision, sample=None):
    with Archive(path) as archive:
        if sample is not None:
            revision = archive.put(product,sample)
        rid = Tasks(archive).run_geometry(product,revision)
        old = Results(archive).read(rid)
        if sample is not None:
            changed = json.loads(json.dumps(sample))
            changed["lattice"][0][0] = 2
            new_revision = archive.put(product,changed)
            newer = Results(archive).read(
                Tasks(archive).run_geometry(product,new_revision))
            assert newer["volume"] == F(4)
            assert old["volume"] == F(2)
    with Archive(path) as archive:
        assert Results(archive).read(rid) == old
        packet = json.loads(Results(archive).export(rid))
    volume, expected = offline_values(packet)
    assert old["volume"] == volume
    for row in old["bonds"]:
        assert row["displacement"] == expected[row["id"]]["displacement"]
        assert row["length_squared"] == expected[row["id"]]["length_squared"]
    if sample is not None:
        assert expected["uv"]["displacement"] == [F(5,4),F(3,2),F(0)]
        assert expected["uv"]["length_squared"] == F(61,16)
    def pair(q):
        return [q.numerator,q.denominator]
    expected_json = {
        key:{"displacement":[pair(x) for x in row["displacement"]],
             "length_squared":pair(row["length_squared"])}
        for key,row in expected.items()}
    actual_json = {
        row["id"]:{"displacement":row["displacement"],
                    "length_squared":row["length_squared"]}
        for row in packet["result"]["bonds"]}
    print("POINTWISE_SCENARIO " + json.dumps({
        "case_id":"N8",
        "input":packet["source"]["cell"],
        "expected":{"volume":pair(volume),"bonds":expected_json},
        "actual":{"volume":packet["result"]["volume"],"bonds":actual_json},
        "downstream_expected":{k:v["length_squared"]
                               for k,v in expected_json.items()},
        "downstream_actual":{k:v["length_squared"]
                             for k,v in actual_json.items()},
        "code_paths":["crystal/archive.py","crystal/tasks.py",
                      "crystal/executor.py","crystal/result_store.py",
                      "crystal/results.py","crystal/export.py"],
        "oracle":"Fraction recomputation from exported source after archive close"}))
    return 0


def main():
    parser = argparse.ArgumentParser(description="Primitive-cell laboratory workflow")
    parser.add_argument("--database",type=Path)
    parser.add_argument("--product")
    parser.add_argument("--revision",type=int)
    args = parser.parse_args()
    if args.database is not None:
        if not args.database.is_file() or not args.product or args.revision is None:
            parser.error("existing --database requires --product and --revision")
        return demonstrate(args.database,args.product,args.revision)
    if args.product is not None or args.revision is not None:
        parser.error("--product and --revision require --database")
    sample = json.loads((Path(__file__).resolve().parent/
                        "examples/day_cell.json").read_text(encoding="utf-8"))
    with tempfile.TemporaryDirectory(prefix="crystal-demo-") as directory:
        status = demonstrate(Path(directory)/"lab.sqlite","day-cell",None,sample)
        if status:
            return status
        return demonstrate_connectivity(Path(directory)/"connectivity.sqlite")


if __name__ == "__main__":
    raise SystemExit(main())
