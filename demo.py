import sys
sys.dont_write_bytecode = True
import argparse
from fractions import Fraction as F
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
        return demonstrate(Path(directory)/"lab.sqlite","day-cell",None,sample)


if __name__ == "__main__":
    raise SystemExit(main())
