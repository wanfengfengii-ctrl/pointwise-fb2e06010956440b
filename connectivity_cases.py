"""Public inputs and independent certificate checks, not a connectivity solver."""
import argparse
import copy
from fractions import Fraction
import hashlib
import json
from pathlib import Path
import sqlite3
import tempfile
from crystal import Archive, Tasks, Results


def cell_with(edges, names=("u",)):
    return {
        "unit": "nm", "lattice": [[1, 0, 0], [0, 1, 0], [0, 0, 1]],
        "sites": [{"id": name, "fractional": [0, 0, 0]} for name in names],
        "bonds": [{"id": key, "u": u, "v": v, "offset": list(t)}
                  for key, u, v, t in edges],
    }


def query(u, v, t):
    return {"u": u, "v": v, "translation": list(t)}


def fixtures():
    def self_queries(vectors):
        return [query("u", "u", t) for t in vectors]
    two = cell_with([("a", "u", "u", (2, 0, 0))])
    both = cell_with([("a", "u", "u", (2, 0, 0)),
                      ("b", "u", "u", (3, 0, 0)),
                      ("zero", "u", "u", (0, 0, 0))])
    yield "index-two", two, self_queries(
        [(0, 0, 0), (1, 0, 0), (2, 0, 0), (-4, 0, 0), (0, 1, 0)]), [
            True, False, True, True, False]
    yield "coprime", both, self_queries(
        [(1, 0, 0), (-7, 0, 0), (10**30, 0, 0), (0, 0, 1)]), [
            True, True, True, False]
    yield "tree", cell_with([("edge", "u", "v", (1, 0, 0))], ("u", "v")), [
        query("u", "v", (1, 0, 0)), query("u", "v", (0, 0, 0)),
        query("v", "u", (-1, 0, 0)), query("v", "u", (0, 0, 0)),
        query("u", "u", (0, 0, 0)), query("u", "u", (2, 0, 0))], [
            True, False, True, False, True, False]
    yield "parallel", cell_with([
        ("a", "u", "v", (1, 0, 0)),
        ("b", "u", "v", (3, 0, 0))], ("u", "v")), [
        query("u", "v", (2, 0, 0)), query("u", "v", (3, 0, 0)),
        query("u", "v", (-1, 0, 0)), query("v", "u", (-3, 0, 0)),
        query("v", "u", (-2, 0, 0))], [False, True, True, True, False]
    yield "separate-components", cell_with([
        ("a", "u", "u", (2, 0, 0)),
        ("b", "v", "v", (3, 0, 0))], ("u", "v")), [
        query("u", "u", (1, 0, 0)), query("v", "v", (1, 0, 0)),
        query("u", "v", (0, 0, 0)), query("v", "v", (3, 0, 0))], [
            False, False, False, True]
    yield "coupled-rank-two", cell_with([
        ("a", "u", "u", (2, 1, 0)),
        ("b", "u", "u", (0, 3, 0))]), self_queries([
            (2, 0, 0), (2, 1, 0), (0, 3, 0), (4, 2, 0),
            (1, 0, 0), (0, 0, 1)]), [False, True, True, True, False, False]
    yield "coupled-rank-three", cell_with([
        ("a", "u", "u", (2, 1, 0)),
        ("b", "u", "u", (0, 3, 1)),
        ("c", "u", "u", (0, 0, 5))]), self_queries([
            (2, 4, 6), (2, 4, 1), (2, 4, 2), (0, 1, 0)]), [
                True, True, False, False]
    m = 10**30
    yield "large-bezout", cell_with([
        ("a", "u", "u", (2*m, 0, 0)),
        ("b", "u", "u", (2*m+1, 0, 0))]), self_queries([
            (1, 0, 0), (-1, 0, 0), (0, 1, 0)]), [True, True, False]
    yield "isolated-identities", cell_with([], ("u", "v")), [
        query("u", "u", (0, 0, 0)), query("u", "u", (1, 0, 0)),
        query("u", "v", (0, 0, 0))], [True, False, False]
    yield "empty-queries", cell_with([], ()), [], []


def audit_answers(cell, queries, result, expected=None):
    assert set(result) == {"source", "answers"}
    assert set(result["source"]) == {
        "product_id", "revision", "digest", "calculation"}
    assert result["source"]["calculation"] == "connectivity-v1"
    answers = result["answers"]
    assert type(answers) is list and len(answers) == len(queries)
    ids = {s["id"] for s in cell["sites"]}
    edges = {b["id"]: b for b in cell["bonds"]}
    adjacency = {u: set() for u in ids}
    for edge in edges.values():
        adjacency[edge["u"]].add(edge["v"])
        adjacency[edge["v"]].add(edge["u"])
    for number, (q, answer) in enumerate(zip(queries, answers)):
        assert set(answer) == {"query", "reachable", "coefficients"}
        assert answer["query"] == {
            "u": q["u"], "v": q["v"], "translation": list(q["translation"])}
        assert type(answer["reachable"]) is bool
        if expected is not None:
            assert answer["reachable"] is expected[number], (number, answer)
        coefficients = answer["coefficients"]
        if not answer["reachable"]:
            assert coefficients is None
            continue
        component, pending = set(), [q["u"]]
        while pending:
            u = pending.pop()
            if u not in component:
                component.add(u)
                pending.extend(adjacency[u] - component)
        assert q["v"] in component
        assert type(coefficients) is dict
        balance = dict.fromkeys(ids, 0)
        displacement = [0, 0, 0]
        for key, coefficient in coefficients.items():
            assert key in edges and type(coefficient) is int and coefficient != 0
            edge = edges[key]
            assert edge["u"] in component and edge["v"] in component
            balance[edge["u"]] -= coefficient
            balance[edge["v"]] += coefficient
            for axis in range(3):
                displacement[axis] += coefficient * edge["offset"][axis]
        assert balance == {u: int(u == q["v"]) - int(u == q["u"]) for u in ids}
        assert displacement == list(q["translation"])


def audit_packet(text, queries, expected):
    packet = json.loads(text)
    assert set(packet) == {"format", "version", "source", "result"}
    assert packet["format"] == "crystal-connectivity-evidence"
    assert type(packet["version"]) is int and packet["version"] == 1
    source, result = packet["source"], packet["result"]
    assert set(source) == {"product_id", "revision", "digest", "cell", "raw_cell"}
    cell = source["cell"]
    canonical = json.dumps(cell, ensure_ascii=True, sort_keys=True,
                           separators=(",", ":")).encode("ascii")
    assert hashlib.sha256(canonical).hexdigest() == source["digest"]
    for key in ("product_id", "revision", "digest"):
        assert result["source"][key] == source[key]
    raw = source["raw_cell"]
    scale = {"nm": Fraction(1), "angstrom": Fraction(1, 10),
             "pm": Fraction(1, 1000)}[raw["unit"]]
    def q(value):
        return Fraction(*value) if isinstance(value, list) else Fraction(value)
    assert cell["unit"] == "nm"
    assert [[q(x)*scale for x in row] for row in raw["lattice"]] == [
        [q(x) for x in row] for row in cell["lattice"]]
    assert {s["id"]: [q(x) for x in s["fractional"]] for s in raw["sites"]} == {
        s["id"]: [q(x) for x in s["fractional"]] for s in cell["sites"]}
    assert sorted(raw["bonds"], key=lambda b: b["id"]) == cell["bonds"]
    audit_answers(cell, queries, result, expected)
    return packet


def stored_rows(connection):
    tables = sorted(row[0] for row in connection.execute(
        "SELECT name FROM sqlite_master WHERE type='table' "
        "AND name NOT LIKE 'sqlite_%'"))
    return {table: connection.execute(
        'SELECT * FROM "' + table.replace('"', '""') + '"').fetchall()
            for table in tables}


def checks_in_memory():
    with Archive() as archive:
        old_revision = archive.put("legacy", cell_with([]))
        legacy_id = Tasks(archive).run_geometry("legacy", old_revision)
        legacy = Results(archive).read(legacy_id)
        for name, cell, queries, expected in fixtures():
            revision = archive.put(name, cell)
            supplied = copy.deepcopy(queries)
            rid = Tasks(archive).run_connectivity(name, revision, supplied)
            value = Results(archive).read(rid)
            assert value["source"]["product_id"] == name
            assert value["source"]["revision"] == revision
            audit_answers(cell, queries, value, expected)
            packet = audit_packet(Results(archive).export(rid), queries, expected)
            assert packet["result"] == value
            saved = copy.deepcopy(value)
            if supplied:
                supplied[0]["translation"][0] += 101
            supplied.clear()
            value["source"]["revision"] = 999
            if value["answers"] and value["answers"][0]["coefficients"]:
                value["answers"][0]["coefficients"].clear()
            value["answers"].clear()
            assert Results(archive).read(rid) == saved
        name, cell, queries, expected = next(fixtures())
        revision = archive.current_revision(name)
        repeated = [{"u": "u", "v": "u", "translation": (2, 0, 0)} for _ in range(16)]
        repeated_id = Tasks(archive).run_connectivity(name, revision, repeated)
        audit_answers(cell, repeated, Results(archive).read(repeated_id), [True]*16)
        for large in (
                cell_with([], tuple("s"+str(i) for i in range(17))),
                cell_with([("e"+str(i), "u", "u", (0, 0, 0)) for i in range(33)])):
            large_revision = archive.put("too-large", large)
            snapshot = stored_rows(archive.connection)
            try:
                Tasks(archive).run_connectivity("too-large", large_revision, [])
            except ValueError:
                pass
            else:
                raise AssertionError("oversize connectivity input was accepted")
            assert stored_rows(archive.connection) == snapshot
        before = stored_rows(archive.connection)
        bad_inputs = [
            None, (), [None], [{}], queries * 4,
            [query("absent", "u", (0, 0, 0))],
            [{"u": "u", "v": "u", "translation": [True, 0, 0]}],
            [{"u": "u", "v": "u", "translation": [Fraction(1), 0, 0]}],
            [{"u": "u", "v": "u", "translation": [0, 0]}],
            [dict(queries[0], unexpected=1)],
            [queries[0], {"u": "u", "v": "u", "translation": [0, False, 0]}],
        ]
        for invalid in bad_inputs:
            try:
                Tasks(archive).run_connectivity(name, revision, invalid)
            except ValueError:
                pass
            else:
                raise AssertionError(("accepted invalid query", invalid))
            assert stored_rows(archive.connection) == before
        for product, number in (("missing", 1), (name, revision+1)):
            try:
                Tasks(archive).run_connectivity(product, number, [])
            except KeyError:
                pass
            else:
                raise AssertionError("missing source was accepted")
        # Reject publication writes irrespective of internal result tables.
        rejections = []
        for number, table in enumerate(before):
            escaped = table.replace('"', '""')
            for operation in ("INSERT", "UPDATE", "DELETE"):
                trigger = "connectivity_reject_" + str(number) + "_" + operation
                rejections.append(trigger)
                archive.connection.execute(
                    "CREATE TEMP TRIGGER " + trigger + " BEFORE " + operation +
                    ' ON "' + escaped + '" BEGIN '
                    "SELECT RAISE(ABORT,'connectivity publication rejected'); END")
        try:
            try:
                fresh_queries = [query("u", "u", (6, 0, 0)),
                                 query("u", "u", (7, 0, 0))]
                Tasks(archive).run_connectivity(name, revision, fresh_queries)
            except sqlite3.IntegrityError as exc:
                assert "connectivity publication rejected" in str(exc)
            else:
                raise AssertionError("publication rejection did not propagate")
            assert stored_rows(archive.connection) == before
            assert Results(archive).read(legacy_id) == legacy
        finally:
            for trigger in rejections:
                archive.connection.execute("DROP TRIGGER " + trigger)
    print("CONNECTIVITY_MEMORY_CHECKS_OK")


def persistent_demo(path):
    q = [query("u", "u", (1, 0, 0))]
    old_cell = cell_with([("a", "u", "u", (2, 0, 0))])
    with Archive(path) as archive:
        revision = archive.put("specimen", old_cell)
        old_id = Tasks(archive).run_connectivity("specimen", revision, q)
        old = Results(archive).read(old_id)
        audit_answers(old_cell, q, old, [False])
        new_cell = copy.deepcopy(old_cell)
        new_cell["bonds"].append(
            {"id": "b", "u": "u", "v": "u", "offset": [3, 0, 0]})
        newer = archive.put("specimen", new_cell)
        new_id = Tasks(archive).run_connectivity("specimen", newer, q)
        new = Results(archive).read(new_id)
        audit_answers(new_cell, q, new, [True])
    with Archive(path) as archive:
        assert Results(archive).read(old_id) == old
        assert Results(archive).read(new_id) == new
        old_packet = Results(archive).export(old_id)
        new_packet = Results(archive).export(new_id)
    # These checks execute after the database is closed.
    audit_packet(old_packet, q, [False])
    audit_packet(new_packet, q, [True])
    print("POINTWISE_SCENARIO " + json.dumps({
        "case_id": "periodic-connectivity",
        "old_revision": revision, "new_revision": newer,
        "expected": [False, True],
        "actual": [old["answers"][0]["reachable"], new["answers"][0]["reachable"]],
        "integer_chain": new["answers"][0]["coefficients"],
        "offline_source_digest": json.loads(new_packet)["source"]["digest"],
    }, sort_keys=True))


def fixture_self_check():
    # Only validates supplied inputs and the positive-certificate checker.
    cases = list(fixtures())
    with Archive() as archive:
        for name, cell, queries, expected in cases:
            archive.put(name, cell)
            assert len(queries) == len(expected) <= 16
    source = {"product_id": "p", "revision": 1, "digest": "0"*64,
              "calculation": "connectivity-v1"}
    cell = cell_with([("a", "u", "u", (2, 0, 0)),
                      ("b", "u", "u", (3, 0, 0))])
    q = query("u", "u", (1, 0, 0))
    result = {"source": source, "answers": [
        {"query": q, "reachable": True, "coefficients": {"a": -1, "b": 1}}]}
    audit_answers(cell, [q], result, [True])
    for coefficients in ({"a": True}, {"a": 0}, {"a": 1}, {"absent": 1}):
        bad = copy.deepcopy(result)
        bad["answers"][0]["coefficients"] = coefficients
        try:
            audit_answers(cell, [q], bad, [True])
        except AssertionError:
            pass
        else:
            raise AssertionError("certificate checker accepted invalid evidence")
    print("CONNECTIVITY_FIXTURE_SELF_CHECK_OK", len(cases))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--memory", action="store_true")
    parser.add_argument("--self-check", action="store_true")
    args = parser.parse_args()
    if args.self_check:
        fixture_self_check()
        return
    checks_in_memory()
    if not args.memory:
        with tempfile.TemporaryDirectory(prefix="connectivity-check-") as directory:
            persistent_demo(Path(directory) / "archive.sqlite")


if __name__ == "__main__":
    main()
