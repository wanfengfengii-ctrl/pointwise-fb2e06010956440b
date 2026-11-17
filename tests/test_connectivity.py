"""Independent regression for the periodic-image connectivity chain.

Hand-computed expectations and certificate audits are written directly in
this file; the exported packet audit uses only json and hashlib, never the
product solver, mirroring the V1 export recalculation test.
"""
import copy
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path
import sqlite3
import tempfile
import unittest
from crystal import Archive, Tasks, Results
from crystal.connectivity_result import (
    ConnectivityAnswer, ConnectivityResult,
    encode_connectivity_result, decode_connectivity_result,
    connectivity_result_to_mapping)
from crystal.executor import execute_connectivity
from crystal.smith_normal_form import smith_normal_form, solve_integer_system
from tests.test_intake import skew_cell

DIGEST_ZERO = "0" * 64


def cell_with(edges, names=("u",)):
    return {
        "unit": "nm", "lattice": [[1, 0, 0], [0, 1, 0], [0, 0, 1]],
        "sites": [{"id": name, "fractional": [0, 0, 0]} for name in names],
        "bonds": [{"id": key, "u": u, "v": v, "offset": list(t)}
                  for key, u, v, t in edges],
    }


def query(u, v, translation):
    return {"u": u, "v": v, "translation": list(translation)}


def audit_certificate(cell, answers, queries):
    """Independent positive-certificate audit written for this suite."""
    ids = {site["id"] for site in cell["sites"]}
    edges = {bond["id"]: bond for bond in cell["bonds"]}
    adjacency = {identity: set() for identity in ids}
    for bond in cell["bonds"]:
        adjacency[bond["u"]].add(bond["v"])
        adjacency[bond["v"]].add(bond["u"])
    assert len(answers) == len(queries)
    for answer, q in zip(answers, queries):
        assert set(answer) == {"query", "reachable", "coefficients"}
        assert answer["query"] == {
            "u": q["u"], "v": q["v"], "translation": list(q["translation"])}
        assert type(answer["reachable"]) is bool
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
        assert displacement == list(q["translation"])


def audit_packet_offline(text):
    """Off-database packet audit using only json and hashlib."""
    packet = json.loads(text)
    assert set(packet) == {"format", "version", "source", "result"}
    assert packet["format"] == "crystal-connectivity-evidence"
    assert packet["version"] == 1
    source, result = packet["source"], packet["result"]
    assert set(source) == {"product_id", "revision", "digest",
                           "cell", "raw_cell"}
    canonical = json.dumps(source["cell"], ensure_ascii=True, sort_keys=True,
                           separators=(",", ":")).encode("ascii")
    assert hashlib.sha256(canonical).hexdigest() == source["digest"]
    assert set(result) == {"source", "answers"}
    assert result["source"]["calculation"] == "connectivity-v1"
    for key in ("product_id", "revision", "digest"):
        assert result["source"][key] == source[key]
    return packet


def table_rows(connection):
    tables = sorted(row[0] for row in connection.execute(
        "SELECT name FROM sqlite_master WHERE type='table' "
        "AND name NOT LIKE 'sqlite_%'"))
    return {table: connection.execute(
        'SELECT * FROM "' + table + '"').fetchall() for table in tables}


class TestHandComputedChains(unittest.TestCase):
    def test_single_bond_index_two(self):
        cell = cell_with([("a", "u", "u", (2, 0, 0))])
        queries = [query("u", "u", t) for t in
                   [(0, 0, 0), (1, 0, 0), (2, 0, 0), (-2, 0, 0), (0, 0, 5)]]
        with Archive() as archive:
            revision = archive.put("p", cell)
            rid = Tasks(archive).run_connectivity("p", revision, queries)
            answers = Results(archive).read(rid)["answers"]
        self.assertEqual([a["reachable"] for a in answers],
                         [True, False, True, True, False])
        self.assertEqual(answers[0]["coefficients"], {})
        self.assertEqual(answers[2]["coefficients"], {"a": 1})
        self.assertEqual(answers[3]["coefficients"], {"a": -1})
        audit_certificate(cell, answers, queries)

    def test_coprime_offsets_extend_reachability(self):
        cell = cell_with([("a", "u", "u", (2, 0, 0)),
                          ("b", "u", "u", (3, 0, 0))])
        queries = [query("u", "u", t) for t in
                   [(1, 0, 0), (-1, 0, 0), (5, 0, 0), (0, 1, 0)]]
        with Archive() as archive:
            revision = archive.put("p", cell)
            rid = Tasks(archive).run_connectivity("p", revision, queries)
            answers = Results(archive).read(rid)["answers"]
        self.assertEqual([a["reachable"] for a in answers],
                         [True, True, True, False])
        audit_certificate(cell, answers, queries)

    def test_parallel_bonds_reach_odd_translations_only(self):
        cell = cell_with([("a", "u", "v", (1, 0, 0)),
                          ("b", "u", "v", (3, 0, 0))], ("u", "v"))
        queries = [query("u", "v", (2, 0, 0)), query("u", "v", (3, 0, 0)),
                   query("u", "v", (-1, 0, 0)), query("v", "u", (-3, 0, 0))]
        with Archive() as archive:
            revision = archive.put("p", cell)
            rid = Tasks(archive).run_connectivity("p", revision, queries)
            answers = Results(archive).read(rid)["answers"]
        self.assertEqual([a["reachable"] for a in answers],
                         [False, True, True, True])
        audit_certificate(cell, answers, queries)

    def test_components_and_isolated_sites(self):
        cell = cell_with([("a", "u", "u", (2, 0, 0)),
                          ("b", "v", "v", (3, 0, 0))], ("u", "v", "w"))
        queries = [query("u", "v", (0, 0, 0)), query("u", "u", (2, 0, 0)),
                   query("v", "v", (1, 0, 0)), query("w", "w", (0, 0, 0)),
                   query("w", "w", (2, 0, 0))]
        with Archive() as archive:
            revision = archive.put("p", cell)
            rid = Tasks(archive).run_connectivity("p", revision, queries)
            answers = Results(archive).read(rid)["answers"]
        self.assertEqual([a["reachable"] for a in answers],
                         [False, True, False, True, False])
        audit_certificate(cell, answers, queries)

    def test_coupled_rank_and_large_bezout(self):
        cell = cell_with([("a", "u", "u", (2, 1, 0)),
                          ("b", "u", "u", (0, 3, 1)),
                          ("c", "u", "u", (0, 0, 5))])
        queries = [query("u", "u", (2, 4, 6)), query("u", "u", (2, 4, 2)),
                   query("u", "u", (0, 1, 0))]
        big = 10 ** 30
        large_cell = cell_with([("a", "u", "u", (2 * big, 0, 0)),
                                ("b", "u", "u", (2 * big + 1, 0, 0))])
        large_queries = [query("u", "u", (1, 0, 0)),
                         query("u", "u", (-1, 0, 0)),
                         query("u", "u", (0, 1, 0))]
        with Archive() as archive:
            revision = archive.put("p", cell)
            answers = Results(archive).read(
                Tasks(archive).run_connectivity("p", revision, queries)
            )["answers"]
            self.assertEqual([a["reachable"] for a in answers],
                             [True, False, False])
            audit_certificate(cell, answers, queries)
            large_revision = archive.put("q", large_cell)
            large_answers = Results(archive).read(
                Tasks(archive).run_connectivity(
                    "q", large_revision, large_queries))["answers"]
            self.assertEqual([a["reachable"] for a in large_answers],
                             [True, True, False])
            audit_certificate(large_cell, large_answers, large_queries)

    def test_tree_path_and_chord_cycle_interact(self):
        cell = cell_with([("e1", "u", "v", (1, 0, 0)),
                          ("e2", "v", "w", (0, 1, 0)),
                          ("c", "u", "w", (2, 2, 0))], ("u", "v", "w"))
        queries = [query("u", "w", (3, 3, 0)), query("u", "w", (1, 1, 0)),
                   query("u", "w", (2, 2, 0)), query("u", "w", (2, 3, 0)),
                   query("w", "u", (-3, -3, 0))]
        with Archive() as archive:
            revision = archive.put("p", cell)
            rid = Tasks(archive).run_connectivity("p", revision, queries)
            answers = Results(archive).read(rid)["answers"]
        self.assertEqual([a["reachable"] for a in answers],
                         [True, True, True, False, True])
        audit_certificate(cell, answers, queries)

    def test_zero_offset_and_duplicate_queries(self):
        cell = cell_with([("z", "u", "u", (0, 0, 0)),
                          ("a", "u", "u", (2, 0, 0))])
        queries = [query("u", "u", (2, 0, 0))] * 4
        with Archive() as archive:
            revision = archive.put("p", cell)
            rid = Tasks(archive).run_connectivity("p", revision, queries)
            answers = Results(archive).read(rid)["answers"]
        self.assertEqual(len(answers), 4)
        self.assertTrue(all(a["reachable"] for a in answers))
        audit_certificate(cell, answers, queries)

    def test_empty_cell_and_empty_queries(self):
        with Archive() as archive:
            revision = archive.put("empty", cell_with([], ()))
            rid = Tasks(archive).run_connectivity("empty", revision, [])
            value = Results(archive).read(rid)
            self.assertEqual(value["answers"], [])
            self.assertEqual(set(value), {"source", "answers"})
            self.assertEqual(value["source"]["calculation"], "connectivity-v1")
            packet = audit_packet_offline(Results(archive).export(rid))
            self.assertEqual(packet["result"], value)


class TestInputBoundaries(unittest.TestCase):
    def test_size_limits_only_restrict_connectivity(self):
        with Archive() as archive:
            many_sites = cell_with(
                [], tuple("s" + str(i) for i in range(17)))
            many_bonds = cell_with(
                [("e" + str(i), "u", "u", (0, 0, 0)) for i in range(33)])
            for cell in (many_sites, many_bonds):
                revision = archive.put("p", cell)
                before = table_rows(archive.connection)
                with self.assertRaises(ValueError):
                    Tasks(archive).run_connectivity("p", revision, [])
                self.assertEqual(table_rows(archive.connection), before)
                # The V1 geometry chain is not limited by these bounds.
                self.assertIsInstance(
                    Tasks(archive).run_geometry("p", revision), int)
            exact_sites = cell_with(
                [], tuple("s" + str(i) for i in range(16)))
            exact_bonds = cell_with(
                [("e" + str(i), "u", "u", (0, 0, 0)) for i in range(32)])
            for cell in (exact_sites, exact_bonds):
                revision = archive.put("ok", cell)
                self.assertIsInstance(
                    Tasks(archive).run_connectivity("ok", revision, []), int)

    def test_invalid_queries_reject_whole_batch(self):
        cell = cell_with([("a", "u", "u", (2, 0, 0))])
        with Archive() as archive:
            revision = archive.put("p", cell)
            before = table_rows(archive.connection)
            bad_inputs = [
                None, (), "queries", [None], [{}],
                [query("u", "u", (0, 0, 0))] * 17,
                [query("absent", "u", (0, 0, 0))],
                [query("u", "absent", (0, 0, 0))],
                [{"u": "u", "v": "u", "translation": [True, 0, 0]}],
                [{"u": "u", "v": "u", "translation": [F(1), 0, 0]}],
                [{"u": "u", "v": "u", "translation": [0.0, 0, 0]}],
                [{"u": "u", "v": "u", "translation": [[2, 1], 0, 0]}],
                [{"u": "u", "v": "u", "translation": [0, 0]}],
                [{"u": "u", "v": "u", "translation": None}],
                [dict(query("u", "u", (0, 0, 0)), extra=1)],
                [{"u": "u", "v": "u"}],
                [query("u", "u", (2, 0, 0)),
                 {"u": "u", "v": "u", "translation": [0, False, 0]}],
            ]
            for invalid in bad_inputs:
                with self.subTest(invalid=invalid), \
                        self.assertRaises(ValueError):
                    Tasks(archive).run_connectivity("p", revision, invalid)
                self.assertEqual(table_rows(archive.connection), before)
            with self.assertRaises(KeyError):
                Tasks(archive).run_connectivity("missing", 1, [])
            with self.assertRaises(KeyError):
                Tasks(archive).run_connectivity("p", revision + 1, [])
            self.assertEqual(table_rows(archive.connection), before)

    def test_query_tuple_translation_and_big_ints(self):
        cell = cell_with([("a", "u", "u", (2, 0, 0))])
        with Archive() as archive:
            revision = archive.put("p", cell)
            rid = Tasks(archive).run_connectivity(
                "p", revision,
                [{"u": "u", "v": "u", "translation": (2, 0, 0)},
                 {"u": "u", "v": "u",
                  "translation": [2 * 10 ** 40, 0, 0]}])
            answers = Results(archive).read(rid)["answers"]
            self.assertEqual(answers[0]["query"]["translation"], [2, 0, 0])
            self.assertEqual(answers[1]["coefficients"], {"a": 10 ** 40})


class TestFrozenResults(unittest.TestCase):
    def test_correction_mutation_and_reopen_do_not_drift(self):
        old_cell = cell_with([("a", "u", "u", (2, 0, 0))])
        queries = [query("u", "u", (1, 0, 0))]
        with tempfile.TemporaryDirectory(prefix="crystal-test-") as directory:
            path = Path(directory) / "archive.sqlite"
            with Archive(path) as archive:
                old_revision = archive.put("p", old_cell)
                old_id = Tasks(archive).run_connectivity(
                    "p", old_revision, queries)
                old = Results(archive).read(old_id)
                self.assertFalse(old["answers"][0]["reachable"])
                new_cell = copy.deepcopy(old_cell)
                new_cell["bonds"].append(
                    {"id": "b", "u": "u", "v": "u", "offset": [3, 0, 0]})
                new_revision = archive.put("p", new_cell)
                new_id = Tasks(archive).run_connectivity(
                    "p", new_revision, queries)
                new = Results(archive).read(new_id)
                self.assertTrue(new["answers"][0]["reachable"])
                # Mutating returned containers and the original queries
                # changes nothing stored.
                queries[0]["translation"][0] = 99
                new["answers"][0]["coefficients"].clear()
                new["answers"].clear()
                new["source"]["revision"] = 999
                self.assertEqual(Results(archive).read(new_id)["answers"],
                                 Results(archive).read(new_id)["answers"])
                self.assertEqual(Results(archive).read(old_id), old)
            with Archive(path) as archive:
                self.assertEqual(Results(archive).read(old_id), old)
                self.assertEqual(Results(archive).read(new_id)["source"]
                                 ["revision"], new_revision)
                old_packet = Results(archive).export(old_id)
                new_packet = Results(archive).export(new_id)
            # Off-database audit after close.
            for text, reachable in ((old_packet, False), (new_packet, True)):
                packet = audit_packet_offline(text)
                self.assertEqual(
                    packet["result"]["answers"][0]["reachable"], reachable)
            chain = json.loads(new_packet)["result"]["answers"][0]
            audit_certificate(new_cell, [chain], queries[:1] and
                              [{"u": "u", "v": "u", "translation": [1, 0, 0]}])
            self.assertEqual(
                archive_closed_digest(new_packet),
                hashlib.sha256(json.dumps(
                    json.loads(new_packet)["source"]["cell"],
                    ensure_ascii=True, sort_keys=True,
                    separators=(",", ":")).encode("ascii")).hexdigest())

    def test_result_ids_share_one_sequence_with_geometry(self):
        with Archive() as archive:
            revision = archive.put("p", skew_cell())
            first = Tasks(archive).run_geometry("p", revision)
            second = Tasks(archive).run_connectivity(
                "p", revision, [query("u", "v", (0, 1, 0))])
            third = Tasks(archive).run_geometry("p", revision)
            self.assertEqual((first, second, third),
                             tuple(sorted((first, second, third))))
            self.assertEqual(len({first, second, third}), 3)
            self.assertEqual(Results(archive).read(second)["source"]
                             ["calculation"], "connectivity-v1")
            self.assertEqual(Results(archive).read(first)["source"]
                             ["calculation"], "primitive-v1")
            self.assertEqual(Results(archive).read(third)["volume"], F(2))
            with self.assertRaises(KeyError):
                Results(archive).read(third + 1)

    def test_publication_failure_leaves_no_trace(self):
        with Archive() as archive:
            revision = archive.put("p", cell_with([("a", "u", "u", (2, 0, 0))]))
            old_id = Tasks(archive).run_connectivity(
                "p", revision, [query("u", "u", (2, 0, 0))])
            old = Results(archive).read(old_id)
            before = table_rows(archive.connection)
            archive.connection.execute(
                "CREATE TEMP TRIGGER reject_connectivity BEFORE INSERT "
                "ON tasks BEGIN SELECT RAISE(ABORT,'connectivity write "
                "rejected'); END")
            with self.assertRaisesRegex(sqlite3.IntegrityError,
                                        "connectivity write rejected"):
                Tasks(archive).run_connectivity(
                    "p", revision, [query("u", "u", (4, 0, 0))])
            self.assertEqual(table_rows(archive.connection), before)
            self.assertEqual(Results(archive).read(old_id), old)
            archive.connection.execute("DROP TRIGGER reject_connectivity")
            new_id = Tasks(archive).run_connectivity(
                "p", revision, [query("u", "u", (4, 0, 0))])
            self.assertNotEqual(new_id, old_id)
            self.assertEqual(archive.connection.execute(
                "PRAGMA foreign_key_check").fetchall(), [])

    def test_missing_result_and_body(self):
        with Archive() as archive:
            revision = archive.put("p", cell_with([]))
            rid = Tasks(archive).run_connectivity("p", revision, [])
            with self.assertRaises(KeyError):
                Results(archive).read(rid + 1)
            with self.assertRaises(KeyError):
                Results(archive).export(rid + 1)
            for argument in (True, 0, 1.5):
                with self.assertRaises(ValueError):
                    Results(archive).read(argument)
            archive.connection.execute(
                "DELETE FROM result_bodies WHERE result_id=?", (rid,))
            with self.assertRaises(ValueError):
                Results(archive).read(rid)


def archive_closed_digest(packet_text):
    return json.loads(packet_text)["source"]["digest"]


V1_TASKS_DDL = """CREATE TABLE tasks (
        task_id INTEGER PRIMARY KEY AUTOINCREMENT,
        product_id TEXT NOT NULL,
        revision INTEGER NOT NULL,
        digest TEXT NOT NULL,
        calculation TEXT NOT NULL CHECK(calculation = 'primitive-v1'),
        FOREIGN KEY(product_id,revision)
            REFERENCES revisions(product_id,revision))"""


class TestV1DatabaseMigration(unittest.TestCase):
    def build_v1_database(self, path):
        """Create a genuine V1-format database file with V1 results."""
        from crystal.intake import accept_cell
        from crystal.revision_codec import (
            encode_revision_cell, encode_raw_cell)
        from crystal.executor import execute_geometry
        from crystal.result_codec import encode_result
        connection = sqlite3.connect(path)
        try:
            connection.execute("PRAGMA foreign_keys=ON")
            connection.execute(
                """CREATE TABLE products (
                    product_id TEXT PRIMARY KEY,
                    current_revision INTEGER,
                    FOREIGN KEY(product_id,current_revision)
                        REFERENCES revisions(product_id,revision)
                        DEFERRABLE INITIALLY DEFERRED)""")
            connection.execute(
                """CREATE TABLE revisions (
                    product_id TEXT NOT NULL REFERENCES products(product_id),
                    revision INTEGER NOT NULL CHECK(revision > 0),
                    digest TEXT NOT NULL,
                    canonical_json TEXT NOT NULL,
                    raw_json TEXT NOT NULL,
                    PRIMARY KEY(product_id,revision))""")
            connection.execute(V1_TASKS_DDL)
            connection.execute(
                """CREATE TABLE results (
                    result_id INTEGER PRIMARY KEY AUTOINCREMENT,
                    task_id INTEGER NOT NULL UNIQUE
                        REFERENCES tasks(task_id))""")
            connection.execute(
                """CREATE TABLE result_bodies (
                    result_id INTEGER PRIMARY KEY
                        REFERENCES results(result_id),
                    payload TEXT NOT NULL)""")
            cell, digest = accept_cell(skew_cell())
            connection.execute(
                "INSERT INTO products(product_id) VALUES ('p')")
            connection.execute(
                "INSERT INTO revisions VALUES ('p',1,?,?,?)",
                (digest, encode_revision_cell(cell),
                 encode_raw_cell(cell.raw_mapping())))
            connection.execute(
                "UPDATE products SET current_revision=1 "
                "WHERE product_id='p'")
            result = execute_geometry(
                cell, ("p", 1, digest, "primitive-v1"))
            connection.execute(
                "INSERT INTO tasks(product_id,revision,digest,calculation) "
                "VALUES ('p',1,?,'primitive-v1')", (digest,))
            connection.execute("INSERT INTO results(task_id) VALUES (1)")
            connection.execute(
                "INSERT INTO result_bodies VALUES (1,?)",
                (encode_result(result),))
            connection.execute("PRAGMA user_version=1")
            connection.commit()
        finally:
            connection.close()
        return digest

    def test_v1_database_keeps_working_and_accepts_connectivity(self):
        with tempfile.TemporaryDirectory(prefix="crystal-test-") as directory:
            path = Path(directory) / "v1.sqlite"
            digest = self.build_v1_database(path)
            with Archive(path) as archive:
                # The old CHECK constraint was widened in place.
                sql = archive.connection.execute(
                    "SELECT sql FROM sqlite_master WHERE name='tasks'"
                ).fetchone()[0]
                self.assertIn("connectivity-v1", sql)
                self.assertEqual(archive.connection.execute(
                    "PRAGMA user_version").fetchone()[0], 1)
                # Old revision, pointer and V1 result survive the migration.
                self.assertEqual(archive.current_revision("p"), 1)
                legacy = Results(archive).read(1)
                self.assertEqual(legacy["source"]["digest"], digest)
                self.assertEqual(legacy["volume"], F(2))
                # New connectivity tasks land in the same id sequence.
                rid = Tasks(archive).run_connectivity(
                    "p", 1, [query("u", "v", (0, 1, 0))])
                self.assertGreater(rid, 1)
                value = Results(archive).read(rid)
                self.assertTrue(value["answers"][0]["reachable"])
                # V1 geometry still publishes after the migration.
                later = Tasks(archive).run_geometry("p", 1)
                self.assertGreater(later, rid)
                self.assertEqual(Results(archive).read(later), legacy)
                self.assertEqual(archive.connection.execute(
                    "PRAGMA foreign_key_check").fetchall(), [])
            with Archive(path) as archive:
                # Reopening is stable: no repeated migration, no drift.
                self.assertEqual(Results(archive).read(1), legacy)
                self.assertEqual(Results(archive).read(rid), value)
                packet = json.loads(Results(archive).export(rid))
                self.assertEqual(packet["format"],
                                 "crystal-connectivity-evidence")
                self.assertEqual(json.loads(Results(archive).export(1))
                                 ["format"], "crystal-primitive-evidence")


class TestConnectivityCodec(unittest.TestCase):
    def sample_result(self):
        source = ("p", 1, DIGEST_ZERO, "connectivity-v1")
        answers = (
            ConnectivityAnswer("u", "u", (2, 0, 0), True, (("a", 1),)),
            ConnectivityAnswer("u", "u", (1, 0, 0), False, None),
            ConnectivityAnswer("u", "u", (0, 0, 0), True, ()),
        )
        return ConnectivityResult(source, answers)

    def test_roundtrip_and_mapping(self):
        result = self.sample_result()
        text = encode_connectivity_result(result)
        self.assertEqual(decode_connectivity_result(text), result)
        mapping = connectivity_result_to_mapping(result)
        self.assertEqual(set(mapping), {"source", "answers"})
        self.assertEqual(mapping["answers"][0]["coefficients"], {"a": 1})
        self.assertIsNone(mapping["answers"][1]["coefficients"])
        self.assertEqual(mapping["answers"][2]["coefficients"], {})
        self.assertEqual(
            decode_connectivity_result(json.dumps(json.loads(text))), result)

    def test_decode_rejects_malformed_payloads(self):
        result = self.sample_result()
        base = json.loads(encode_connectivity_result(result))
        variants = []
        value = copy.deepcopy(base)
        value["version"] = True
        variants.append(value)
        value = copy.deepcopy(base)
        value["result"]["source"]["calculation"] = "primitive-v1"
        variants.append(value)
        value = copy.deepcopy(base)
        value["result"]["answers"][0]["coefficients"] = {"a": 0}
        variants.append(value)
        value = copy.deepcopy(base)
        value["result"]["answers"][0]["coefficients"] = {"a": True}
        variants.append(value)
        value = copy.deepcopy(base)
        value["result"]["answers"][1]["coefficients"] = {}
        variants.append(value)
        value = copy.deepcopy(base)
        value["result"]["answers"][0]["reachable"] = 1
        variants.append(value)
        value = copy.deepcopy(base)
        value["result"]["answers"][0]["query"]["translation"] = [1, 0]
        variants.append(value)
        value = copy.deepcopy(base)
        del value["result"]["answers"][0]["query"]["v"]
        variants.append(value)
        for variant in variants:
            with self.subTest(variant=variant), self.assertRaises(ValueError):
                decode_connectivity_result(json.dumps(variant))

    def test_executor_rejects_forged_evidence(self):
        cell, digest = None, None
        with Archive() as archive:
            revision = archive.put("p", cell_with([("a", "u", "u", (2, 0, 0))]))
            fixed = archive.revisions.read_revision("p", revision)
            cell, digest = fixed.cell, fixed.digest
        source = ("p", 1, digest, "connectivity-v1")
        good = ConnectivityResult(source, (
            ConnectivityAnswer("u", "u", (2, 0, 0), True, (("a", 1),)),))
        from crystal.executor import validate_connectivity_result
        self.assertIs(validate_connectivity_result(cell, good), good)
        forged = [
            ConnectivityResult(source, (
                ConnectivityAnswer("u", "u", (2, 0, 0), True, (("a", 2),)),)),
            ConnectivityResult(source, (
                ConnectivityAnswer("u", "u", (2, 0, 0), True,
                                   (("absent", 1),)),)),
            ConnectivityResult(source, (
                ConnectivityAnswer("u", "u", (1, 0, 0), True, (("a", 1),)),)),
            ConnectivityResult(source, (
                ConnectivityAnswer("u", "absent", (0, 0, 0), True, ()),)),
        ]
        for result in forged:
            with self.subTest(result=result), self.assertRaises(ValueError):
                validate_connectivity_result(cell, result)
        with self.assertRaises(ValueError):
            execute_connectivity(cell, ("p", 1, DIGEST_ZERO,
                                        "connectivity-v1"), ())


class TestSmithNormalForm(unittest.TestCase):
    def multiply(self, a, b):
        return [[sum(x * y for x, y in zip(row, col))
                 for col in zip(*b)] for row in a]

    def test_diagonalization_and_solve(self):
        cases = [
            [[2, 3, 0], [0, 0, 0], [0, 0, 0]],
            [[2, 0, 0], [1, 3, 0], [0, 1, 5]],
            [[6, 10, 15], [0, 0, 0], [0, 0, 0]],
            [[2 * 10 ** 30, 2 * 10 ** 30 + 1], [0, 0], [0, 0]],
            [[0, 0], [0, 0], [0, 0]],
            [[4], [6], [10]],
        ]
        for matrix in cases:
            with self.subTest(matrix=matrix):
                left, diagonal, right = smith_normal_form(matrix)
                self.assertEqual(
                    self.multiply(self.multiply(left, matrix), right),
                    diagonal)
                size = min(len(diagonal), len(diagonal[0]))
                previous = 1
                for i in range(size):
                    entry = diagonal[i][i]
                    self.assertGreaterEqual(entry, 0)
                    if entry:
                        self.assertEqual(entry % previous, 0)
                        previous = entry
                    for j in range(len(diagonal[0])):
                        if j != i:
                            self.assertEqual(diagonal[i][j], 0)
                # Solvability is decided exactly, and solutions verify.
                for target in ([1, 0, 0], [2, 0, 0], [0, 1, 0], [6, 0, 0]):
                    solution = solve_integer_system(matrix, target)
                    if solution is None:
                        continue
                    self.assertEqual(
                        [sum(row[j] * solution[j]
                             for j in range(len(solution)))
                         for row in matrix], target)


if __name__ == "__main__":
    unittest.main()
