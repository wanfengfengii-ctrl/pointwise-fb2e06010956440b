from dataclasses import replace
import sqlite3
import unittest
from crystal import Archive, Tasks, Results
from crystal.result_codec import encode_result
from crystal.result_store import ResultStore
from crystal.tasks import GeometryTask
from tests.test_intake import skew_cell


def table_counts(connection):
    return tuple(connection.execute("SELECT COUNT(*) FROM " + table).fetchone()[0]
                 for table in ("tasks", "results", "result_bodies"))


class TestTasks(unittest.TestCase):
    def test_publication_failure(self):
        with Archive() as archive:
            revision = archive.put("p",skew_cell())
            tasks, results = Tasks(archive), Results(archive)
            old_id = tasks.run_geometry("p",revision)
            old = results.read(old_id)
            before = table_counts(archive.connection)
            archive.connection.execute(
                "CREATE TEMP TRIGGER reject_body BEFORE INSERT ON result_bodies "
                "BEGIN SELECT RAISE(ABORT,'body write rejected'); END")
            with self.assertRaisesRegex(sqlite3.IntegrityError,"body write rejected"):
                tasks.run_geometry("p",revision)
            self.assertEqual(table_counts(archive.connection),before)
            self.assertEqual(results.read(old_id),old)
            self.assertEqual(archive.read_revision("p",revision)["unit"],"nm")
            archive.connection.execute("DROP TRIGGER reject_body")
            new_id = tasks.run_geometry("p",revision)
            self.assertNotEqual(new_id,old_id)
            self.assertEqual(results.read(new_id),old)
            self.assertEqual(table_counts(archive.connection),(2,2,2))
            self.assertEqual(archive.connection.execute(
                "PRAGMA foreign_key_check").fetchall(),[])

    def test_publication_success_and_missing(self):
        with Archive() as archive:
            tasks, results = Tasks(archive), Results(archive)
            for product, revision in (("missing",1),("p",1)):
                with self.assertRaises(KeyError):
                    tasks.run_geometry(product,revision)
            self.assertEqual(table_counts(archive.connection),(0,0,0))
            r = archive.put("p",skew_cell())
            rid = tasks.run_geometry("p",r)
            self.assertEqual(table_counts(archive.connection),(1,1,1))
            self.assertEqual(results.read(rid)["source"]["revision"],r)
            with self.assertRaises(KeyError):
                results.read(rid+1)
            with self.assertRaises(KeyError):
                results.export(rid+1)
            for argument in (True,0,1.5):
                with self.assertRaises(ValueError):
                    results.read(argument)
            archive.connection.execute("DELETE FROM result_bodies WHERE result_id=?",
                                       (rid,))
            with self.assertRaises(ValueError):
                results.read(rid)

    def test_binding_mismatch_is_rejected(self):
        with Archive() as archive:
            r = archive.put("p",skew_cell())
            rid = Tasks(archive).run_geometry("p",r)
            store = ResultStore(archive.connection)
            result = store.read_bundle(rid)
            task = GeometryTask("other",r,result.source[2])
            with self.assertRaises(ValueError):
                store.publish_bundle(task,result)
            changed = replace(result,source=("other",r,result.source[2],"primitive-v1"))
            archive.connection.execute(
                "UPDATE result_bodies SET payload=? WHERE result_id=?",
                (encode_result(changed),rid))
            with self.assertRaises(ValueError):
                Results(archive).read(rid)
            self.assertEqual(table_counts(archive.connection),(1,1,1))
