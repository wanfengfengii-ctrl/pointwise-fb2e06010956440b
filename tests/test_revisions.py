from fractions import Fraction as F
from pathlib import Path
import sqlite3
import tempfile
import unittest
from crystal import Archive, Tasks, Results
from tests.test_intake import skew_cell


class TestRevisions(unittest.TestCase):
    def test_frozen_result_reopen(self):
        with tempfile.TemporaryDirectory(prefix="crystal-test-") as directory:
            path = Path(directory)/"archive.sqlite"
            with Archive(path) as archive:
                r1 = archive.put("p",skew_cell())
                id1 = Tasks(archive).run_geometry("p",r1)
                original = Results(archive).read(id1)
                changed = skew_cell(); changed["lattice"][0][0] = 2
                r2 = archive.put("p",changed)
                id2 = Tasks(archive).run_geometry("p",r2)
                again = Tasks(archive).run_geometry("p",r1)
                self.assertEqual(Results(archive).read(again),original)
                self.assertEqual(Results(archive).read(id2)["volume"],F(4))
                self.assertEqual(Results(archive).read(id2)["bonds"][0]
                                 ["length_squared"],F(85,16))
                returned = Results(archive).read(id1)
                returned["source"]["revision"] = r2
                returned["positions"][0]["cartesian"][0] = F(99)
                returned["bonds"][0]["offset"][0] = 99
                source = archive.read_revision("p",r1)
                source["lattice"][0][0][0] = 99
                source["sites"].clear()
                self.assertEqual(Results(archive).read(id1),original)
            with Archive(path) as archive:
                self.assertEqual(archive.current_revision("p"),r2)
                self.assertEqual(Results(archive).read(id1),original)
                self.assertEqual(archive.read_revision("p",r1)["lattice"][0][0],
                                 [1,1])
                self.assertEqual(Results(archive).read(id2)["source"]["revision"],r2)
                self.assertEqual(archive.connection.execute(
                    "PRAGMA foreign_key_check").fetchall(),[])

    def test_correction_rollback(self):
        with Archive() as archive:
            old = archive.put("p",skew_cell())
            result_id = Tasks(archive).run_geometry("p",old)
            result = Results(archive).read(result_id)
            archive.connection.execute(
                "CREATE TEMP TRIGGER reject_pointer BEFORE UPDATE ON products "
                "BEGIN SELECT RAISE(ABORT,'pointer write rejected'); END")
            changed = skew_cell(); changed["lattice"][0][0] = 3
            with self.assertRaisesRegex(sqlite3.IntegrityError,"pointer write rejected"):
                archive.put("p",changed)
            self.assertEqual(archive.current_revision("p"),old)
            self.assertEqual(archive.connection.execute(
                "SELECT COUNT(*) FROM revisions").fetchone()[0],1)
            with self.assertRaises(KeyError):
                archive.read_revision("p",old+1)
            self.assertEqual(Results(archive).read(result_id),result)
            archive.connection.execute("DROP TRIGGER reject_pointer")
            self.assertEqual(archive.put("p",changed),old+1)

    def test_first_revision_rollback(self):
        with Archive() as archive:
            archive.connection.execute(
                "CREATE TEMP TRIGGER reject_revision BEFORE INSERT ON revisions "
                "BEGIN SELECT RAISE(ABORT,'revision write rejected'); END")
            with self.assertRaises(sqlite3.IntegrityError):
                archive.put("new",skew_cell())
            for table in ("products","revisions"):
                self.assertEqual(archive.connection.execute(
                    "SELECT COUNT(*) FROM "+table).fetchone()[0],0)
            archive.connection.execute("DROP TRIGGER reject_revision")
            self.assertEqual(archive.put("new",skew_cell()),1)
