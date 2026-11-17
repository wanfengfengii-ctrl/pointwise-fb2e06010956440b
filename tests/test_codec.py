from dataclasses import replace
from fractions import Fraction as F
import json
import unittest
from crystal import Archive, Tasks, Results
from crystal.exact_numbers import read_rational, decode_rational, json_load
from crystal.revision_codec import (
    encode_revision_cell, decode_revision_cell, encode_raw_cell, decode_raw_cell)
from crystal.result_codec import encode_result, decode_result, result_to_mapping
from crystal.schema import initialize_schema
from tests.test_intake import skew_cell


class TestCodec(unittest.TestCase):
    def test_exact_numeric_boundaries(self):
        for value in (True,1.0,None,[1,0],[True,2],[1,False],"1/2"):
            with self.subTest(value=value), self.assertRaises(ValueError):
                read_rational(value)
        self.assertEqual(read_rational([2,-6]),F(-1,3))
        for value in ([2,6],[1,-3],[0,2],[1,True],1):
            with self.assertRaises(ValueError):
                decode_rational(value)
        for text in ('{"a":1,"a":2}','{"a":0.5}','{"a":NaN}'):
            with self.assertRaises(ValueError):
                json_load(text)

    def test_revision_and_result_roundtrip(self):
        with Archive() as archive:
            r = archive.put("p",skew_cell())
            fixed = archive.revisions.read_revision("p",r)
            text = encode_revision_cell(fixed.cell)
            self.assertEqual(encode_revision_cell(decode_revision_cell(text)),text)
            raw = encode_raw_cell(fixed.cell.raw_mapping())
            self.assertEqual(encode_raw_cell(decode_raw_cell(raw)),raw)
            rid = Tasks(archive).run_geometry("p",r)
            result = Results(archive).store.read_bundle(rid)
            encoded = encode_result(result)
            self.assertEqual(decode_result(encoded),result)
            self.assertEqual(result_to_mapping(result),Results(archive).read(rid))
            value = json.loads(encoded)
            del value["result"]["bonds"][0]["offset"]
            with self.assertRaises(ValueError):
                decode_result(json.dumps(value))
            value = json.loads(encoded); value["version"] = True
            with self.assertRaises(ValueError):
                decode_result(json.dumps(value))
            value = json.loads(encoded)
            value["result"]["volume"] = [4,2]
            with self.assertRaises(ValueError):
                decode_result(json.dumps(value))
            value = json.loads(text); del value["cell"]["sites"]
            with self.assertRaises(ValueError):
                decode_revision_cell(json.dumps(value))

    def test_incompatible_schema_preserves_data(self):
        with Archive() as archive:
            r = archive.put("p",skew_cell())
            before = archive.read_revision("p",r)
            archive.connection.execute("PRAGMA user_version=2")
            with self.assertRaises(ValueError):
                initialize_schema(archive.connection)
            self.assertEqual(archive.read_revision("p",r),before)
            self.assertEqual(archive.connection.execute(
                "PRAGMA user_version").fetchone()[0],2)
