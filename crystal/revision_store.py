from dataclasses import dataclass, replace
from .database import transaction
from .exact_numbers import read_integer
from .intake import validate_product_id
from .input_digest import verify_cell_digest
from .revision_codec import (
    encode_revision_cell, decode_revision_cell, encode_raw_cell, decode_raw_cell)


@dataclass(frozen=True)
class CellRevision:
    product_id: str
    revision: int
    digest: str
    cell: object


class RevisionStore:
    def __init__(self, connection):
        self.connection = connection

    def append_revision(self, product_id, cell, digest):
        validate_product_id(product_id)
        verify_cell_digest(cell, digest)
        canonical = encode_revision_cell(cell)
        raw = encode_raw_cell(cell.raw_mapping())
        with transaction(self.connection):
            self.connection.execute(
                "INSERT OR IGNORE INTO products(product_id) VALUES (?)",
                (product_id,))
            row = self.connection.execute(
                "SELECT COALESCE(MAX(revision),0)+1 FROM revisions "
                "WHERE product_id=?", (product_id,)).fetchone()
            revision = row[0]
            self.connection.execute(
                "INSERT INTO revisions VALUES (?,?,?,?,?)",
                (product_id, revision, digest, canonical, raw))
            self.connection.execute(
                "UPDATE products SET current_revision=? WHERE product_id=?",
                (revision, product_id))
        return revision

    def read_revision(self, product_id, revision):
        validate_product_id(product_id)
        read_integer(revision, "revision")
        if revision < 1:
            raise ValueError("revision must be positive")
        row = self.connection.execute(
            "SELECT digest,canonical_json,raw_json FROM revisions "
            "WHERE product_id=? AND revision=?", (product_id, revision)).fetchone()
        if row is None:
            raise KeyError((product_id, revision))
        cell = decode_revision_cell(row[1])
        verify_cell_digest(cell, row[0])
        raw = decode_raw_cell(row[2])
        cell = replace(cell, raw_json=encode_raw_cell(raw))
        return CellRevision(product_id, revision, row[0], cell)

    def current_revision(self, product_id):
        validate_product_id(product_id)
        row = self.connection.execute(
            "SELECT current_revision FROM products WHERE product_id=?",
            (product_id,)).fetchone()
        if row is None or row[0] is None:
            raise KeyError(product_id)
        return row[0]
