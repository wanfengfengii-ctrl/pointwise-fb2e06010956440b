from .database import open_database, close_database
from .intake import accept_cell, validate_product_id
from .input_model import cell_to_mapping
from .revision_store import RevisionStore


class Archive:
    def __init__(self, path=":memory:"):
        self.connection = open_database(path)
        self.revisions = RevisionStore(self.connection)
        self._closed = False

    def put(self, product_id, cell):
        validate_product_id(product_id)
        accepted, digest = accept_cell(cell)
        return self.revisions.append_revision(product_id, accepted, digest)

    def read_revision(self, product_id, revision):
        return cell_to_mapping(
            self.revisions.read_revision(product_id, revision).cell)

    def current_revision(self, product_id):
        return self.revisions.current_revision(product_id)

    def close(self):
        if not self._closed:
            close_database(self.connection)
            self._closed = True

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        self.close()
