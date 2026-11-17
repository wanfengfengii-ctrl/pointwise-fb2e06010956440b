from dataclasses import dataclass
from .result_codec import CALCULATION
from .executor import execute_geometry
from .result_store import ResultStore


@dataclass(frozen=True)
class GeometryTask:
    product_id: str
    revision: int
    digest: str
    calculation: str = CALCULATION

    @property
    def source(self):
        return (self.product_id, self.revision, self.digest, self.calculation)


class Tasks:
    def __init__(self, archive):
        self.revisions = archive.revisions
        self.store = ResultStore(archive.connection)

    def run_geometry(self, product_id, revision):
        fixed = self.revisions.read_revision(product_id, revision)
        task = GeometryTask(fixed.product_id, fixed.revision, fixed.digest)
        result = execute_geometry(fixed.cell, task.source)
        return self.store.publish_bundle(task, result)
