from dataclasses import dataclass
from .result_codec import CALCULATION
from .connectivity_result import CONNECTIVITY_CALCULATION
from .connectivity_queries import (
    check_connectivity_limits, read_connectivity_queries)
from .executor import execute_geometry, execute_connectivity
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


@dataclass(frozen=True)
class ConnectivityTask:
    product_id: str
    revision: int
    digest: str
    calculation: str = CONNECTIVITY_CALCULATION

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

    def run_connectivity(self, product_id, revision, queries):
        fixed = self.revisions.read_revision(product_id, revision)
        canonical = read_connectivity_queries(
            queries, {site.id for site in fixed.cell.sites})
        check_connectivity_limits(fixed.cell)
        task = ConnectivityTask(fixed.product_id, fixed.revision, fixed.digest)
        result = execute_connectivity(fixed.cell, task.source, canonical)
        return self.store.publish_bundle(task, result)
