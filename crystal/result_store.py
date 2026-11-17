from .database import transaction
from .exact_numbers import read_integer
from .revision_store import RevisionStore
from .result_codec import (
    CALCULATION, encode_result, decode_result, source_to_mapping)
from .connectivity_result import (
    CONNECTIVITY_CALCULATION, connectivity_source_to_mapping,
    encode_connectivity_result, decode_connectivity_result)


def validate_any_source(source):
    if len(source) == 4 and source[3] == CONNECTIVITY_CALCULATION:
        return connectivity_source_to_mapping(source)
    return source_to_mapping(source)


def encode_any_result(result):
    if result.source[3] == CONNECTIVITY_CALCULATION:
        return encode_connectivity_result(result)
    return encode_result(result)


def decode_any_result(calculation, payload):
    if calculation == CALCULATION:
        return decode_result(payload)
    if calculation == CONNECTIVITY_CALCULATION:
        return decode_connectivity_result(payload)
    raise ValueError("unsupported calculation")


class ResultStore:
    def __init__(self, connection):
        self.connection = connection
        self.revisions = RevisionStore(connection)

    def publish_bundle(self, task, result):
        source = task.source
        validate_any_source(source)
        if result.source != source:
            raise ValueError("task/result source mismatch")
        payload = encode_any_result(result)
        with transaction(self.connection):
            revision = self.revisions.read_revision(source[0], source[1])
            if revision.digest != source[2]:
                raise ValueError("task/revision source mismatch")
            cursor = self.connection.execute(
                "INSERT INTO tasks(product_id,revision,digest,calculation) "
                "VALUES (?,?,?,?)", source)
            task_id = cursor.lastrowid
            cursor = self.connection.execute(
                "INSERT INTO results(task_id) VALUES (?)", (task_id,))
            result_id = cursor.lastrowid
            self.connection.execute(
                "INSERT INTO result_bodies(result_id,payload) VALUES (?,?)",
                (result_id, payload))
        return result_id

    def read_bundle(self, result_id):
        read_integer(result_id, "result_id")
        if result_id < 1:
            raise ValueError("result_id must be positive")
        row = self.connection.execute(
            "SELECT t.product_id,t.revision,t.digest,t.calculation,b.payload "
            "FROM results r JOIN tasks t ON t.task_id=r.task_id "
            "LEFT JOIN result_bodies b ON b.result_id=r.result_id "
            "WHERE r.result_id=?", (result_id,)).fetchone()
        if row is None:
            raise KeyError(result_id)
        if row[4] is None:
            raise ValueError("result body is missing")
        result = decode_any_result(row[3], row[4])
        if result.source != tuple(row[:4]):
            raise ValueError("result header/body binding mismatch")
        revision = self.revisions.read_revision(row[0], row[1])
        if revision.digest != row[2]:
            raise ValueError("result/revision binding mismatch")
        return result
