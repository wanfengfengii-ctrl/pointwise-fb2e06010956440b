from .result_store import ResultStore
from .result_codec import result_to_mapping
from .export import build_export_packet, encode_export_packet


class Results:
    def __init__(self, archive):
        self.store = ResultStore(archive.connection)
        self.revisions = archive.revisions

    def read(self, result_id):
        return result_to_mapping(self.store.read_bundle(result_id))

    def export(self, result_id):
        result = self.store.read_bundle(result_id)
        revision = self.revisions.read_revision(
            result.source[0], result.source[1])
        return encode_export_packet(build_export_packet(revision, result))
