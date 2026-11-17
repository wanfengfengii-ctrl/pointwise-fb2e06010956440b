from .exact_numbers import exact_keys, json_dump, json_load
from .input_model import snapshot_raw_input, cell_to_mapping
from .intake import accept_cell

REVISION_FORMAT_VERSION = 1


def encode_revision_cell(cell):
    return json_dump({"version": REVISION_FORMAT_VERSION,
                      "cell": cell_to_mapping(cell)})


def decode_revision_cell(text):
    payload = json_load(text)
    exact_keys(payload, ("version", "cell"))
    if type(payload["version"]) is not int or payload["version"] != 1:
        raise ValueError("unsupported revision format")
    cell, _ = accept_cell(payload["cell"])
    if cell_to_mapping(cell) != payload["cell"]:
        raise ValueError("revision cell is not canonical")
    return cell


def encode_raw_cell(value):
    return json_dump(snapshot_raw_input(value))


def decode_raw_cell(text):
    raw = snapshot_raw_input(json_load(text))
    accept_cell(raw)
    return raw
