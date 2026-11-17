"""Connectivity result value objects and their exact JSON codec.

The stored payload is {"version": 1, "result": {...}}; the result mirrors
the read mapping exactly: source plus answers in input order. Coefficients
are the net forward uses of each original bond id, kept as a sorted tuple
of (bond id, nonzero strict integer) pairs, or None when unreachable.
"""
from dataclasses import dataclass
from .exact_numbers import (
    exact_keys, read_identity, read_integer, json_dump, json_load)
from .bonds import read_integer_offset

CONNECTIVITY_CALCULATION = "connectivity-v1"
CONNECTIVITY_RESULT_FORMAT_VERSION = 1


def connectivity_source_to_mapping(source):
    product, revision, digest, calculation = source
    read_identity(product)
    read_integer(revision, "revision")
    if revision < 1 or type(digest) is not str or len(digest) != 64:
        raise ValueError("invalid source binding")
    if any(c not in "0123456789abcdef" for c in digest):
        raise ValueError("invalid source digest")
    if calculation != CONNECTIVITY_CALCULATION:
        raise ValueError("unsupported calculation")
    return {"product_id": product, "revision": revision,
            "digest": digest, "calculation": calculation}


def read_connectivity_source(value):
    exact_keys(value, ("product_id", "revision", "digest", "calculation"))
    source = tuple(value[k] for k in
                   ("product_id", "revision", "digest", "calculation"))
    connectivity_source_to_mapping(source)
    return source


def read_coefficients(value):
    if type(value) is not dict:
        raise ValueError("coefficients must be a bond id mapping")
    pairs = []
    for key, coefficient in value.items():
        read_identity(key)
        read_integer(coefficient, "coefficient")
        if coefficient == 0:
            raise ValueError("coefficients omit zero terms")
        pairs.append((key, coefficient))
    return tuple(sorted(pairs))


def canonical_coefficients(value):
    if type(value) is not dict:
        try:
            value = dict(value)
        except (TypeError, ValueError) as exc:
            raise ValueError(
                "coefficients must be a bond id mapping") from exc
    return read_coefficients(value)


@dataclass(frozen=True)
class ConnectivityAnswer:
    u: str
    v: str
    translation: tuple
    reachable: bool
    coefficients: object  # sorted (bond id, int) pairs, or None

    def __post_init__(self):
        read_identity(self.u)
        read_identity(self.v)
        object.__setattr__(self, "translation",
                           read_integer_offset(self.translation))
        if type(self.reachable) is not bool:
            raise ValueError("reachable must be a bool")
        if self.reachable:
            if self.coefficients is None:
                raise ValueError("reachable answer requires coefficients")
            object.__setattr__(self, "coefficients",
                               canonical_coefficients(self.coefficients))
        elif self.coefficients is not None:
            raise ValueError("unreachable answer carries no coefficients")


@dataclass(frozen=True)
class ConnectivityResult:
    source: tuple
    answers: tuple

    def __post_init__(self):
        connectivity_source_to_mapping(self.source)
        object.__setattr__(self, "answers", tuple(self.answers))
        for answer in self.answers:
            if not isinstance(answer, ConnectivityAnswer):
                raise ValueError("answers must be connectivity answers")


def connectivity_result_to_mapping(result):
    return {
        "source": connectivity_source_to_mapping(result.source),
        "answers": [
            {"query": {"u": answer.u, "v": answer.v,
                       "translation": list(answer.translation)},
             "reachable": answer.reachable,
             "coefficients": (None if answer.coefficients is None
                              else dict(answer.coefficients))}
            for answer in result.answers],
    }


def encode_connectivity_result(result):
    text = json_dump({"version": CONNECTIVITY_RESULT_FORMAT_VERSION,
                      "result": connectivity_result_to_mapping(result)})
    decode_connectivity_result(text)
    return text


def decode_connectivity_result(text):
    payload = json_load(text)
    exact_keys(payload, ("version", "result"))
    if (type(payload["version"]) is not int
            or payload["version"] != CONNECTIVITY_RESULT_FORMAT_VERSION):
        raise ValueError("unsupported result format")
    value = exact_keys(payload["result"], ("source", "answers"))
    source = read_connectivity_source(value["source"])
    if type(value["answers"]) is not list:
        raise ValueError("incomplete answers ledger")
    answers = []
    for item in value["answers"]:
        exact_keys(item, ("query", "reachable", "coefficients"))
        query = exact_keys(item["query"], ("u", "v", "translation"))
        reachable = item["reachable"]
        if type(reachable) is not bool:
            raise ValueError("reachable must be a bool")
        raw = item["coefficients"]
        if not reachable and raw is not None:
            raise ValueError("unreachable answer carries no coefficients")
        coefficients = None if raw is None else read_coefficients(raw)
        answers.append(ConnectivityAnswer(
            query["u"], query["v"],
            read_integer_offset(query["translation"]),
            reachable, coefficients))
    return ConnectivityResult(source, tuple(answers))
