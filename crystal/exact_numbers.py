from fractions import Fraction
import json


def read_integer(value, name="integer"):
    if type(value) is not int:
        raise ValueError(f"{name} must be an integer, excluding bool")
    return value


def read_rational(value):
    if type(value) is int:
        return Fraction(value)
    if isinstance(value, Fraction):
        return value
    if type(value) in (list, tuple) and len(value) == 2:
        n, d = (read_integer(x) for x in value)
        if d:
            return Fraction(n, d)
    raise ValueError("expected integer, Fraction or [numerator, denominator]")


def encode_rational(value):
    q = read_rational(value)
    return [q.numerator, q.denominator]


def decode_rational(value):
    if type(value) is not list or len(value) != 2:
        raise ValueError("encoded rational must be an integer pair")
    q = read_rational(value)
    if encode_rational(q) != value:
        raise ValueError("noncanonical rational encoding")
    return q


def exact_keys(value, keys):
    if type(value) is not dict or set(value) != set(keys):
        raise ValueError("missing or unexpected fields")
    return value


def read_identity(value):
    if type(value) is not str or not value:
        raise ValueError("identity must be a nonempty string")
    return value


def json_dump(value):
    return json.dumps(value, ensure_ascii=True, sort_keys=True,
                      separators=(",", ":"), allow_nan=False)


def json_load(text):
    def pairs(items):
        result = {}
        for key, value in items:
            if key in result:
                raise ValueError("duplicate JSON member")
            result[key] = value
        return result

    def reject(value):
        raise ValueError("noninteger JSON number: " + value)

    try:
        return json.loads(text, object_pairs_hook=pairs, parse_float=reject,
                          parse_constant=reject)
    except (TypeError, json.JSONDecodeError) as exc:
        raise ValueError("invalid exact JSON") from exc
