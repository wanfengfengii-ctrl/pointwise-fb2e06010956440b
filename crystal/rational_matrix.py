from fractions import Fraction
from .exact_numbers import read_rational


def read_vector3(value):
    if type(value) not in (list, tuple) or len(value) != 3:
        raise ValueError("expected vector of length three")
    return tuple(read_rational(x) for x in value)


def read_matrix3(value):
    if type(value) not in (list, tuple) or len(value) != 3:
        raise ValueError("expected three matrix rows")
    return tuple(read_vector3(row) for row in value)


def determinant3(a):
    return (a[0][0] * (a[1][1] * a[2][2] - a[1][2] * a[2][1])
            - a[0][1] * (a[1][0] * a[2][2] - a[1][2] * a[2][0])
            + a[0][2] * (a[1][0] * a[2][1] - a[1][1] * a[2][0]))


def dot3(a, b):
    return sum((a[i] * b[i] for i in range(3)), Fraction(0))


def matvec3(a, vector):
    return tuple(dot3(row, vector) for row in a)


def add3(a, b):
    return tuple(a[i] + b[i] for i in range(3))


def subtract3(a, b):
    return tuple(a[i] - b[i] for i in range(3))
