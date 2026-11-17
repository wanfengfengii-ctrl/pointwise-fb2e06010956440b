"""Exact integer linear algebra over the integers (stdlib only).

smith_normal_form computes unimodular U, V with U @ matrix @ V diagonal.
solve_integer_system decides whether matrix @ n == target has an integer
solution and returns one such vector, without floats, rationals or windows.
"""


def _identity(size):
    return [[1 if i == j else 0 for j in range(size)] for i in range(size)]


def smith_normal_form(matrix):
    """Return (left, diagonal, right) with left @ matrix @ right == diagonal.

    left and right are unimodular integer square matrices; diagonal has the
    same shape as matrix, is zero off the diagonal, and its diagonal entries
    are nonnegative with each dividing the next.
    """
    diagonal = [list(row) for row in matrix]
    rows = len(diagonal)
    cols = len(diagonal[0]) if rows else 0
    left = _identity(rows)
    right = _identity(cols)

    def swap_rows(i, j):
        diagonal[i], diagonal[j] = diagonal[j], diagonal[i]
        left[i], left[j] = left[j], left[i]

    def swap_cols(i, j):
        for row in diagonal:
            row[i], row[j] = row[j], row[i]
        for row in right:
            row[i], row[j] = row[j], row[i]

    def add_row(source, target, scale):
        diagonal[target] = [x + scale * y
                            for x, y in zip(diagonal[target], diagonal[source])]
        left[target] = [x + scale * y
                        for x, y in zip(left[target], left[source])]

    def add_col(source, target, scale):
        for row in diagonal:
            row[target] += scale * row[source]
        for row in right:
            row[target] += scale * row[source]

    def scale_row(index, factor):
        diagonal[index] = [factor * x for x in diagonal[index]]
        left[index] = [factor * x for x in left[index]]

    k = 0
    while k < rows and k < cols:
        pivot = None
        for i in range(k, rows):
            for j in range(k, cols):
                if diagonal[i][j] and (pivot is None or abs(diagonal[i][j])
                                       < abs(diagonal[pivot[0]][pivot[1]])):
                    pivot = (i, j)
        if pivot is None:
            break
        if pivot[0] != k:
            swap_rows(k, pivot[0])
        if pivot[1] != k:
            swap_cols(k, pivot[1])
        while True:
            settled = True
            for i in range(k + 1, rows):
                if diagonal[i][k]:
                    settled = False
                    quotient = diagonal[i][k] // diagonal[k][k]
                    add_row(k, i, -quotient)
                    if diagonal[i][k]:
                        swap_rows(i, k)
            for j in range(k + 1, cols):
                if diagonal[k][j]:
                    settled = False
                    quotient = diagonal[k][j] // diagonal[k][k]
                    add_col(k, j, -quotient)
                    if diagonal[k][j]:
                        swap_cols(j, k)
            if settled:
                divisible = True
                for i in range(k + 1, rows):
                    for j in range(k + 1, cols):
                        if diagonal[i][j] % diagonal[k][k]:
                            divisible = False
                            add_row(i, k, 1)
                            break
                    if not divisible:
                        break
                if divisible:
                    break
        if diagonal[k][k] < 0:
            scale_row(k, -1)
        k += 1
    return left, diagonal, right


def prepare_integer_system(matrix):
    """Precompute the Smith normal form of an integer system matrix.

    matrix is a list of equal-length integer rows. Returns a prepared
    system for solve_prepared.
    """
    rows = len(matrix)
    if any(len(row) != len(matrix[0]) for row in matrix):
        raise ValueError("ragged integer system")
    cols = len(matrix[0]) if rows else 0
    left, diagonal, right = smith_normal_form(matrix)
    return (rows, cols, left, diagonal, right)


def solve_prepared(prepared, target):
    """Return an integer vector n with matrix @ n == target, or None.

    Exact divisibility decides solvability; no floats or rationals.
    """
    rows, cols, left, diagonal, right = prepared
    if len(target) != rows:
        raise ValueError("integer system target has wrong length")
    transformed = [sum(left[i][j] * target[j] for j in range(rows))
                   for i in range(rows)]
    particular = [0] * cols
    for i in range(rows):
        pivot = diagonal[i][i] if i < cols else 0
        if pivot:
            if transformed[i] % pivot:
                return None
            particular[i] = transformed[i] // pivot
        elif transformed[i]:
            return None
    return [sum(right[i][j] * particular[j] for j in range(cols))
            for i in range(cols)]


def solve_integer_system(matrix, target):
    """Return an integer vector n with matrix @ n == target, or None.

    matrix is a list of equal-length integer rows; target is an integer
    vector of the matching length. Exact divisibility decides solvability.
    """
    return solve_prepared(prepare_integer_system(matrix), target)
