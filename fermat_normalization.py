# fermat_normalization.py
#
# This script normalizes the two Fermat sixer representatives produced by
# fermat_sixers.py.
#
# Input:
#   orbit 54 representative: (0, 4, 10, 12, 20, 25)
#   orbit 18 representative: (0, 4, 8, 10, 14, 15)
#
# For each representative we choose:
#
#   L0    = first line,
#   Linf  = second line,
#   L1    = third line.
#
# We then put L0 and Linf in standard position, write the remaining lines as
# graphs of linear maps U -> W, and normalize by M_i -> M_1^{-1} M_i.
#
# Output:
#   the matrices M2, M3, M4 for each of the two Fermat orbits.


from fractions import Fraction


class Eisenstein:
    """
    Element a + b*omega of Q(omega), where omega^2 + omega + 1 = 0.
    """

    def __init__(self, a=0, b=0):
        self.a = Fraction(a)
        self.b = Fraction(b)

    def __add__(self, other):
        other = to_eisenstein(other)
        return Eisenstein(self.a + other.a, self.b + other.b)

    def __radd__(self, other):
        return self + other

    def __neg__(self):
        return Eisenstein(-self.a, -self.b)

    def __sub__(self, other):
        return self + (-to_eisenstein(other))

    def __rsub__(self, other):
        return to_eisenstein(other) - self

    def __mul__(self, other):
        other = to_eisenstein(other)

        # omega^2 = -1 - omega
        return Eisenstein(
            self.a * other.a - self.b * other.b,
            self.a * other.b + self.b * other.a - self.b * other.b
        )

    def __rmul__(self, other):
        return self * other

    def conjugate(self):
        # omega maps to omega^2 = -1 - omega.
        return Eisenstein(self.a - self.b, -self.b)

    def norm(self):
        value = self * self.conjugate()
        if value.b != 0:
            raise ValueError("Norm should be rational.")
        return value.a

    def inverse(self):
        nrm = self.norm()
        if nrm == 0:
            raise ZeroDivisionError("division by zero")
        conjugate = self.conjugate()
        return Eisenstein(conjugate.a / nrm, conjugate.b / nrm)

    def __truediv__(self, other):
        return self * to_eisenstein(other).inverse()

    def __eq__(self, other):
        other = to_eisenstein(other)
        return self.a == other.a and self.b == other.b

    def is_zero(self):
        return self.a == 0 and self.b == 0

    def __repr__(self):
        def show_fraction(q):
            if q.denominator == 1:
                return str(q.numerator)
            return f"{q.numerator}/{q.denominator}"

        a = self.a
        b = self.b

        if b == 0:
            return show_fraction(a)

        if a == 0:
            if b == 1:
                return "omega"
            if b == -1:
                return "-omega"
            return f"{show_fraction(b)}*omega"

        sign = "+" if b > 0 else "-"
        coeff = abs(b)

        if coeff == 1:
            omega_part = "omega"
        else:
            omega_part = f"{show_fraction(coeff)}*omega"

        return f"{show_fraction(a)}{sign}{omega_part}"


def to_eisenstein(value):
    if isinstance(value, Eisenstein):
        return value
    return Eisenstein(value, 0)


zero = Eisenstein(0)
one = Eisenstein(1)
omega = Eisenstein(0, 1)
omega2 = omega * omega


def rref(matrix):
    """
    Reduced row echelon form over Q(omega).
    """
    mat = [[to_eisenstein(x) for x in row] for row in matrix]
    rows = len(mat)
    cols = len(mat[0])
    pivot_columns = []
    pivot_row = 0

    for col in range(cols):
        pivot = None

        for row in range(pivot_row, rows):
            if not mat[row][col].is_zero():
                pivot = row
                break

        if pivot is None:
            continue

        mat[pivot_row], mat[pivot] = mat[pivot], mat[pivot_row]

        pivot_inverse = mat[pivot_row][col].inverse()
        mat[pivot_row] = [x * pivot_inverse for x in mat[pivot_row]]

        for row in range(rows):
            if row != pivot_row and not mat[row][col].is_zero():
                factor = mat[row][col]
                mat[row] = [
                    mat[row][j] - factor * mat[pivot_row][j]
                    for j in range(cols)
                ]

        pivot_columns.append(col)
        pivot_row += 1

        if pivot_row == rows:
            break

    return mat, pivot_columns


def nullspace(matrix):
    """
    Basis for the nullspace of a matrix over Q(omega).
    """
    reduced, pivot_columns = rref(matrix)
    cols = len(reduced[0])
    free_columns = [j for j in range(cols) if j not in pivot_columns]

    basis = []

    for free_col in free_columns:
        vector = [zero for _ in range(cols)]
        vector[free_col] = one

        for row, pivot_col in enumerate(pivot_columns):
            vector[pivot_col] = -reduced[row][free_col]

        basis.append(vector)

    return basis


def identity_matrix(n):
    return [
        [one if i == j else zero for j in range(n)]
        for i in range(n)
    ]


def matrix_multiply(A, B):
    return [
        [
            sum((A[i][k] * B[k][j] for k in range(len(B))), zero)
            for j in range(len(B[0]))
        ]
        for i in range(len(A))
    ]


def matrix_inverse(A):
    """
    Inverse of a square matrix over Q(omega).
    """
    n = len(A)
    augmented = [
        A[i][:] + identity_matrix(n)[i]
        for i in range(n)
    ]

    reduced, pivot_columns = rref(augmented)

    if pivot_columns[:n] != list(range(n)):
        raise ValueError("Matrix is not invertible.")

    return [row[n:] for row in reduced]


def matrix_from_columns(columns):
    """
    Convert a list of column vectors into a matrix.
    """
    return [
        [columns[j][i] for j in range(len(columns))]
        for i in range(len(columns[0]))
    ]


# The same ordering of the 27 Fermat lines as in fermat_sixers.py.

roots = [one, omega, omega2]

pairings = [
    ((0, 1), (2, 3)),
    ((0, 2), (1, 3)),
    ((0, 3), (1, 2)),
]


def line_equations(line_index):
    """
    Return the two equations defining the Fermat line with the given index.
    """
    pairing_index = line_index // 9
    remainder = line_index % 9

    root_index_a = remainder // 3
    root_index_b = remainder % 3

    a = roots[root_index_a]
    b = roots[root_index_b]

    (i, j), (k, ell) = pairings[pairing_index]

    row1 = [zero, zero, zero, zero]
    row2 = [zero, zero, zero, zero]

    row1[i] = one
    row1[j] = a

    row2[k] = one
    row2[ell] = b

    return [row1, row2]


def line_basis(line_index):
    """
    Return a basis for the two-dimensional vector subspace defining the line.
    """
    return nullspace(line_equations(line_index))


def graph_matrix(L0, Linf, line_index):
    """
    Write a line as a graph U -> W after choosing L0 = P(U) and Linf = P(W).
    """
    basis_U = line_basis(L0)
    basis_W = line_basis(Linf)
    basis_line = line_basis(line_index)

    change_of_basis = matrix_from_columns(basis_U + basis_W)
    line_matrix = matrix_from_columns(basis_line)

    coordinates = matrix_multiply(
        matrix_inverse(change_of_basis),
        line_matrix
    )

    upper = coordinates[:2]
    lower = coordinates[2:]

    return matrix_multiply(lower, matrix_inverse(upper))


def normalize_sixer(representative):
    """
    Normalize a sixer in the form

        0, infinity, I_2, M2, M3, M4.

    The first three lines of the representative are used as L0, Linf, L1.
    """
    L0 = representative[0]
    Linf = representative[1]

    graph_matrices = [
        graph_matrix(L0, Linf, line_index)
        for line_index in representative[2:]
    ]

    M1_inverse = matrix_inverse(graph_matrices[0])

    return [
        matrix_multiply(M1_inverse, M)
        for M in graph_matrices[1:]
    ]


def print_matrix(M):
    for row in M:
        print(row)


representative_54 = (0, 4, 10, 12, 20, 25)
representative_18 = (0, 4, 8, 10, 14, 15)


print("orbit 54 representative:")
print(representative_54)
print("chosen normalization:")
print("L0 =", representative_54[0],
      "Linf =", representative_54[1],
      "L1 =", representative_54[2])

M54 = normalize_sixer(representative_54)

for index, matrix in enumerate(M54, start=2):
    print()
    print(f"M{index} =")
    print_matrix(matrix)


print()
print("orbit 18 representative:")
print(representative_18)
print("chosen normalization:")
print("L0 =", representative_18[0],
      "Linf =", representative_18[1],
      "L1 =", representative_18[2])

M18 = normalize_sixer(representative_18)

for index, matrix in enumerate(M18, start=2):
    print()
    print(f"M{index} =")
    print_matrix(matrix)
