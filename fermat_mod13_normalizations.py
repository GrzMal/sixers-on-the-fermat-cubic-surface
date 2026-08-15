# fermat_mod13_normalizations.py
#
# This script checks the modulo 13 determinant behavior for all ordered
# normalizations of the two Fermat sixer representatives.
#
# Input:
#   orbit 54 representative: (0, 4, 10, 12, 20, 25)
#   orbit 18 representative: (0, 4, 8, 10, 14, 15)
#
# For each representative, we consider every ordered triple
#
#       (L0, Linf, L1)
#
# of distinct lines in the sixer.  There are 6*5*4 = 120 such triples.
#
# We work over F_13 and reduce omega by
#
#       omega -> 3,
#
# since 3^2 + 3 + 1 = 0 mod 13.
#
# For each normalization, we compute the determinants of the ten pairwise
# differences among
#
#       0, I_2, M_2, M_3, M_4.
#
# The output records, for each orbit, how many normalizations have a given
# number of nonsquare determinants modulo 13.


from itertools import permutations
from collections import Counter


p = 13
omega = 3
roots = [1, omega, omega * omega % p]

pairings = [
    ((0, 1), (2, 3)),
    ((0, 2), (1, 3)),
    ((0, 3), (1, 2)),
]


def inv(a):
    return pow(a % p, p - 2, p)


def rref(matrix):
    """
    Reduced row echelon form over F_13.
    """
    mat = [[x % p for x in row] for row in matrix]
    rows = len(mat)
    cols = len(mat[0])
    pivot_columns = []
    pivot_row = 0

    for col in range(cols):
        pivot = None

        for row in range(pivot_row, rows):
            if mat[row][col] != 0:
                pivot = row
                break

        if pivot is None:
            continue

        mat[pivot_row], mat[pivot] = mat[pivot], mat[pivot_row]

        pivot_inverse = inv(mat[pivot_row][col])
        mat[pivot_row] = [
            pivot_inverse * x % p
            for x in mat[pivot_row]
        ]

        for row in range(rows):
            if row != pivot_row and mat[row][col] != 0:
                factor = mat[row][col]
                mat[row] = [
                    (mat[row][j] - factor * mat[pivot_row][j]) % p
                    for j in range(cols)
                ]

        pivot_columns.append(col)
        pivot_row += 1

        if pivot_row == rows:
            break

    return mat, pivot_columns


def nullspace(matrix):
    """
    Basis for the nullspace of a matrix over F_13.
    """
    reduced, pivot_columns = rref(matrix)
    cols = len(reduced[0])
    free_columns = [j for j in range(cols) if j not in pivot_columns]

    basis = []

    for free_col in free_columns:
        vector = [0 for _ in range(cols)]
        vector[free_col] = 1

        for row, pivot_col in enumerate(pivot_columns):
            vector[pivot_col] = -reduced[row][free_col] % p

        basis.append(vector)

    return basis


def matrix_from_columns(columns):
    """
    Convert a list of column vectors into a matrix.
    """
    return [
        [columns[j][i] % p for j in range(len(columns))]
        for i in range(len(columns[0]))
    ]


def matrix_multiply(A, B):
    return [
        [
            sum(A[i][k] * B[k][j] for k in range(len(B))) % p
            for j in range(len(B[0]))
        ]
        for i in range(len(A))
    ]


def identity_matrix(n):
    return [
        [1 if i == j else 0 for j in range(n)]
        for i in range(n)
    ]


def matrix_inverse(A):
    """
    Inverse of a square matrix over F_13.
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

    row1 = [0, 0, 0, 0]
    row2 = [0, 0, 0, 0]

    row1[i] = 1
    row1[j] = a

    row2[k] = 1
    row2[ell] = b

    return [row1, row2]


def line_basis(line_index):
    """
    Return a basis for the two-dimensional vector space defining the line.
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


def normalize_sixer(representative, triple):
    """
    Normalize a sixer using the ordered triple

        (L0, Linf, L1).

    The normalization is M_i -> M_1^{-1} M_i.
    """
    L0, Linf, L1 = triple

    remaining = [
        line_index
        for line_index in representative
        if line_index not in triple
    ]

    M1 = graph_matrix(L0, Linf, L1)
    M1_inverse = matrix_inverse(M1)

    return [
        matrix_multiply(M1_inverse, graph_matrix(L0, Linf, line_index))
        for line_index in remaining
    ]


def matrix_subtract(A, B):
    return [
        [
            (A[i][j] - B[i][j]) % p
            for j in range(2)
        ]
        for i in range(2)
    ]


def determinant(A):
    return (A[0][0] * A[1][1] - A[0][1] * A[1][0]) % p


Z = [
    [0, 0],
    [0, 0]
]

Id = [
    [1, 0],
    [0, 1]
]

squares = {
    a * a % p
    for a in range(1, p)
}


def nonsquare_count_for_normalization(representative, triple):
    """
    Return the number of nonsquare determinants among the ten pairwise
    differences for the normalization determined by triple.
    """
    normalized_matrices = normalize_sixer(representative, triple)
    all_matrices = [Z, Id] + normalized_matrices

    values = []

    for i in range(len(all_matrices)):
        for j in range(i + 1, len(all_matrices)):
            delta = matrix_subtract(all_matrices[i], all_matrices[j])
            det_delta = determinant(delta)

            if det_delta == 0:
                raise ValueError("A determinant vanished.")

            values.append(det_delta)

    return sum(
        1 for value in values
        if value not in squares
    )


def normalization_profile(representative):
    """
    For all 120 ordered choices of (L0, Linf, L1), record the number of
    nonsquare determinants.
    """
    counter = Counter()

    for triple in permutations(representative, 3):
        count = nonsquare_count_for_normalization(representative, triple)
        counter[count] += 1

    return dict(sorted(counter.items()))


representative_54 = (0, 4, 10, 12, 20, 25)
representative_18 = (0, 4, 8, 10, 14, 15)


print("orbit 54 nonsquare determinant counts:")
print(normalization_profile(representative_54))

print()

print("orbit 18 nonsquare determinant counts:")
print(normalization_profile(representative_18))

